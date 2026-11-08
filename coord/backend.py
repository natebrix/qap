"""Storage, queue and ledger for the coordinator (docs/distributed/plan.md
sections 3-4), behind one small interface with two implementations:

  local:DIR          a directory: files for the store, one SQLite database
                     (WAL, shared by the coordinator and the worker processes on
                     one machine) for the queue and the ledger
  aws:BUCKET[/PREFIX][?region=R&name=N]
                     S3 for the store, two SQS queues (units-large,
                     units-small) and two DynamoDB tables (runs, units) named
                     after N (default "qap"); `qapctl aws-setup` creates them

The three parts:

  Store   put/get of objects by key (unit states, results, the incumbent,
          control flags, totals, status, snapshots)
  Queue   at-least-once messages with a visibility timeout, the lease: a
          received message is invisible until its visibility runs out, the
          worker extends it while it works and deletes it when the unit is
          recorded; two lanes polled in order (large first)
  Ledger  the runs and units tables. Every unit update is conditional
          (expected status, attempt), so a late or duplicate worker cannot
          overwrite a newer attempt and nothing is merged twice; merge
          sequence numbers come from an atomic counter.
"""
import json
import os
import sqlite3
import threading
import time
import urllib.parse


def open_backend(spec):
    if spec.startswith("local:"):
        return LocalBackend(spec[len("local:"):])
    if spec.startswith("aws:"):
        return AwsBackend(spec[len("aws:"):])
    raise ValueError(f"unknown backend {spec!r} (local:DIR or aws:BUCKET[/PREFIX][?region=R&name=N])")


class Message:
    __slots__ = ("handle", "body", "lane", "receive_count")

    def __init__(self, handle, body, lane, receive_count):
        self.handle, self.body, self.lane, self.receive_count = handle, body, lane, receive_count


# columns kept outside the JSON document (queried or conditioned on)
COLUMNS = ("status", "attempt", "merged_seq")


# ------------------------------------------------------------------- local

class LocalBackend:
    def __init__(self, root):
        self.root = os.path.abspath(root)
        os.makedirs(os.path.join(self.root, "store"), exist_ok=True)
        self.db_path = os.path.join(self.root, "ledger.db")
        self._local = threading.local()
        with self._conn() as c:
            c.executescript("""
              CREATE TABLE IF NOT EXISTS runs (run_id TEXT PRIMARY KEY, doc TEXT, ver INTEGER);
              CREATE TABLE IF NOT EXISTS units (run_id TEXT, uid TEXT, status TEXT, attempt INTEGER,
                                                merged_seq INTEGER, doc TEXT, ver INTEGER,
                                                PRIMARY KEY (run_id, uid));
              CREATE INDEX IF NOT EXISTS units_status ON units (run_id, status);
              CREATE INDEX IF NOT EXISTS units_seq ON units (run_id, merged_seq);
              CREATE TABLE IF NOT EXISTS counters (run_id TEXT, name TEXT, value INTEGER,
                                                   PRIMARY KEY (run_id, name));
              CREATE TABLE IF NOT EXISTS messages (id INTEGER PRIMARY KEY AUTOINCREMENT, lane TEXT,
                                                   body TEXT, visible_at REAL, receive_count INTEGER);
              CREATE INDEX IF NOT EXISTS messages_vis ON messages (lane, visible_at);
            """)

    def _conn(self):
        # one connection per thread and process (worker slots are threads; a
        # forked process must not share its parent's connection)
        loc = self._local
        if getattr(loc, "db", None) is None or loc.pid != os.getpid():
            loc.db = sqlite3.connect(self.db_path, timeout=60, isolation_level=None)
            loc.db.execute("PRAGMA journal_mode=WAL")
            loc.db.execute("PRAGMA synchronous=NORMAL")
            loc.pid = os.getpid()
        return loc.db

    def _tx(self):
        c = self._conn()
        c.execute("BEGIN IMMEDIATE")
        return c

    # store
    def _path(self, key):
        p = os.path.normpath(os.path.join(self.root, "store", key))
        assert p.startswith(os.path.join(self.root, "store")), key
        return p

    def put_bytes(self, key, data):
        p = self._path(key)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        tmp = p + f".tmp{os.getpid()}"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, p)

    def put_file(self, key, path):
        with open(path, "rb") as f:
            self.put_bytes(key, f.read())

    def get_bytes(self, key):
        try:
            with open(self._path(key), "rb") as f:
                return f.read()
        except FileNotFoundError:
            return None

    def get_file(self, key, path):
        data = self.get_bytes(key)
        if data is None:
            return False
        with open(path, "wb") as f:
            f.write(data)
        return True

    def get_if_changed(self, key, tag):
        """(data, tag): data None when unchanged since tag or missing"""
        try:
            st = os.stat(self._path(key))
        except FileNotFoundError:
            return None, tag
        t = f"{st.st_mtime_ns}-{st.st_size}"
        if t == tag:
            return None, tag
        return self.get_bytes(key), t

    def exists(self, key):
        return os.path.exists(self._path(key))

    def delete(self, key):
        try:
            os.remove(self._path(key))
        except FileNotFoundError:
            pass

    def list(self, prefix):
        base = os.path.join(self.root, "store")
        out = []
        top = self._path(prefix) if prefix else base
        d = top if os.path.isdir(top) else os.path.dirname(top)
        for dirpath, _, files in os.walk(d):
            for f in files:
                k = os.path.relpath(os.path.join(dirpath, f), base)
                if k.startswith(prefix) and ".tmp" not in f:
                    out.append(k)
        return sorted(out)

    # queue
    def send(self, lane, body):
        with self._tx() as c:
            c.execute("INSERT INTO messages (lane, body, visible_at, receive_count) VALUES (?, ?, ?, 0)",
                      (lane, json.dumps(body), time.time()))

    def receive(self, lanes, visibility):
        now = time.time()
        c = self._tx()
        try:
            for lane in lanes:
                row = c.execute("SELECT id, body, receive_count FROM messages WHERE lane = ? AND visible_at <= ? "
                                "ORDER BY id LIMIT 1", (lane, now)).fetchone()
                if row:
                    c.execute("UPDATE messages SET visible_at = ?, receive_count = receive_count + 1 WHERE id = ?",
                              (now + visibility, row[0]))
                    c.execute("COMMIT")
                    return Message(row[0], json.loads(row[1]), lane, row[2] + 1)
            c.execute("COMMIT")
        except Exception:
            c.execute("ROLLBACK")
            raise
        return None

    def extend(self, msg, visibility):
        with self._tx() as c:
            c.execute("UPDATE messages SET visible_at = ? WHERE id = ?", (time.time() + visibility, msg.handle))

    def release(self, msg):
        self.extend(msg, 0)

    def delete_message(self, msg):
        with self._tx() as c:
            c.execute("DELETE FROM messages WHERE id = ?", (msg.handle,))

    def depth(self, lane):
        return self._conn().execute("SELECT COUNT(*) FROM messages WHERE lane = ? AND visible_at <= ?",
                                    (lane, time.time())).fetchone()[0]

    # ledger: runs
    def put_run(self, run_id, doc):
        with self._tx() as c:
            c.execute("INSERT OR REPLACE INTO runs (run_id, doc, ver) VALUES (?, ?, 0)", (run_id, json.dumps(doc)))

    def get_run(self, run_id):
        row = self._conn().execute("SELECT doc FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        return json.loads(row[0]) if row else None

    def update_run(self, run_id, fields):
        c = self._tx()
        try:
            row = c.execute("SELECT doc FROM runs WHERE run_id = ?", (run_id,)).fetchone()
            doc = json.loads(row[0])
            doc.update(fields)
            c.execute("UPDATE runs SET doc = ?, ver = ver + 1 WHERE run_id = ?", (json.dumps(doc), run_id))
            c.execute("COMMIT")
            return doc
        except Exception:
            c.execute("ROLLBACK")
            raise

    def next_seq(self, run_id, name="merge"):
        c = self._tx()
        try:
            c.execute("INSERT OR IGNORE INTO counters (run_id, name, value) VALUES (?, ?, 0)", (run_id, name))
            c.execute("UPDATE counters SET value = value + 1 WHERE run_id = ? AND name = ?", (run_id, name))
            v = c.execute("SELECT value FROM counters WHERE run_id = ? AND name = ?", (run_id, name)).fetchone()[0]
            c.execute("COMMIT")
            return v
        except Exception:
            c.execute("ROLLBACK")
            raise

    # ledger: units
    @staticmethod
    def _row(doc):
        return (doc.get("status"), doc.get("attempt", 0), doc.get("merged_seq"))

    def put_unit(self, run_id, doc):
        """insert if absent; True when inserted"""
        c = self._tx()
        try:
            cur = c.execute("INSERT OR IGNORE INTO units (run_id, uid, status, attempt, merged_seq, doc, ver) "
                            "VALUES (?, ?, ?, ?, ?, ?, 0)", (run_id, doc["uid"]) + self._row(doc) + (json.dumps(doc),))
            c.execute("COMMIT")
            return cur.rowcount == 1
        except Exception:
            c.execute("ROLLBACK")
            raise

    def get_unit(self, run_id, uid):
        row = self._conn().execute("SELECT doc FROM units WHERE run_id = ? AND uid = ?", (run_id, uid)).fetchone()
        return json.loads(row[0]) if row else None

    def update_unit(self, run_id, uid, fields, expect=None, check=None):
        """Set fields when every expect[k] matches (a value, or a tuple/list of
        allowed values) and check(doc) is true; returns the new doc or None."""
        c = self._tx()
        try:
            row = c.execute("SELECT doc FROM units WHERE run_id = ? AND uid = ?", (run_id, uid)).fetchone()
            if not row:
                c.execute("COMMIT")
                return None
            doc = json.loads(row[0])
            if not _matches(doc, expect) or (check and not check(doc)):
                c.execute("COMMIT")
                return None
            doc.update(fields)
            c.execute("UPDATE units SET status = ?, attempt = ?, merged_seq = ?, doc = ?, ver = ver + 1 "
                      "WHERE run_id = ? AND uid = ?", self._row(doc) + (json.dumps(doc), run_id, uid))
            c.execute("COMMIT")
            return doc
        except Exception:
            c.execute("ROLLBACK")
            raise

    def units_by_status(self, run_id, status):
        rows = self._conn().execute("SELECT doc FROM units WHERE run_id = ? AND status = ?", (run_id, status))
        return [json.loads(r[0]) for r in rows]

    def count_by_status(self, run_id):
        rows = self._conn().execute("SELECT status, COUNT(*) FROM units WHERE run_id = ? GROUP BY status", (run_id,))
        return {s: n for s, n in rows}

    def units_merged_after(self, run_id, seq):
        rows = self._conn().execute("SELECT doc FROM units WHERE run_id = ? AND merged_seq > ? ORDER BY merged_seq",
                                    (run_id, seq))
        return [json.loads(r[0]) for r in rows]

    def all_units(self, run_id):
        return [json.loads(r[0]) for r in self._conn().execute("SELECT doc FROM units WHERE run_id = ?", (run_id,))]


def _matches(doc, expect):
    for k, v in (expect or {}).items():
        if isinstance(v, (tuple, list, set)):
            if doc.get(k) not in v:
                return False
        elif doc.get(k) != v:
            return False
    return True


# --------------------------------------------------------------------- aws

class AwsBackend:
    """S3 + SQS + DynamoDB. Units: table <name>-units, key (run_id, uid), the
    document as a JSON string plus status, attempt, merged_seq and run_status
    (= run_id#status, partition key of the by_status index) as attributes;
    updates are optimistic (a version attribute) with the expectations checked
    on the read item. Runs and counters: table <name>-runs (key run_id;
    counters as items run_id#counter#name)."""

    def __init__(self, spec):
        import boto3
        u = urllib.parse.urlparse("aws://" + spec)
        q = dict(urllib.parse.parse_qsl(u.query))
        self.bucket = u.netloc
        self.prefix = u.path.strip("/")
        self.name = q.get("name", "qap")
        region = q.get("region") or os.environ.get("AWS_REGION") or os.environ.get("AWS_DEFAULT_REGION") or "us-east-1"
        self.region = region
        kw = {"region_name": region}
        if q.get("endpoint"):
            kw["endpoint_url"] = q["endpoint"]
        self.s3 = boto3.client("s3", **kw)
        self.sqs = boto3.client("sqs", **kw)
        self.ddb = boto3.client("dynamodb", **kw)
        self._qurl = {}

    # setup (also what Terraform will declare)
    def setup(self):
        try:
            self.s3.head_bucket(Bucket=self.bucket)
        except Exception:
            # outside us-east-1, S3 refuses a bucket without its region
            # (IllegalLocationConstraintException)
            loc = {} if self.region == "us-east-1" else \
                {"CreateBucketConfiguration": {"LocationConstraint": self.region}}
            self.s3.create_bucket(Bucket=self.bucket, **loc)
        for lane in ("large", "small"):
            self.sqs.create_queue(QueueName=f"{self.name}-units-{lane}",
                                  Attributes={"VisibilityTimeout": "1800", "MessageRetentionPeriod": "1209600"})
        existing = self.ddb.list_tables()["TableNames"]
        if f"{self.name}-runs" not in existing:
            self.ddb.create_table(TableName=f"{self.name}-runs", BillingMode="PAY_PER_REQUEST",
                                  KeySchema=[{"AttributeName": "run_id", "KeyType": "HASH"}],
                                  AttributeDefinitions=[{"AttributeName": "run_id", "AttributeType": "S"}])
        if f"{self.name}-units" not in existing:
            self.ddb.create_table(
                TableName=f"{self.name}-units", BillingMode="PAY_PER_REQUEST",
                KeySchema=[{"AttributeName": "run_id", "KeyType": "HASH"}, {"AttributeName": "uid", "KeyType": "RANGE"}],
                AttributeDefinitions=[{"AttributeName": a, "AttributeType": t} for a, t in
                                      (("run_id", "S"), ("uid", "S"), ("run_status", "S"), ("merged_seq", "N"))],
                GlobalSecondaryIndexes=[
                    {"IndexName": "by_status", "Projection": {"ProjectionType": "ALL"},
                     "KeySchema": [{"AttributeName": "run_status", "KeyType": "HASH"},
                                   {"AttributeName": "uid", "KeyType": "RANGE"}]},
                    {"IndexName": "by_seq", "Projection": {"ProjectionType": "ALL"},
                     "KeySchema": [{"AttributeName": "run_id", "KeyType": "HASH"},
                                   {"AttributeName": "merged_seq", "KeyType": "RANGE"}]}])

    # store
    def _key(self, key):
        return f"{self.prefix}/{key}" if self.prefix else key

    def put_bytes(self, key, data):
        self.s3.put_object(Bucket=self.bucket, Key=self._key(key), Body=data)

    def put_file(self, key, path):
        self.s3.upload_file(path, self.bucket, self._key(key))

    def get_bytes(self, key):
        try:
            return self.s3.get_object(Bucket=self.bucket, Key=self._key(key))["Body"].read()
        except self.s3.exceptions.NoSuchKey:
            return None

    def get_file(self, key, path):
        data = self.get_bytes(key)
        if data is None:
            return False
        with open(path, "wb") as f:
            f.write(data)
        return True

    def get_if_changed(self, key, tag):
        try:
            kw = {"Bucket": self.bucket, "Key": self._key(key)}
            if tag:
                kw["IfNoneMatch"] = tag
            r = self.s3.get_object(**kw)
            return r["Body"].read(), r["ETag"]
        except self.s3.exceptions.NoSuchKey:
            return None, tag
        except Exception as e:     # 304 Not Modified arrives as a ClientError
            if "304" in str(e) or "Not Modified" in str(e):
                return None, tag
            raise

    def exists(self, key):
        try:
            self.s3.head_object(Bucket=self.bucket, Key=self._key(key))
            return True
        except Exception:
            return False

    def delete(self, key):
        self.s3.delete_object(Bucket=self.bucket, Key=self._key(key))

    def list(self, prefix):
        out = []
        pfx = self._key(prefix)
        for page in self.s3.get_paginator("list_objects_v2").paginate(Bucket=self.bucket, Prefix=pfx):
            for o in page.get("Contents", []):
                k = o["Key"]
                out.append(k[len(self.prefix) + 1:] if self.prefix else k)
        return sorted(out)

    # queue
    def _url(self, lane):
        if lane not in self._qurl:
            self._qurl[lane] = self.sqs.get_queue_url(QueueName=f"{self.name}-units-{lane}")["QueueUrl"]
        return self._qurl[lane]

    def send(self, lane, body):
        self.sqs.send_message(QueueUrl=self._url(lane), MessageBody=json.dumps(body))

    def receive(self, lanes, visibility):
        for lane in lanes:
            r = self.sqs.receive_message(QueueUrl=self._url(lane), MaxNumberOfMessages=1,
                                         VisibilityTimeout=int(visibility), WaitTimeSeconds=1,
                                         AttributeNames=["ApproximateReceiveCount"])
            for m in r.get("Messages", []):
                return Message((lane, m["ReceiptHandle"]), json.loads(m["Body"]), lane,
                               int(m["Attributes"].get("ApproximateReceiveCount", "1")))
        return None

    def extend(self, msg, visibility):
        lane, h = msg.handle
        self.sqs.change_message_visibility(QueueUrl=self._url(lane), ReceiptHandle=h, VisibilityTimeout=int(visibility))

    def release(self, msg):
        self.extend(msg, 0)

    def delete_message(self, msg):
        lane, h = msg.handle
        self.sqs.delete_message(QueueUrl=self._url(lane), ReceiptHandle=h)

    def depth(self, lane):
        a = self.sqs.get_queue_attributes(QueueUrl=self._url(lane), AttributeNames=["ApproximateNumberOfMessages"])
        return int(a["Attributes"]["ApproximateNumberOfMessages"])

    # ledger: runs
    def put_run(self, run_id, doc):
        self.ddb.put_item(TableName=f"{self.name}-runs",
                          Item={"run_id": {"S": run_id}, "doc": {"S": json.dumps(doc)}, "ver": {"N": "0"}})

    def get_run(self, run_id):
        r = self.ddb.get_item(TableName=f"{self.name}-runs", Key={"run_id": {"S": run_id}}, ConsistentRead=True)
        return json.loads(r["Item"]["doc"]["S"]) if "Item" in r else None

    def update_run(self, run_id, fields):
        while True:
            r = self.ddb.get_item(TableName=f"{self.name}-runs", Key={"run_id": {"S": run_id}}, ConsistentRead=True)
            doc, ver = json.loads(r["Item"]["doc"]["S"]), r["Item"]["ver"]["N"]
            doc.update(fields)
            try:
                self.ddb.put_item(TableName=f"{self.name}-runs",
                                  Item={"run_id": {"S": run_id}, "doc": {"S": json.dumps(doc)},
                                        "ver": {"N": str(int(ver) + 1)}},
                                  ConditionExpression="ver = :v", ExpressionAttributeValues={":v": {"N": ver}})
                return doc
            except self.ddb.exceptions.ConditionalCheckFailedException:
                continue

    def next_seq(self, run_id, name="merge"):
        r = self.ddb.update_item(TableName=f"{self.name}-runs", Key={"run_id": {"S": f"{run_id}#counter#{name}"}},
                                 UpdateExpression="ADD #v :one", ExpressionAttributeNames={"#v": "value"},
                                 ExpressionAttributeValues={":one": {"N": "1"}}, ReturnValues="UPDATED_NEW")
        return int(r["Attributes"]["value"]["N"])

    # ledger: units
    def _item(self, run_id, doc, ver):
        it = {"run_id": {"S": run_id}, "uid": {"S": doc["uid"]}, "doc": {"S": json.dumps(doc)},
              "ver": {"N": str(ver)}, "status": {"S": doc.get("status") or ""},
              "run_status": {"S": f"{run_id}#{doc.get('status')}"}, "attempt": {"N": str(doc.get("attempt", 0))}}
        if doc.get("merged_seq") is not None:
            it["merged_seq"] = {"N": str(doc["merged_seq"])}
        return it

    def put_unit(self, run_id, doc):
        try:
            self.ddb.put_item(TableName=f"{self.name}-units", Item=self._item(run_id, doc, 0),
                              ConditionExpression="attribute_not_exists(uid)")
            return True
        except self.ddb.exceptions.ConditionalCheckFailedException:
            return False

    def get_unit(self, run_id, uid):
        r = self.ddb.get_item(TableName=f"{self.name}-units", Key={"run_id": {"S": run_id}, "uid": {"S": uid}},
                              ConsistentRead=True)
        return json.loads(r["Item"]["doc"]["S"]) if "Item" in r else None

    def update_unit(self, run_id, uid, fields, expect=None, check=None):
        while True:
            r = self.ddb.get_item(TableName=f"{self.name}-units", Key={"run_id": {"S": run_id}, "uid": {"S": uid}},
                                  ConsistentRead=True)
            if "Item" not in r:
                return None
            doc, ver = json.loads(r["Item"]["doc"]["S"]), int(r["Item"]["ver"]["N"])
            if not _matches(doc, expect) or (check and not check(doc)):
                return None
            doc.update(fields)
            try:
                self.ddb.put_item(TableName=f"{self.name}-units", Item=self._item(run_id, doc, ver + 1),
                                  ConditionExpression="ver = :v", ExpressionAttributeValues={":v": {"N": str(ver)}})
                return doc
            except self.ddb.exceptions.ConditionalCheckFailedException:
                continue

    def _query(self, **kw):
        out = []
        for page in self.ddb.get_paginator("query").paginate(TableName=f"{self.name}-units", **kw):
            out += [json.loads(i["doc"]["S"]) for i in page.get("Items", [])]
        return out

    def units_by_status(self, run_id, status):
        return self._query(IndexName="by_status", KeyConditionExpression="run_status = :s",
                           ExpressionAttributeValues={":s": {"S": f"{run_id}#{status}"}})

    def count_by_status(self, run_id):
        out = {}
        for d in self.all_units(run_id):
            out[d.get("status")] = out.get(d.get("status"), 0) + 1
        return out

    def units_merged_after(self, run_id, seq):
        return self._query(IndexName="by_seq", KeyConditionExpression="run_id = :r AND merged_seq > :q",
                           ExpressionAttributeValues={":r": {"S": run_id}, ":q": {"N": str(seq)}})

    def all_units(self, run_id):
        return self._query(KeyConditionExpression="run_id = :r", ExpressionAttributeValues={":r": {"S": run_id}})
