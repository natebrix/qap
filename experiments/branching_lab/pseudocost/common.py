"""Shared by describe.py, predict.py and truth.py: the tables, the features, the models, the metrics."""
import os, pickle, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
SCALE = 1e4   # ylog = log1p(SCALE * shortfall / |cutoff|): 0 for a fathomed child


def load(run):
    d = pd.DataFrame(pickle.load(open(os.path.join(HERE, "data", run + ".pkl"), "rb")))
    return derive(d)


def derive(d):
    a = d.cutoff.abs()
    d["open"] = 1 - d.fath1
    d["srel"] = d.splus / a
    d["ylog"] = np.log1p(SCALE * d.srel)
    d["pshort"] = (d.cutoff - d.pdnn) / a           # the parent's own first-solve shortfall
    d["cqpb_s"] = (d.cutoff - d.cqpb) / a           # the child's QPB shortfall (the baseline)
    d["gainrel"] = d.gaincap / a                    # min(child dnn, cutoff) - parent bound, over |cutoff|
    d["pair"] = d.row * 1000 + d.col
    return d


BASE = ["depth", "nfree", "pshort", "prelgap", "pcut_ran", "rsum", "rnnz", "csum", "cnnz",
        "rclu", "cclu", "dfix", "gl", "agree", "par_agree", "row", "col"]
QPB = ["cqpb_s", "pdnn_minus_qpb"]


def feats(d, with_qpb):
    f = BASE + (QPB if with_qpb else [])
    return [c for c in f if c in d and d[c].notna().any()]


# ---------------- models: each returns a function test_df -> score (larger: more likely open / short)

def m_qpb(tr):
    return lambda te: te.cqpb_s.values


def m_random(tr, seed=0):
    rng = np.random.default_rng(seed)
    return lambda te: rng.random(len(te))


def _onehot(d, cols, vocab=None):
    if vocab is None:
        vocab = {c: {v: i for i, v in enumerate(sorted(d[c].unique()))} for c in cols}
    blocks = []
    for c in cols:
        idx = d[c].map(vocab[c]).fillna(-1).astype(int).values
        ok = idx >= 0
        blocks.append(sparse.csr_matrix((np.ones(ok.sum()), (np.nonzero(ok)[0], idx[ok])),
                                        shape=(len(d), len(vocab[c]))))
    return sparse.hstack(blocks).tocsr(), vocab


def m_pseudocost(tr, cols=("row", "col", "pair"), alpha=3.0, target="ylog"):
    """a pseudocost table with shrinkage: additive row, column and pair effects, ridge-shrunk to the mean
    (alpha = the number of observations at which an effect gets half its raw weight, roughly)"""
    X, vocab = _onehot(tr, list(cols))
    mu = tr[target].mean()
    mdl = Ridge(alpha=alpha, fit_intercept=False).fit(X, tr[target].values - mu)
    return lambda te: mdl.predict(_onehot(te, list(cols), vocab)[0]) + mu


def _prep(d, fs, med=None):
    X = d[fs].astype(float).copy()
    if med is None: med = X.median()
    X = X.fillna(med)
    return X, med


def m_linear(tr, with_qpb=False, target="ylog"):
    fs = [f for f in feats(tr, with_qpb) if f not in ("row", "col", "rclu", "cclu")]
    X, med = _prep(tr, fs)
    mu, sd = X.mean(), X.std().replace(0, 1)
    mdl = Ridge(alpha=1.0).fit((X - mu) / sd, tr[target].values)
    return lambda te: mdl.predict((_prep(te, fs, med)[0] - mu) / sd)


def m_gbt(tr, with_qpb=False, target="ylog", fs=None, seed=0, within=False):
    """within: fit the target less its parent's mean (what ranking a parent's children needs), with the
    child's QPB shortfall less its parent's mean as an extra feature when with_qpb"""
    fs = fs or feats(tr, with_qpb)
    if within:
        tr = tr.copy()
        tr["_y"] = tr[target] - tr.groupby("pkey")[target].transform("mean"); target = "_y"
        if with_qpb:
            fs = fs + ["cqpb_w"]
            tr["cqpb_w"] = tr.cqpb_s - tr.groupby("pkey").cqpb_s.transform("mean")
    cat = [f in ("row", "col", "rclu", "cclu") for f in fs]
    mdl = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                        min_samples_leaf=20, l2_regularization=1.0,
                                        categorical_features=cat, random_state=seed)
    mdl.fit(tr[fs].astype(float).values, tr[target].values)
    def pred(te):
        if within and with_qpb:
            te = te.assign(cqpb_w=te.cqpb_s - te.groupby("pkey").cqpb_s.transform("mean"))
        return mdl.predict(te[fs].astype(float).values)
    return pred


# ---------------- metrics within each parent

def within_parent(te, score, min_kids=3):
    """te: a test frame; score: array aligned with te. Mean over parents of: AUC of score for open against
    fathomed (parents with both), Spearman of score with the shortfall (parents with an open child),
    R-precision (the k children ranked highest, k = the parent's open children: the share that are open),
    and recall/precision of a shortlist of the top quarter of each parent's children."""
    t = te.assign(s=score)
    auc, sp, rp, rec_n, rec_d, pre_n, pre_d = [], [], [], 0, 0, 0, 0
    npar = 0
    for _, g in t.groupby("pkey", sort=False):
        if len(g) < min_kids: continue
        y = g.open.values; s = g.s.values; k = int(y.sum())
        npar += 1
        if k == 0: continue
        if np.ptp(s) > 0 and k < len(g):
            r = rankdata(s)
            auc.append((r[y == 1].sum() - k * (k + 1) / 2) / (k * (len(g) - k)))
        else:
            auc.append(0.5) if k < len(g) else None
        if np.ptp(s) > 0 and np.ptp(g.splus.values) > 0:
            sp.append(spearmanr(s, g.splus.values).statistic)
        o = np.argsort(-s, kind="stable")
        rp.append(y[o[:k]].mean())
        q = max(1, int(np.ceil(len(g) / 4)))
        rec_n += y[o[:q]].sum(); rec_d += k; pre_n += y[o[:q]].sum(); pre_d += q
    return dict(parents=npar, scored=len(auc), auc=np.mean(auc) if auc else np.nan,
                spearman=np.nanmean(sp) if sp else np.nan, rprec=np.mean(rp) if rp else np.nan,
                recall25=rec_n / rec_d if rec_d else np.nan, prec25=pre_n / pre_d if pre_d else np.nan)


def fmt(x, nd=2):
    return "-" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))
