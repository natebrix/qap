/*
============================================================

 The settings registry (docs/settings_plan.md, Phases 1 and 2;
 CONTRIBUTING.md, "Where settings live").

 Every setting the solver reads is an entry in one table (settings.C):
 its par key and/or environment name, type, default, range, class and
 one line of description. The table is read once, the first time
 anything asks for a setting (settings_load is idempotent), so every
 reader sees the same value and nothing is re-read mid-run.

 Where a setting lives:
   policy       changes the search: a key in the params file's
                [settings] section (keys dnn.*, cut.*, pool.*, ship.*,
                branch.*, est.*, fw.*). The flags -x and -w are aliases
                of dnn.depths/maxgap/giveup/warm and fw.patience and
                override the file. Recorded in state files.
   experiment   an A/B switch in flight: an environment variable
                (DNN_*), also recorded; it graduates to a key or a
                default, or goes.
   run          this run's identity: environment (QAP_EST_SEED,
                QAP_ALLOW_UNKNOWN) or the -p path; not recorded.
   diagnostic   output only: environment (QAP_DNN_TRACE, QAP_NODE_TRACE, ...).

 Unknown names are errors (exit 2): a QAP_DNN_/QAP_EST_/QAP_XBAR_/QAP_NODE_/
 DNN_ variable not in the table, a [settings] key not in the table, or a
 value that does not parse or is out of range. The environment name a
 policy setting had before Phase 2 (QAP_DNN_BETA, ...) is an error whose
 message names its key. QAP_ALLOW_UNKNOWN=1 downgrades unknown and old
 names to warnings, for bisecting old commits. Other QAP_* names belong
 to the launchers and are not checked.

 The record: settings_record() is "<hash> NAME=VALUE ..." where the
 hash covers every policy and experiment entry (defaults included) and
 the list is the entries that differ from their defaults. write_state
 puts it in the state file (QAP-STATE 8) and read_state compares the
 hash on -r (--allow-settings-change overrides).

============================================================
*/
#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <iosfwd>

enum SettingClass { SET_POLICY, SET_EXPERIMENT, SET_RUN, SET_DIAGNOSTIC };
enum SettingSource { SRC_DEFAULT, SRC_ENV, SRC_FLAG, SRC_PAR };

struct Setting {
  const char *key;       // the params-file key, or NULL
  const char *env;       // the environment name, or NULL
  const char *old_env;   // the environment name it had before Phase 2, or NULL
  char type;             // 'i' integer, 'd' double, 's' string
  const char *def;       // the default, as text ("" for an unset string)
  double lo, hi;         // inclusive range for 'i' and 'd'
  SettingClass cls;
  const char *desc;
  // resolved at load
  std::string text;
  SettingSource src;
  long i;
  double d;
  const char *name() const { return key ? key : env; }
};

// Read the environment into the table; exits 2 on an unknown name or a
// bad value. Called by the first setting_* lookup if not before.
void settings_load();

// The [settings] section of a params file (keys and values, "key = value"
// per line, # comments). Exits 2 on an unknown key or a bad value. A
// value from the file does not override one given by a flag.
void settings_read_par(const char *path);

// Lookups by key or environment name. An unknown name is a programming
// error (abort).
const Setting &setting(const char *name);
long setting_i(const char *name);
double setting_d(const char *name);
const char *setting_s(const char *name);      // the default when unset
bool setting_given(const char *name);         // source is not the default

// Set an entry from a flag or the par file (validated like the
// environment; exits 2 on a bad value). SRC_PAR yields to SRC_FLAG.
void settings_set(const char *name, const std::string &value, SettingSource src);

// Replace an entry's default (dnntest's DNNParams defaults differ from the
// search's). Only for an entry still at its default; exits 2 otherwise.
void settings_override_default(const char *name, const std::string &value);

// The startup dump: one line per entry,
//   setting <name> <value> <source> <class>
// ("" for an empty value), after a line "settings <count> hash <hash>".
void settings_dump(std::ostream &os);

// The reference table as Markdown (RUNNING.md includes it; tests/settings_doc.py
// checks that the two agree), and the old-name -> key map, one "OLD KEY" per
// line (bench/envpar.py turns old environment recipes into [settings] sections).
void settings_table(std::ostream &os);
void settings_env_map(std::ostream &os);

// "<hash> NAME=VALUE ..." (see above), and its check against a record read
// from a state file: true when the hashes agree; `why` explains otherwise.
std::string settings_record();
bool settings_check_record(const std::string &record, std::string &why);

// 64-bit FNV-1a of a file's bytes as 16 hex digits, "" if it cannot be read.
std::string settings_hash_file(const char *path);

#endif
