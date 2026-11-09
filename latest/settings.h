/*
============================================================

 The settings registry (docs/settings_plan.md, Phase 1).

 Every setting the solver reads from the environment, and the policy
 flags -x and -w, is an entry in one table (settings.C): name, type,
 default, range, class, one-line description. The table is read once,
 from the environment, the first time anything asks for a setting
 (settings_load is idempotent), so every reader sees the same value
 and nothing is re-read mid-run.

 Classes:
   policy       changes the search (goes into the state-file record)
   experiment   an A/B switch in flight: also changes the search, also
                recorded; graduates to a par key or a default, or goes
   run          this run's identity or I/O (seed, allow-unknown)
   diagnostic   output only (trace files, debug prints)

 Unknown names are errors: a variable in the environment whose name
 starts with one of the solver's prefixes (QAP_DNN_, QAP_EST_,
 QAP_XBAR_, DNN_) and is not in the table stops the run with exit code
 2, as does a value that does not parse or is out of range.
 QAP_ALLOW_UNKNOWN=1 downgrades unknown names to a warning (for
 bisecting old commits). Other QAP_* names belong to the launchers
 (QAP_ROLE, QAP_BUCKET, ...) and are not checked.

 The record: settings_record() is "<hash> NAME=VALUE ..." where the
 hash covers every policy and experiment entry (defaults included) and
 the list is the entries that differ from their defaults. write_state
 puts it in the state file (QAP-STATE 8) and read_state compares the
 hash on -r (driver.C: --allow-settings-change overrides).

============================================================
*/
#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <iosfwd>

enum SettingClass { SET_POLICY, SET_EXPERIMENT, SET_RUN, SET_DIAGNOSTIC };
enum SettingSource { SRC_DEFAULT, SRC_ENV, SRC_FLAG, SRC_PAR };

struct Setting {
  const char *name;
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
};

// Read the environment into the table; exits 2 on an unknown name or a
// bad value. Called by the first setting_* lookup if not before.
void settings_load();

// Lookups. An unknown name is a programming error (abort).
const Setting &setting(const char *name);
long setting_i(const char *name);
double setting_d(const char *name);
const char *setting_s(const char *name);      // the default when unset
bool setting_given(const char *name);         // source is not the default

// Set an entry from a flag or the par file (validated like the
// environment; exits 2 on a bad value).
void settings_set(const char *name, const std::string &value, SettingSource src);

// The startup dump: one line per entry,
//   setting <name> <value> <source> <class>
// ("" for an empty value), after a line "settings <count> hash <hash>".
void settings_dump(std::ostream &os);

// "<hash> NAME=VALUE ..." (see above), and its check against a record read
// from a state file: true when the hashes agree; `why` explains otherwise.
std::string settings_record();
bool settings_check_record(const std::string &record, std::string &why);

// 64-bit FNV-1a of a file's bytes as 16 hex digits, "" if it cannot be read.
std::string settings_hash_file(const char *path);

#endif
