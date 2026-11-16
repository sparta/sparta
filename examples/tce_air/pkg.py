"""Helpers shared by the steps of the tce_air package.

Every step defines its runs (Run objects) in runs.py.  The run list is
used both to run and to analyze, so the analysis never depends on how or
where a run happened.  Each step has one size, "quick", chosen so that the
whole package runs on a 4-core machine in about an hour; set_size()
registers other sizes (more particles, steps or seeds for the same decks)
from outside the package.

Executables come from environment variables, or from bin/ as written by
build.sh:
  SPARTA_EXE        executable of this branch (default bin/spa_release)
  SPARTA_EXE_ARGS   extra arguments of that executable
  SPARTA_CMP_EXE    comparison build with the earlier TCE treatments
                    (react_modify vib_energy dof, elec_energy yes|no|micro);
                    optional (default bin/spa_compare)
  SPARTA_CMP_ARGS   extra arguments of the comparison build
  SPARTA_NP         MPI ranks per run (default 4)
  SPARTA_MPIRUN     MPI launcher (default "mpirun"), e.g. "mpiexec"
  TCE_AIR_RESULTS   results directory (default ./results next to this file)
"""
import json
import math
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "data"))


def results_dir():
    return os.environ.get("TCE_AIR_RESULTS") or os.path.join(HERE, "results")


# ------------------------------------------------------------------ sizes
_EXTRA_SIZES = {}


def set_size(step, name, params):
    """register size 'name' of a step (a dict like the entries of the
    step's SIZES) from outside the package"""
    _EXTRA_SIZES.setdefault(step, {})[name] = dict(params)


def size_params(step, sizes, name):
    if name in _EXTRA_SIZES.get(step, {}):
        return _EXTRA_SIZES[step][name]
    if name in sizes:
        return sizes[name]
    raise KeyError("step %s has no size '%s' (it has: %s)"
                   % (step, name, ", ".join(sorted(set(sizes) |
                                            set(_EXTRA_SIZES.get(step, {}))))))


# ------------------------------------------------------------------ runs
class Run:
    """one SPARTA run: deck (path relative to the step directory), -var
    settings, executable kind ('release' or 'compare'), and a unique tag.
    files: extra files the deck reads, copied into the run directory."""

    def __init__(self, step, tag, deck, variables, kind="release",
                 files=(), stepdir=None):
        self.step, self.tag, self.deck = step, tag, deck
        self.variables, self.kind, self.files = dict(variables), kind, files
        self._stepdir = stepdir

    @property
    def stepdir(self):
        """directory of the deck and its files: the step's directory in
        this package, unless the run was defined elsewhere"""
        return self._stepdir or os.path.join(HERE, self.step)

    @property
    def workdir(self):
        return os.path.join(results_dir(), self.step, self.tag)

    @property
    def log(self):
        return os.path.join(self.workdir, "log.sparta")

    def inputs(self):
        """the files the deck reads: every word of a command line (after
        substituting the run's variables) that names a file, as a path or
        relative to the step directory (whose files are copied into the
        run directory); e.g. the data files ${DATA}/air.*"""
        text = open(os.path.join(self.stepdir, self.deck)).read()
        for k, v in self.variables.items():
            text = text.replace("${%s}" % k, str(v))
        found = set()
        for line in text.split("\n"):
            for w in line.split("#")[0].split():
                for p in (w, os.path.join(self.stepdir, w)):
                    if os.path.isfile(p):
                        found.add(os.path.realpath(p))
                        break
        found.update(os.path.realpath(os.path.join(self.stepdir, f))
                     for f in self.files)
        return sorted(found)

    def record(self):
        """what defines the run: deck (name and content), the content of
        the files it reads, variables, executable kind"""
        import hashlib

        def sha1(path):
            with open(path, "rb") as fh:
                return hashlib.sha1(fh.read()).hexdigest()
        return {"deck": self.deck,
                "deck_sha1": sha1(os.path.join(self.stepdir, self.deck)),
                "inputs_sha1": {os.path.basename(p): sha1(p)
                                for p in self.inputs()},
                "kind": self.kind,
                "variables": {k: str(v) for k, v in self.variables.items()}}

    def done(self):
        """the log finished, and run.json shows the run was made with the
        current deck (name and content), input files, variables and
        executable kind: a run whose definition changed is run again"""
        rec = os.path.join(self.workdir, "run.json")
        if not (os.path.exists(self.log) and os.path.exists(rec)):
            return False
        try:
            same = json.load(open(rec)) == self.record()
        except ValueError:
            return False
        return same and "Loop time" in open(self.log).read()

    def command(self, exe=None, args=None):
        """argument list of the SPARTA call (without the MPI launcher);
        exe and args replace the executable of the run's kind"""
        e, extra = exe_and_args(self.kind)
        if exe:
            e, extra = exe, list(args or [])
        cmd = [e] + extra + ["-in", os.path.join(self.stepdir, self.deck),
                             "-log", self.log, "-screen", "none"]
        for k, v in self.variables.items():
            cmd += ["-var", k, str(v)]
        return cmd


def exe_and_args(kind="release"):
    if kind == "release":
        e = os.environ.get("SPARTA_EXE") or os.path.join(BIN, "spa_release")
        a = os.environ.get("SPARTA_EXE_ARGS", "")
    elif kind == "compare":
        e = os.environ.get("SPARTA_CMP_EXE") or os.path.join(BIN, "spa_compare")
        a = os.environ.get("SPARTA_CMP_ARGS", "")
    elif kind == "precap":
        e = os.environ.get("SPARTA_PRECAP_EXE") or os.path.join(BIN,
                                                                "spa_precap")
        a = os.environ.get("SPARTA_PRECAP_ARGS", "")
    else:
        raise ValueError(kind)
    return e, shlex.split(a)


def have(kind):
    return os.path.exists(exe_and_args(kind)[0])


def prepare(run):
    os.makedirs(run.workdir, exist_ok=True)
    with open(os.path.join(run.workdir, "run.json"), "w") as fh:
        json.dump(run.record(), fh, indent=1)
    if os.path.exists(run.log):
        os.remove(run.log)
    for f in run.files:
        src = f if os.path.isabs(f) else os.path.join(run.stepdir, f)
        dst = os.path.join(run.workdir, os.path.basename(f))
        if not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)


def start_cpu(run, np=None):
    """start one run on the CPU with np MPI ranks, without waiting"""
    if not have(run.kind):
        raise FileNotFoundError("no '%s' executable: run ./build.sh first"
                                % run.kind)
    prepare(run)
    cmd = shlex.split(os.environ.get("SPARTA_MPIRUN", "mpirun")) + \
        ["-np", str(np or int(os.environ.get("SPARTA_NP", "4")))] + \
        run.command()
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1",
               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    out = open(os.path.join(run.workdir, "stdout"), "w")
    return subprocess.Popen(cmd, cwd=run.workdir, env=env, stdout=out,
                            stderr=subprocess.STDOUT)


# ------------------------------------------------------------ log parsing
def tallies(log):
    """reaction tallies printed at the end of a run: {label: count}.
    SPARTA prints counts of a million and more in %g format (1.0273e+06,
    six significant digits), so the count is parsed as a float."""
    out = {}
    for line in open(log):
        m = re.match(r"\s*reaction (.*): ([0-9.eE+-]+)\s*$", line)
        if m:
            out[" ".join(m.group(1).split())] = int(round(float(m.group(2))))
    return out


def stats(log):
    """last stats table of a log: (header list, rows of floats)"""
    head, rows, cur = None, [], None
    for ln in open(log).read().splitlines():
        w = ln.split()
        if w and w[0] == "Step":
            cur, rows, head = w, [], w
            continue
        if cur is not None:
            if ln.startswith("Loop time"):
                cur = None
                continue
            try:
                rows.append([float(x) for x in w])
            except ValueError:
                pass
    return head, rows


def warnings(log):
    return [ln.strip() for ln in open(log) if ln.startswith("WARNING")]


# ------------------------------------------------------------- reporting
class Report:
    """PASS/FAIL lines of one step, written to results/<step>/results.txt"""

    def __init__(self, step):
        """step: the step's directory name, optionally followed by a
        label such as "(quick)" for the report's first line"""
        self.step, self.lines, self.nfail = step, [], 0
        self.path = os.path.join(results_dir(), step.split()[0],
                                 "results.txt")

    def info(self, text):
        self.lines.append("  " + text)
        print("  " + text)

    def check(self, ok, text):
        line = "  [%s] %s" % ("PASS" if ok else "FAIL", text)
        self.lines.append(line)
        print(line)
        self.nfail += 0 if ok else 1

    def write(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with open(self.path, "w") as fh:
            fh.write("%s\n%s\nfailed checks: %d\n"
                     % (self.step, "\n".join(self.lines), self.nfail))
        return self.nfail


def mean_err_corr(x):
    """mean of a correlated time series and its standard error, from the
    integrated autocorrelation time (Sokal's window: sum the
    autocorrelation up to the first lag M >= 5 tau)"""
    import numpy as np
    x = np.asarray(x, float)
    n = len(x)
    if n < 4:
        return float(x.mean()), float("inf")
    d = x - x.mean()
    var = d.dot(d) / n
    if var == 0.0:
        return float(x.mean()), 0.0
    tau = 0.5
    for m in range(1, n // 2):
        tau += d[:-m].dot(d[m:]) / ((n - m) * var)
        if m >= 5 * tau:
            break
    tau = max(tau, 0.5)
    return float(x.mean()), float(np.sqrt(2.0 * tau * var / n))


def mean_sd(xs):
    n = len(xs)
    m = sum(xs) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1)) if n > 1 else 0.0
    return m, sd
