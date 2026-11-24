#!/usr/bin/env python3
"""Driver of the tce_air reproduction package.

  ./tce_air.py list    STEP [-v]                      runs and their status
  ./tce_air.py run     STEP [--np N] [--jobs J]       run the pending runs, J at
                            [--group G] [--shard i/n]  a time with N MPI ranks
                                                       each (N x J <= cores)
  ./tce_air.py analyze STEP                           figures, tables, checks

STEP is a step directory name (e.g. step1_rates) or 'all'.  Runs whose log
already finished are skipped, so 'run' can be repeated after an
interruption.  Results go to results/<step>/ (or $TCE_AIR_RESULTS).
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pkg  # noqa: E402

STEPS = sorted(d for d in os.listdir(HERE)
               if d.startswith("step") and os.path.isdir(os.path.join(HERE, d)))


def load(step, name):
    path = os.path.join(HERE, step, name + ".py")
    spec = importlib.util.spec_from_file_location("%s_%s" % (step, name), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def steps(arg):
    return STEPS if arg == "all" else [arg]


def selected(step, size, group):
    rl = load(step, "runs").runs(size)
    return [r for r in rl if group is None or r.tag.startswith(group + "_")]


def cmd_list(a):
    for st in steps(a.step):
        rl = selected(st, a.size, a.group)
        ndone = sum(r.done() for r in rl)
        print("%s (%s): %d runs, %d done" % (st, a.size, len(rl), ndone))
        if a.verbose:
            for r in rl:
                print("  %-4s %-8s %s" % ("done" if r.done() else "-",
                                          r.kind, r.tag))


def shard(rl, spec):
    """keep every n-th run starting at i, for spec 'i/n'"""
    if not spec:
        return rl
    i, n = (int(x) for x in spec.split("/"))
    return [r for k, r in enumerate(rl) if k % n == i]


def cmd_run(a):
    """run the pending runs of the step(s), --jobs at a time with --np MPI
    ranks each"""
    import time
    pending, skip = [], set()
    for st in steps(a.step):
        for r in shard(selected(st, a.size, a.group), a.shard):
            if r.done():
                continue
            if pkg.have(r.kind):
                pending.append(r)
            else:
                skip.add(r.kind)
    if skip:
        print("skipping the runs that need the '%s' executable (not built; "
              "see README)" % "/".join(sorted(skip)))
    print("%d runs, %d at a time" % (len(pending), a.jobs), flush=True)
    busy, failed, n0 = {}, [], len(pending)
    while pending or busy:
        for k in [k for k in range(a.jobs) if k not in busy]:
            if not pending:
                break
            r = pending.pop(0)
            busy[k] = (r, pkg.start_cpu(r, a.np), time.time())
            print("[%d/%d] start %s/%s" % (n0 - len(pending), n0, r.step,
                                            r.tag), flush=True)
        time.sleep(0.5)
        for k, (r, proc, t0) in list(busy.items()):
            if proc.poll() is not None:
                ok = proc.returncode == 0 and r.done()
                if not ok:
                    failed.append("%s/%s" % (r.step, r.tag))
                print("        %s %s/%s (%.0f s)"
                      % ("done" if ok else "FAILED", r.step, r.tag,
                         time.time() - t0), flush=True)
                del busy[k]
    if failed:
        print("%d runs FAILED (see log.sparta and stdout in their run "
              "directories): %s" % (len(failed), " ".join(failed)))
        sys.exit(1)


# figure and table names written by earlier versions of the package (numbered after
# an old draft of the paper); removed when a step is analyzed again
OLD_FIGS = ("fig03_progression.pdf", "fig05_eq_discretevib.pdf",
            "fig06_noneq.pdf", "fig07_reservoir.pdf", "fig04_qk_compare.pdf",
            "fig02_zeldovich.pdf", "fig09_closedbox.pdf", "fig10_keq_cea.pdf",
            "table6.tex", "table6.txt", "table5_qct.tex", "table5_qct.txt",
            "table7_recomb.tex", "table7_recomb.txt", "table8_energy.tex",
            "table8_energy.txt", "table11_cellT.tex", "table11_cellT.txt",
            "table10_captest.tex", "table10_captest.txt")


def cmd_analyze(a):
    nfail = 0
    for st in steps(a.step):
        print("== %s" % st)
        for f in OLD_FIGS:
            old = os.path.join(pkg.results_dir(), st, f)
            if os.path.exists(old):
                os.remove(old)
        nfail += load(st, "analyze").analyze(a.size)
    sys.exit(min(nfail, 255))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("list", "run", "analyze"):
        s = sub.add_parser(name)
        s.add_argument("step")
        s.add_argument("--size", default="quick")
        s.add_argument("--group", default=None)
        if name == "list":
            s.add_argument("-v", "--verbose", action="store_true")
        if name == "run":
            s.add_argument("--np", type=int, default=None,
                           help="MPI ranks per CPU run (default "
                           "$SPARTA_NP or 4)")
            s.add_argument("--jobs", type=int, default=1,
                           help="CPU runs at a time (e.g. --jobs 6 --np 2 "
                           "on 12 cores)")
            s.add_argument("--shard", default=None,
                           help="i/n: only every n-th run starting at i")
    a = p.parse_args()
    {"list": cmd_list, "run": cmd_run, "analyze": cmd_analyze}[a.cmd](a)


if __name__ == "__main__":
    main()
