#!/usr/bin/env python3
"""Regression tests for SPARTA rigid-body surface objects (fix rigid).

Runs a set of small input decks and checks the results against analytic
expectations with tolerances.  Deterministic tests (no particles) must
produce identical results on any number of procs; statistical tests use
loose tolerances.

Usage:
  python3 run_tests.py --exe /path/to/spa_serial
  python3 run_tests.py --exe /path/to/spa_mpi --mpi "mpirun -np 4"

Exit code = number of failed tests.
"""

import argparse
import os
import shlex
import subprocess
import sys

THISDIR = os.path.dirname(os.path.abspath(__file__))

# mass of N in air.species, fnum from the decks

MASS_N = 2.325e-26
FNUM = 0.001


def run_deck(exe_cmd, deck, extra=None, expect_error=False):
    """Run one deck, return (returncode, stdout+stderr)."""
    cmd = exe_cmd + ["-in", deck] + (extra or [])
    proc = subprocess.run(cmd, cwd=THISDIR, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=600)
    return proc.returncode, proc.stdout


def parse_stats(output):
    """Parse the last stats table in the output into a list of dicts."""
    header = None
    rows = []
    for line in output.splitlines():
        toks = line.split()
        if not toks:
            continue
        if toks[0] == "Step":
            header = toks
            rows = []
            continue
        if header is None:
            continue
        if toks[0].startswith("Loop"):
            header = None
            continue
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            continue
        if len(vals) == len(header):
            rows.append(dict(zip(header, vals)))
    return rows


def approx(a, b, rel=0.0, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


# ----------------------------------------------------------------------
# individual tests: each returns a list of failure strings (empty = pass)
# ----------------------------------------------------------------------

def test_ballistic(exe_cmd):
    rc, out = run_deck(exe_cmd, "in.test.ballistic")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    last = rows[-1]
    fails = []
    # tolerances limited by the ~8 significant digits of stats output
    # com = com0 + v * t,  t = 1000 * 1e-4 = 0.1
    if not approx(last["f_1[1]"], 3.0 + 12.0 * 0.1, rel=1e-7):
        fails.append("xcm = %.12g, expected 4.2" % last["f_1[1]"])
    if not approx(last["f_1[2]"], 4.0 + 7.0 * 0.1, rel=1e-7):
        fails.append("ycm = %.12g, expected 4.7" % last["f_1[2]"])
    if not approx(last["f_1[4]"], 12.0, rel=1e-7):
        fails.append("vx = %.12g, expected 12" % last["f_1[4]"])
    if not approx(last["f_1[5]"], 7.0, rel=1e-7):
        fails.append("vy = %.12g, expected 7" % last["f_1[5]"])
    # omega_z = Lz / Izz = 1e-22 / 1.6666667e-23
    womega = 1.0e-22 / 1.6666667e-23
    if not approx(last["f_1[15]"], womega, rel=1e-7):
        fails.append("omega = %.12g, expected %.12g" % (last["f_1[15]"], womega))
    return fails


def test_bounce(exe_cmd):
    fails = []
    # elastic cases: both force laws must rebound at -50 within 1%
    # damped case: DEM spring-dashpot must rebound slower, with a
    #   coefficient of restitution in a loose band around the analytic
    #   value for the effective 2-contact linear spring-dashpot
    cases = (
        ("linear", "1.0e-18", "0.0", (-50.5, -49.5)),
        ("hertz", "4.0e-18", "0.0", (-50.5, -49.5)),
        ("linear-damped", "1.0e-18", "3.0e-21", (-32.5, -17.5)),
    )
    for label, pk, pdamp, vband in cases:
        pstyle = label.split("-")[0]
        rc, out = run_deck(exe_cmd, "in.test.bounce",
                           extra=["-var", "pstyle", pstyle,
                                  "-var", "pk", pk,
                                  "-var", "pdamp", pdamp])
        if rc:
            fails.append("%s: run failed with exit code %d" % (label, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("%s: no stats output" % label)
            continue
        # force engages when body corner is cutoff=0.5 from wall at x=8,
        # i.e. xcm <= 7.5; no tunneling means xcm never reaches the wall
        xmax = max(r["f_1[1]"] for r in rows)
        if xmax > 7.55:
            fails.append("%s: max xcm = %.6g, tunneled past push-off zone"
                         % (label, xmax))
        vfinal = rows[-1]["f_1[4]"]
        if not vband[0] <= vfinal <= vband[1]:
            fails.append("%s: final vx = %.6g, expected within [%g,%g]"
                         % (label, vfinal, vband[0], vband[1]))
    return fails


def test_momentum(exe_cmd):
    rc, out = run_deck(exe_cmd, "in.test.momentum")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    mass_body = 1.0e-22

    def total_px(r):
        return MASS_N * FNUM * r["c_r"] + mass_body * r["f_1[4]"]

    p0 = total_px(rows[0])
    fails = []
    for r in rows:
        if not approx(total_px(r), p0, rel=0.015):
            fails.append("step %d: total px = %.6g vs initial %.6g "
                         "(> 1.5%% drift)" % (int(r["Step"]), total_px(r), p0))
    # body must have absorbed a significant momentum fraction by the end
    pbody = mass_body * rows[-1]["f_1[4]"]
    if pbody < 0.25 * p0:
        fails.append("body momentum %.3g < 25%% of gas momentum %.3g, "
                     "coupling too weak" % (pbody, p0))
    return fails


def test_overrun(exe_cmd):
    rc, out = run_deck(exe_cmd, "in.test.overrun")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    fails = []
    # f_1 (scalar) = cumulative particles deleted inside the body
    # the body sweeps 40% of a grid cell per step through nearly
    # stationary gas.  swept collision coverage adds the body surfs to
    # the collision lists of every cell they sweep into during the step,
    # so no particle is overtaken undetected: every particle in the path
    # is reflected off the moving surf rather than deleted.  the deletion
    # count must therefore not grow at all after the initial setup.  any
    # growth means a particle tunneled into the body and was deleted.
    ndel = rows[-1]["f_1"] - rows[0]["f_1"]
    nptotal = rows[0]["Np"]
    if ndel != 0:
        fails.append("deleted %g of %g particles after setup: a particle "
                     "was overtaken by the moving body instead of being "
                     "reflected (swept collision coverage failed)"
                     % (ndel, nptotal))
    return fails


def test_remap(exe_cmd):
    results = {}
    fails = []
    for mode in ("cutcell", "incremental"):
        rc, out = run_deck(exe_cmd, "in.test.remap",
                           extra=["-var", "mode", mode])
        if rc:
            fails.append("mode %s: run failed with exit code %d" % (mode, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("mode %s: no stats output" % mode)
            continue
        last = rows[-1]
        results[mode] = (last["f_1[1]"], last["f_1[2]"], last["f_1[15]"])
    if fails:
        return fails
    ref = results["cutcell"]
    for i, name in enumerate(("xcm", "ycm", "omega")):
        if not approx(results["incremental"][i], ref[i], rel=1e-10, abs_=1e-13):
            fails.append("incremental: %s = %.15g differs from cutcell %.15g"
                         % (name, results["incremental"][i], ref[i]))
    return fails


def test_multiremap(exe_cmd):
    # two gas-driven bodies: cutcell and incremental must give identical
    # trajectories, verifying multi-body incremental re-cut
    results = {}
    fails = []
    for mode in ("cutcell", "incremental"):
        rc, out = run_deck(exe_cmd, "in.test.multiremap",
                           extra=["-var", "mode", mode])
        if rc:
            fails.append("mode %s: run failed with exit code %d" % (mode, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("mode %s: no stats output" % mode)
            continue
        last = rows[-1]
        results[mode] = tuple(last[k] for k in
                              ("f_1[1]", "f_1[2]", "f_1[15]",
                               "f_2[1]", "f_2[2]", "f_2[15]"))
    if fails:
        return fails
    labels = ("b1 xcm", "b1 ycm", "b1 omega", "b2 xcm", "b2 ycm", "b2 omega")
    for i, name in enumerate(labels):
        if not approx(results["incremental"][i], results["cutcell"][i],
                      rel=1e-10, abs_=1e-13):
            fails.append("%s: incremental %.15g differs from cutcell %.15g"
                         % (name, results["incremental"][i],
                            results["cutcell"][i]))
    return fails


def test_pushpair(exe_cmd):
    # ASYMMETRIC body-body contact: heavy large body overtakes a light
    # small one, corner-vs-face contact. Total momentum of the pair must
    # be conserved (contact forces are equal-and-opposite on both
    # bodies); tolerance is set by the 8-digit stats output, not physics
    rc, out = run_deck(exe_cmd, "in.test.pushpair")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    m1, m2 = 4.0e-22, 1.0e-22
    px0 = m1 * rows[0]["f_1[4]"] + m2 * rows[0]["f_2[4]"]
    py0 = m1 * rows[0]["f_1[5]"] + m2 * rows[0]["f_2[5]"]
    fails = []
    for r in rows:
        px = m1 * r["f_1[4]"] + m2 * r["f_2[4]"]
        py = m1 * r["f_1[5]"] + m2 * r["f_2[5]"]
        if not approx(px, px0, rel=1e-6):
            fails.append("step %d: px = %.10e vs initial %.10e, body-body "
                         "contact violates momentum conservation"
                         % (int(r["Step"]), px, px0))
        if abs(py - py0) > 1e-6 * abs(px0):
            fails.append("step %d: py = %.3e drifted from %.3e"
                         % (int(r["Step"]), py, py0))
    # the collision must actually have happened
    if rows[-1]["f_2[4]"] < 20.0:
        fails.append("final body2 vx = %.6g, no significant collision "
                     "occurred; test geometry is broken"
                     % rows[-1]["f_2[4]"])
    return fails


def test_twobody(exe_cmd):
    rc, out = run_deck(exe_cmd, "in.test.twobody")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    last = rows[-1]
    v1, v2 = last["f_1[4]"], last["f_2[4]"]
    fails = []
    # head-on symmetric collision: velocities reverse, ~elastic
    if not approx(v1, -30.0, rel=0.02):
        fails.append("body1 final vx = %.6g, expected -30 within 2%%" % v1)
    if not approx(v2, 30.0, rel=0.02):
        fails.append("body2 final vx = %.6g, expected +30 within 2%%" % v2)
    # total momentum ~ 0
    if abs(v1 + v2) > 0.05:
        fails.append("momentum asymmetry |v1+v2| = %.4g > 0.05" % abs(v1 + v2))
    return fails


def negative_test(exe_cmd, deck, message):
    rc, out = run_deck(exe_cmd, deck, expect_error=True)
    fails = []
    if rc == 0:
        fails.append("run succeeded but an error was expected")
    if message not in out:
        fails.append("expected error message not found: '%s'" % message)
    return fails


def test_badmoi(exe_cmd):
    return negative_test(exe_cmd, "in.test.badmoi", "triangle inequality")


def test_notwatertight(exe_cmd):
    return negative_test(exe_cmd, "in.test.notwatertight", "not watertight")


TESTS = [
    ("ballistic", test_ballistic),
    ("bounce", test_bounce),
    ("momentum", test_momentum),
    ("overrun", test_overrun),
    ("remap", test_remap),
    ("multiremap", test_multiremap),
    ("twobody", test_twobody),
    ("pushpair", test_pushpair),
    ("badmoi", test_badmoi),
    ("notwatertight", test_notwatertight),
]

# tests whose decks support -var dist 1 (global surfs explicit/distributed)
# remap/multiremap are excluded: incremental re-cut requires
# non-distributed surfs

DIST_TESTS = {"ballistic", "bounce", "momentum", "overrun",
              "twobody", "pushpair"}


def main():
    global run_deck
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True,
                        help="path to SPARTA executable")
    parser.add_argument("--mpi", default="",
                        help='MPI launcher prefix, e.g. "mpirun -np 4"')
    parser.add_argument("--tests", default="",
                        help="comma-separated subset of tests to run")
    parser.add_argument("--dist", action="store_true",
                        help="run with distributed surfs "
                             "(global surfs explicit/distributed)")
    args = parser.parse_args()

    exe_cmd = shlex.split(args.mpi) + [os.path.abspath(args.exe)]

    subset = None
    if args.tests:
        subset = set(args.tests.split(","))

    if args.dist:
        if subset is None:
            subset = set(DIST_TESTS)
        else:
            subset &= DIST_TESTS
        base_run_deck = run_deck

        def dist_run_deck(exe_cmd, deck, extra=None, expect_error=False):
            extra = (extra or []) + ["-var", "dist", "1"]
            return base_run_deck(exe_cmd, deck, extra, expect_error)

        run_deck = dist_run_deck

    nfail = 0
    for name, func in TESTS:
        if subset and name not in subset:
            continue
        fails = func(exe_cmd)
        if fails:
            nfail += 1
            print("FAIL %s" % name)
            for f in fails:
                print("     %s" % f)
        else:
            print("PASS %s" % name)

    print("%d test(s) failed" % nfail)
    return nfail


if __name__ == "__main__":
    sys.exit(main())
