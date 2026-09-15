#!/usr/bin/env python3
"""Regression tests for SPARTA rigid-body surface objects (fix rigid).

Runs a set of small input decks and checks the results against analytic
expectations with tolerances.  Deterministic tests (no particles) must
produce identical results on any number of procs; statistical tests use
loose tolerances.

Usage:
  python3 run_tests.py --exe /path/to/spa_serial
  python3 run_tests.py --exe /path/to/spa_mpi --mpi "mpirun -np 4"
  python3 run_tests.py --exe /path/to/spa_kokkos --args "-k on -sf kk"

Exit code = number of failed tests.
"""

import argparse
import glob
import os
import shlex
import subprocess
import sys

THISDIR = os.path.dirname(os.path.abspath(__file__))

# mass of N in air.species, fnum from the decks

MASS_N = 2.325e-26
FNUM = 0.001


def run_deck(exe_cmd, deck, extra=None):
    """Run one deck, return (returncode, stdout+stderr)."""
    cmd = exe_cmd + ["-in", deck] + (extra or [])
    try:
        proc = subprocess.run(cmd, cwd=THISDIR, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True,
                              timeout=600)
    except subprocess.TimeoutExpired as e:
        # a hang (e.g. a collective entered by only some ranks) is a
        # failure, not a reason to stall the suite
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else \
            (e.stdout or "")
        return -1, out + "\nTIMEOUT after 600 s\n"
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


def test_rotation(exe_cmd):
    # torque-free tumbling of an asymmetric body: angular momentum is
    # exactly conserved, so the rotational kinetic energy T = L.w/2 can
    # be formed from the reported angular velocity and the known L
    Lx = Ly = Lz = 1.0e-23
    fails = []
    drift = {}
    for rot in ("euler", "richardson"):
        rc, out = run_deck(exe_cmd, "in.test.rotation",
                           extra=["-var", "rot", rot])
        if rc:
            fails.append("%s: run failed with exit code %d" % (rot, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("%s: no stats output" % rot)
            continue
        ke = [0.5 * (Lx * r["f_1[13]"] + Ly * r["f_1[14]"] + Lz * r["f_1[15]"])
              for r in rows]
        ke0 = ke[0]
        drift[rot] = max(abs(k - ke0) for k in ke) / abs(ke0)

        # the body must actually tumble: for an asymmetric free body the
        # angular velocity is not constant.  a constant w would mean the
        # rotational dynamics are not being integrated at all
        wx = [r["f_1[13]"] for r in rows]
        if max(wx) - min(wx) < 0.01 * abs(max(wx)):
            fails.append("%s: angular velocity is nearly constant; an "
                         "asymmetric torque-free body must tumble" % rot)
    if fails:
        return fails

    # euler is first order, richardson second: at this timestep the
    # energy drift should differ by orders of magnitude
    if drift["euler"] > 1.0e-2:
        fails.append("euler: rotational energy drift %.3e is too large"
                     % drift["euler"])
    if drift["richardson"] > 1.0e-5:
        fails.append("richardson: rotational energy drift %.3e, expected "
                     "second-order accuracy" % drift["richardson"])
    if drift["richardson"] > drift["euler"]:
        fails.append("richardson drift %.3e exceeds euler %.3e; the "
                     "higher-order scheme is not better"
                     % (drift["richardson"], drift["euler"]))
    return fails


def test_force(exe_cmd):
    # constant external force on the COM: velocity Verlet is exact,
    # x = x0 + v0*t + a*t^2/2, v = v0 + a*t
    fx = 1.0e-21
    rc, out = run_deck(exe_cmd, "in.test.ballistic",
                       extra=["-var", "fx", repr(fx)])
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    last = rows[-1]
    fails = []
    a = fx / 1.0e-22
    dt = 1.0e-4
    n = 1000
    t = n * dt
    xexp = 3.0 + 12.0 * t + 0.5 * a * t * t
    vexp = 12.0 + a * t
    if not approx(last["f_1[1]"], xexp, rel=1e-7):
        fails.append("xcm = %.12g, expected %.12g" % (last["f_1[1]"], xexp))
    if not approx(last["f_1[4]"], vexp, rel=1e-7):
        fails.append("vx = %.12g, expected %.12g" % (last["f_1[4]"], vexp))
    if not approx(last["f_1[2]"], 4.7, rel=1e-7):
        fails.append("ycm = %.12g, expected 4.7 (force is x only)"
                     % last["f_1[2]"])
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


def test_restitution(exe_cmd):
    # analytic properties of the push-off contact model, measured as the
    # rebound speed of a body launched at a static wall with no gas
    def bounce(pstyle, pk, pdamp, v0):
        rc, out = run_deck(exe_cmd, "in.test.restitution",
                           extra=["-var", "pstyle", pstyle, "-var", "pk", pk,
                                  "-var", "pdamp", pdamp, "-var", "v0",
                                  repr(v0)])
        if rc:
            return None
        rows = parse_stats(out)
        if not rows:
            return None
        vout = rows[-1]["f_1[4]"]
        if vout >= 0.0:          # never rebounded: passed through the wall
            return None
        return abs(vout) / v0

    fails = []

    # elastic contact conserves energy exactly: e = 1 for both force laws
    for pstyle, pk in (("linear", "1.0e-18"), ("hertz", "4.0e-18")):
        for v0 in (10.0, 25.0):
            e = bounce(pstyle, pk, "0.0", v0)
            if e is None:
                fails.append("%s elastic v0=%g: no rebound" % (pstyle, v0))
            elif not approx(e, 1.0, abs_=1e-3):
                fails.append("%s elastic v0=%g: restitution %.6f, expected 1"
                             % (pstyle, v0, e))

    # the wall is 0.2 thick, less than 2*cutoff, so a corner pt which
    # crosses the near face is within range of the far face too; contact
    # is one-sided, so the far face must not push the body on through
    # (v0=90 is within the 4-contact capacity of the linear spring)
    e = bounce("linear", "1.0e-18", "0.0", 90.0)
    if e is None:
        fails.append("linear elastic v0=90: no rebound, thin wall pushed "
                     "the body through instead of repelling it")
    elif not approx(e, 1.0, abs_=1e-3):
        fails.append("linear elastic v0=90: restitution %.6f, expected 1" % e)

    # a linear spring-dashpot has a restitution independent of impact
    # speed; this is the property which distinguishes it from Hertzian
    lin = [bounce("linear", "1.0e-18", "3.0e-21", v0) for v0 in (10.0, 25.0)]
    if None in lin:
        fails.append("linear damped: no rebound")
    elif not approx(lin[0], lin[1], rel=1e-4):
        fails.append("linear damped: restitution %.6f at v0=10 vs %.6f at "
                     "v0=25, should not depend on impact speed"
                     % (lin[0], lin[1]))

    # a Hertzian spring-dashpot restitution does depend on impact speed
    her = [bounce("hertz", "4.0e-18", "3.0e-21", v0) for v0 in (10.0, 25.0)]
    if None in her:
        fails.append("hertz damped: no rebound")
    elif approx(her[0], her[1], rel=1e-3):
        fails.append("hertz damped: restitution %.6f at v0=10 and %.6f at "
                     "v0=25 are equal; a Hertzian contact must depend on "
                     "impact speed" % (her[0], her[1]))

    # damping must actually remove energy, else the checks above are vacuous
    if lin[0] is not None and lin[0] > 0.9:
        fails.append("linear damped restitution %.4f is too close to elastic "
                     "for the test to be meaningful" % lin[0])
    return fails


def test_momentum(exe_cmd):
    rc, out = run_deck(exe_cmd, "in.test.momentum")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    mass_body = 1.0e-22

    dt = 1.0e-4                      # timestep of the deck

    def gas_plus_body(r):
        return MASS_N * FNUM * r["c_r"] + mass_body * r["f_1[4]"]

    # compute surf tallies the momentum the gas gave up on a step; with
    # velocity Verlet the body receives half of it at the end of that step
    # and half at the start of the next, so half a step's impulse
    # 0.5*dt*fcm is always in flight.  Including it makes the invariant
    # exact rather than accurate to the statistical size of one step's
    # transfer.

    def total_px(r):
        return gas_plus_body(r) + 0.5 * dt * r["f_1[7]"]

    p0 = total_px(rows[0])
    fails = []
    for r in rows:
        if not approx(total_px(r), p0, rel=1.0e-11):
            fails.append("step %d: lag-corrected total px = %.15g vs initial "
                         "%.15g" % (int(r["Step"]), total_px(r), p0))
            break

    # the uncorrected sum should drift by roughly one step's transfer and
    # no more: a bounded offset, not a leak
    raw = [abs(gas_plus_body(r) - gas_plus_body(rows[0])) for r in rows]
    if max(raw) > 0.05 * abs(gas_plus_body(rows[0])):
        fails.append("uncorrected momentum drifted by %.3g, far more than "
                     "one step's impulse; this looks like a leak"
                     % (max(raw) / abs(gas_plus_body(rows[0]))))
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


def compare_remap_modes(exe_cmd, deck, keys, labels):
    """Run deck in both remap modes, require identical final values."""
    results = {}
    fails = []
    for mode in ("cutcell", "incremental"):
        rc, out = run_deck(exe_cmd, deck, extra=["-var", "mode", mode])
        if rc:
            fails.append("mode %s: run failed with exit code %d" % (mode, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("mode %s: no stats output" % mode)
            continue
        results[mode] = rows
    if fails:
        return fails, results
    last_c = results["cutcell"][-1]
    last_i = results["incremental"][-1]
    for key, name in zip(keys, labels):
        if not approx(last_i[key], last_c[key], rel=1e-10, abs_=1e-13):
            fails.append("%s: incremental %.15g differs from cutcell %.15g"
                         % (name, last_i[key], last_c[key]))
    return fails, results


def test_remap(exe_cmd):
    # one gas-driven body: cutcell and incremental must give identical
    # trajectories, verifying the incremental re-cut against the full
    # rebuild; the body must actually have moved
    fails, results = compare_remap_modes(
        exe_cmd, "in.test.remap",
        ("f_1[1]", "f_1[2]", "f_1[15]"), ("xcm", "ycm", "omega"))
    if fails:
        return fails
    if abs(results["cutcell"][-1]["f_1[1]"] - 5.0) < 0.01:
        fails.append("body barely moved (xcm = %.6g), test is too weak"
                     % results["cutcell"][-1]["f_1[1]"])
    return fails


def test_staticdist(exe_cmd):
    # body next to a 200-segment static circle: the incremental re-cut
    # must handle cells holding many static surfs; with --dist on several
    # procs most static surfs are ghost surfs on any one proc, exercising
    # the mover's ghost-cell collision tests and the local body copies
    fails, results = compare_remap_modes(
        exe_cmd, "in.test.staticdist",
        ("f_1[1]", "f_1[2]", "f_1[15]"), ("xcm", "ycm", "omega"))
    if fails:
        return fails
    for mode in ("cutcell", "incremental"):
        ndel = results[mode][-1]["f_1"] - results[mode][0]["f_1"]
        if ndel != 0:
            fails.append("mode %s: %g particles deleted inside the body "
                         "during the run" % (mode, ndel))
    # the gas must actually have moved the body (it drifts ~0.006 in x
    # over the run), else the two modes agree trivially
    if abs(results["cutcell"][-1]["f_1[1]"] - 6.0) < 0.001:
        fails.append("body barely moved (xcm = %.6g), test is too weak"
                     % results["cutcell"][-1]["f_1[1]"])
    return fails


def test_staticdist3d(exe_cmd):
    # 3d: cube body drifting past a 1200-triangle static sphere, both
    # remap modes must agree and no particle may be deleted; with --dist
    # on several procs most sphere triangles are ghost surfs
    fails, results = compare_remap_modes(
        exe_cmd, "in.test.staticdist3d",
        ("f_1[1]", "f_1[2]", "f_1[3]", "f_1[13]", "f_1[14]", "f_1[15]"),
        ("xcm", "ycm", "zcm", "wx", "wy", "wz"))
    if fails:
        return fails
    for mode in ("cutcell", "incremental"):
        ndel = results[mode][-1]["f_1"] - results[mode][0]["f_1"]
        if ndel != 0:
            fails.append("mode %s: %g particles deleted inside the body "
                         "during the run" % (mode, ndel))
    if abs(results["cutcell"][-1]["f_1[1]"] - 4.0) < 0.001:
        fails.append("body barely moved (xcm = %.6g), test is too weak"
                     % results["cutcell"][-1]["f_1[1]"])
    return fails


def test_restart(exe_cmd):
    # restart continuation: a deterministic push-off run split across a
    # write_restart/read_restart must reproduce the one-shot trajectory.
    # the split points straddle the contact with the wall, which is the
    # sensitive case: the body moves on a step under the force and torque
    # accumulated on the previous one, so a continuation which resumed
    # with zero force would lose that impulse.  two bodies, each with its
    # own outfile, so that the state written for a body other than the
    # last-defined one is also checked to be the complete end-of-step state
    total = 1500
    fails = []
    rc, out = run_deck(exe_cmd, "in.test.restart.oneshot",
                       extra=["-var", "nrun", str(total)])
    if rc:
        return ["one-shot run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["one-shot run produced no stats output"]
    ref = rows[-1]

    # the body must actually reach the wall and rebound, else the test
    # never exercises a restart under load

    if ref["f_1[4]"] > -40.0:
        return ["one-shot final vx = %.6g, body did not rebound; test "
                "geometry is broken" % ref["f_1[4]"]]

    # 480 falls inside body 1's contact with the wall (steps ~405-565),
    # so its stored force and torque are nonzero at that split
    for split in (300, 480, 700, 1100):
        rc, out1 = run_deck(exe_cmd, "in.test.restart.part1",
                            extra=["-var", "nrun", str(split)])
        if rc:
            fails.append("split %d: first half failed with exit code %d"
                         % (split, rc))
            continue
        rows1 = parse_stats(out1)
        rc, out2 = run_deck(exe_cmd, "in.test.restart.part2",
                            extra=["-var", "nrun", str(total - split)])
        if rc:
            fails.append("split %d: continuation failed with exit code %d"
                         % (split, rc))
            continue
        rows2 = parse_stats(out2)
        if not rows1 or not rows2:
            fails.append("split %d: a half produced no stats output"
                         % split)
            continue
        keys = (("f_1[1]", "xcm"), ("f_1[4]", "vx"), ("f_1[15]", "omega"),
                ("f_2[1]", "xcm2"), ("f_2[4]", "vx2"), ("f_2[15]", "omega2"))

        # the state read back from the outfile at the start of the
        # continuation must equal the state at the end of the first half
        # exactly: the outfile stores 17 digits, which round-trip a double
        # (15 digits, the previous format, lost the last few ulp)

        first = rows2[0]
        end1 = rows1[-1]
        for key, name in keys:
            if name.startswith("omega"):
                # omega is not stored: it is re-derived from the angular
                # momentum through the inertia eigensolver, to round-off
                if not approx(first[key], end1[key], rel=1e-12, abs_=1e-300):
                    fails.append("split %d: %s re-derived from the outfile "
                                 "as %.17g, was %.17g"
                                 % (split, name, first[key], end1[key]))
            elif first[key] != end1[key]:
                fails.append("split %d: %s read back from the outfile as "
                             "%.17g, written from %.17g"
                             % (split, name, first[key], end1[key]))

        # the rest of the continuation re-derives the body-frame geometry
        # from the restarted surfs, so it tracks the one-shot run to
        # round-off rather than exactly

        last = rows2[-1]
        for key, name in keys:
            if not approx(last[key], ref[key], rel=1e-12, abs_=1e-300):
                fails.append("split %d: %s = %.12g differs from one-shot "
                             "%.12g" % (split, name, last[key], ref[key]))

    for f in glob.glob(os.path.join(THISDIR, "tmp.rigid*")):
        os.remove(f)
    return fails


def test_gridchange(exe_cmd):
    # the grid changing underneath the body must not perturb it: a
    # no-particle push-off trajectory is identical with and without
    # fix balance (random style, full rebuild every 25 steps) and
    # fix adapt (refine/coarsen on the body surfs every 100 steps)
    results = {}
    fails = []
    for pert in ("none", "balance", "balancecell", "adapt"):
        rc, out = run_deck(exe_cmd, "in.test.gridchange",
                           extra=["-var", "pert", pert])
        if rc:
            fails.append("pert %s: run failed with exit code %d" % (pert, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("pert %s: no stats output" % pert)
            continue
        results[pert] = rows
    if fails:
        return fails
    ref = results["none"]
    for pert in ("balance", "balancecell", "adapt"):
        if len(results[pert]) != len(ref):
            fails.append("pert %s: %d stats rows vs %d unperturbed"
                         % (pert, len(results[pert]), len(ref)))
            continue
        for r, r0 in zip(results[pert], ref):
            for key in ("f_1[1]", "f_1[4]", "f_1[15]", "f_1[20]"):
                if not approx(r[key], r0[key], rel=1e-12, abs_=1e-30):
                    fails.append("pert %s, step %d: %s = %.15g differs from "
                                 "unperturbed %.15g"
                                 % (pert, int(r["Step"]), key, r[key], r0[key]))
                    break
            if fails:
                break
    # the body must actually have bounced, else the test is vacuous
    if ref[-1]["f_1[4]"] > -40.0:
        fails.append("final vx = %.6g, body did not rebound; test geometry "
                     "is broken" % ref[-1]["f_1[4]"])
    return fails


def test_splitcell(exe_cmd):
    # body sweeping alongside a diagonal wall that creates split cells:
    # particles entering swept split cells must be reflected, not
    # overrun, so the deletion count must not grow; both remap modes
    # must agree (incremental falls back to a full re-map near split
    # cells, which must not change the result)
    fails, results = compare_remap_modes(
        exe_cmd, "in.test.splitcell",
        ("f_1[1]", "f_1[2]"), ("xcm", "ycm"))
    if fails:
        return fails
    for mode in ("cutcell", "incremental"):
        ndel = results[mode][-1]["f_1"] - results[mode][0]["f_1"]
        if ndel != 0:
            fails.append("mode %s: %g particles overrun by the body in "
                         "split cells" % (mode, ndel))
        if results[mode][-1]["Nscoll"] == 0:
            fails.append("mode %s: no surface collisions, test geometry "
                         "is broken" % mode)
    return fails


def test_splitbalance(exe_cmd):
    # split cells plus a fix balance in the same run: after a full re-map
    # fix rigid reassigns split-cell particles to sub cells, which leaves
    # them unsorted, and the balance later in the step must re-sort
    # before migrating cells.  the box is periodic with no emission and
    # no deletion, so the particle count must stay at its initial value;
    # with stale lists a random balance quadrupled it on 4 ranks
    fails = []
    for mode in ("cutcell", "incremental"):
        rc, out = run_deck(exe_cmd, "in.test.splitbalance",
                           extra=["-var", "mode", mode])
        if rc:
            fails.append("mode %s: run failed with exit code %d" % (mode, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("mode %s: no stats output" % mode)
            continue
        np0 = rows[0]["Np"]
        for row in rows[1:]:
            if row["Np"] != np0:
                fails.append("mode %s: step %d has %d particles, started "
                             "with %d" % (mode, row["Step"], row["Np"], np0))
                break
        if rows[-1]["f_1"] != 0:
            fails.append("mode %s: %g particles deleted inside the body"
                         % (mode, rows[-1]["f_1"]))
    return fails


def test_transplane(exe_cmd):
    # a fast body sweeps over and then vacates cells cut only by a
    # transparent plane; the total flow volume per step must agree
    # between the incremental and cutcell remap modes
    vols = {}
    fails = []
    for mode in ("cutcell", "incremental"):
        rc, out = run_deck(exe_cmd, "in.test.transplane",
                           extra=["-var", "mode", mode])
        if rc:
            fails.append("mode %s: run failed with exit code %d" % (mode, rc))
            continue
        rows = parse_stats(out)
        if len(rows) < 15:
            fails.append("mode %s: expected 15 stats rows, got %d"
                         % (mode, len(rows)))
            continue
        vols[mode] = [r["c_tvol"] for r in rows]
    if fails:
        return fails
    for i, (vi, vc) in enumerate(zip(vols["incremental"], vols["cutcell"])):
        if not approx(vi, vc, rel=1e-12):
            fails.append("step %d: incremental flow volume %.17g differs "
                         "from cutcell %.17g" % (i, vi, vc))
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


def test_recoil(exe_cmd):
    """single specular collision of a particle with a body only 10x heavier:
    exact two-body elastic result, frictionless normal impulse"""
    m = 1.0e10 * 2.325e-26          # fnum * m_N
    u = 1000.0
    fails = []

    # 2d: square, I = M/6, hit at r = (-0.5, yhit-5, 0), n = (-1,0,0)
    # mratio = 1 is the equal-mass exchange: particle stops, body takes u
    for yhit, mratio in ((5.0, 10.0), (5.4, 10.0), (5.0, 1.0)):
        M = mratio * m
        I2 = M / 6.0
        ry = yhit - 5.0
        J = 2.0 * m * u / (1.0 + m * (1.0 / M + ry * ry / I2))
        vx = u - J / m
        vcm = J / M
        omega = -ry * J / I2
        label = "2d yhit=%g M/m=%g" % (yhit, mratio)
        rc, out = run_deck(exe_cmd, "in.test.recoil",
                           ["-var", "yhit", str(yhit),
                            "-var", "mratio", str(mratio)])
        if rc != 0:
            fails.append("%s: run failed with exit code %d" % (label, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("%s: no stats output" % label)
            continue
        row = rows[-1]
        if row["Np"] != 1:
            fails.append("%s: particle lost" % label)
        for key, want in (("c_rvx", vx), ("c_rvy", 0.0), ("f_1[4]", vcm),
                          ("f_1[5]", 0.0), ("f_1[15]", omega)):
            if not approx(row[key], want, rel=1.0e-10, abs_=1.0e-10 * u):
                fails.append("%s: %s = %.12g, expected %.12g"
                             % (label, key, row[key], want))
        e0 = 0.5 * m * u * u
        e1 = 0.5 * m * (row["c_rvx"] ** 2 + row["c_rvy"] ** 2) \
            + 0.5 * M * (row["f_1[4]"] ** 2 + row["f_1[5]"] ** 2) \
            + 0.5 * I2 * row["f_1[15]"] ** 2
        if not approx(e1, e0, rel=1.0e-10):
            fails.append("%s: energy %.12g vs %.12g" % (label, e1, e0))

    # 3d: unit cube, I = M/6, hit at r = (-0.5, 0.3, 0.2), n = (-1,0,0)
    M = 10.0 * m
    I3 = M / 6.0
    ry, rz = 0.3, 0.2
    J = 2.0 * m * u / (1.0 + m * (1.0 / M + (ry * ry + rz * rz) / I3))
    vx = u - J / m
    vcm = J / M
    wy, wz = rz * J / I3, -ry * J / I3
    rc, out = run_deck(exe_cmd, "in.test.recoil3d")
    if rc != 0:
        fails.append("3d: run failed with exit code %d" % rc)
        return fails
    rows = parse_stats(out)
    if not rows:
        fails.append("3d: no stats output")
        return fails
    row = rows[-1]
    for key, want in (("c_rvx", vx), ("c_rvy", 0.0), ("c_rvz", 0.0),
                      ("f_1[4]", vcm), ("f_1[5]", 0.0), ("f_1[6]", 0.0),
                      ("f_1[13]", 0.0), ("f_1[14]", wy), ("f_1[15]", wz)):
        if not approx(row[key], want, rel=1.0e-10, abs_=1.0e-10 * u):
            fails.append("3d: %s = %.12g, expected %.12g"
                         % (key, row[key], want))
    e0 = 0.5 * m * u * u
    e1 = 0.5 * m * (row["c_rvx"] ** 2 + row["c_rvy"] ** 2 + row["c_rvz"] ** 2) \
        + 0.5 * M * (row["f_1[4]"] ** 2 + row["f_1[5]"] ** 2 + row["f_1[6]"] ** 2) \
        + 0.5 * I3 * (row["f_1[13]"] ** 2 + row["f_1[14]"] ** 2
                      + row["f_1[15]"] ** 2)
    if not approx(e1, e0, rel=1.0e-10):
        fails.append("3d: energy %.12g vs %.12g" % (e1, e0))
    return fails


def test_exitbox(exe_cmd):
    # two bodies leaving through opposite faces: the run must complete
    # (no false "surfs may enclose the box" error) and the motion is
    # ballistic, ycm = y0 + vy * t with t = 200 * 1e-4
    rc, out = run_deck(exe_cmd, "in.test.exitbox")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    last = rows[-1]
    fails = []
    if not approx(last["f_1[2]"], 3.0 - 400.0 * 0.02, rel=1e-10):
        fails.append("body 1 ycm = %.12g, expected %.12g"
                     % (last["f_1[2]"], 3.0 - 400.0 * 0.02))
    if not approx(last["f_2[2]"], 7.0 + 400.0 * 0.02, rel=1e-10):
        fails.append("body 2 ycm = %.12g, expected %.12g"
                     % (last["f_2[2]"], 7.0 + 400.0 * 0.02))
    return fails


def negative_test(exe_cmd, deck, message):
    rc, out = run_deck(exe_cmd, deck)
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


def test_facetbounce(exe_cmd):
    # a square rebounding from a 200-segment static circle: kinetic
    # energy in free flight after the contact must equal the launch
    # energy (contacts with a faceted surface must be conservative)
    rc, out = run_deck(exe_cmd, "in.test.facetbounce")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if len(rows) < 21:
        return ["expected 21 stats rows, got %d" % len(rows)]
    free = [r for r in rows if r["f_1[20]"] == 0.0 and r["f_1[21]"] == 0.0]
    if len(free) == len(rows):
        return ["the body never touched the circle, test geometry is broken"]
    if rows[-1]["f_1[4]"] >= 0.0:
        return ["body did not rebound (vx = %.6g)" % rows[-1]["f_1[4]"]]
    e0 = rows[0]["v_ke"]
    e1 = free[-1]["v_ke"]
    if not approx(e1, e0, rel=1e-4):
        return ["kinetic energy after the rebound %.10g vs %.10g before "
                "(%.3g%%)" % (e1, e0, 100.0 * (e1 / e0 - 1.0))]
    return []


def test_vacate(exe_cmd):
    # gas collisions with a fast body spanning whole interior cells: the
    # cells the body partly vacates within a step hold particles at zero
    # flow volume until the re-cut; the run must complete with the
    # particle count constant and nothing deleted after step 0
    rc, out = run_deck(exe_cmd, "in.test.vacate")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if len(rows) < 11:
        return ["expected 11 stats rows, got %d" % len(rows)]
    fails = []
    np0 = rows[0]["Np"]
    for r in rows[1:]:
        if r["Np"] != np0:
            fails.append("step %d: Np %g != %g" % (r["Step"], r["Np"], np0))
            break
    if rows[-1]["f_1"] != rows[0]["f_1"]:
        fails.append("%g particles deleted during the run"
                     % (rows[-1]["f_1"] - rows[0]["f_1"]))
    if rows[-1]["f_1[1]"] < 5.5:
        fails.append("body barely moved (xcm = %.6g), test is too weak"
                     % rows[-1]["f_1[1]"])
    return fails


def test_prenofix(exe_cmd):
    return negative_test(exe_cmd, "in.test.prenofix",
                         "not initialized before the run")


def test_refix(exe_cmd):
    # a fix rigid re-defined between runs, then balance_grid before the
    # next run: the rigid map rebuild must not reach the deleted fix;
    # the second run continues the ballistic body from where the
    # re-definition placed it (xcm 5.08 + 20*0.0001*20)
    rc, out = run_deck(exe_cmd, "in.test.refix")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows or rows[-1]["Step"] != 60:
        return ["run did not reach step 60"]
    x = rows[-1]["f_1[1]"]
    if not approx(x, 5.08 + 20.0 * 1.0e-4 * 20, rel=1e-12):
        return ["xcm after the second run %.17g, expected %.17g"
                % (x, 5.08 + 20.0 * 1.0e-4 * 20)]
    return []


def test_inward(exe_cmd):
    # a body traversed the wrong way round (normals pointing into its
    # interior, a container) is rejected in 2d and 3d
    fails = negative_test(exe_cmd, "in.test.inward", "normals point inward")
    fails += negative_test(exe_cmd, "in.test.inward3d",
                           "normals point inward")
    return fails


def test_modifyafter(exe_cmd):
    return negative_test(exe_cmd, "in.test.modifyafter",
                         "attributes were changed")


def test_wallmotion(exe_cmd):
    return negative_test(exe_cmd, "in.test.wallmotion", "own wall motion")


def timestep_independence(exe_cmd, deck, keys, coarse, fine, tol=1.0e-9,
                          hitkey="c_rvx"):
    """run the same physical problem at two timesteps and require the same
    answer: the deck is built so that the only timestep-dependent piece is
    the moving-surf collision test; the particle is launched in +x and
    must have been turned around by the body (hitkey < 0), else the
    moving-surf test was never exercised"""
    fails = []
    last = {}
    for label, (dt, nsteps) in (("coarse", coarse), ("fine", fine)):
        rc, out = run_deck(exe_cmd, deck, ["-var", "dt", dt,
                                           "-var", "nsteps", nsteps])
        if rc:
            fails.append("%s: run failed with exit code %d" % (label, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("%s: no stats output" % label)
            continue
        if rows[-1]["Np"] != 1:
            fails.append("%s: the particle was lost" % label)
            continue
        if rows[-1][hitkey] >= 0.0:
            fails.append("%s: %s = %.6g, the particle never hit the body"
                         % (label, hitkey, rows[-1][hitkey]))
            continue
        last[label] = rows[-1]
    if fails:
        return fails
    for key in keys:
        c, f = last["coarse"][key], last["fine"][key]
        if not approx(c, f, rel=tol, abs_=tol):
            fails.append("%s: %.14g at the coarse timestep, %.14g at the "
                         "fine one" % (key, c, f))
    return fails


def test_rotwall(exe_cmd):
    # the body is heavy enough that the collision does not perturb it and
    # spins at a constant rate, so its pose is exact at any timestep, and
    # the particle is ballistic on either side of the collision.  the hit
    # fraction of the moving-surf test is then the only thing that can
    # make the two runs differ
    return timestep_independence(
        exe_cmd, "in.test.rotwall",
        ("c_rx", "c_ry", "c_rvx", "c_rvy"),
        ("1.0e-3", "5"), ("2.0e-5", "250"))


def test_rotwall3d(exe_cmd):
    return timestep_independence(
        exe_cmd, "in.test.rotwall3d",
        ("c_rx", "c_ry", "c_rz", "c_rvx", "c_rvy", "c_rvz"),
        ("1.0e-3", "5"), ("2.0e-5", "250"))


def test_customemit(exe_cmd):
    # a fix emit/surf which spreads a custom per-surf attribute, defined
    # before fix rigid, across two runs with distributed surfs: the
    # per-surf status flags fix rigid resets after a grid rebuild gate a
    # collective re-spread in the emit fix's init, so they must be reset
    # on every rank or none, else the second run's init hangs
    rc, out = run_deck(exe_cmd, "in.test.customemit")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if len(rows) < 2:
        return ["fewer than two stats rows, the second run did not start"]
    return []


def test_emitsurf(exe_cmd):
    return negative_test(exe_cmd, "in.test.emitsurf",
                         "cannot emit from fix rigid body surfs")


def test_renumber(exe_cmd):
    fails = negative_test(exe_cmd, "in.test.renumber", "were renumbered")
    fails += negative_test(exe_cmd, "in.test.renumber2", "were renumbered")
    return fails


def test_axistuck(exe_cmd):
    # no rigid body: the mover's moving-surf re-hit rule must leave the
    # legitimate repeated hits on one static line alone (axisymmetric)
    rc, out = run_deck(exe_cmd, "in.test.axistuck")
    if rc:
        return ["run failed with exit code %d" % rc]
    rows = parse_stats(out)
    if not rows:
        return ["no stats output"]
    fails = []
    for row in rows:
        if row["Np"] != 2000:
            fails.append("step %d: np = %d, particles were deleted"
                         % (row["Step"], row["Np"]))
            break
    if rows[-1]["Nscoll"] < 1000:
        fails.append("too few surf collisions (%d) for the test to be "
                     "meaningful" % rows[-1]["Nscoll"])
    return fails


def test_tallyorder(exe_cmd):
    # a fix ave/surf on a static wall must see the same per-step hit
    # counts whether it is defined before or after fix rigid
    results = {}
    fails = []
    for order in ("0", "1"):
        rc, out = run_deck(exe_cmd, "in.test.tallyorder",
                           extra=["-var", "order", order])
        if rc:
            fails.append("order %s: run failed with exit code %d"
                         % (order, rc))
            continue
        rows = parse_stats(out)
        if not rows:
            fails.append("order %s: no stats output" % order)
            continue
        results[order] = rows
    if fails:
        return fails
    for r0, r1 in zip(results["0"], results["1"]):
        if r0["c_red"] != r1["c_red"]:
            fails.append("step %d: wall hits %g (ave/surf before fix rigid)"
                         " vs %g (after)" % (r0["Step"], r0["c_red"],
                                             r1["c_red"]))
            break
    if results["0"][-1]["c_red"] < 10:
        fails.append("too few wall hits (%g) for the test to be meaningful"
                     % results["0"][-1]["c_red"])
    return fails


def test_badinfile(exe_cmd):
    return negative_test(exe_cmd, "in.test.badinfile",
                         "Invalid floating point number")


def test_mixture(exe_cmd):
    return negative_test(exe_cmd, "in.test.mixture",
                         "mixture must contain all species")


def test_zerothick(exe_cmd):
    fails = negative_test(exe_cmd, "in.test.zerothick", "encloses zero area")
    fails += negative_test(exe_cmd, "in.test.zerothick3d",
                           "encloses zero volume")
    return fails


def read_state(name):
    """Read the 16 (or 22) body parameters from a fix rigid outfile.
       Returns the list of floats, or None if the file is missing."""
    path = os.path.join(THISDIR, name)
    if not os.path.exists(path):
        return None
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line:
            continue
        return [float(w) for w in line.split()]
    return None


def test_density(exe_cmd):
    """dstyle = density on a unit cube: mass, COM and moi from geometry.

    A cube is exactly representable by flat triangles, so every value is
    analytic and the only error is round-off.  The tolerance is set by the
    outfile's 17-digit format, not by any property of the shape."""
    fails = []
    out_name = "tmp.density.state"
    path = os.path.join(THISDIR, out_name)
    if os.path.exists(path):
        os.remove(path)

    rc, out = run_deck(exe_cmd, "in.test.density")
    if rc:
        return ["run failed with exit code %d" % rc]

    v = read_state(out_name)
    if v is None:
        return ["no state file written"]
    if len(v) < 16:
        return ["state file has %d values, expected at least 16" % len(v)]

    # unit cube, density 1
    TOL = 1.0e-12
    if abs(v[0] - 1.0) > TOL:
        fails.append("mass %.17g != 1 (unit cube, density 1)" % v[0])
    for k, nm in enumerate("xyz"):
        if abs(v[1+k]) > TOL:
            fails.append("com %s %.3e != 0 (cube is centred on the origin)"
                         % (nm, v[1+k]))
    for k, nm in enumerate(("ixx", "iyy", "izz")):
        if abs(v[4+k] - 1.0/6.0) > TOL:
            fails.append("%s %.17g != M/6 = %.17g" % (nm, v[4+k], 1.0/6.0))
    for k, nm in enumerate(("ixy", "ixz", "iyz")):
        if abs(v[7+k]) > TOL:
            fails.append("%s %.3e != 0 (cube has no products of inertia)"
                         % (nm, v[7+k]))

    os.remove(path)
    return fails


def test_density2d(exe_cmd):
    """dstyle = density in 2d, plus the vcom/angmom keywords.

    Unit square plate: izz = M/6, ixx = iyy = izz/2, all products zero,
    and ixz = iyz = 0 exactly as 2d requires.  vcom and angmom are user
    input even under density style, so they must round-trip unchanged."""
    fails = []
    out_name = "tmp.density2d.state"
    path = os.path.join(THISDIR, out_name)
    if os.path.exists(path):
        os.remove(path)

    rc, out = run_deck(exe_cmd, "in.test.density2d")
    if rc:
        return ["run failed with exit code %d" % rc]

    v = read_state(out_name)
    if v is None:
        return ["no state file written"]
    if len(v) < 16:
        return ["state file has %d values, expected at least 16" % len(v)]

    TOL = 1.0e-12
    if abs(v[0] - 1.0) > TOL:
        fails.append("mass %.17g != 1 (unit square, density 1)" % v[0])

    # the body moves during the one step it is run, so the COM is compared
    # against its start-of-step value plus the drift from vcom
    if abs(v[3]) > TOL:
        fails.append("com z %.3e != 0 for a 2d body" % v[3])

    if abs(v[4] - 1.0/12.0) > TOL:
        fails.append("ixx %.17g != M/12 = %.17g" % (v[4], 1.0/12.0))
    if abs(v[5] - 1.0/12.0) > TOL:
        fails.append("iyy %.17g != M/12 = %.17g" % (v[5], 1.0/12.0))
    if abs(v[6] - 1.0/6.0) > TOL:
        fails.append("izz %.17g != M/6 = %.17g" % (v[6], 1.0/6.0))
    if abs(v[7]) > TOL:
        fails.append("ixy %.3e != 0 (square has no product of inertia)" % v[7])

    # a 2d body must have these exactly zero, not merely small: setup_body()
    # requires a principal axis along z
    if v[8] != 0.0 or v[9] != 0.0:
        fails.append("ixz,iyz = %.3e,%.3e must be exactly 0 for 2d"
                     % (v[8], v[9]))

    # vcom and angmom are user input, unchanged by the geometry integration
    for k, (got, want, nm) in enumerate(
            ((v[10], 12.0, "vxcm"), (v[11], -3.0, "vycm"),
             (v[12], 0.0, "vzcm"), (v[15], 5.0e-3, "lz"))):
        if abs(got - want) > 1.0e-12 * max(abs(want), 1.0):
            fails.append("%s %.17g != %.17g as given" % (nm, got, want))
    if v[13] != 0.0 or v[14] != 0.0:
        fails.append("lx,ly = %.3e,%.3e must be exactly 0 for 2d"
                     % (v[13], v[14]))

    os.remove(path)
    return fails


def test_baddensity(exe_cmd):
    return negative_test(exe_cmd, "in.test.baddensity",
                         "body density must be positive")


def test_badvcom(exe_cmd):
    return negative_test(exe_cmd, "in.test.badvcom",
                         "vcom keyword requires density style")


TESTS = [
    ("ballistic", test_ballistic),
    ("force", test_force),
    ("rotation", test_rotation),
    ("bounce", test_bounce),
    ("restitution", test_restitution),
    ("momentum", test_momentum),
    ("recoil", test_recoil),
    ("overrun", test_overrun),
    ("remap", test_remap),
    ("multiremap", test_multiremap),
    ("transplane", test_transplane),
    ("splitbalance", test_splitbalance),
    ("staticdist", test_staticdist),
    ("staticdist3d", test_staticdist3d),
    ("splitcell", test_splitcell),
    ("gridchange", test_gridchange),
    ("exitbox", test_exitbox),
    ("restart", test_restart),
    ("twobody", test_twobody),
    ("pushpair", test_pushpair),
    ("badmoi", test_badmoi),
    ("notwatertight", test_notwatertight),
    ("zerothick", test_zerothick),
    ("inward", test_inward),
    ("refix", test_refix),
    ("prenofix", test_prenofix),
    ("vacate", test_vacate),
    ("facetbounce", test_facetbounce),
    ("modifyafter", test_modifyafter),
    ("wallmotion", test_wallmotion),
    ("customemit", test_customemit),
    ("emitsurf", test_emitsurf),
    ("renumber", test_renumber),
    ("mixture", test_mixture),
    ("badinfile", test_badinfile),
    ("rotwall", test_rotwall),
    ("rotwall3d", test_rotwall3d),
    ("axistuck", test_axistuck),
    ("tallyorder", test_tallyorder),
    ("density", test_density),
    ("density2d", test_density2d),
    ("baddensity", test_baddensity),
    ("badvcom", test_badvcom),
]

# tests whose decks support -var dist 1 (global surfs explicit/distributed)
# remap, multiremap, staticdist, and splitcell verify the incremental
# re-cut against the full rebuild in distributed mode; staticdist is the
# one whose static surfs are not local on every proc when run on
# several procs

DIST_TESTS = {"ballistic", "force", "rotation", "bounce", "restitution",
              "momentum",
              "overrun",
              "remap", "multiremap", "transplane", "staticdist",
              "staticdist3d",
              "splitcell", "gridchange", "exitbox", "twobody", "pushpair",
              "tallyorder", "rotwall", "rotwall3d", "customemit",
              "splitbalance",
              "vacate", "facetbounce"}


def main():
    global run_deck
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True,
                        help="path to SPARTA executable")
    parser.add_argument("--mpi", default="",
                        help='MPI launcher prefix, e.g. "mpirun -np 4"')
    parser.add_argument("--args", default="",
                        help='extra SPARTA command-line args placed after '
                             'the executable, e.g. "-k on -sf kk"')
    parser.add_argument("--tests", default="",
                        help="comma-separated subset of tests to run")
    parser.add_argument("--dist", action="store_true",
                        help="run with distributed surfs "
                             "(global surfs explicit/distributed)")
    args = parser.parse_args()

    exe_cmd = (shlex.split(args.mpi) + [os.path.abspath(args.exe)] +
               shlex.split(args.args))

    subset = None
    if args.tests:
        subset = set(args.tests.split(","))

    if args.dist:
        if subset is None:
            subset = set(DIST_TESTS)
        else:
            subset &= DIST_TESTS
        base_run_deck = run_deck

        def dist_run_deck(exe_cmd, deck, extra=None):
            extra = (extra or []) + ["-var", "dist", "1"]
            return base_run_deck(exe_cmd, deck, extra)

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
