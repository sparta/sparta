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
    # constant external force on the COM: the semi-implicit Euler
    # trajectory is x_n = x0 + n*v0*dt + a*dt^2*n*(n+1)/2, v_n = v0 + n*a*dt
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
    xexp = 3.0 + n * 12.0 * dt + a * dt * dt * n * (n + 1) / 2.0
    vexp = 12.0 + n * a * dt
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

    # compute surf tallies the momentum the gas gave up on a step and fix
    # rigid applies it to the body on the next one, so one step's impulse
    # dt*fcm is always in flight.  Including it makes the invariant exact
    # rather than accurate to the statistical size of one step's transfer.

    def total_px(r):
        return gas_plus_body(r) + dt * r["f_1[7]"]

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
    # with zero force would lose that impulse
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

    for split in (300, 700, 1100):
        rc, _ = run_deck(exe_cmd, "in.test.restart.part1",
                         extra=["-var", "nrun", str(split)])
        if rc:
            fails.append("split %d: first half failed with exit code %d"
                         % (split, rc))
            continue
        rc, out2 = run_deck(exe_cmd, "in.test.restart.part2",
                            extra=["-var", "nrun", str(total - split)])
        if rc:
            fails.append("split %d: continuation failed with exit code %d"
                         % (split, rc))
            continue
        rows2 = parse_stats(out2)
        if not rows2:
            fails.append("split %d: continuation produced no stats output"
                         % split)
            continue
        last = rows2[-1]
        for key, name in (("f_1[1]", "xcm"), ("f_1[4]", "vx"),
                          ("f_1[15]", "omega")):
            if not approx(last[key], ref[key], rel=1e-7, abs_=1e-12):
                fails.append("split %d: %s = %.12g differs from one-shot "
                             "%.12g" % (split, name, last[key], ref[key]))

    for f in glob.glob(os.path.join(THISDIR, "tmp.rigid.*")):
        os.remove(f)
    return fails


def test_gridchange(exe_cmd):
    # the grid changing underneath the body must not perturb it: a
    # no-particle push-off trajectory is identical with and without
    # fix balance (random style, full rebuild every 25 steps) and
    # fix adapt (refine/coarsen on the body surfs every 100 steps)
    results = {}
    fails = []
    for pert in ("none", "balance", "adapt"):
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
    for pert in ("balance", "adapt"):
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


def test_zerothick(exe_cmd):
    fails = negative_test(exe_cmd, "in.test.zerothick", "encloses zero area")
    fails += negative_test(exe_cmd, "in.test.zerothick3d",
                           "encloses zero volume")
    return fails


TESTS = [
    ("ballistic", test_ballistic),
    ("force", test_force),
    ("rotation", test_rotation),
    ("bounce", test_bounce),
    ("restitution", test_restitution),
    ("momentum", test_momentum),
    ("overrun", test_overrun),
    ("remap", test_remap),
    ("multiremap", test_multiremap),
    ("staticdist", test_staticdist),
    ("staticdist3d", test_staticdist3d),
    ("splitcell", test_splitcell),
    ("gridchange", test_gridchange),
    ("restart", test_restart),
    ("twobody", test_twobody),
    ("pushpair", test_pushpair),
    ("badmoi", test_badmoi),
    ("notwatertight", test_notwatertight),
    ("zerothick", test_zerothick),
]

# tests whose decks support -var dist 1 (global surfs explicit/distributed)
# remap, multiremap, staticdist, and splitcell verify the incremental
# re-cut against the full rebuild in distributed mode; staticdist is the
# one whose static surfs are not local on every proc when run on
# several procs

DIST_TESTS = {"ballistic", "force", "rotation", "bounce", "restitution",
              "momentum",
              "overrun",
              "remap", "multiremap", "staticdist", "staticdist3d",
              "splitcell", "gridchange", "twobody", "pushpair"}


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
