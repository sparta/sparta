#!/usr/bin/env python3
"""
audit_conservation.py - independent, machine-precision audit of SWPM
                        conservation from raw particle dumps.

Motivation: the regression decks and verify_swpm.py check conservation
through SPARTA's weighted computes.  This audit removes that dependency:
it reads raw per-particle data (velocity, erot, evib, and the stochastic
weight custom attribute) from dump files and evaluates every moment in
Python, so an error in the SWPM dynamics cannot be masked by a matching
error in the diagnostics.

Part A isolates the REDUCTION operator: collisions are disabled with a
near-zero VSS diameter (no collisions => no splits either), the reduction
group is the whole cell (Ngmax >> N), and dumps are taken before and after
the single merge event.  For each reduction scheme it asserts, at machine
precision:

  scheme    conserved exactly                       provably NOT conserved
  ------    -----------------------------------     ----------------------
  energy    W, P(momentum), E_tr, E_rot, E_vib      stress deviator, q
  heat      + heat-flux vector q                    stress deviator
  stress    + full stress tensor P_ij (and q)       4th moment (kurtosis)

The "NOT conserved" assertions matter: they prove the audit is sensitive
enough to see a violation if one existed, and they mark the exact boundary
of each scheme's guarantee.

Part B audits the FULL method (splitting + equal-weight collisions +
reduction, translational and polyatomic) over a long run: total weight,
momentum, and total energy (translational + rotational + vibrational)
must show only floating-point-accumulation drift (<1e-9 relative).
Collisions legitimately relax higher moments, so only the true collision
invariants are asserted here.

Usage:
    python3 audit_conservation.py /path/to/spa_binary
    python3 audit_conservation.py "mpirun -np 4 /path/to/spa_binary"

Exit code 0 = all checks passed.
"""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

KB = 1.380658e-23
MASS_N2 = 4.65e-26

# N2 species with vibration; VibRel raised so Part B exchanges vib energy
SPECIES_FILE = """# N2 test species (fast vibrational relaxation, audit only)
N2  28.016   4.65E-26  2    0.2   2    0.2    3371.0    1.0      0.0
"""

# near-zero diameter: collision cross-section ~1e-60 m^2 => nattempt == 0,
# so the only weight/velocity-changing operation left is reduction
VSS_NOCOLL = """# N2 VSS with near-zero diameter: disables collisions (audit only)
N2      1.0E-30    0.74    273.15    1.0
"""

VSS_REAL = """# N2 VSS (Gopalan 2016)
N2      4.110E-10    0.740    273.15    1.360
"""

# ---------------------------------------------------------------------------
# input decks
# ---------------------------------------------------------------------------

PART_A_DECK = """
seed               12137
dimension          3

species            n2.species N2
mixture            air N2 vstream 30.0 -20.0 50.0 temp 5000.0
mixture            air N2 frac 1.0

boundary           p p p
create_box         0.0 1.0e-5 0.0 1.0e-5 0.0 1.0e-5
create_grid        1 1 1
balance_grid       rcb part

global             nrho 4.247e22
global             fnum {fnum}

collide            vss air n2.nocoll.vss
collide_modify     remain no
collide_modify     rotate smooth vibrate smooth

fix                fweight stochastic_weight
collide_modify     stochastic_weight yes
collide_modify     split 8 1.0
collide_modify     reduce {reduce_spec}

create_particles   air n 0

dump               d particle all 1 {dumpfile} id vx vy vz erot evib p_stochastic_wt
dump_modify        d format float %.17g

timestep           1.0e-9
run                2
"""

PART_B_DECK = """
seed               12137
dimension          3

species            n2.species N2
mixture            air N2 vstream 30.0 -20.0 50.0 temp 5000.0
mixture            air N2 frac 1.0

boundary           p p p
create_box         0.0 1.0e-5 0.0 1.0e-5 0.0 1.0e-5
create_grid        {grid}
balance_grid       rcb part

global             nrho 4.247e22
global             fnum {fnum}
global             comm/sort yes

collide            vss air n2.vss
collide_modify     remain no
collide_modify     rotate smooth vibrate smooth

fix                fweight stochastic_weight
collide_modify     stochastic_weight yes
collide_modify     split 8 1.0
collide_modify     reduce {reduce_spec}

create_particles   air n 0

dump               d particle all {ndump} {dumpfile} id vx vy vz erot evib p_stochastic_wt
dump_modify        d format float %.17g

timestep           1.0e-9
run                {steps}
"""

# ---------------------------------------------------------------------------
# dump parsing and moment evaluation (pure python: no numpy dependency)
# ---------------------------------------------------------------------------

def read_dump(path):
    """return list of (timestep, particles) where particles is a list of
    dicts with vx,vy,vz,erot,evib,w"""
    snaps = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i+1])
            assert lines[i+2].startswith("ITEM: NUMBER")
            n = int(lines[i+3])
            j = i + 4
            while not lines[j].startswith("ITEM: ATOMS"):
                j += 1
            cols = lines[j].split()[2:]
            idx = {c: k for k, c in enumerate(cols)}
            parts = []
            for k in range(n):
                t = lines[j+1+k].split()
                parts.append({
                    "vx": float(t[idx["vx"]]),
                    "vy": float(t[idx["vy"]]),
                    "vz": float(t[idx["vz"]]),
                    "erot": float(t[idx["erot"]]),
                    "evib": float(t[idx["evib"]]),
                    "w": float(t[idx["p_stochastic_wt"]]),
                })
            snaps.append((step, parts))
            i = j + 1 + n
        else:
            i += 1
    return snaps


def moments(parts, mass=MASS_N2):
    """weighted moments of a particle set, computed independently of SPARTA:
    W       total weight
    P[3]    momentum                  sum w m v
    Etr     translational energy      sum w 1/2 m |v|^2
    Erot    rotational energy         sum w erot
    Evib    vibrational energy        sum w evib
    Pij[6]  central stress            sum w m c_i c_j   (xx,yy,zz,xy,xz,yz)
    q[3]    central heat flux         sum w 1/2 m c_i |c|^2
    M4      central 4th moment        sum w m |c|^4
    where c = v - V and V is the weighted mean velocity."""
    W = sum(p["w"] for p in parts)
    P = [mass*sum(p["w"]*p["v"+d] for p in parts) for d in "xyz"]
    Etr = 0.5*mass*sum(p["w"]*(p["vx"]**2+p["vy"]**2+p["vz"]**2)
                       for p in parts)
    Erot = sum(p["w"]*p["erot"] for p in parts)
    Evib = sum(p["w"]*p["evib"] for p in parts)
    V = [P[d]/(mass*W) for d in range(3)]
    Pij = [0.0]*6
    q = [0.0, 0.0, 0.0]
    M4 = 0.0
    for p in parts:
        c = [p["vx"]-V[0], p["vy"]-V[1], p["vz"]-V[2]]
        c2 = c[0]*c[0]+c[1]*c[1]+c[2]*c[2]
        wm = p["w"]*mass
        Pij[0] += wm*c[0]*c[0]; Pij[1] += wm*c[1]*c[1]
        Pij[2] += wm*c[2]*c[2]; Pij[3] += wm*c[0]*c[1]
        Pij[4] += wm*c[0]*c[2]; Pij[5] += wm*c[1]*c[2]
        for d in range(3):
            q[d] += 0.5*wm*c[d]*c2
        M4 += wm*c2*c2
    return {"W": W, "P": P, "Etr": Etr, "Erot": Erot, "Evib": Evib,
            "Pij": Pij, "q": q, "M4": M4}


# ---------------------------------------------------------------------------

def run_sparta(cmd, deck_text, workdir, tag):
    infile = os.path.join(workdir, "in." + tag)
    with open(infile, "w") as f:
        f.write(deck_text)
    argv = [os.path.abspath(t) if os.path.exists(t) else t
            for t in shlex.split(cmd)]
    argv += ["-in", infile, "-log", os.path.join(workdir, "log." + tag)]
    proc = subprocess.run(argv, cwd=workdir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=1200)
    return proc.returncode, proc.stdout


def check(cond, msg, failures):
    print("  [%s] %s" % ("ok" if cond else "FAIL", msg))
    if not cond:
        failures.append(msg)


def rel(a, b, scale):
    return abs(a - b) / max(abs(scale), 1e-300)


def vrel(a, b, scale):
    return max(rel(x, y, scale) for x, y in zip(a, b))


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    failures = []

    workdir = tempfile.mkdtemp(prefix="swpm_audit_")
    with open(os.path.join(workdir, "n2.species"), "w") as f:
        f.write(SPECIES_FILE)
    with open(os.path.join(workdir, "n2.nocoll.vss"), "w") as f:
        f.write(VSS_NOCOLL)
    with open(os.path.join(workdir, "n2.vss"), "w") as f:
        f.write(VSS_REAL)

    NP = 200
    FNUM = 4.247e22 * 1.0e-5**3 / NP
    TIGHT = 1e-11      # "conserved exactly": FP roundoff of ~200-term sums
    LOOSE = 1e-3       # "changed": scheme does not conserve this moment

    # ------------------------------------------------------------------
    # Part A: reduction operator in isolation (no collisions, no splits;
    # whole cell merged as a single group in one event at step 1)
    # ------------------------------------------------------------------

    part_a = [
        # scheme, extra exact quantities, quantities that must change
        ("energy binary 50 4 1000000", "energy", [], ["q", "Pdev"]),
        ("heat binary 50 4 1000000", "heat", ["q"], ["Pdev"]),
        ("stress binary 50 8 1000000", "stress", ["q", "Pij"], ["M4"]),
    ]

    for spec, name, extra_exact, must_change in part_a:
        print("== Part A (reduction only): %s" % name)
        dumpfile = "a_%s.dump" % name
        deck = PART_A_DECK.format(fnum=FNUM, reduce_spec=spec,
                                  dumpfile=dumpfile)
        rc, out = run_sparta(cmd, deck, workdir, "a_" + name)
        if rc != 0:
            check(False, "run completed (rc=%d)" % rc, failures)
            print(out[-1500:])
            continue
        snaps = read_dump(os.path.join(workdir, dumpfile))
        pre = moments(snaps[0][1])
        post = moments(snaps[-1][1])
        npre, npost = len(snaps[0][1]), len(snaps[-1][1])
        check(npost < npre,
              "reduction happened (np %d -> %d)" % (npre, npost), failures)

        escale = pre["Etr"] + pre["Erot"] + pre["Evib"]
        pscale = MASS_N2 * pre["W"] * 400.0   # thermal momentum scale
        check(rel(post["W"], pre["W"], pre["W"]) < TIGHT,
              "total weight conserved", failures)
        check(vrel(post["P"], pre["P"], pscale) < TIGHT,
              "momentum vector conserved", failures)
        check(rel(post["Etr"], pre["Etr"], escale) < TIGHT,
              "translational energy conserved", failures)
        check(rel(post["Erot"], pre["Erot"], escale) < TIGHT,
              "rotational energy conserved", failures)
        check(rel(post["Evib"], pre["Evib"], escale) < TIGHT,
              "vibrational energy conserved", failures)
        # trace of stress == 2*Etr in the drift frame: conserved by all
        tr_pre = pre["Pij"][0]+pre["Pij"][1]+pre["Pij"][2]
        tr_post = post["Pij"][0]+post["Pij"][1]+post["Pij"][2]
        check(rel(tr_post, tr_pre, tr_pre) < TIGHT,
              "stress trace (thermal energy) conserved", failures)

        if "q" in extra_exact:
            qscale = max(abs(x) for x in pre["q"])
            check(vrel(post["q"], pre["q"], qscale) < 1e-9,
                  "heat-flux vector conserved (scheme guarantee)", failures)
        if "Pij" in extra_exact:
            check(vrel(post["Pij"], pre["Pij"], tr_pre) < 1e-9,
                  "full stress tensor conserved (scheme guarantee)", failures)

        for qty in must_change:
            if qty == "q":
                qscale = max(abs(x) for x in pre["q"])
                changed = vrel(post["q"], pre["q"], qscale) > LOOSE
                label = "heat flux NOT conserved (as expected for %s)" % name
            elif qty == "Pdev":
                dev_pre = [pre["Pij"][0]-tr_pre/3, pre["Pij"][3],
                           pre["Pij"][4], pre["Pij"][5]]
                dev_post = [post["Pij"][0]-tr_post/3, post["Pij"][3],
                            post["Pij"][4], post["Pij"][5]]
                changed = vrel(dev_post, dev_pre, tr_pre) > LOOSE/1e3
                label = ("stress deviator NOT conserved "
                         "(as expected for %s)" % name)
            elif qty == "M4":
                changed = rel(post["M4"], pre["M4"], pre["M4"]) > LOOSE
                label = ("4th moment NOT conserved "
                         "(as expected: beyond %s guarantee)" % name)
            check(changed, label + " - audit is sensitive", failures)

    # ------------------------------------------------------------------
    # Part B: full method (splits + collisions + reductions), long run.
    # only the true collision invariants are asserted; drift must be at
    # the floating-point accumulation level.
    # ------------------------------------------------------------------

    part_b = [
        ("full-heat-weight-multicell", "2 2 1", 4,
         "heat weight 40 6 16", 5000),
        ("full-energy-binary", "1 1 1", 1,
         "energy binary 100 4 16", 5000),
        ("full-stress-binary", "1 1 1", 1,
         "stress binary 100 8 24", 5000),
    ]

    for name, grid, ncells, spec, steps in part_b:
        print("== Part B (full method, %d steps): %s" % (steps, name))
        dumpfile = "b_%s.dump" % name
        deck = PART_B_DECK.format(fnum=FNUM/ncells, grid=grid,
                                  reduce_spec=spec, steps=steps,
                                  ndump=steps//10, dumpfile=dumpfile)
        rc, out = run_sparta(cmd, deck, workdir, "b_" + name)
        if rc != 0:
            check(False, "run completed (rc=%d)" % rc, failures)
            print(out[-1500:])
            continue
        snaps = read_dump(os.path.join(workdir, dumpfile))
        ref = moments(snaps[0][1])
        escale = ref["Etr"] + ref["Erot"] + ref["Evib"]
        pscale = MASS_N2 * ref["W"] * 400.0
        wd = ed = pd = 0.0
        nred = 0
        for step, parts in snaps[1:]:
            m = moments(parts)
            wd = max(wd, rel(m["W"], ref["W"], ref["W"]))
            ed = max(ed, rel(m["Etr"]+m["Erot"]+m["Evib"], escale, escale))
            pd = max(pd, vrel(m["P"], ref["P"], pscale))
        check(wd < 1e-9, "weight drift %.2e < 1e-9 over run" % wd, failures)
        check(pd < 1e-9, "momentum drift %.2e < 1e-9 over run" % pd, failures)
        check(ed < 1e-9, "total energy drift %.2e < 1e-9 over run" % ed,
              failures)
        counts = [len(p) for _, p in snaps]
        check(min(counts[1:]) < counts[0] and max(counts) > 0,
              "splits/reductions active (np range %d..%d)" %
              (min(counts), max(counts)), failures)

    print()
    if failures:
        print("%d AUDIT CHECK(S) FAILED:" % len(failures))
        for f in failures:
            print("  - " + f)
        sys.exit(1)
    print("All SWPM conservation audit checks passed.")
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
