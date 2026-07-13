#!/usr/bin/env python3
"""
Compare serial vs Kokkos SWPM output for a few representative cases.

Runs the same 1D Couette/Fourier SWPM cases with:
  - spa_serial  (non-Kokkos path)
  - spa_kokkos  (Kokkos path: binary + '-k on -sf kk')

then checks that the averaged profiles match within tolerance.

With SPARTA_KOKKOS_EXACT the Kokkos serial backend uses the same RNG
trajectory as the non-Kokkos code, so profiles should be bit-identical
(same seed, same cell order, same collisions).  We use a generous
tolerance (1%) to allow for any floating-point ordering differences.
"""

import os, math, subprocess, sys
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
SPA_S  = os.path.normpath(os.path.join(HERE, "..", "..", "..", "src", "spa_serial"))
SPA_K  = os.path.normpath(os.path.join(
    HERE, "..", "..", "..", "build_kokkos_serial", "src", "spa_serial"))
SPECIES= os.path.normpath(os.path.join(HERE, "..", "air.species"))
VSS    = os.path.normpath(os.path.join(HERE, "..", "air.vss"))
OUT    = os.path.join(HERE, "out_compare")
os.makedirs(OUT, exist_ok=True)

KB    = 1.380658e-23
MASS  = 6.63e-26
DREF  = 4.11e-10
NRHO  = 1.0e22
NY    = 40        # coarser for speed
PPC   = 20
T0    = 273.0
UW    = 50.0
TCOLD = 200.0
THOT  = 400.0
NSTEADY = 2000
NEVERY, NREPEAT = 5, 4000
NSAMPLE = NEVERY * NREPEAT
TOL   = 0.02     # 2% relative tolerance on mean absolute profile value

# Small set of representative cases: (flow, reduce_mode, kn)
CASES = [
    ("couette", "dsmc",   0, 0.1),
    ("couette", "energy", 1, 0.1),
    ("couette", "heat",   2, 0.1),
    ("couette", "stress", 3, 0.1),
    ("fourier", "energy", 1, 0.1),
    ("fourier", "stress", 3, 0.1),
]

def mfp():
    return 1.0 / (math.sqrt(2.0) * math.pi * DREF * DREF * NRHO)

def make_input(case, modename, modeid, kn, tag):
    ly = mfp() / kn
    dy = ly / NY
    lx = lz = dy
    fnum = NRHO * (lx * lz * dy) / PPC
    tref = T0 if case == "couette" else THOT
    vth  = math.sqrt(2.0 * KB * tref / MASS)
    dt   = 0.05 * dy / vth
    nsteady = NSTEADY if case == "couette" else int(NSTEADY * 1.5)
    outfile = os.path.join(OUT, f"{case}_{modename}_kn{kn:g}_{tag}.grid")

    if modeid == 0:
        swpm = ""
    else:
        choice = {1: "energy", 2: "heat", 3: "stress"}[modeid]
        swpm = (f"fix             fw stochastic_weight\n"
                f"collide_modify  stochastic_weight yes\n"
                f"collide_modify  split {PPC} 1.0\n"
                f"collide_modify  reduce {choice} binary {2*PPC} 8 16\n")

    if case == "couette":
        walls = (f"surf_collide    bot diffuse {T0} 1.0 translate {-UW} 0 0\n"
                 f"surf_collide    top diffuse {T0} 1.0 translate {UW} 0 0\n"
                 f"bound_modify    ylo collide bot\n"
                 f"bound_modify    yhi collide top\n")
        seed  = (f'variable        yv internal 0.0\n'
                 f'variable        vxp equal "{UW}*(2.0*v_yv/{ly}-1.0)"\n'
                 f"create_particles air n 0 vstream vxp NULL NULL NULL yv NULL\n")
        comps = ("compute         c1 grid all all u\n"
                 "compute         c2 pflux/grid all all momxy\n")
        mixtemp = T0
    else:
        walls = (f"surf_collide    bot diffuse {TCOLD} 1.0\n"
                 f"surf_collide    top diffuse {THOT} 1.0\n"
                 f"bound_modify    ylo collide bot\n"
                 f"bound_modify    yhi collide top\n")
        w1 = 1.81
        seed  = (f'variable        yv internal 0.0\n'
                 f'variable        Tp equal "({TCOLD}^{w1}+({THOT}^{w1}-{TCOLD}^{w1})*v_yv/{ly})^(1.0/{w1})"\n'
                 f"create_particles air n 0 temperature Tp NULL yv NULL\n")
        comps = ("compute         c1 thermal/grid all all temp\n"
                 "compute         c2 eflux/grid all all heaty\n")
        mixtemp = 0.5 * (TCOLD + THOT)

    txt = f"""\
seed            12345
dimension       3
global          gridcut 0.0 comm/sort yes

species         {SPECIES} Ar
mixture         air Ar vstream 0 0 0 temp {mixtemp}
mixture         air Ar frac 1.0

boundary        p s p
create_box      0 {lx} 0 {ly} 0 {lz}
create_grid     1 {NY} 1
balance_grid    rcb part

global          nrho {NRHO:.6e} fnum {fnum:.6e}

{walls}
collide         vss air {VSS}
collide_modify  remain no
{swpm}
{seed}
{comps}
timestep        {dt:.6e}
run             {nsteady}

reset_timestep  0
fix             favg ave/grid all {NEVERY} {NREPEAT} {NSAMPLE} c_c1[1] c_c2[1]
dump            dg grid all {NSAMPLE} {outfile} yc f_favg[1] f_favg[2]
run             {NSAMPLE}
"""
    return txt, outfile

def parse_grid(path):
    with open(path) as f:
        lines = f.readlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: CELLS")]
    i = starts[-1]
    n = None
    for j in range(i, -1, -1):
        if lines[j].startswith("ITEM: NUMBER OF CELLS"):
            n = int(lines[j + 1]); break
    rows = [list(map(float, lines[i + 1 + k].split())) for k in range(n)]
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 1], a[:, 2]

def run_case(binary, args, case, modename, modeid, kn, tag):
    txt, outfile = make_input(case, modename, modeid, kn, tag)
    if os.path.exists(outfile):
        os.remove(outfile)
    infile = os.path.join(OUT, f"in.{case}_{modename}_kn{kn:g}_{tag}")
    with open(infile, "w") as f:
        f.write(txt)
    cmd = [binary] + args + ["-in", infile, "-log", "none"]
    p = subprocess.run(cmd, cwd=OUT, capture_output=True, text=True, timeout=300)
    if p.returncode != 0 or not os.path.exists(outfile):
        err = (p.stdout + p.stderr).strip().splitlines()
        return None, (err[-1] if err else "no output")
    return parse_grid(outfile), "ok"

def relerr(a, b):
    scale = max(np.mean(np.abs(a)), np.mean(np.abs(b)), 1e-30)
    return np.max(np.abs(a - b)) / scale

def main():
    if not os.path.exists(SPA_K):
        print(f"ERROR: Kokkos binary not found at {SPA_K}")
        print("Build it first:  cd build_kokkos_serial && cmake ... && make")
        sys.exit(1)

    print(f"serial:  {SPA_S}")
    print(f"kokkos:  {SPA_K} -k on -sf kk")
    print(f"running {len(CASES)} cases x 2 builds = {2*len(CASES)} simulations\n")

    all_pass = True
    for (case, mn, mid, kn) in CASES:
        (v1s, v2s), st_s = run_case(SPA_S, [], case, mn, mid, kn, "serial")
        (v1k, v2k), st_k = run_case(SPA_K, ["-k", "on", "-sf", "kk"], case, mn, mid, kn, "kokkos")

        if v1s is None:
            print(f"FAIL  {case:8s} {mn:7s} Kn={kn}  serial:  {st_s}")
            all_pass = False; continue
        if v1k is None:
            print(f"FAIL  {case:8s} {mn:7s} Kn={kn}  kokkos:  {st_k}")
            all_pass = False; continue

        e1 = relerr(v1s, v1k)
        e2 = relerr(v2s, v2k)
        ok = e1 <= TOL and e2 <= TOL
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"{status}  {case:8s} {mn:7s} Kn={kn}  rel_err_v1={e1:.4f}  rel_err_v2={e2:.4f}  (tol={TOL})")

    print()
    if all_pass:
        print("ALL PASS — Kokkos SWPM matches serial within tolerance.")
    else:
        print("SOME CASES FAILED.")
    sys.exit(0 if all_pass else 1)

if __name__ == "__main__":
    main()
