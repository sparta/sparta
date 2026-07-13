#!/usr/bin/env python3
"""
SWPM validation suite: 1D Couette and 1D Fourier flows.

For each flow and a sweep of Knudsen numbers (continuum -> free molecular) it
runs pure DSMC and the stochastic weighted particle method (SWPM) with each of
the three particle-reduction schemes (energy, heat, stress), then plots:

  Couette: x-velocity profile u(y) and shear stress Pxy(y)
  Fourier: temperature profile T(y) and heat flux qy(y)

The flow is 1D: the domain is discretized only in y (the gradient direction);
x and z are periodic.  The bottom/top y boundaries are no-slip diffuse walls.
Each run seeds the expected steady-state profile in create_particles (linear
velocity for Couette, linear temperature for Fourier) so it reaches steady
state quickly.  Runs are small and executed up to MAXPROC at a time.

Usage:  python3 run_validation.py
Outputs: couette.png, fourier.png (and per-run grid dumps in out/).
"""

import os, math, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- paths
HERE   = os.path.dirname(os.path.abspath(__file__))
SPA    = os.path.normpath(os.path.join(HERE, "..", "..", "..", "src", "spa_serial"))
SPECIES= os.path.normpath(os.path.join(HERE, "..", "air.species"))
VSS    = os.path.normpath(os.path.join(HERE, "..", "air.vss"))
OUT    = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- physics / numerics
KB    = 1.380658e-23
MASS  = 6.63e-26       # argon
DREF  = 4.11e-10       # argon VSS reference diameter

NRHO  = 1.0e22         # fixed number density -> fixed mean free path
NY    = 200            # cells across the gap (>= 40, resolves the mfp)
PPC   = 40             # target particles per cell

T0    = 273.0          # Couette wall temperature
UW    = 50.0           # Couette wall speed (+/- in x)
TCOLD = 200.0          # Fourier cold wall
THOT  = 400.0          # Fourier hot wall

# Knudsen number is set by varying the channel height: Kn = mfp / Ly, with the
# density (mfp) held fixed.  cell size dy = Ly/NY stays <= mfp for Kn >= 0.005.
KN_LIST = [0.01, 0.1, 1.0]            # continuum -> slip -> transition
MODES   = [("dsmc", 0), ("energy", 1), ("heat", 2), ("stress", 3)]

NSTEADY = 8000                         # warm-up steps to equilibrate before measuring
NEVERY, NREPEAT = 5, 10000
NSAMPLE = NEVERY * NREPEAT             # long averaging window for the flux signal
MAXPROC = 8
RUN_TIMEOUT = 600                      # per-run wall-clock cap (s); guards against runaway

def mfp():
    return 1.0 / (math.sqrt(2.0) * math.pi * DREF * DREF * NRHO)

# ---------------------------------------------------------------- helpers
def make_input(case, modename, modeid, kn):
    ly = mfp() / kn                 # channel height set by Kn (fixed density)
    dy = ly / NY
    lx = lz = dy                    # cubic cells (x,z periodic, one cell wide)
    nrho = NRHO
    fnum = nrho * (lx * lz * dy) / PPC
    tref = T0 if case == "couette" else THOT
    vth  = math.sqrt(2.0 * KB * tref / MASS)
    # small CFL: SWPM needs a small timestep so per-step splitting stays modest
    # relative to the once-per-step reduce (large dt -> weight-concentration runaway)
    dt   = 0.05 * dy / vth
    # Couette's linear velocity seed is the exact continuum steady state, so it
    # equilibrates immediately.  Fourier conducts heat on a diffusion timescale
    # ~L^2/alpha (steps ~ 1/Kn), so it needs a much longer warm-up.
    nsteady = NSTEADY if case == "couette" else min(120000, int(NSTEADY * 1.5 / kn))
    outfile = os.path.join(OUT, f"{case}_{modename}_kn{kn:g}.grid")

    # SWPM collision block (split keeps >=PPC, reduce caps at 2*PPC)
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
    else:  # fourier
        walls = (f"surf_collide    bot diffuse {TCOLD} 1.0\n"
                 f"surf_collide    top diffuse {THOT} 1.0\n"
                 f"bound_modify    ylo collide bot\n"
                 f"bound_modify    yhi collide top\n")
        # steady Fourier profile for conductivity k ~ T^omega (Ar VSS omega=0.81):
        # heat flux const => T^(omega+1) is linear in y.  seeding this curved
        # profile (not a straight line) starts the run near equilibrium.
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

global          nrho {nrho:.6e} fnum {fnum:.6e}

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

def run_one(args):
    case, modename, modeid, kn = args
    txt, outfile = make_input(case, modename, modeid, kn)
    if os.path.exists(outfile):
        os.remove(outfile)
    infile = os.path.join(OUT, f"in.{case}_{modename}_kn{kn:g}")
    with open(infile, "w") as f:
        f.write(txt)
    try:
        p = subprocess.run([SPA, "-in", infile, "-log", "none"],
                           cwd=OUT, capture_output=True, text=True, timeout=RUN_TIMEOUT)
    except subprocess.TimeoutExpired:
        return (case, modename, kn, None, "FAIL: timeout (likely SWPM attempt runaway)")
    if p.returncode != 0 or not os.path.exists(outfile):
        err = (p.stdout + p.stderr).strip().splitlines()
        msg = err[-1] if err else "no output"
        return (case, modename, kn, None, f"FAIL: {msg}")
    y, v1, v2 = parse_grid(outfile)
    return (case, modename, kn, (y, v1, v2), "ok")

def parse_grid(path):
    """read the LAST snapshot of a dump-grid file: columns yc val1 val2."""
    with open(path) as f:
        lines = f.readlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: CELLS")]
    i = starts[-1]
    # number of cells from the preceding NUMBER OF CELLS item
    n = None
    for j in range(i, -1, -1):
        if lines[j].startswith("ITEM: NUMBER OF CELLS"):
            n = int(lines[j + 1]); break
    rows = [list(map(float, lines[i + 1 + k].split())) for k in range(n)]
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]          # sort by yc
    return a[:, 0], a[:, 1], a[:, 2]

# ---------------------------------------------------------------- plot
COLORS = {"dsmc": "k", "energy": "tab:blue", "heat": "tab:green", "stress": "tab:red"}
STYLE  = {"dsmc": "o-", "energy": "s--", "heat": "^--", "stress": "v--"}

def plot_case(results, case, ylabels, fname):
    fig, axes = plt.subplots(2, len(KN_LIST), figsize=(4.2 * len(KN_LIST), 7),
                             squeeze=False)
    for col, kn in enumerate(KN_LIST):
        for row in range(2):
            ax = axes[row][col]
            for mn, _ in MODES:
                d = results.get((case, mn, kn))
                if d is None:
                    continue
                y, v1, v2 = d
                yn = (y - y.min()) / (y.max() - y.min())   # normalize to [0,1]
                val = v1 if row == 0 else v2
                ax.plot(yn, val, STYLE[mn], color=COLORS[mn], ms=3,
                        lw=1.2, mfc="none", label=mn)
            if row == 0:
                ax.set_title(f"Kn = {kn:g}")
            if col == 0:
                ax.set_ylabel(ylabels[row])
            if row == 1:
                ax.set_xlabel("y / L")
            ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=8, loc="best")
    fig.suptitle(f"1D {case.capitalize()} flow: DSMC vs SWPM (energy/heat/stress)")
    fig.tight_layout()
    out = os.path.join(HERE, fname)
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")

def main():
    tasks = [(case, mn, mid, kn)
             for case in ("couette", "fourier")
             for kn in KN_LIST
             for (mn, mid) in MODES]

    print(f"running {len(tasks)} simulations ({MAXPROC} at a time)...")
    results = {}
    with ThreadPoolExecutor(max_workers=MAXPROC) as ex:
        for case, mn, kn, data, status in ex.map(run_one, tasks):
            print(f"  {case:8s} {mn:7s} Kn={kn:<5g} -> {status}")
            results[(case, mn, kn)] = data

    plot_case(results, "couette", ["u  [m/s]", "Pxy  (momxy)"], "couette.png")
    plot_case(results, "fourier", ["T  [K]", "qy  (heaty)"], "fourier.png")
    print("done.")

if __name__ == "__main__":
    main()
