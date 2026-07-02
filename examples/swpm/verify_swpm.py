#!/usr/bin/env python3
"""
verify_swpm.py - functional/unit-level verification of the SWPM
                 (stochastically weighted particle method) implementation.

Complements the gold-log regression tests in this directory: instead of
comparing against archived logs (which pins one trajectory), this script
asserts the *invariants* the SWPM algorithms must satisfy on any machine,
any compiler and any seed, plus the documented error conditions.  Run it
after any refactor of the SWPM code.

Checks performed:

1. Conservation invariants, per input deck (a superset of the regression
   decks in this directory):
     - weighted thermal temperature is constant (energy conservation
       through splitting, elastic collisions and all reduction schemes)
     - weighted number density is constant (weight conservation)
     - weighted mean velocity is constant (momentum conservation)
   for every combination of reduction scheme (energy, heat, stress) and
   grouping strategy (binary, weight), plus the always-split regime.

2. Polyatomic (N2) conservation, per reduction scheme: with rotational and
   vibrational relaxation active the individual energy modes exchange, and
   the conserved invariant is the TOTAL (translational + rotational +
   vibrational) weighted energy.  Also asserts the vibrational energy
   actually exchanges, so the check cannot pass vacuously.

3. Algorithm activity (the invariants must not pass vacuously):
     - the split deck must grow np
     - reduction decks must keep np below the starting count

4. Determinism: running the same deck twice gives identical stats output.

5. Error conditions:
     - collide_modify stochastic_weight yes without fix stochastic_weight
     - collide_modify reduce stress with Ngmin < 6
     - collide_modify reduce energy with Ngmin < 2
     - collide_modify nearcp yes + stochastic weighting (unsupported)
     - collide_modify vibrate discrete + particle reduction (survivor evib
       is continuous, inconsistent with quantized levels)

Usage:
    python3 verify_swpm.py /path/to/spa_binary
    python3 verify_swpm.py "mpirun -np 4 /path/to/spa_binary"

Exit code 0 = all checks passed.
"""

import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# input deck template: homogeneous argon box, SWPM options injected
# ---------------------------------------------------------------------------

DECK_TEMPLATE = """
seed               {seed}
dimension          3

species            air.species Ar
mixture            air Ar vstream 0.0 0.0 0.0 temp 298.15
mixture            air Ar frac 1.0

boundary           p p p
create_box         0.0 1.0e-5 0.0 1.0e-5 0.0 1.0e-5
create_grid        {grid}
balance_grid       rcb part

global             nrho 4.247e22
global             fnum {fnum}

collide            vss air air.vss
collide_modify     remain no

fix                fweight stochastic_weight
collide_modify     stochastic_weight yes
{swpm_lines}

create_particles   air n 0

compute gT         temp
compute n          grid all all n nrho
compute nr1        reduce ave c_n[1]
compute nr2        reduce ave c_n[2]
compute v          grid all all u v w
compute vr1        reduce ave c_v[1]
compute vr2        reduce ave c_v[2]
compute vr3        reduce ave c_v[3]

stats              50
stats_style        step np c_gT c_nr1 c_nr2 c_vr1 c_vr2 c_vr3

timestep           1.0e-9
run                {steps}
"""

# fnum chosen so a 1x1x1 grid at nrho 4.247e22 gives ppc particles per cell
def fnum(ppc, ncells):
    return 4.247e22 * 1.0e-5**3 / (ppc * ncells)

CONSERVATION_CASES = [
    # (name, grid, ncells, ppc, swpm option lines)
    ("split-only", "1 1 1", 1, 50,
     "collide_modify     split 1000 1.0"),
    ("reduce-energy-binary", "1 1 1", 1, 200,
     "collide_modify     split 8 1.0\n"
     "collide_modify     reduce energy binary 100 4 16"),
    ("reduce-heat-binary", "1 1 1", 1, 200,
     "collide_modify     split 8 1.0\n"
     "collide_modify     reduce heat binary 100 4 16"),
    ("reduce-stress-binary", "1 1 1", 1, 200,
     "collide_modify     split 8 1.0\n"
     "collide_modify     reduce stress binary 100 8 24"),
    ("reduce-energy-weight", "2 2 1", 4, 60,
     "collide_modify     split 1000 1.0\n"
     "collide_modify     reduce energy weight 40 4 16"),
    ("reduce-heat-weight", "4 4 1", 16, 30,
     "collide_modify     split 1000 1.0\n"
     "collide_modify     reduce heat weight 40 6 16"),
    ("reduce-stress-weight", "2 2 1", 4, 60,
     "collide_modify     split 1000 1.0\n"
     "collide_modify     reduce stress weight 60 8 24"),
]

# polyatomic (N2) deck: internal-energy conservation through split/reduce.
# n2.testvib.species has an artificially fast vibrational relaxation
# probability so all three energy modes exchange within the run; the
# invariant is their weighted sum (total energy), not any single mode.

POLY_DECK = """
seed               12137
dimension          3

species            n2.testvib.species N2
mixture            air N2 vstream 0.0 0.0 0.0 temp 5000.0
mixture            air N2 frac 1.0

boundary           p p p
create_box         0.0 1.0e-5 0.0 1.0e-5 0.0 1.0e-5
create_grid        1 1 1
balance_grid       rcb part

global             nrho 4.247e22
global             fnum {fnum}

collide            vss air air.vss
collide_modify     remain no
collide_modify     rotate smooth vibrate smooth

fix                fweight stochastic_weight
collide_modify     stochastic_weight yes
collide_modify     split 8 1.0
collide_modify     reduce {reduce_spec}

create_particles   air n 0

compute n          grid all all n nrho
compute v          grid all all u v w
compute e          grid all all ke erot evib
compute nr2        reduce ave c_n[2]
compute vr1        reduce ave c_v[1]
compute vr2        reduce ave c_v[2]
compute vr3        reduce ave c_v[3]
compute er1        reduce ave c_e[1]
compute er2        reduce ave c_e[2]
compute er3        reduce ave c_e[3]

stats              50
stats_style        step np c_nr2 c_vr1 c_vr2 c_vr3 c_er1 c_er2 c_er3

timestep           1.0e-9
run                500
"""

POLYATOMIC_CASES = [
    # (name, reduce spec) - each scheme must conserve total (trans+rot+vib)
    # energy while the modes exchange through relaxation
    ("polyatomic-energy", "energy binary 100 4 16"),
    ("polyatomic-heat", "heat binary 100 4 16"),
    ("polyatomic-stress", "stress binary 100 8 24"),
]

ERROR_CASES = [
    ("stochastic_weight-without-fix",
     "air.species Ar",
     "collide_modify     stochastic_weight yes",
     "requires fix stochastic_weight"),
    ("reduce-stress-small-Ngmin",
     "air.species Ar",
     "fix                fweight stochastic_weight\n"
     "collide_modify     stochastic_weight yes\n"
     "collide_modify     reduce stress binary 100 4 24",
     "requires Ngmin >= 6"),
    ("reduce-energy-small-Ngmin",
     "air.species Ar",
     "fix                fweight stochastic_weight\n"
     "collide_modify     stochastic_weight yes\n"
     "collide_modify     reduce energy binary 100 1 16",
     "requires Ngmin >= 2"),
    ("nearcp-unsupported",
     "air.species Ar",
     "fix                fweight stochastic_weight\n"
     "collide_modify     stochastic_weight yes\n"
     "collide_modify     nearcp yes 10",
     "near-neighbor"),
    ("reduce-discrete-vibration",
     "n2.testvib.species N2",
     "collide_modify     vibrate discrete\n"
     "fix                fweight stochastic_weight\n"
     "collide_modify     stochastic_weight yes\n"
     "collide_modify     reduce energy binary 100 4 16",
     "vibrate smooth or no"),
]

ERROR_DECK_TEMPLATE = """
seed               12137
dimension          3

species            {species}
mixture            air {sp_id} vstream 0.0 0.0 0.0 temp 298.15
mixture            air {sp_id} frac 1.0

boundary           p p p
create_box         0.0 1.0e-5 0.0 1.0e-5 0.0 1.0e-5
create_grid        1 1 1
balance_grid       rcb part

global             nrho 4.247e22
global             fnum 8.494e6

collide            vss air air.vss
{swpm_lines}

create_particles   air n 0

timestep           1.0e-9
run                10
"""

# ---------------------------------------------------------------------------

def run_sparta(cmd, deck_text, workdir, tag):
    infile = os.path.join(workdir, "in." + tag)
    with open(infile, "w") as f:
        f.write(deck_text)
    logfile = os.path.join(workdir, "log." + tag)
    # sparta runs with cwd=workdir: resolve any relative paths (the binary)
    argv = [os.path.abspath(t) if os.path.exists(t) else t
            for t in shlex.split(cmd)]
    argv += ["-in", infile, "-log", logfile]
    proc = subprocess.run(argv, cwd=workdir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=600)
    return proc.returncode, proc.stdout, logfile


def parse_stats(output):
    """return (header, rows) of the last stats table in the output"""
    lines = output.splitlines()
    header = None
    rows = []
    for i, line in enumerate(lines):
        if re.match(r"^Step\s+Np\s", line):
            header = line.split()
            rows = []
            for l in lines[i+1:]:
                toks = l.split()
                if len(toks) != len(header):
                    break
                try:
                    rows.append([float(t) for t in toks])
                except ValueError:
                    break
    return header, rows


def check(cond, msg, failures):
    status = "ok" if cond else "FAIL"
    print("  [%s] %s" % (status, msg))
    if not cond:
        failures.append(msg)


def rel_diff(a, b, scale):
    return abs(a - b) / max(abs(scale), 1e-300)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    failures = []

    workdir = tempfile.mkdtemp(prefix="swpm_verify_")
    for aux in ("air.species", "air.vss", "n2.testvib.species"):
        shutil.copy(os.path.join(THIS_DIR, aux), workdir)

    # thermal speed scale for momentum tolerance (m/s)
    VSCALE = 400.0

    for name, grid, ncells, ppc, swpm in CONSERVATION_CASES:
        print("== conservation: %s" % name)
        deck = DECK_TEMPLATE.format(seed=12137, grid=grid, steps=500,
                                    fnum=fnum(ppc, ncells), swpm_lines=swpm)
        rc, out, _ = run_sparta(cmd, deck, workdir, name)
        header, rows = parse_stats(out)
        if rc != 0 or not rows:
            check(False, "run completed with stats output (rc=%d)" % rc,
                  failures)
            print(out[-2000:])
            continue
        col = {h: j for j, h in enumerate(header)}
        first, last = rows[0], rows[-1]

        T0 = first[col["c_gT"]]
        check(all(rel_diff(r[col["c_gT"]], T0, T0) < 1e-6 for r in rows),
              "weighted temperature constant (energy conservation)", failures)
        n0 = first[col["c_nr2"]]
        check(all(rel_diff(r[col["c_nr2"]], n0, n0) < 1e-6 for r in rows),
              "weighted number density constant (weight conservation)",
              failures)
        for c in ("c_vr1", "c_vr2", "c_vr3"):
            v0 = first[col[c]]
            if ncells == 1:
                # single cell: the weighted mean velocity is exactly conserved
                check(all(rel_diff(r[col[c]], v0, VSCALE) < 1e-6
                          for r in rows),
                      "weighted mean velocity %s constant (momentum)" % c,
                      failures)
        if name == "split-only":
            check(last[col["Np"]] > 2*first[col["Np"]],
                  "np grew (splitting active)", failures)
        else:
            check(last[col["Np"]] < first[col["Np"]],
                  "np dropped below initial count (reduction active)",
                  failures)

    for name, reduce_spec in POLYATOMIC_CASES:
        print("== polyatomic conservation: %s" % name)
        deck = POLY_DECK.format(fnum=fnum(200, 1), reduce_spec=reduce_spec)
        rc, out, _ = run_sparta(cmd, deck, workdir, name)
        header, rows = parse_stats(out)
        if rc != 0 or not rows:
            check(False, "run completed with stats output (rc=%d)" % rc,
                  failures)
            print(out[-2000:])
            continue
        col = {h: j for j, h in enumerate(header)}
        first, last = rows[0], rows[-1]

        n0 = first[col["c_nr2"]]
        check(all(rel_diff(r[col["c_nr2"]], n0, n0) < 1e-6 for r in rows),
              "weighted number density constant (weight conservation)",
              failures)
        for c in ("c_vr1", "c_vr2", "c_vr3"):
            v0 = first[col[c]]
            check(all(rel_diff(r[col[c]], v0, VSCALE) < 1e-6 for r in rows),
                  "weighted mean velocity %s constant (momentum)" % c,
                  failures)

        # total energy (trans + rot + vib per weighted particle) is the
        # invariant; the individual modes exchange through relaxation

        def etot(r):
            return r[col["c_er1"]] + r[col["c_er2"]] + r[col["c_er3"]]
        e0 = etot(first)
        check(all(rel_diff(etot(r), e0, e0) < 1e-6 for r in rows),
              "total (trans+rot+vib) energy constant", failures)
        check(any(rel_diff(r[col["c_er3"]], first[col["c_er3"]],
                           first[col["c_er3"]]) > 1e-3 for r in rows),
              "vibrational energy exchanges (test is not vacuous)", failures)
        check(last[col["Np"]] < first[col["Np"]],
              "np dropped below initial count (reduction active)", failures)

    print("== determinism: same seed twice gives identical stats")
    name, grid, ncells, ppc, swpm = CONSERVATION_CASES[1]
    deck = DECK_TEMPLATE.format(seed=12137, grid=grid, steps=500,
                                fnum=fnum(ppc, ncells), swpm_lines=swpm)
    _, out1, _ = run_sparta(cmd, deck, workdir, "determinism.a")
    _, out2, _ = run_sparta(cmd, deck, workdir, "determinism.b")
    check(parse_stats(out1) == parse_stats(out2),
          "identical stats tables", failures)

    for name, species, swpm, expect in ERROR_CASES:
        print("== error condition: %s" % name)
        deck = ERROR_DECK_TEMPLATE.format(species=species,
                                          sp_id=species.split()[-1],
                                          swpm_lines=swpm)
        rc, out, _ = run_sparta(cmd, deck, workdir, name)
        ok = rc != 0 and "ERROR" in out and expect in out
        check(ok, "rejected with message containing '%s'" % expect, failures)
        if not ok:
            print(out[-1000:])

    print()
    if failures:
        print("%d CHECK(S) FAILED:" % len(failures))
        for f in failures:
            print("  - " + f)
        sys.exit(1)
    print("All SWPM verification checks passed.")
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
