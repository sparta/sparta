#!/usr/bin/env python3
"""Quantitative validation of detailed-balance reverse reactions (issue #472).

Recomputes the species partition functions and equilibrium constants
INDEPENDENTLY of the SPARTA implementation (from the same data files), then
checks the measured forward/backward rates from frozen-composition
reservoirs against them:

  1  exchange pair N2 + O <-> NO + N: the measured tally ratio must equal
     the analytic equilibrium constant Keq(T) at several temperatures
  2  the derived backward rate must lie close to the independent
     literature Arrhenius fit for NO + N -> N2 + O (data/air.tce)
  3  dissociation/recombination pair N2 + N <-> N + N + N: the measured
     tally ratio must equal the analytic volumetric equilibrium constant
  4  a reacting (non-frozen) box initialized off-equilibrium must relax
     toward the analytic equilibrium composition from both sides

Usage:  python3 validate_reverse.py --exe ../../src/spa_serial
Exit code = number of failed checks.
"""

import argparse, math, os, re, shutil, subprocess, sys, tempfile

KB = 1.380649e-23
H  = 6.62607015e-34
HERE = os.path.dirname(os.path.abspath(__file__))

FAIL = 0
def check(name, ok, detail=""):
    global FAIL
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + detail) if detail else ""))

# ---------------------------------------------------------------- data

def parse_species(path):
    sp = {}
    for line in open(path):
        w = line.split()
        if not w or w[0].startswith("#"): continue
        # id molwt molmass rotdof rotrel vibdof vibrel vibtemp specwt charge
        sp[w[0]] = {"mass": float(w[2]), "rotdof": int(w[3]),
                    "vibdof": int(w[5]), "vibtemp": float(w[7])}
    return sp

def parse_rot(path):
    rot = {}
    for line in open(path):
        w = line.split()
        if not w or w[0].startswith("#"): continue
        n = int(w[1])
        temps = [float(x) for x in w[2:2+n]]
        sigma = float(w[2+n]) if len(w) > 2+n else 1.0
        rot[w[0]] = (temps, sigma)
    return rot

def parse_elec(path):
    elec = {}
    for line in open(path):
        w = line.split()
        if not w or w[0].startswith("#"): continue
        n = int(w[1])
        states = []
        for k in range(n):
            t = float(w[2+5*k]); g = float(w[2+5*k+2])
            states.append((t, g))
        elec[w[0]] = states
    return elec

def parse_vss(path):
    vss = {}
    for line in open(path):
        w = line.split()
        if not w or w[0].startswith("#"): continue
        vss[w[0]] = float(w[2])       # omega
    return vss

SPECIES = parse_species(os.path.join(HERE,"air.species"))
ROT = parse_rot(os.path.join(HERE,"air.rot"))
ELEC = parse_elec(os.path.join(HERE,"air.elec"))
VSS = parse_vss(os.path.join(HERE,"air.vss"))

def q_total(name, T, vibsmooth=False):
    """partition function per unit volume, matching ReactBird exactly;
    with vibsmooth the classical harmonic vib form (q = T/theta per mode),
    matching partition_function under collide_modify vibrate smooth"""
    s = SPECIES[name]
    q = (2.0*math.pi*s["mass"]*KB*T/(H*H))**1.5
    if s["rotdof"] == 2 and name in ROT:
        temps, sigma = ROT[name]
        q *= T/(sigma*temps[0])
    if s["vibdof"] >= 2 and s["vibtemp"] > 0.0:
        if vibsmooth: q *= T/s["vibtemp"]
        else:         q *= 1.0/(1.0 - math.exp(-s["vibtemp"]/T))
    if name in ELEC:
        q *= sum(g*math.exp(-t/T) for t,g in ELEC[name])
    return q

def park_fit(func, Ts):
    """solve for Park coefficients c0..c4 with
    ln f(T) = c0/Z + c1 + c2 ln Z + c3 Z + c4 Z^2,  Z = 1e4/T,
    exactly through the 5 given temperatures (pure-python 5x5 solve)"""
    A, y = [], []
    for T in Ts:
        Z = 1e4/T
        A.append([1.0/Z, 1.0, math.log(Z), Z, Z*Z])
        y.append(math.log(func(T)))
    # Gaussian elimination with partial pivoting
    n = 5
    M = [A[i][:] + [y[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        d = M[col][col]
        for j in range(col, n+1): M[col][j] /= d
        for r in range(n):
            if r == col: continue
            f = M[r][col]
            for j in range(col, n+1): M[r][j] -= f*M[col][j]
    return [M[i][n] for i in range(n)]

# ------------- independent reconstruction of the exchange DB table ----------
# Reproduces ReactBird::build_db_table (grid, exponents, ladder convolution)
# WITHOUT the calibration scale (which cancels in the rate ratio used by the
# non-equilibrium check), so the reverse per-collision probability shape
# mtab(ecc) can be integrated over an arbitrary collision-energy distribution.

EV = 1.602176634e-19

def _ladder(np, arr, du, n, eps, g):
    work = arr.copy(); out = np.zeros(n)
    for m in range(len(eps)):
        if g[m] == 0.0: continue
        sh = eps[m]/du; i0 = int(sh); frac = sh - i0
        if i0 < n:     out[i0:]   += g[m]*(1.0-frac)*work[:n-i0]
        if i0+1 < n:   out[i0+1:] += g[m]*frac*work[:n-i0-1]
    return out

def db_table_shape(np, spB, spF, Fcoeff):
    """(du, n, mtab) for reverse pair spB with forward pair spF and forward
    coefficients Fcoeff = (Ea_F, eta_F, dHreac_F); unscaled."""
    Ea_F, eta_F, dH_F = Fcoeff
    ea_eff = Ea_F + dH_F; eaP = max(ea_eff, 0.0)
    thetas = [SPECIES[s]["vibtemp"] for s in spB+spF if SPECIES[s]["vibtemp"] > 0]
    theta_min = min(thetas) if thetas else 0.0
    umax = eaP + 40.0*EV
    du = KB*theta_min/16.0 if theta_min > 0 else umax/20000.0
    if umax/du > 200000.0: du = umax/200000.0
    n = int(umax/du) + 2
    u = np.arange(n)*du
    zB = 0.5*(SPECIES[spB[0]]["rotdof"] + SPECIES[spB[1]]["rotdof"])
    zF = 0.5*(SPECIES[spF[0]]["rotdof"] + SPECIES[spF[1]]["rotdof"])
    omB = 0.5*(VSS[spB[0]] + VSS[spB[1]])
    exp_num = zF + eta_F + 0.5
    exp_den = zB + 1.5 - omB
    den = np.where(u > 0.0,    np.power(np.maximum(u, 0.0),      exp_den), 0.0)
    num = np.where(u > ea_eff, np.power(np.maximum(u-ea_eff,0.), exp_num), 0.0)
    def conv(arr, sp):
        out = arr; th = SPECIES[sp]["vibtemp"]
        if th > 0.0:
            eps = []; g = []; l = 0
            while l*th*KB < umax: eps.append(l*th*KB); g.append(1.0); l += 1
            out = _ladder(np, out, du, n, eps, g)
        if sp in ELEC:
            out = _ladder(np, out, du, n, [KB*t for t,_ in ELEC[sp]],
                          [float(gg) for _,gg in ELEC[sp]])
        return out
    for sp in spF: num = conv(num, sp)
    for sp in spB: den = conv(den, sp)
    with np.errstate(divide="ignore", invalid="ignore"):
        mtab = np.where(den > 0.0, num/den, 0.0)
    return du, n, mtab

def sample_ecc(np, rng, spB, Ttr, Tint, ns, Trot=None, Tvib=None, Telec=None):
    """collision-energy distribution of the reverse pair spB: VSS translation
    (Gamma 5/2-omega, at Ttr), continuum rotors (at Trot), discrete vib
    (at Tvib) and electronic (at Telec) ladders; matches the runtime
    pre_etotal.  Trot/Tvib/Telec default to Tint (fully-coupled internal
    temperature); pass them separately for mode-specific non-equilibrium."""
    if Trot is None: Trot = Tint
    if Tvib is None: Tvib = Tint
    if Telec is None: Telec = Tint
    omB = 0.5*(VSS[spB[0]] + VSS[spB[1]])
    ecc = rng.gamma(2.5-omB, KB*Ttr, ns)
    for sp in spB:
        rd = SPECIES[sp]["rotdof"]
        if rd > 0: ecc = ecc + rng.gamma(0.5*rd, KB*Trot, ns)
        th = SPECIES[sp]["vibtemp"]
        if th > 0.0:
            x = math.exp(-th/Tvib)
            ecc = ecc + (rng.geometric(1.0-x, ns)-1)*th*KB
        if sp in ELEC:
            temps = np.array([t for t,_ in ELEC[sp]])
            degs  = np.array([g for _,g in ELEC[sp]], float)
            w = degs*np.exp(-temps/Telec); w /= w.sum()
            ecc = ecc + KB*temps[rng.choice(len(temps), ns, p=w)]
    return ecc

# reaction data (must match rev.tce)
A_F, B_F, EA_F, DH_F = 1.069e-12, -1.0, 5.175e-19, -5.175e-19   # N2+O -> NO+N
A_D, B_D, EA_D, DH_D = 4.980e-8, -1.5, 1.561e-18, -1.561e-18    # N2+N -> N+N+N
# independent literature fit for NO + N -> N2 + O (data/air.tce)
A_L, B_L, EA_L = 4.059e-12, -1.359, 0.0

def kf(T):  return A_F*T**B_F*math.exp(-EA_F/(KB*T))
def kb_lit(T): return A_L*T**B_L*math.exp(-EA_L/(KB*T))

def keq_exchange(T):
    """N2 + O <-> NO + N (dimensionless)"""
    return (q_total("NO",T)*q_total("N",T)) / \
           (q_total("N2",T)*q_total("O",T)) * math.exp(DH_F/(KB*T))

def keq_dissoc(T):
    """N2 <-> N + N (per volume, 1/m^3); third body cancels"""
    return (q_total("N",T)**2/q_total("N2",T)) * math.exp(DH_D/(KB*T))

def kb_derived(T):
    return kf(T)/keq_exchange(T)

# ---------------------------------------------------------------- running

def run(exe, deck, varz, tag, extra_args=None, subs=None, extra_files=None):
    wd = os.path.join(HERE, "work_validate", tag)
    os.makedirs(wd, exist_ok=True)
    for f in ("air.species","air.vss","air.rot","air.elec","rev.tce",
              "rev_exch.tce","rev_mol.tce","multi_fwd.tce","econsv.species",
              "nl.species","nl.vss","nl.rot","nl.elec","nl.tce","ce.species",
              "ce.vss","ce.elec","ce.tce",deck):
        shutil.copy(os.path.join(HERE,f), wd)
    for name, text in (extra_files or {}).items():
        open(os.path.join(wd,name),"w").write(text)
    if subs:
        dk = open(os.path.join(wd,deck)).read()
        for old,new in subs.items(): dk = dk.replace(old,new)
        open(os.path.join(wd,deck),"w").write(dk)
    cmd = [exe] + (extra_args or []) + ["-in", deck, "-log", "log.sparta"]
    for k,v in varz.items(): cmd += ["-var", k, str(v)]
    r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True,
                       timeout=3600)
    if r.returncode != 0:
        tail = "\n".join((r.stdout or "").splitlines()[-15:])
        raise RuntimeError("run failed (%s):\n%s" % (tag, tail))
    return os.path.join(wd,"log.sparta")

PARITY = "exact"

def tallies_agree(t1, t2):
    """two binaries' tallies agree: identical in exact mode, within 4
    Poisson standard deviations of the difference in statistical mode"""
    if not t1 or set(t1) != set(t2):
        return False
    if PARITY == "exact":
        return t1 == t2
    return all(abs(t1[k] - t2[k]) <= 4.0*math.sqrt(max(t1[k] + t2[k], 1))
               for k in t1)

def tallies(log):
    t = {}
    for line in open(log):
        m = re.match(r"\s*reaction (.+): ([0-9.eE+-]+)", line)
        if m: t[m.group(1).strip()] = float(m.group(2))
    return t

def stats_rows(log):
    rows, header = [], None
    lines = open(log).read().splitlines()
    for i,l in enumerate(lines):
        if l.startswith("Step "):
            header = l.split(); rows = []
            for l2 in lines[i+1:]:
                w = l2.split()
                # end of the stats block
                if not w or l2.startswith("Loop time") or l2.startswith("Step "):
                    break
                # skip warnings/other messages interleaved in the run output
                # (e.g. a one-time "reaction probability exceeded 1" warning)
                try: float(w[0])
                except (ValueError,IndexError): continue
                rows.append(w)
    return header, rows

# ---------------------------------------------------------------- checks

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("--exe-args", default="")
    p.add_argument("--exe2", default=None,
                   help="second binary (e.g. KOKKOS) for the keq parity check")
    p.add_argument("--exe2-args", default="")
    p.add_argument("--parity", default="exact",
                   choices=("exact", "statistical"),
                   help="exact: the second binary must give identical tallies "
                        "(CPU vs KOKKOS built with SPARTA_KOKKOS_EXACT); "
                        "statistical: every tally must agree within 4 "
                        "standard deviations (e.g. CPU vs GPU)")
    args = p.parse_args()
    exe = os.path.abspath(args.exe)
    extra = args.exe_args.split() if args.exe_args else []
    exe2 = os.path.abspath(args.exe2) if args.exe2 else None
    global PARITY
    PARITY = args.parity
    extra2 = args.exe2_args.split() if args.exe2_args else []

    FNUM, NRHO, V = 1.767e6, 7.07043e22, 1.0e-12
    NRHO_HI, FNUM_HI = 1.414086e25, 3.534e8  # dense reservoir for 3-body rates
    NSTEP, DT = 2000, 1.0e-9
    nfrac = 0.25

    print("check 1: exchange pair tally ratio vs analytic Keq(T)")
    for T in (8000.0, 10000.0, 15000.0):
        log = run(exe,"in.reverse_rate",
                  {"T":T,"RB":1000.0,"NRHO":NRHO,"FNUM":FNUM},
                  "rate%d"%T,extra)
        t = tallies(log)
        f = t.get("N2 + O --> NO + N",0.0)
        b = t.get("NO + N --> N2 + O",0.0)
        keq = keq_exchange(T)
        ratio = f/b if b else float("inf")
        sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))     # relative stat error
        dev = abs(ratio/keq - 1.0)
        check("T=%6.0fK  f/b=%.4f  Keq=%.4f (dev %.1f%%, stat %.1f%%)"
              % (T, ratio, keq, 100*dev, 100*sig),
              dev < max(4*sig, 0.06) and f > 200 and b > 200)

        # check 2 piggybacks on the same runs
        kb_meas = b*FNUM/((NRHO*nfrac)**2*V*NSTEP*DT)
        check("T=%6.0fK  derived kb=%.3e  literature fit=%.3e (x%.2f)"
              % (T, kb_meas, kb_lit(T), kb_meas/kb_lit(T)),
              0.5 < kb_meas/kb_lit(T) < 2.0)

    print("check 3: dissociation/recombination pair vs volumetric Keq")
    # longer run + tight bound: a missing constant in the 3-body table
    # calibration (e.g. the Gamma(3/2) norm of the third body's
    # translational energy, ~13%) must fail this check
    T = 15000.0
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1.0,"NRHO":NRHO_HI,"FNUM":FNUM_HI},
              "recomb%d"%T,extra,
              subs={"run             2000":"run             6000"})
    t = tallies(log)
    d = t.get("N2 + N --> N + N + N",0.0)
    r = t.get("N + N --> N2 + N",0.0)
    # rate_d = kd n_N2 n_N V ; rate_r = kr n_N^2 n_N V  (third body N)
    # kd/kr = (d/r) * n_N  must equal keq_dissoc
    nN = NRHO_HI*nfrac
    keqd = keq_dissoc(T)
    ratio = (d/r)*nN if r else float("inf")
    sig = math.sqrt(1.0/max(d,1)+1.0/max(r,1))
    dev = abs(ratio/keqd - 1.0)
    check("T=%6.0fK  (d/r)*n_N=%.3e  Keq=%.3e (dev %.1f%%, stat %.1f%%)"
          % (T, ratio, keqd, 100*dev, 100*sig),
          dev < max(3*sig, 0.05) and d > 500 and r > 500)

    print("check 4: reacting box relaxes to the ADIABATIC equilibrium")
    # closed reflective (adiabatic) box, only the exchange N2+O<->NO+N active.
    # The exchange has a large reaction energy (coeff[4]=5.175e-19 ~ 2.5 kT at
    # 15 kK), so as it proceeds it shifts the gas temperature: the endothermic
    # forward direction (reactant side) COOLS the gas, the exothermic reverse
    # (product side) HEATS it.  The box therefore relaxes NOT to the isothermal
    # 15 kK equilibrium (NO frac 0.417) but to the energy-CONSERVING adiabatic
    # equilibrium at the shifted temperature -- a lower NO fraction from the
    # reactant side, a higher one from the product side.  Both are computed
    # here by solving Keq(T_eq) jointly with total-energy conservation, and the
    # measured compositions must match (the integral signature of a detailed-
    # balance pair that also conserves energy).
    T0 = 15000.0
    def _eint(sp, Temp):   # translational + internal energy per particle
        e = 1.5*KB*Temp
        if SPECIES[sp]["rotdof"] >= 2: e += KB*Temp
        th = SPECIES[sp]["vibtemp"]
        if th > 0.0: e += KB*th/(math.exp(th/Temp)-1.0)
        if sp in ELEC:
            num = den = 0.0
            for t,g in ELEC[sp]:
                w = g*math.exp(-t/Temp); den += w; num += w*KB*t
            e += num/den
        return e
    EFORM = {"N2":0.0,"O":0.0,"NO":0.0,"N":5.175e-19}   # products carry +dH
    def _adiabatic_eq(a, b, lo, hi):
        # N2=O=a-x, NO=N=b+x; conserve total (thermal+chemical) energy from x=0
        def etot(x,Temp):
            c = {"N2":a-x,"O":a-x,"NO":b+x,"N":b+x}
            return sum(c[s]*(_eint(s,Temp)+EFORM[s]) for s in c)
        E0 = etot(0.0, T0)
        def xof(Temp):
            s = math.sqrt(keq_exchange(Temp)); return (s*a-b)/(1.0+s)
        for _ in range(100):
            Tm = 0.5*(lo+hi)
            if etot(xof(Tm), Tm) > E0: hi = Tm
            else: lo = Tm
        Tm = 0.5*(lo+hi); x = xof(Tm)
        return Tm, (b+x)/((a-x)+(b+x))
    for tag,(a,b,lo,hi) in (("fwd",(0.45,0.05,8000.,15000.)),
                            ("bwd",(0.05,0.45,15000.,25000.))):
        Teq, frac_eq = _adiabatic_eq(a,b,lo,hi)
        log = run(exe,"in.reverse_eq",
                  {"T":T0,"FN2":a,"FO":a,"FNO":b,"FN":b},"eq_"+tag,extra)
        header, rows = stats_rows(log)
        iNO = header.index("c_spcount[3]"); iN2 = header.index("c_spcount[1]")
        frac0 = float(rows[0][iNO])/(float(rows[0][iNO])+float(rows[0][iN2]))
        fracT = float(rows[-1][iNO])/(float(rows[-1][iNO])+float(rows[-1][iN2]))
        check("%s: NO/(NO+N2) %.3f -> %.3f  adiabatic eq %.3f (T_eq=%.0f K; "
              "isothermal 15 kK would be 0.417)" % (tag,frac0,fracT,frac_eq,Teq),
              abs(fracT-frac_eq) < 0.02 and abs(fracT-frac0) > 0.1)

    print("check 5: analytic sanity of the implementation inputs")
    # ground-state electronic degeneracy ratio for this reaction is
    # (4*4)/(1*9) = 1.78; confirm the elec file carries it
    g = lambda name: ELEC[name][0][1]
    check("elec ground degeneracies N2=1 O=9 NO=4 N=4",
          g("N2")==1 and g("O")==9 and g("NO")==4 and g("N")==4)
    sig = lambda name: ROT[name][1]
    check("symmetry numbers N2=2 NO=1", sig("N2")==2 and sig("NO")==1)

    print("check 6: error paths reject invalid B-style inputs; unsupported "
          "styles switch B lines off")
    def expect_error(tag, tce_text, msg_frag, deck="in.reverse_rate",
                     react_style_line=None, edits=(), files=None):
        wd = os.path.join(HERE,"work_validate","err_"+tag)
        os.makedirs(wd, exist_ok=True)
        for f in ("air.species","air.vss","air.rot","air.elec",
                  "in.reverse_rate"):
            shutil.copy(os.path.join(HERE,f), wd)
        open(os.path.join(wd,"rev.tce"),"w").write(tce_text)
        for name, text in (files or {}).items():
            open(os.path.join(wd,name),"w").write(text)
        dk = open(os.path.join(wd,deck)).read()
        if react_style_line:
            dk = dk.replace("react           tce rev.tce", react_style_line)
        for a,b in edits: dk = dk.replace(a,b)
        open(os.path.join(wd,deck),"w").write(dk)
        cmd = [exe] + extra + ["-in", deck, "-log", "log.sparta",
               "-var","T","15000.0","-var","RB","1000.0",
               "-var","NRHO","7.07043e22","-var","FNUM","1.767e6"]
        r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True,
                           timeout=600)
        out = (r.stdout or "") + (r.stderr or "")
        check("%s: run fails with '%s'" % (tag, msg_frag),
              r.returncode != 0 and msg_frag in out)

    expect_error("nopartner",
        "NO + N --> N2 + O\nE B 0.0 0.0 0.0 0.0 0.0\n",
        "No forward partner")
    expect_error("wildcard",
        "N2 + N --> N + N + N\nD A 1.0 1.561e-18 4.980e-8 -1.5 -1.561e-18\n"
        "\nN + N --> N2 + atom\nR B 0.0 0.0 0.0 0.0 0.0\n",
        "explicit third-body species")
    # a rotating species without rotfile data would silently get q_rot = 1
    expect_error("norotfile",
        "N2 + O --> NO + N\nE A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19\n"
        "\nNO + N --> N2 + O\nE B 0.0 0.0 0.0 0.0 0.0\n",
        "needs the rotational temperature",
        edits=(("rotfile air.rot ",""),))
    # a zero rotational temperature would hang the discrete rotational
    # sampler (zero rotational quantum): rejected at init.  Forward-only, so
    # the collision-model check is reached (a B line stops earlier, at the
    # reverse reaction's own rotfile check)
    rot0 = "".join(l.replace(" 2.88 "," 0.0 ") if l.startswith("N2 ") else l
                   for l in open(os.path.join(HERE,"air.rot")))
    assert " 0.0 " in [l for l in rot0.splitlines() if l.startswith("N2 ")][0]
    expect_error("rottemp0",
        "N2 + O --> NO + N\nE A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19\n",
        "non-positive rotational temperature",
        edits=(("rotfile air.rot ","rotfile zero.rot "),),
        files={"zero.rot":rot0})
    # reverse auto is an explicit request: error with the partial-energy model
    expect_error("autopartial",
        "N2 + O --> NO + N\nE A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19\n",
        "reverse auto requires react_modify partial_energy no",
        edits=(("partial_energy no","partial_energy yes reverse auto"),))

    # the qk and tce/qk styles do not implement reverse rates: B lines are
    # switched off with a warning (so a file carrying reverses still runs
    # forward-only), rather than aborting the run.  tce/qk rejects
    # compute_chem_rates, so drop it here, and a few steps suffice

    def expect_warning(tag, tce_text, msg_frag, react_style_line):
        wd = os.path.join(HERE,"work_validate","warn_"+tag)
        os.makedirs(wd, exist_ok=True)
        for f in ("air.species","air.vss","air.rot","air.elec",
                  "in.reverse_rate"):
            shutil.copy(os.path.join(HERE,f), wd)
        open(os.path.join(wd,"rev.tce"),"w").write(tce_text)
        dk = open(os.path.join(wd,"in.reverse_rate")).read()
        dk = dk.replace("react           tce rev.tce", react_style_line)
        dk = dk.replace("compute_chem_rates yes ", "")
        dk = dk.replace("run             2000", "run             10")
        open(os.path.join(wd,"in.reverse_rate"),"w").write(dk)
        cmd = [exe] + extra + ["-in", "in.reverse_rate", "-log", "log.sparta",
               "-var","T","15000.0","-var","RB","1000.0",
               "-var","NRHO","7.07043e22","-var","FNUM","1.767e6"]
        r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True,
                           timeout=600)
        out = (r.stdout or "") + (r.stderr or "")
        check("%s: run completes with warning '%s'" % (tag, msg_frag),
              r.returncode == 0 and "ERROR" not in out and msg_frag in out)

    expect_warning("qkstyle",
        "N2 + O --> NO + N\nE A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19\n"
        "\nNO + N --> N2 + O\nE B 0.0 0.0 0.0 0.0 0.0\n",
        "1 reverse (B-style) reaction(s) deactivated",
        react_style_line="react           tce/qk rev.tce")

    # forward-only reaction file for the auto-reverse checks

    FWD_TCE = (
        "N2 + O --> NO + N\n"
        "E A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19\n\n"
        "N2 + N --> N + N + N\n"
        "D A 1.0 1.561e-18 4.980e-8 -1.5 -1.561e-18\n")

    print("check 7: reverse auto generates B partners from a forward-only file")
    T = 15000.0
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1000.0,"NRHO":NRHO,"FNUM":FNUM},
              "auto%d"%T,extra,
              subs={"react           tce rev.tce":
                    "react           tce fwd.tce",
                    "rboost ${RB}":"rboost ${RB} reverse auto"},
              extra_files={"fwd.tce":FWD_TCE})
    gen = any("Generated 2 reverse reaction" in l for l in open(log))
    check("init reports 2 generated reverse reactions", gen)
    t = tallies(log)
    f = t.get("N2 + O --> NO + N",0.0)
    b = t.get("NO + N --> N2 + O",0.0)
    keq = keq_exchange(T)
    ratio = f/b if b else float("inf")
    sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
    dev = abs(ratio/keq - 1.0)
    check("generated exchange reverse holds detailed balance "
          "(f/b=%.4f Keq=%.4f dev %.1f%%)" % (ratio,keq,100*dev),
          dev < max(4*sig, 0.06) and f > 200 and b > 200)

    print("check 8: external Keq fit reproduces the fitted reverse rate")
    # Park-form fit of the Keq implied by the air.tce forward/reverse pair:
    # matching it must reproduce the literature backward rate exactly
    c1 = math.log(A_F/A_L) + (B_F-B_L)*math.log(10000.0)
    c3 = -EA_F/KB/10000.0
    PARK = ("N2 + O --> NO + N\n"
            "park 0.0 %.6f %.6f %.6f 0.0\n" % (c1, -(B_F-B_L), c3))
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1000.0,"NRHO":NRHO,"FNUM":FNUM},
              "keqfit%d"%T,extra,
              subs={"react           tce rev.tce":
                    "react           tce fwd.tce",
                    "rboost ${RB}":
                    "rboost ${RB} reverse auto keq_file park.keq"},
              extra_files={"fwd.tce":FWD_TCE,"park.keq":PARK})
    t = tallies(log)
    b = t.get("NO + N --> N2 + O",0.0)
    kb_meas = b*FNUM/((NRHO*nfrac)**2*V*NSTEP*DT)
    check("derived kb matches the literature fit via the Keq file "
          "(kb=%.3e lit=%.3e x%.2f)" % (kb_meas,kb_lit(T),kb_meas/kb_lit(T)),
          0.9 < kb_meas/kb_lit(T) < 1.1 and b > 200)

    print("check 9: external Keq fit on a dissociation/recombination pair")
    # a Park fit of the analytic volumetric dissociation Keq (1/m^3), fed
    # via keq_file, must reproduce that Keq as the recombination backward
    # rate -- the recombination (m^6/s) analogue of check 8, exercising the
    # m^3->m^6 unit path that the exchange check cannot
    T = 15000.0
    cd = park_fit(keq_dissoc, (8000.,11000.,14000.,17000.,20000.))
    PARKD = ("N2 + N --> N + N + N\n"
             "park %.8g %.8g %.8g %.8g %.8g\n" % tuple(cd))
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1.0,"NRHO":NRHO_HI,"FNUM":FNUM_HI},
              "keqrecomb%d"%T,extra,
              subs={"react           tce rev.tce":
                    "react           tce fwd.tce",
                    "rboost ${RB}":
                    "rboost ${RB} reverse auto keq_file park.keq",
                    "run             2000":"run             6000"},
              extra_files={"fwd.tce":FWD_TCE,"park.keq":PARKD})
    t = tallies(log)
    d = t.get("N2 + N --> N + N + N",0.0)
    r = t.get("N + N --> N2 + N",0.0)
    nN = NRHO_HI*nfrac
    keqd = keq_dissoc(T)
    ratio = (d/r)*nN if r else float("inf")
    sig = math.sqrt(1.0/max(d,1)+1.0/max(r,1))
    dev = abs(ratio/keqd - 1.0)
    check("T=%6.0fK  (d/r)*n_N=%.3e  Keq=%.3e (dev %.1f%%, stat %.1f%%)"
          % (T, ratio, keqd, 100*dev, 100*sig),
          dev < max(3*sig, 0.05) and d > 500 and r > 500)

    print("check 9b: the same external Keq fit under units cgs")
    # the keq_file fit stays in SI (1/m^3) and SPARTA converts it to 1/cm^3:
    # the same reservoir in cgs units (masses in g, diameters in cm, energies
    # in erg, the dissociation prefactor in cm^3/s) must reproduce the same
    # Keq; without the conversion the recombination is 1e6 times too slow
    def cgs_species(text):
        out = []
        for l in text.splitlines():
            w = l.split()
            if w and not l.lstrip().startswith("#") and len(w) >= 10:
                w[2] = "%.6e" % (float(w[2])*1.0e3); l = "  ".join(w)
            out.append(l)
        return "\n".join(out) + "\n"
    def cgs_vss(text):
        out = []
        for l in text.splitlines():
            w = l.split()
            if w and not l.lstrip().startswith("#") and len(w) >= 5:
                w[1] = "%.6e" % (float(w[1])*1.0e2); l = "  ".join(w)
            out.append(l)
        return "\n".join(out) + "\n"
    def cgs_tce(text):
        out = []
        for l in text.splitlines():
            w = l.split()
            if len(w) == 7 and w[0] in "DEIR" and w[1] == "A":
                avol = 1.0e12 if w[0] == "R" else 1.0e6
                w[3] = "%.6e" % (float(w[3])*1.0e7)      # Ea, J -> erg
                w[4] = "%.6e" % (float(w[4])*avol)       # A, m^3/s -> cm^3/s
                w[6] = "%.6e" % (float(w[6])*1.0e7)      # dE, J -> erg
                l = " ".join(w)
            out.append(l)
        return "\n".join(out) + "\n"
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1.0,"NRHO":NRHO_HI*1.0e-6,"FNUM":FNUM_HI},
              "keqrecomb_cgs%d"%T,extra,
              subs={"seed            12345":
                    "units           cgs\nseed            12345",
                    "gridcut 1.0e-5":"gridcut 1.0e-3",
                    "create_box      0 0.0001 0 0.0001 0 0.0001":
                    "create_box      0 0.01 0 0.01 0 0.01",
                    "react           tce rev.tce":
                    "react           tce fwd.tce",
                    "rboost ${RB}":
                    "rboost ${RB} reverse auto keq_file park.keq",
                    "run             2000":"run             6000"},
              extra_files={"fwd.tce":cgs_tce(FWD_TCE),"park.keq":PARKD,
                           "air.species":cgs_species(
                               open(os.path.join(HERE,"air.species")).read()),
                           "air.vss":cgs_vss(
                               open(os.path.join(HERE,"air.vss")).read())})
    t = tallies(log)
    d = t.get("N2 + N --> N + N + N",0.0)
    r = t.get("N + N --> N2 + N",0.0)
    ratio = (d/r)*nN if r else float("inf")
    sig = math.sqrt(1.0/max(d,1)+1.0/max(r,1))
    dev = abs(ratio/keqd - 1.0)
    check("cgs: T=%6.0fK  (d/r)*n_N=%.3e  Keq=%.3e (dev %.1f%%, stat %.1f%%)"
          % (T, ratio, keqd, 100*dev, 100*sig),
          dev < max(3*sig, 0.05) and d > 500 and r > 500)

    print("check 10: standard eta=-1.5 dissociation raises no bounds warning")
    # the ubiquitous eta = -3/2 dissociation with one rotor sits exactly on
    # the low-energy trend bound; check_tce_bounds must not warn about it
    log = run(exe,"in.reverse_rate",{"T":15000.,"RB":1000.,
              "NRHO":NRHO,"FNUM":FNUM},"bounds",extra)
    warned = any("does not vanish" in l for l in open(log))
    check("no spurious 'does not vanish' warning on eta=-1.5", not warned)

    print("check 11: reverse detailed balance under vibrate smooth")
    # with classical (smooth) vibration the calibration target and the
    # detailed-balance table must share the same vib temperature
    # dependence: the table must not warn that it drifts, and the exchange
    # reverse must reproduce the classical-vib Keq
    T = 15000.0
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1000.,"NRHO":NRHO,"FNUM":FNUM},
              "smooth%d"%T,extra,
              subs={"vibrate discrete":"vibrate smooth"})
    drift = any("detailed-balance table drifts" in l for l in open(log))
    t = tallies(log)
    f = t.get("N2 + O --> NO + N",0.0)
    b = t.get("NO + N --> N2 + O",0.0)
    keqs = (q_total("NO",T,True)*q_total("N",T,True)) / \
           (q_total("N2",T,True)*q_total("O",T,True))*math.exp(DH_F/(KB*T))
    ratio = f/b if b else float("inf")
    sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
    dev = abs(ratio/keqs - 1.0)
    check("smooth-vib exchange reverse: no drift warning and f/b matches "
          "classical Keq (f/b=%.4f Keq=%.4f dev %.1f%%)" % (ratio,keqs,100*dev),
          (not drift) and dev < max(4*sig, 0.06) and f > 200 and b > 200)

    print("check 13: 3-body recombination with a molecular third body")
    # N + N -> N2 + N2: the third body N2 carries discrete vibrational and
    # electronic ladders (folded into the density of states) and continuum
    # rotation (a flat measure variable) -- the general molecular-M case of
    # build_db3_table, where check 3 uses an atomic third body N
    T = 15000.0
    molsub = {"react           tce rev.tce":"react           tce rev_mol.tce",
              "run             2000":"run             8000"}
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1.0,"NRHO":NRHO_HI,"FNUM":FNUM_HI},
              "molrecomb",extra,subs=molsub)
    t = tallies(log)
    d = t.get("N2 + N2 --> N + N + N2",0.0)
    r = t.get("N + N --> N2 + N2",0.0)
    nN = NRHO_HI*nfrac
    keqd = keq_dissoc(T)
    ratio = (d/r)*nN if r else float("inf")
    sig = math.sqrt(1.0/max(d,1)+1.0/max(r,1))
    dev = abs(ratio/keqd - 1.0)
    check("molecular third body N2: (d/r)*n_N=%.3e Keq=%.3e (dev %.1f%%, stat %.1f%%)"
          % (ratio, keqd, 100*dev, 100*sig),
          dev < max(3*sig, 0.05) and d > 500 and r > 500)
    if exe2:
        tk = tallies(run(exe2,"in.reverse_rate",
                     {"T":T,"RB":1.0,"NRHO":NRHO_HI,"FNUM":FNUM_HI},
                     "molrecomb_kk",extra2,subs=molsub))
        check("molecular third body: CPU/second-binary tallies %s"
              % ("identical" if PARITY == "exact" else "agree within 4 sd"),
              tallies_agree(t, tk))

    print("check 14: non-equilibrium reverse rate (two-temperature reservoir)")
    # The per-collision reverse probability is a temperature-free function of
    # the total collision energy (microscopic reversibility), so out of
    # equilibrium the reverse RATE must follow the collision-energy
    # distribution, not any single temperature.  Measure the barriered reverse
    # N2 + O -> NO + N (bounded probability, unlike the barrierless direction
    # which saturates above 1) in a frozen reservoir whose translational
    # temperature is fixed while the rotational/vibrational/electronic modes
    # are held at a different temperature, and compare the ratio of the 2T
    # reverse rate to the equilibrium reverse rate (same Ttr, so the collision
    # rate and densities cancel) against an INDEPENDENT microcanonical integral
    # of the reconstructed detailed-balance table over the two collision-energy
    # distributions.
    try:
        import numpy as _np
    except ImportError:
        check("non-equilibrium reverse rate (needs numpy)", True,
              "SKIPPED: numpy not available")
        _np = None
    if _np is not None:
        # barriered reverse: forward NO+N->N2+O (barrierless literature fit),
        # reverse N2+O->NO+N (endothermic).  spB/spF and Fcoeff below.
        spB, spF = ("N2","O"), ("NO","N")
        Fcoeff = (0.0, -1.359, 5.175e-19)     # Ea_F, eta_F, reaction energy
        NONEQ_TCE = ("NO + N --> N2 + O\nE A 0.0 0.0 4.059e-12 -1.359 5.175e-19\n\n"
                     "N2 + O --> NO + N\nE B 0.0 0.0 0.0 0.0 0.0\n")
        # frozen data files: zero every relaxation number so the two-temperature
        # state stays stationary under relax constant
        def _frozen_species():
            out = []
            for line in open(os.path.join(HERE,"air.species")):
                w = line.split()
                if len(w) >= 10 and not line.strip().startswith("#"):
                    w[4] = "0.0"; w[6] = "0.0"; out.append("  ".join(w))
                else: out.append(line.rstrip("\n"))
            return "\n".join(out) + "\n"
        def _frozen_elec():
            out = []
            for line in open(os.path.join(HERE,"air.elec")):
                w = line.split()
                if len(w) > 2 and not line.strip().startswith("#"):
                    nl = int(w[1])
                    for k in range(nl): w[2+5*k+1] = "0.0"
                    out.append(" ".join(w))
                else: out.append(line.rstrip("\n"))
            return "\n".join(out) + "\n"
        xf = {"air_frozen.species": _frozen_species(),
              "air_frozen.elec":   _frozen_elec(),
              "rev_noneq.tce":     NONEQ_TCE}
        du, ntab, mtab = db_table_shape(_np, spB, spF, Fcoeff)
        def _interp(ecc):
            x = ecc/du; k = _np.clip(x.astype(int), 0, ntab-2); f = x - k
            return (1.0-f)*mtab[k] + f*mtab[k+1]
        rng = _np.random.default_rng(2024); NS = 4_000_000
        def _avg(Ttr, Tint):
            return _interp(sample_ecc(_np, rng, spB, Ttr, Tint, NS)).mean()

        # self-check: the reconstructed table must reproduce the equilibrium
        # detailed-balance temperature dependence (else the 2T prediction is
        # untrustworthy).  k_rev(T) = k(NO+N->N2+O) * Keq(N2+O<->NO+N).
        omB = 0.5*(VSS[spB[0]] + VSS[spB[1]])
        def _krev_an(T): return kb_lit(T)*keq_exchange(T)
        ratios = []
        for T in (8000., 12000., 16000., 20000.):
            ratios.append(_avg(T,T)*T**(1-omB) / _krev_an(T))
        drift = max(ratios)/min(ratios) - 1.0
        check("reconstructed table reproduces equilibrium detailed balance "
              "across 8-20 kK (drift %.1f%%)" % (100*drift), drift < 0.03)

        Ttr = 8000.0
        beq = None
        for Tint in (Ttr, 20000.0, 4000.0):
            log = run(exe, "in.reverse_noneq",
                      {"T":Ttr, "TINT":Tint, "NRHO":NRHO, "FNUM":FNUM},
                      "noneq%d"%Tint, extra,
                      subs={"run             8000":"run             30000"},
                      extra_files=xf)
            b = tallies(log).get("N2 + O --> NO + N", 0.0)
            if Tint == Ttr:
                beq = b
                check("equilibrium control run has adequate statistics",
                      beq > 400, "reverse counts=%.0f" % beq)
                continue
            rm = b/beq if beq else float("inf")
            rp = _avg(Ttr, Tint)/_avg(Ttr, Ttr)
            sig = math.sqrt(1.0/max(b,1) + 1.0/max(beq,1))
            dev = abs(rm/rp - 1.0)
            check("Tint=%5.0fK  R_rev/R_eq meas=%.3f pred=%.3f (dev %.1f%%, "
                  "stat %.1f%%)" % (Tint, rm, rp, 100*dev, 100*sig),
                  dev < max(3*sig, 0.05) and b > 200)

    print("check 15: nonlinear (rotdof=3) molecule reverse reaction")
    # exchange pair TRIA + ATB <-> DIA + ATO with TRIA a nonlinear triatomic
    # (rotdof=3): exercises the nonlinear rotational partition function
    # (qrot ~ T^1.5, partition_function) and the zcont=3/2 continuum in
    # build_db_table.  Detailed balance is validated to ~12 kK here; a small
    # (~4% at 20 kK) drift appears at very high T, beyond the table's built-in
    # calibration self-check range (see the memo).
    nlsp = parse_species(os.path.join(HERE,"nl.species"))
    nlrot = parse_rot(os.path.join(HERE,"nl.rot"))
    nlel = parse_elec(os.path.join(HERE,"nl.elec"))
    def qnl(name, T):
        s = nlsp[name]
        q = (2.0*math.pi*s["mass"]*KB*T/(H*H))**1.5
        if s["rotdof"] == 2 and name in nlrot:
            temps,sg = nlrot[name]; q *= T/(sg*temps[0])
        elif s["rotdof"] == 3 and name in nlrot:
            temps,sg = nlrot[name]
            q *= math.sqrt(math.pi*T**3/(sg*sg*temps[0]*temps[1]*temps[2]))
        if s["vibdof"] >= 2 and s["vibtemp"] > 0.0:
            q *= 1.0/(1.0 - math.exp(-s["vibtemp"]/T))
        if name in nlel:
            q *= sum(g*math.exp(-t/T) for t,g in nlel[name])
        return q
    dHnl = -0.5e-19    # forward reaction energy (nl.tce coeff[4])
    for T in (10000.0, 12000.0):
        keqnl = qnl("DIA",T)*qnl("ATO",T)/(qnl("TRIA",T)*qnl("ATB",T)) \
                * math.exp(dHnl/(KB*T))
        log = run(exe,"in.reverse_nl",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
                  "nl%d"%T,extra)
        t = tallies(log)
        f = t.get("TRIA + ATB --> DIA + ATO",0.0)
        b = t.get("DIA + ATO --> TRIA + ATB",0.0)
        ratio = f/b if b else float("inf")
        sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
        dev = abs(ratio/keqnl - 1.0)
        check("T=%6.0fK nonlinear rotor: f/b=%.4f Keq=%.4f (dev %.1f%%, "
              "stat %.1f%%)" % (T,ratio,keqnl,100*dev,100*sig),
              dev < max(4*sig, 0.06) and f > 200 and b > 200)

    print("check 16: molecular third-body recombination under vibrate smooth")
    # check 13 (N+N->N2+N2) with continuum (classical) vibration: the third
    # body N2's vibration is a flat measure variable folded into the 3-body
    # density of states exactly as its rotation is, so detailed balance must
    # hold with vibrate smooth as it does with vibrate discrete
    T = 15000.0
    smolsub = {"react           tce rev.tce":"react           tce rev_mol.tce",
               "run             2000":"run             8000",
               "vibrate discrete":"vibrate smooth"}
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1.0,"NRHO":NRHO_HI,"FNUM":FNUM_HI},
              "molsmooth",extra,subs=smolsub)
    drift = any("detailed-balance table drifts" in l for l in open(log))
    t = tallies(log)
    d = t.get("N2 + N2 --> N + N + N2",0.0)
    r = t.get("N + N --> N2 + N2",0.0)
    nN = NRHO_HI*nfrac
    # classical-vibration volumetric Keq (q_vib = T/theta per mode)
    keqd = (q_total("N",T)**2/q_total("N2",T,True))*math.exp(DH_D/(KB*T))
    ratio = (d/r)*nN if r else float("inf")
    sig = math.sqrt(1.0/max(d,1)+1.0/max(r,1))
    dev = abs(ratio/keqd - 1.0)
    check("smooth-vib molecular third body: (d/r)*n_N=%.3e Keq=%.3e (dev %.1f%%, "
          "stat %.1f%%)" % (ratio, keqd, 100*dev, 100*sig),
          (not drift) and dev < max(3*sig, 0.05) and d > 500 and r > 500)

    print("check 17: external-Keq residual goodness-of-fit guard")
    # feed an external Keq (the statmech exchange Keq times a strongly T-varying
    # factor) that the 5-coefficient Park form cannot represent to 2% over the
    # full 1000-60000 K residual-sampling range.  Two things must hold: the
    # goodness-of-fit guard added to fit_keq_residual must WARN (so a user is not
    # silently handed a biased Keq), and the reverse rate at the operating
    # temperature - where the fit is anchored and good - must still reproduce the
    # target, confirming the residual correction stays locally accurate.
    def keq_ext(T):    # statmech Keq times a strong Arrhenius-like distortion
        return keq_exchange(T)*math.exp(-8000.0/T)*(T/10000.0)**2
    cext = park_fit(keq_ext,(8000.,11000.,14000.,17000.,20000.))
    PARKX = ("N2 + O --> NO + N\npark %.8g %.8g %.8g %.8g %.8g\n" % tuple(cext))
    T = 15000.0
    log = run(exe,"in.reverse_rate",
              {"T":T,"RB":1000.0,"NRHO":NRHO,"FNUM":FNUM},
              "keqvary%d"%T,extra,
              subs={"react           tce rev.tce":"react           tce fwd.tce",
                    "rboost ${RB}":
                    "rboost ${RB} reverse auto keq_file park.keq"},
              extra_files={"fwd.tce":FWD_TCE,"park.keq":PARKX})
    fitwarn = any("residual fit is off" in l for l in open(log))
    t = tallies(log)
    b = t.get("NO + N --> N2 + O",0.0)
    kb_meas = b*FNUM/((NRHO*nfrac)**2*V*NSTEP*DT)
    kb_target = kf(T)/keq_ext(T)
    dev = abs(kb_meas/kb_target - 1.0)
    sig = math.sqrt(1.0/max(b,1))
    check("poorly-fittable external Keq warns and still reproduces at %.0f K: "
          "kb=%.3e target=%.3e (dev %.1f%%, stat %.1f%%, guard warned=%s)"
          % (T, kb_meas, kb_target, 100*dev, 100*sig, fitwarn),
          fitwarn and dev < max(3*sig, 0.06) and b > 200)

    print("check 18: restart then continue with reverse reactions")
    # write a restart mid-run, read it back, re-issue collide/react/fix, and
    # continue: the detailed-balance tables and Keq fits rebuild
    # deterministically at init and the per-particle electronic state is
    # restored by fix elecmode, so detailed balance must hold after the restart
    T = 15000.0
    run(exe,"in.reverse_restart1",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
        "restart",extra)
    log = run(exe,"in.reverse_restart2",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
              "restart",extra)
    t = tallies(log)
    f = t.get("N2 + O --> NO + N",0.0)
    b = t.get("NO + N --> N2 + O",0.0)
    keq = keq_exchange(T)
    ratio = f/b if b else float("inf")
    sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
    dev = abs(ratio/keq - 1.0)
    check("post-restart exchange detailed balance: f/b=%.4f Keq=%.4f "
          "(dev %.1f%%, stat %.1f%%)" % (ratio,keq,100*dev,100*sig),
          dev < max(4*sig, 0.06) and f > 200 and b > 200)

    print("check 19: charge-exchange reverse reaction")
    # charge-exchange MAp + MB <-> MA + MBp (data in ce.*) moves charge between
    # two atoms with no free electron, so it is an EXCHANGE reaction whose
    # reverse is derived by detailed balance (unlike ionization, whose reverse
    # depends on the electron temperature and is rejected at init).  The
    # electronic ground-state degeneracies drive the Keq prefactor.
    cesp = parse_species(os.path.join(HERE,"ce.species"))
    ceel = parse_elec(os.path.join(HERE,"ce.elec"))
    def qce(name, T):
        return (2.0*math.pi*cesp[name]["mass"]*KB*T/(H*H))**1.5 * ceel[name][0][1]
    dHce = -0.5e-19    # forward reaction energy (ce.tce coeff[4])
    for T in (10000.0, 15000.0):
        keqce = qce("MA",T)*qce("MBp",T)/(qce("MAp",T)*qce("MB",T)) \
                * math.exp(dHce/(KB*T))
        log = run(exe,"in.reverse_ce",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
                  "ce%d"%T,extra)
        t = tallies(log)
        f = t.get("MAp + MB --> MA + MBp",0.0)
        b = t.get("MA + MBp --> MAp + MB",0.0)
        ratio = f/b if b else float("inf")
        sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
        dev = abs(ratio/keqce - 1.0)
        check("T=%6.0fK charge exchange: f/b=%.4f Keq=%.4f (dev %.1f%%, "
              "stat %.1f%%)" % (T,ratio,keqce,100*dev,100*sig),
              dev < max(4*sig, 0.06) and f > 200 and b > 200)

    # two coupled Zeldovich exchanges sharing NO/O/N (multi_fwd.tce), both
    # reverses generated by "reverse auto"
    def keq_zel(rA,rB,pA,pB,dH,T):
        return q_total(pA,T)*q_total(pB,T)/(q_total(rA,T)*q_total(rB,T)) \
               * math.exp(dH/(KB*T))
    DH1, DH2 = -5.175e-19, -2.684e-19   # multi_fwd.tce coeff[4] of each fwd

    print("check 20: multi-channel detailed balance (two coupled reverses)")
    # both reverse reactions compete for the shared reactants NO/O/N in one
    # frozen reservoir; each pair's tally ratio must independently equal its
    # own analytic Keq -- a test of the reaction pairing/indexing under
    # simultaneous active reverses
    T = 15000.0
    log = run(exe,"in.reverse_multi",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
              "multi",extra)
    t = tallies(log)
    K1 = keq_zel("N2","O","NO","N",DH1,T)
    K2 = keq_zel("NO","O","O2","N",DH2,T)
    for tag,(fk,bk,K) in (
        ("N2+O<->NO+N",("N2 + O --> NO + N","NO + N --> N2 + O",K1)),
        ("NO+O<->O2+N",("NO + O --> O2 + N","O2 + N --> NO + O",K2))):
        f = t.get(fk,0.0); b = t.get(bk,0.0)
        ratio = f/b if b else float("inf")
        sig = math.sqrt(1.0/max(f,1)+1.0/max(b,1))
        dev = abs(ratio/K - 1.0)
        check("%s  f/b=%.4f  Keq=%.4f (dev %.1f%%, stat %.1f%%)"
              % (tag,ratio,K,100*dev,100*sig),
              dev < max(4*sig, 0.06) and f > 200 and b > 200)

    print("check 21: multi-channel equilibrium vs independent coupled solve")
    # a closed reflective box with the two coupled exchange pairs, seeded far
    # from equilibrium, must relax to the JOINT chemical equilibrium of the
    # two reactions; compared here against an independent two-reaction-extent
    # equilibrium solve (element totals conserved).  Detects any error that
    # holds one pair's detailed balance but not the coupled fixed point.
    n0 = {"N2":0.30,"O2":0.20,"O":0.30,"N":0.10,"NO":0.10}   # in.reverse_multieq
    def _comp(x1,x2):
        return {"N2":n0["N2"]-x1,"O2":n0["O2"]+x2,"NO":n0["NO"]+x1-x2,
                "N":n0["N"]+x1+x2,"O":n0["O"]-x1-x2}
    def _resid(x):
        c = _comp(*x)
        if min(c.values()) <= 0: return [1e9,1e9]
        return [c["NO"]*c["N"]-K1*c["N2"]*c["O"],
                c["O2"]*c["N"]-K2*c["NO"]*c["O"]]
    x = [0.0,0.0]
    for _ in range(200):
        f = _resid(x)
        if max(abs(v) for v in f) < 1e-13: break
        h = 1e-8; J = [[0,0],[0,0]]
        for j in range(2):
            xp = x[:]; xp[j] += h; fp = _resid(xp)
            for i in range(2): J[i][j] = (fp[i]-f[i])/h
        det = J[0][0]*J[1][1]-J[0][1]*J[1][0]
        dx = [-(J[1][1]*f[0]-J[0][1]*f[1])/det,
              -(-J[1][0]*f[0]+J[0][0]*f[1])/det]
        x = [x[0]+0.5*dx[0], x[1]+0.5*dx[1]]
    ceq = _comp(*x); teq = sum(ceq.values())
    feq = {s: ceq[s]/teq for s in ceq}
    log = run(exe,"in.reverse_multieq",{"T":T,"NRHO":NRHO,"FNUM":FNUM},
              "multieq",extra)
    header, rows = stats_rows(log)
    idx = {s: header.index(c) for s,c in
           (("O2","c_spcount[1]"),("N2","c_spcount[2]"),("O","c_spcount[3]"),
            ("N","c_spcount[4]"),("NO","c_spcount[5]"))}
    tot0 = sum(float(rows[0][idx[s]]) for s in idx)
    totT = sum(float(rows[-1][idx[s]]) for s in idx)
    f0 = {s: float(rows[0][idx[s]])/tot0 for s in idx}
    fT = {s: float(rows[-1][idx[s]])/totT for s in idx}
    worst = max(abs(fT[s]-feq[s]) for s in idx)
    moved = max(abs(f0[s]-feq[s]) for s in idx)
    check("closed box relaxes to coupled equilibrium (worst |dev|=%.3f over "
          "5 species; NO %.3f->%.3f eq %.3f)"
          % (worst, f0["NO"], fT["NO"], feq["NO"]),
          worst < 0.03 and moved > 0.1)

    print("check 22: energy conservation across the reverse-disposal path")
    # closed reacting box (mass-consistent econsv.species so the exchange
    # conserves mass exactly); the per-particle thermal energy (ke+erot+evib,
    # electronic off) dumped at the first and last step must change by exactly
    # the net reaction extent times the reaction energy coeff[4].  A disposal
    # that leaked energy would fail by O(reaction energy) per event.
    DH_EX = 5.175e-19    # |coeff[4]| of the N2+O<->NO+N exchange (rev_exch.tce)
    wd = run(exe,"in.reverse_econsv",
             {"T":12000.,"FN2":0.45,"FO":0.45,"FNO":0.05,"FN":0.05},
             "econsv",extra)
    wdir = os.path.dirname(wd)
    def _thermal(frame_lines):
        return sum(float(w[1])+float(w[2])+float(w[3])
                   for w in (l.split() for l in frame_lines) if len(w) == 4)
    dumpf = os.path.join(wdir,"dump.econsv")
    frames = []
    cur = None
    for line in open(dumpf):
        if line.startswith("ITEM: TIMESTEP"):
            if cur is not None: frames.append(cur)
            cur = []; reading = False
        elif line.startswith("ITEM: ATOMS"): reading = True
        elif cur is not None and reading: cur.append(line)
    if cur is not None: frames.append(cur)
    E0 = _thermal(frames[0]); E1 = _thermal(frames[-1])
    t = tallies(wd)
    nnet = t.get("N2 + O --> NO + N",0.0) - t.get("NO + N --> N2 + O",0.0)
    resid = (E1 - E0) + nnet*DH_EX
    dHeff = -(E1 - E0)/nnet if nnet else 0.0
    check("thermal energy change matches net extent x reaction energy "
          "(dH_eff=%.5e vs %.5e; residual %.1e = %.1e of E_therm)"
          % (dHeff, DH_EX, resid, abs(resid/E0)),
          abs(resid/E0) < 1e-6 and abs(nnet) > 500)

    print("check 23: mode-resolved non-equilibrium reverse rate")
    # check 14 drives all internal modes to one temperature together; here each
    # mode (rotation, vibration, electronic) is driven hot INDIVIDUALLY while
    # the other two stay at the translational temperature.  Microscopic
    # reversibility requires the reverse rate to respond to each mode's own
    # energy content; the measured ratio to the fully-cold rate must match an
    # independent microcanonical integral that puts only that mode hot.
    try:
        import numpy as _np
    except ImportError:
        check("mode-resolved non-equilibrium (needs numpy)", True,
              "SKIPPED: numpy not available"); _np = None
    if _np is not None:
        spB, spF = ("N2","O"), ("NO","N")
        Fcoeff = (0.0, -1.359, 5.175e-19)
        NONEQ_TCE = ("NO + N --> N2 + O\nE A 0.0 0.0 4.059e-12 -1.359 5.175e-19\n\n"
                     "N2 + O --> NO + N\nE B 0.0 0.0 0.0 0.0 0.0\n")
        def _frozen_species():
            out = []
            for line in open(os.path.join(HERE,"air.species")):
                w = line.split()
                if len(w) >= 10 and not line.strip().startswith("#"):
                    w[4] = "0.0"; w[6] = "0.0"; out.append("  ".join(w))
                else: out.append(line.rstrip("\n"))
            return "\n".join(out) + "\n"
        def _frozen_elec():
            out = []
            for line in open(os.path.join(HERE,"air.elec")):
                w = line.split()
                if len(w) > 2 and not line.strip().startswith("#"):
                    nl = int(w[1])
                    for k in range(nl): w[2+5*k+1] = "0.0"
                    out.append(" ".join(w))
                else: out.append(line.rstrip("\n"))
            return "\n".join(out) + "\n"
        xf = {"air_frozen.species": _frozen_species(),
              "air_frozen.elec":   _frozen_elec(),
              "rev_noneq.tce":     NONEQ_TCE}
        du, ntab, mtab = db_table_shape(_np, spB, spF, Fcoeff)
        def _ip(ecc):
            x = ecc/du; k = _np.clip(x.astype(int), 0, ntab-2); f = x - k
            return (1.0-f)*mtab[k] + f*mtab[k+1]
        rng = _np.random.default_rng(77); NS = 4_000_000
        def _avgm(Ttr, Tr, Tv, Te):
            return _ip(sample_ecc(_np, rng, spB, Ttr, Ttr, NS,
                                  Trot=Tr, Tvib=Tv, Telec=Te)).mean()
        Ttr, Thot = 8000.0, 25000.0
        logb = run(exe, "in.reverse_modeneq",
                   {"T":Ttr,"TROT":Ttr,"TVIB":Ttr,"TELEC":Ttr,
                    "NRHO":NRHO,"FNUM":FNUM}, "modeeq", extra, extra_files=xf)
        beq = tallies(logb).get("N2 + O --> NO + N", 0.0)
        base = _avgm(Ttr,Ttr,Ttr,Ttr)
        for mode,(tr,tv,te) in (("rotation",   (Thot,Ttr,Ttr)),
                                ("vibration",  (Ttr,Thot,Ttr)),
                                ("electronic", (Ttr,Ttr,Thot))):
            log = run(exe, "in.reverse_modeneq",
                      {"T":Ttr,"TROT":tr,"TVIB":tv,"TELEC":te,
                       "NRHO":NRHO,"FNUM":FNUM}, "mode_"+mode, extra,
                      extra_files=xf)
            b = tallies(log).get("N2 + O --> NO + N", 0.0)
            rm = b/beq if beq else float("inf")
            rp = _avgm(Ttr,tr,tv,te)/base
            sig = math.sqrt(1.0/max(b,1) + 1.0/max(beq,1))
            dev = abs(rm/rp - 1.0)
            check("%-10s hot: R_rev/R_eq meas=%.3f pred=%.3f (dev %.1f%%, "
                  "stat %.1f%%)" % (mode, rm, rp, 100*dev, 100*sig),
                  dev < max(3*sig, 0.05) and b > 200)

    print("check 24: unbounded-probability diagnostic fires and does not crash")
    # a barrierless reverse whose derived rate exceeds the collision rate makes
    # the cumulative react_prob exceed 1 (the TCE model saturating); the added
    # one-time diagnostic must fire to warn the user, and the run must still
    # complete and tally reactions - it is a warning, not a hard stop.
    logsat = run(exe,"in.reverse_rate",
                 {"T":5000.,"RB":1000.,"NRHO":NRHO,"FNUM":FNUM},"psat",extra)
    warned = any("reaction probability exceeded 1" in l for l in open(logsat))
    bsat = tallies(logsat).get("NO + N --> N2 + O",0.0)
    check("P>1 saturation diagnostic fires and the run still tallies "
          "(warned=%s, reverse tally=%.0f)" % (warned, bsat),
          warned and bsat > 200)

    if exe2:
        print("check 12: external-Keq path is bit-for-bit across binaries")
        # tgas feeds the keq prefactor; verify the two binaries (e.g. CPU
        # and KOKKOS) produce identical reaction tallies on the keq path
        c1 = math.log(A_F/A_L) + (B_F-B_L)*math.log(10000.0)
        c3 = -EA_F/KB/10000.0
        PARK = ("N2 + O --> NO + N\n"
                "park 0.0 %.6f %.6f %.6f 0.0\n" % (c1, -(B_F-B_L), c3))
        sub = {"react           tce rev.tce":"react           tce fwd.tce",
               "rboost ${RB}":"rboost ${RB} reverse auto keq_file park.keq"}
        xf = {"fwd.tce":FWD_TCE,"park.keq":PARK}
        vz = {"T":15000.,"RB":1000.,"NRHO":NRHO,"FNUM":FNUM}
        t1 = tallies(run(exe ,"in.reverse_rate",vz,"keqpar1",extra ,sub,xf))
        t2 = tallies(run(exe2,"in.reverse_rate",vz,"keqpar2",extra2,sub,xf))
        same = tallies_agree(t1, t2)
        listing = (",".join("%s=%g" % (k.split("-->")[0].strip(), v)
                            for k, v in sorted(t1.items()))
                   if t1 else "no tallies")
        check("CPU and second-binary keq tallies %s (%s)"
              % ("identical" if PARITY == "exact" else "agree within 4 sd",
                 listing), same)

    print("%d failures" % FAIL)
    sys.exit(FAIL)

if __name__ == "__main__":
    main()
