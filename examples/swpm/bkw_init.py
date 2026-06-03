#!/usr/bin/env python
"""
bkw_init.py - generate a Bobylev-Krook-Wu (BKW) distributed set of particles
              for a spatially homogeneous relaxation benchmark.

The BKW distribution is an exact, time-dependent solution of the spatially
homogeneous Boltzmann equation for Maxwell molecules.  It is a standard
benchmark for collision operators (DSMC, SWPM, ...): initialize the velocities
from the BKW distribution at a chosen non-equilibrium parameter beta0, run the
collision operator, and verify that the distribution relaxes toward the
Maxwellian as predicted (beta(t) = beta0 * exp(-t/6 * nu), with nu the Maxwell
molecule collision frequency).

The 3D BKW velocity distribution (isotropic, zero mean) is

    f(v) = A^(3/2) / pi^(3/2) * exp(-A v^2) * [1 + beta*(A v^2 - 3/2)]

with A = (1+beta) * m / (2 k T).  beta=0 recovers the Maxwellian at temperature
T.  The bracket is non-negative only for beta <= 2/3, so 0 < beta0 <= 0.65 is a
typical choice.

The output is written in the dump format read by SPARTA's `read_particles`
command (fields: id type x y z vx vy vz at ITEM: TIMESTEP 0).

Usage:
    bkw_init.py Np Lx Ly Lz outfile [beta0] [T] [mass]

    Np      number of particles to generate
    Lx,Ly,Lz  box dimensions (m); positions are uniform in the box
    outfile output dump file (read with: read_particles outfile 0)
    beta0   initial BKW parameter      (default 0.65)
    T       temperature in K           (default 273.15)
    mass    molecular mass in kg       (default 6.63e-26, argon)
"""

from __future__ import print_function
import sys
import numpy as np

KB = 1.380658e-23

def usage():
    print(__doc__)
    sys.exit(1)

if len(sys.argv) < 6:
    usage()

Np   = int(sys.argv[1])
Lx   = float(sys.argv[2])
Ly   = float(sys.argv[3])
Lz   = float(sys.argv[4])
outfile = sys.argv[5]
beta0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.65
T     = float(sys.argv[7]) if len(sys.argv) > 7 else 273.15
mass  = float(sys.argv[8]) if len(sys.argv) > 8 else 6.63e-26

if not (0.0 < beta0 <= 2.0/3.0):
    print("ERROR: beta0 must be in (0, 2/3] for f(v) >= 0")
    sys.exit(1)

kTm = KB * T / mass                 # thermal speed squared scale
A = (1.0 + beta0) / (2.0 * kTm)     # exponent coefficient

def fBKW(v):
    """Unnormalized 3D BKW speed density (per unit velocity volume)."""
    Av2 = A * v * v
    return pow(A, 1.5) * np.exp(-Av2) * (1.0 + beta0 * (Av2 - 1.5)) * pow(np.pi, -1.5)

# rejection sampling of the speed v against the radial density v^2 * f(v).
# sample candidate speeds with density proportional to v^2 (uniform in the
# velocity-space sphere) and accept with probability f(v)/fmax.

vmp = np.sqrt(2.0 * kTm)            # most probable speed of the Maxwellian
vmax = 4.0 * vmp                    # sampling cutoff (f is negligible beyond)

# bound the density on [0,vmax] for rejection sampling
vgrid = np.linspace(0.0, vmax, 2000)
fmax = 1.1 * np.max(fBKW(vgrid))

with open(outfile, "w") as fp:
    fp.write("ITEM: TIMESTEP\n0\n")
    fp.write("ITEM: NUMBER OF ATOMS\n%d\n" % Np)
    fp.write("ITEM: BOX BOUNDS pp pp pp\n")
    fp.write("0 %.6e\n0 %.6e\n0 %.6e\n" % (Lx, Ly, Lz))
    fp.write("ITEM: ATOMS id type x y z vx vy vz\n")

    ip = 0
    while ip < Np:
        v = vmax * pow(np.random.rand(), 1.0/3.0)   # density ~ v^2
        if fBKW(v) > fmax * np.random.rand():
            # isotropic direction
            cz = 2.0 * np.random.rand() - 1.0
            s = np.sqrt(1.0 - cz*cz)
            phi = 2.0 * np.pi * np.random.rand()
            vx = v * s * np.cos(phi)
            vy = v * s * np.sin(phi)
            vz = v * cz

            x = np.random.rand() * Lx
            y = np.random.rand() * Ly
            z = np.random.rand() * Lz

            # id type x y z vx vy vz
            fp.write("%d 1 %.6e %.6e %.6e %.6e %.6e %.6e\n" %
                     (ip+1, x, y, z, vx, vy, vz))
            ip += 1

print("wrote %d BKW particles (beta0=%g, T=%g K, mass=%g kg) to %s" %
      (Np, beta0, T, mass, outfile))
