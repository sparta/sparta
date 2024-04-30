from __future__ import print_function

# error message

def error(str=None):
  if str: print("ERROR:",str)
  else: print("Syntax: bkw_init.py N xl yl zl outfile")
  sys.exit()

# ----------------------------------------------------------------------
# main program

import sys,re
import numpy as np

Np = int(sys.argv[1])
xl = float(sys.argv[2])
yl = float(sys.argv[3])
zl = float(sys.argv[4])
outfile = sys.argv[5]

fp = open(outfile,"w")

# Header
print('ITEM: TIMESTEP',file = fp)
print('0',file = fp)
print('ITEM: NUMBER OF ATOMS',file = fp)
print(str(Np),file = fp)
print('ITEM: BOX BOUNDS pp pp pp',file = fp)
print('0 {0:.2e}'.format(xl),file=fp)
print('0 {0:.2e}'.format(yl),file=fp)
print('0 {0:.2e}'.format(zl),file=fp)
print('ITEM: ATOMS id type x y z vx vy vz',file = fp)

# from paper
'''
m = 1
kb = 1
T = 1

beta0 = 0.65
beta = beta0/(1+beta0)
A = m * (1+beta) / (2*kb*T)
vmp = np.sqrt(2.0*kb*T/m)
vmax = 2.5*vmp

fmax = 3.55;

ip = 0;
vout = np.zeros(3,)
while (ip < Np):
  iv = np.random.rand()*vmax;
  piv = pow(A,1.5)*np.exp(-A*iv*iv)*(1+beta*(A*iv*iv - 1.5))*pow(np.pi,-1.5)/fmax
  if(piv > 1):
    print(str(piv))
    print("over one")
  if(piv > np.random.rand()):
    # print to file
    x = np.random.rand()*xl
    y = np.random.rand()*yl
    z = np.random.rand()*zl

    cz = 2*np.random.rand()-1
    rrr = np.sqrt(1-cz*cz)
    theta = 2.*np.pi*np.random.rand()
    cx = rrr*np.cos(theta)
    cy = rrr*np.sin(theta)
    vx = iv*cx
    vy = iv*cy
    vz = iv*cz

    # id type x y z vx vy vz
    print('{0:1d} {1:2d} {2:.5e} {3:.5e} {4:.5e} {5:.5e} {6:.5e} {7:.5e}'.format(ip,1,x,y,z,vx,vy,vz),file=fp)

    ip = ip + 1

sys.exit()
'''

# BKW form:
# f0 * (f1*v^2) exp(-v^2/C0)

K = 0.6;
f0 = 1.0/pow(np.pi*K,1.5)/K;
f1 = (1.0-K)/K

fmax = 0.04;
ip = 0;
vmax = 4.0

cq = 1.0
cqsq = 1.0

while (ip < Np):
  iv = vmax*np.random.rand();
  piv = 5.0/9.0*(iv*iv/cqsq)*pow(np.sqrt(6*np.pi/5.0)*cq,-3.0)*np.exp(-(iv*iv)/(6/5*cqsq))
  if(piv > fmax):
    print("greater than 1")
  if(piv/fmax > np.random.rand()):
    # print to file
    x = np.random.rand()*xl
    y = np.random.rand()*yl
    z = np.random.rand()*zl

    theta = 2*np.pi*np.random.rand();
    phi = np.arccos(1.0-2.0*np.random.rand());
    ex = np.sin(phi)*np.cos(theta);
    ey = np.sin(phi)*np.sin(theta);
    ez = np.cos(phi);

    vx = iv*ex;
    vy = iv*ey;
    vz = iv*ez;
  

    # id type x y z vx vy vz
    print('{0:1d} {1:2d} {2:.5e} {3:.5e} {4:.5e} {5:.5e} {6:.5e} {7:.5e}'.format(ip,1,x,y,z,vx,vy,vz),file=fp)

    ip = ip + 1

sys.exit()









