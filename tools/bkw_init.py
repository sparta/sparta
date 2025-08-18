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

def fBKW(iv,beta):
    T = 273
    m = 6.64e-26
    k = 1.38e-23
    kTm = k*T/m
    
    A = (1+beta)/kTm*0.5;
    ivsq = iv*iv
    prob = pow(A,1.5)*np.exp(-A*ivsq)*(1+beta*(A*ivsq-1.5))*pow(np.pi,-1.5);
    return prob

# from paper
pmax = 3.5e-9
vmp = np.sqrt(2*1.38e-23*273/6.63e-26);
vmax = 3*vmp;
beta0 = 0.65

ip = 0;
vout = np.zeros(3,)
while (ip < Np):
  vv = vmax*pow(np.random.rand(),1/3)
  probvv = fBKW(vv,beta0);
  if(probvv > pmax*np.random.rand()):
    cz = 2*np.random.rand()-1
    rrr = np.sqrt(1-cz*cz)
    theta = 2*np.pi*np.random.rand()
    cx = rrr*np.cos(theta)
    cy = rrr*np.sin(theta)

    # print to file
    x = np.random.rand()*xl
    y = np.random.rand()*yl
    z = np.random.rand()*zl

    vx = vv*cx;
    vy = vv*cy;
    vz = vv*cz;

    # id type x y z vx vy vz
    print('{0:1d} {1:2d} {2:.5e} {3:.5e} {4:.5e} {5:.5e} {6:.5e} {7:.5e}'.format(ip,1,x,y,z,vx,vy,vz),file=fp)

    ip = ip + 1

sys.exit()










