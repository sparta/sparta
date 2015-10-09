# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# sdata tool

oneline = "Read, create, manipulate SPARTA surf files"

docstr = """
s = sdata()			   create a surf data object
s = sdata(ID,"mem.surf")           read in one or more SPARTA surf files
s = sdata(ID,"mem.part.gz mem.surf")  can be gzipped
s = sdata(ID,"mem.*")		   wildcard expands to multiple files
s.read(ID,"mem.surf")		   read in one or more data files

  all surf data in files becomes one surf with ID
  surf files contain the following kinds of entries in SPARTA format
    points and lines (for 2d), points and triangles (for 3d)
  read() has same argument options as constructor

s.seed = 48379                     set random # seed (def = 12345)
s.circle(ID,x,y,r,n)               create a 2d circle with N lines and ID
s.rect(ID,x1,y1,x2,y2,nx,ny)       create a 2d rect, 2 corner pts, Nx,Ny segs 
s.tri(ID,x1,y1,x2,y2,x3,y3,n1,n2,n3)  create a 2d tri, 3 pts, N1,N2,N3 segs 
s.sphere(ID,x,y,z,r,n)		   create a 3d sphere with NxN sqs per face
s.box(ID,x1,y1,z1,x2,y2,z2,nx,ny,nz)  3d box, 2 corner pts, Nx,Ny,Nz per face
s.spikycircle(ID,x,y,rmin,rmax,n)  2d circle, N lines, Rmin <= rad <= Rmax
s.spikysphere(ID,x,y,z,rmin,rmax,n)  3d sphere, NxN sqs, Rmin <= rad <= Rmax

  tri should be ordered so that (0,0,1) x (pt2-pt1) = outward normal
  spikycircle is same as circle, with each of N pts at random Rmin < rad < Rmax
  spikysphere is same as sphere, with each surf pt at random Rmin < rad < Rmax

s.surf2d(ID,plist,llist)           create a custom 2d surf
s.surf3d(ID,plist,tlist)           create a custom 3d surf

  each point in plist is (x,y) for 2d or (x,y,z) for 3d
  each line in llist is (i,j) for C-style indices into plist
  each triangle in tlist is (i,j,k) for C-style indices into plist

s.center(ID,x,y,z)                 set center point of surf
s.trans(ID,dx,dy,dz)   	 	   translate surf and its center point
s.rotate(ID,theta,Rx,Ry,Rz)        rotate surf by theta around R vector
s.scale(ID,sx,sy,sz)		   scale a surf around center point
s.invert(ID) 	                   invert normal direction of surf

  default center for created surfs is the x,y,z or geometric center
  default center for read-in surf is center of bounding box of all points
  rotation and scaling of surf are relative to its center point

s.join(ID,id1,id2,...)		   combine id1,id2,etc into new surf with ID
s.delete(id1,id2,...)              delete one or more surfs
s.rename(ID,IDnew)                 rename a surf
s.copy(ID,IDnew) 	           create a new surf as copy of old surf

  join does not delete id1,id2,etc
  center for joined surf becomes center of bounding box of all points

s.select(id1,id2,...)              select one or more surfs
s.select()                         select all surfs
s.unselect(id1,id2,...)            unselect one or more surfs
s.unselect()                       unselect all surfs

  selection applies to write() and viz()
  surfs are selected by default when read or created
  
s.write("file")			   write all selected surfs to SPARTA file
s.write("file",id1,id2,...)	   write only listed & selected surfs to file

s.grid(xlo,xhi,ylo,yhi,ny,ny)      bounding box and Nx by Ny grid
s.grid(xlo,xhi,ylo,yhi,zlo,zhi,ny,ny,nz)   ditto for 3d
s.gridfile(xlo,xhi,ylo,yhi,file)   bbox and SPARTA-formatted parent grid file
s.gridfile(xlo,xhi,ylo,yhi,zlo,zhi,file)   ditto for 3d

  grid command superpose a grid, for viz only
  also changes bounding box to xyzlo to xyzhi

index,time,flag = s.iterator(0/1)          loop over single snapshot
time,box,atoms,bonds,tris,lines = s.viz(index)   return list of viz objects

  iterator() and viz() are compatible with equivalent dump calls
  iterator() called with arg = 0 first time, with arg = 1 on subsequent calls
    index = timestep index within dump object (only 0 for data file)
    time = timestep value (only 0 for data file)
    flag = -1 when iteration is done, 1 otherwise
  viz() returns info for selected objs for specified timestep index (must be 0)
    time = 0
    box = [xlo,ylo,zlo,xhi,yhi,zhi] = bounding box
    atoms = NULL
    bonds = NULL
    tris = id,type,x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz for each tri as 2d array
      NULL if triangles do not exist
    lines = id,type,x1,y1,z1,x2,y2,z2 for each line as 2d array
      NULL if lines do not exist
    types are assigned to each surf in ascending order
"""

# History
#   10/12, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   dim = 2 or 3, all surfs must be the same
#   ids = dictionary of IDs that points to surfs index
#   surfs = list of surfs

# Imports and external programs

import sys,glob
from os import popen
from math import pi,cos,sin,sqrt
from copy import deepcopy

try: from DEFAULTS import PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

BIG = 1.0e20

# Class definition

class sdata:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.dim = 0
    self.nselect = 1
    self.seed = 12345
    self.gridflag = 0
    self.ids = {}
    self.surfs = []
    self.plist = []
    self.llist = []
    self.tlist = []

    if len(list) == 1: raise StandardError,"surf ID and surf file required"
    if len(list) > 1: self.read(list[0],*list[1:])

  # --------------------------------------------------------------------

  def read(self,id,*list):

    # flist = list of all surf file names

    words = list[0].split()
    flist = []
    for word in words: flist += glob.glob(word)
    if len(flist) == 0 and len(list) == 1:
      raise StandardError,"no surf file specified"

    # read all surf file as one surf with ID
    
    points = []
    lines = []
    triangles = []

    for file in flist:
      npoints_prev = len(points)

      # test for gzipped file

      if file[-3:] == ".gz":
        f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r')
      else: f = open(file)

      # read file

      pflag = lflag = tflag = 0
      npoints = nlines = ntriangles = 0

      line = f.readline()
      while 1:
        line = f.readline()
        if not line: break
        line = line.strip()
        if not line: continue
        if '#' in line: line = line[:line.index['#']]
        if not line: continue
        if "points" in line: npoints = int(line.split()[0])
        elif "lines" in line: nlines = int(line.split()[0])
        elif "triangles" in line: ntriangles = int(line.split()[0])

        elif "Points" in line:
          if npoints == 0:
            raise StandardError, "invalid SPARTA surf file"
          if nlines == 0 and ntriangles == 0:
            raise StandardError, "invalid SPARTA surf file"
          if nlines and ntriangles:
            raise StandardError, "invalid SPARTA surf file"
          if nlines:
            if self.dim == 3: 
              raise StandardError, "cannot have both 2d/3d surfs"
            self.dim = 2
          if ntriangles:
            if self.dim == 2: 
              raise StandardError, "cannot have both 2d/3d surfs"
            self.dim = 3

          line = f.readline()
          for i in xrange(npoints):
            list = f.readline().split()
            pt = [float(value) for value in list[1:]]
            if self.dim == 2: pt.append(0.0)
            points.append(pt)
          
        elif "Lines" in line:
          if npoints == 0 or nlines == 0:
            raise StandardError, "invalid SPARTA surf file"

          line = f.readline()
          for i in xrange(nlines):
            list = f.readline().split()
            lines.append([int(value)-npoints_prev-1 for value in list[1:]])

        elif "Triangles" in line:
          if npoints == 0 or ntriangles == 0:
            raise StandardError, "invalid SPARTA surf file"

          line = f.readline()
          for i in xrange(ntriangles):
            list = f.readline().split()
            triangles.append([int(value)-npoints_prev-1 for value in list[1:]])

      f.close()

    if self.dim == 2: print "surf %s with %d points, %d lines" % \
        (id,npoints,nlines)
    if self.dim == 3: print "surf %s with %d points, %d triangles" % \
        (id,npoints,ntriangles)

    # create surf
            
    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.triangles = triangles
    box = bbox(points)
    surf.center = [0.5*(box[0]+box[3]),0.5*(box[1]+box[4]),0.5*(box[2]+box[5])]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)
    
  # --------------------------------------------------------------------
  # create a circle

  def circle(self,id,x,y,r,n):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 3: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 2
    
    points = []
    lines = []
    for i in range(n):
      theta = i*2.0*pi/n 
      points.append([x+r*cos(theta),y+r*sin(theta),0.0])
      
    for i in range(n-1):
      lines.append([i,i+1])
    lines.append([n-1,0])

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.center = [x,y,0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a rectangle

  def rect(self,id,x0,y0,x1,y1,nx,ny):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 3: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 2
      
    points = []
    lines = []
    for i in range(ny):
      points.append([x0,y0+i*(y1-y0)/ny,0.0])
    for i in range(nx):
      points.append([x0+i*(x1-x0)/nx,y1,0.0])
    for i in range(ny):
      points.append([x1,y1+i*(y0-y1)/ny,0.0])
    for i in range(nx):
      points.append([x1+i*(x0-x1)/nx,y0,0.0])

    for i in range(2*(nx+ny)-1):
      lines.append([i,i+1])
    lines.append([2*(nx+ny)-1,0])

    print points
    print lines
    
    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.center = [0.5*(x0+x1),0.5*(y0+y1),0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a triangle

  def tri(self,id,x0,y0,x1,y1,x2,y2,n1,n2,n3):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 3: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 2
      
    points = []
    lines = []
    for i in range(n1):
      points.append([x0+i*(x1-x0)/n1,y0+i*(y1-y0)/n1,0.0])
    for i in range(n2):
      points.append([x1+i*(x2-x1)/n2,y1+i*(y2-y1)/n2,0.0])
    for i in range(n3):
      points.append([x2+i*(x0-x2)/n3,y2+i*(y0-y2)/n3,0.0])

    for i in range(n1*n2*n3-1):
      lines.append([i,i+1])
    lines.append([n1*n2*n3-1,0])
    
    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.center = [(x0+x1+x2)/3.0,(y0+y1+y2)/3.0,0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a sphere

  def sphere(self,id,x,y,z,r,n):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    if self.dim == 2: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 3

    pts,triangles = box_triangulate(n,n,n)
    points = []
    for pt in pts:
      ptnew = [pt[0]-0.5,pt[1]-0.5,pt[2]-0.5]
      normalize(ptnew)
      ptnew[0] = x + r*ptnew[0]
      ptnew[1] = y + r*ptnew[1]
      ptnew[2] = z + r*ptnew[2]
      points.append(ptnew)

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.triangles = triangles
    surf.center = [x,y,z]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a 3d box

  def box(self,id,x0,y0,z0,x1,y1,z1,nx,ny,nz):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 2: raise StandardError, "cannot have both 2d/3d surfs"
    if nx <= 0 or ny <= 0 or nz <= 0:
      raise StandardError, "invalid box nx,ny,nz values"
    self.dim = 3

    pts,triangles = box_triangulate(nx,ny,nz)
    points = []
    for pt in pts:
      xnew = x0 + pt[0]*(x1-x0)
      ynew = y0 + pt[1]*(y1-y0)
      znew = z0 + pt[2]*(z1-z0)
      points.append([xnew,ynew,znew])

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.triangles = triangles
    surf.center = [0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1)]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a spiky circle

  def spikycircle(self,id,x,y,rmin,rmax,n):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 3: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 2
    
    points = []
    lines = []
    for i in range(n):
      theta = i*2.0*pi/n
      r = rmin + self.random()*(rmax-rmin)
      points.append([x+r*cos(theta),y+r*sin(theta),0.0])
      
    for i in range(n-1):
      lines.append([i,i+1])
    lines.append([n-1,0])

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.center = [x,y,0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a spiky sphere

  def spikysphere(self,id,x,y,z,rmin,rmax,n):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    if self.dim == 2: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 3

    pts,triangles = box_triangulate(n,n,n)
    points = []
    for pt in pts:
      ptnew = [pt[0]-0.5,pt[1]-0.5,pt[2]-0.5]
      normalize(ptnew)
      r = rmin + self.random()*(rmax-rmin)
      ptnew[0] = x + r*ptnew[0]
      ptnew[1] = y + r*ptnew[1]
      ptnew[2] = z + r*ptnew[2]
      points.append(ptnew)

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.triangles = triangles
    surf.center = [x,y,z]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a custom 2d surf from list of points and lines

  def surf2d(self,id,plist,llist):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 3: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 2
    
    surf = Surface()
    surf.select = 1
    surf.points = plist
    surf.lines = llist
    surf.center = [0.0,0.0,0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # create a custom 3d surf from list of points and lines

  def surf3d(self,id,plist,tlist):
    if self.ids.has_key(id): raise StandardError,"ID %s is already in use" % id
    if self.dim == 2: raise StandardError, "cannot have both 2d/3d surfs"
    self.dim = 3
    
    surf = Surface()
    surf.select = 1
    surf.points = plist
    surf.triangles = tlist
    surf.center = [0.0,0.0,0.0]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # set center pt of a surf

  def center(self,id,x,y,z):
    if not self.ids.has_key(id):
      raise StandardError,"ID %s is not defined" % id
    if self.dim == 2 and z != 0.0:
      raise StandardError,"z center of 2d surf must be 0.0"
    surf = self.surfs[self.ids[id]]
    surf.center = [x,y,z]

  # --------------------------------------------------------------------
  # translate a surf by dx,dy,dz displacement
  # add displacement to its vertices and center pt

  def trans(self,id,dx,dy,dz):
    if not self.ids.has_key(id):
      raise StandardError,"ID %s is not defined" % id
    if self.dim == 2 and dz != 0.0:
      raise StandardError,"dz translation of 2d surf must be 0.0"

    surf = self.surfs[self.ids[id]]
    surf.xc += dx
    surf.yc += dy
    surf.zc += dz

    for i,pt in enumerate(surf.points):
      pt[0] += dx
      pt[1] += dy
      pt[2] += dz

  # --------------------------------------------------------------------
  # rotate a surf by theta around (Rx,Ry,Rz) and center pt
  # convert (Rx,Ry,Rz) to unit vector
  # quat = (cos(theta/2),Rx*sin(theta/2),Ry*sin(theta/2),Rz*sin(theta/2))
  # rotation matrix P =
  #   2 * ( (q0^2+q1^2-q2^2-q3^2)/2   q1q2-q0q3    q1q3+q0q2 )
  #       ( q1q2+q0q3  (q0^2-q1^2+q2^2-q3^2)/2     q2q3-q0q1 )
  #       ( q1q3-q0q2     q2q3+q0q1   q0^2-q1^2-q2^2+q3^2)/2 )	
  # for each point x: xnew = P (x - center) + center

  def rotate(self,id,theta,rx,ry,rz):
    if not self.ids.has_key(id):
      raise StandardError,"ID %s is not defined" % id
    if self.dim == 2 and (rx != 0.0 or ry != 0.0):
      raise StandardError,"rx,ry rotation of 2d surf must be 0.0"

    r = [rx,ry,rz]
    normalize(r)
    angle = pi * theta/2/180.0
    q = (cos(angle),rx*sin(angle),ry*sin(angle),rz*sin(angle))
    p00 = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
    p01 = 2 * (q[1]*q[2] - q[0]*q[3])
    p02 = 2 * (q[1]*q[3] + q[0]*q[2])
    p10 = 2 * (q[1]*q[2] + q[0]*q[3])
    p11 = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
    p12 = 2 * (q[2]*q[3] - q[0]*q[1])
    p20 = 2 * (q[1]*q[3] - q[0]*q[2])
    p21 = 2 * (q[2]*q[3] + q[0]*q[1])
    p22 = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]

    surf = self.surfs[self.ids[id]]
    center = surf.center

    for i,pt in enumerate(surf.points):
      ptnew = 3*[0]
      xc = [pt[0]-center[0],pt[1]-center[1],pt[2]-center[2]]
      ptnew[0] = p00*xc[0] + p01*xc[1] + p02*xc[2] + center[0]
      ptnew[1] = p10*xc[0] + p11*xc[1] + p12*xc[2] + center[1]
      ptnew[2] = p20*xc[0] + p21*xc[1] + p22*xc[2] + center[2]
      surf.points[i] = ptnew
    
  # --------------------------------------------------------------------
  # scale an object by sx,sy,sz factors
  # scale its vertices relative to center pt

  def scale(self,id,sx,sy,sz):
    if not self.ids.has_key(id):
      raise StandardError,"ID %s is not defined" % id
    if self.dim == 2 and sz != 1.0:
      raise StandardError,"sz scale of 2d surf must be 1.0"

    surf = self.surfs[self.ids[id]]
    center = surf.center
    for pt in surf.points:
      pt[0] = center[0] + sx*(pt[0]-center[0])
      pt[1] = center[1] + sy*(pt[1]-center[1])
      pt[2] = center[2] + sz*(pt[2]-center[2])
    
  # --------------------------------------------------------------------
  # invert direction of surf normals by swapping line or triangle indices
  
  def invert(self,id):
    if not self.ids.has_key(id):
      raise StandardError,"ID %s is not defined" % id

    surf = self.surfs[self.ids[id]]
    if self.dim == 2:
      for i,line in enumerate(surf.lines):
        line = [line[1],line[0]]
        surf.lines[i] = line
    if self.dim == 3:
      for i,tri in enumerate(surf.triangles):
        tri = [tri[0],tri[2],tri[1]]
        surf.triangles[i] = tri
  
  # --------------------------------------------------------------------
  # join surfs in list to form a new surf
  
  def join(self,id,*list):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    if len(list) == 0:
      raise StandardError,"list of surfs to join is empty"

    points = []
    lines = []
    triangles = []

    for id in list:
      surf = self.surfs[self.ids[id]]
      npoints_prev = len(points)
      points += surf.points
      if self.dim == 2:
        for line in self.lines:
          line[0] += npoints_prev
          line[1] += npoints_prev
          lines.append(line)
      if self.dim == 3:
        for tri in self.triangles:
          tri[0] += npoints_prev
          tri[1] += npoints_prev
          tri[2] += npoints_prev
          triangles.append(tri)

    surf = Surface()
    surf.select = 1
    surf.points = points
    surf.lines = lines
    surf.triangles = triangles
    box = bbox(points)
    surf.center = [0.5*(box[0]+box[3]),0.5*(box[1]+box[4]),0.5*(box[2]+box[5])]
    self.ids[id] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # delete each surf in list
  # reset values in ids since some indices are decremented
  
  def delete(self,*list):
    for id in list:
      i = self.ids[id]
      del self.ids[id]
      del self.surfs[i]
      for key in self.ids.keys():
        j = self.ids[key]
        if j > i: self.ids[key] = j-1
        
  # --------------------------------------------------------------------
  # rename the ID of an object
  # check that new ID doesn't already exist

  def rename(self,idold,idnew):
    if self.ids.has_key(idnew):
      raise StandardError,"ID %s is already in use" % idnew
    i = self.ids[idold]
    self.ids[idnew] = i
    del self.ids[idold]
    
  # --------------------------------------------------------------------
  # create a deep copy of an object and assign it a new ID
  # check that new name doesn't already exist

  def copy(self,idold,idnew):
    if self.ids.has_key(idnew):
      raise StandardError,"ID %s is already in use" % idnew
    surf = deepcopy(self.surfs[self.ids[idold]])
    surf.select = 1
    self.ids[idnew] = len(self.surfs)
    self.surfs.append(surf)

  # --------------------------------------------------------------------
  # set selection flag for each surf in list
  # if list is empty, select all
  
  def select(self,*list):
    if len(list) == 0: list = self.ids.keys()
    for id in list:
      surf = self.surfs[self.ids[id]]
      surf.select = 1

  # --------------------------------------------------------------------
  # unset selection flag for each surf in list
  # if list is empty, unselect all
  
  def unselect(self,*list):
    if len(list) == 0: list = self.ids.keys()
    for id in list:
      surf = self.surfs[self.ids[id]]
      surf.select = 0

  # --------------------------------------------------------------------
  # write out surfs in list to surf file
  # if list is empty, write all surfs
  
  def write(self,file,*list):
    if not len(list): vlist = range(len(self.surfs))
    else:
      vlist = []
      for id in list: vlist.append(self.ids[id])

    points = []
    lines = []
    triangles = []
    
    for index in vlist:
      obj = self.surfs[index]
      if not obj.select: continue
      npoints_prev = len(points)
      points += obj.points
      if self.dim == 2:
        for line in obj.lines:
          line[0] += npoints_prev
          line[1] += npoints_prev
          lines.append(line)
      if self.dim == 3:
        for tri in obj.triangles:
          tri[0] += npoints_prev
          tri[1] += npoints_prev
          tri[2] += npoints_prev
          triangles.append(tri)
      
    fp = open(file,'w')
    print >>fp,"surf file from Pizza.py"
    print >>fp
    print >>fp,len(points),"points"
    if self.dim == 2: print >>fp,len(lines),"lines"
    if self.dim == 3: print >>fp,len(triangles),"triangles"
    print >>fp
    print >>fp,"Points\n"
    if self.dim == 2:
      for i,point in enumerate(points):
        print >>fp,i+1,point[0],point[1]
    if self.dim == 3:
      for i,point in enumerate(points):
        print >>fp,i+1,point[0],point[1],point[2]
    print >>fp
    if self.dim == 2:
      print >>fp,"Lines\n"
      for i,line in enumerate(lines):
        print >>fp,i+1,line[0]+1,line[1]+1
    if self.dim == 3:
      print >>fp,"Triangles\n"
      for i,tri in enumerate(triangles):
        print >>fp,i+1,tri[0]+1,tri[1]+1,tri[2]+1
    
    fp.close()

  # --------------------------------------------------------------------
  # overlay a top-level grid over surf for viz only

  def grid(self,*args):
    if self.dim == 0: raise StandardError, "dimension must be defined for grid"
    if self.dim == 2 and len(args) != 6:
      raise StandardError, "bad arguments for sdata.grid()"
    if self.dim == 3 and len(args) != 9:
      raise StandardError, "bad arguments for sdata.grid()"

    self.gridflag = 1
    self.idparents = ["0"]
    self.parents = {}
    if self.dim == 2:
      self.parents["0"] = [(args[0],args[1],args[2],args[3]),
                           (args[4],args[5],1)]
    if self.dim == 3:
      self.parents["0"] = [(args[0],args[1],args[2],args[3],args[4],args[5]),
                           (args[6],args[7],args[8])]
        
  # --------------------------------------------------------------------
  # overlay a hierarchical grid over surf for viz only
  # grid comes from SPARTA parent grid file

  def gridfile(self,*args):
    if self.dim == 0: raise StandardError, \
          "dimension must be defined for gridfile"
    if self.dim == 2 and len(args) != 5:
      raise StandardError, "bad arguments for sdata.gridfile()"
    if self.dim == 3 and len(args) != 7:
      raise StandardError, "bad arguments for sdata.gridfile()"

    # read parent file, parent entries should start on line 7

    lines = open(args[4],"r").readlines()
    lines = lines[6:]
    
    self.gridflag = 1
    self.idparents = []
    self.parents = {}
    if self.dim == 2:
      for line in lines:
        words = line.split()
        self.idparents.append(words[1])
        if words[1] == "0":
          self.parents[words[1]] = [(args[0],args[1],args[2],args[3]),
                                    (int(words[2]),int(words[3]),1)]
        else:
          self.parents[words[1]] = [(),(int(words[2]),int(words[3]),1)]
    if self.dim == 3:
      for line in lines:
        words = line.split()
        self.idparents.append(words[1])
        if words[1] == "0":
          self.parents[words[1]] = [(args[0],args[1],args[2],args[3],
                                     args[4],args[5]),
                                    (int(words[2]),int(words[3]),int(words[4]))]
        else:
          self.parents[words[1]] = [(args[0],args[1],args[2],args[3],
                                     args[4],args[5]),
                                    (int(words[2]),int(words[3]),int(words[4]))]

  # --------------------------------------------------------------------
  # iterator called from other tools

  def iterator(self,flag):
    if flag == 0: return 0,0,1
    return 0,0,-1

  # --------------------------------------------------------------------
  # return list of atoms and triangles to viz for cdata object

  def viz(self,isnap):
    if isnap:
      raise StandardError, "cannot call cdata.viz() with isnap != 0"

    # no atoms or bonds
    
    atoms = []
    bonds = []

    # create triangle list from sum of all surfaces and regions
    # id = running count
    # type = type of set of tris

    id = itype = 0
    tris = []
    if self.dim == 3:
      for surf in self.surfs:
        if not surf.select: continue
        itype += 1
        points = surf.points
        triangles = surf.triangles
        for tri in triangles:
          id += 1
          pt1 = points[tri[0]]
          pt2 = points[tri[1]]
          pt3 = points[tri[2]]
          n = normal(pt1,pt2,pt3)
          tris.append([id,itype] + pt1 + pt2 + pt3 + n)

    # create line list from sum of all line objects

    id = itype = 0
    lines = []
    if self.dim == 2:
      itype = 1
      for surf in self.surfs:
        if not surf.select: continue
        points = surf.points
        segments = surf.lines
        for segment in segments:
          id += 1
          lines.append([id,itype] + points[segment[0]] + points[segment[1]])

    # add overlayed grid with new type
    # for each parent, draw its Nx by Ny by Nz sub-lines in 2d or 3d
    # use box stored with parents or compute box from ID
    # grandparent box will always exist due to loop over idparents
    #   which requires a parent cell's grandparent to be earlier in list
          
    if self.gridflag and self.dim == 2:
      for idparent in self.idparents:
        box,subgrid = self.parents[idparent]
        levels = idparent.split('-')
        nlevel = len(levels)
        if idparent == "0": nlevel = 0
        
        if not box:
          if nlevel <= 1:
            idchild = int(idparent)
            idgrandparent = "0"
          else:
            idchild = int(levels[-1])
            idgrandparent = "-".join(levels[:-1])
          gbox,gsubgrid = self.parents[idgrandparent]
          # compute parent box from grandparent box and store in parents hash
          xlo = gbox[0]; xhi = gbox[1]; 
          ylo = gbox[2]; yhi = gbox[3]; 
          nx = gsubgrid[0]; ny = gsubgrid[1]
          ix = (idchild-1) % nx
          iy = (idchild-1) / nx
          box = (xlo + float(ix)*(xhi-xlo)/nx, xlo + float(ix+1)*(xhi-xlo)/nx,
                 ylo + float(iy)*(yhi-ylo)/ny, ylo + float(iy+1)*(yhi-ylo)/ny)
          self.parents[idparent] = [box,subgrid]
          
        xlo = box[0]; xhi = box[1]; 
        ylo = box[2]; yhi = box[3]; 
        nx = subgrid[0]; ny = subgrid[1]

        for i in range(nx+1):
          x = xlo + float(i)*(xhi-xlo)/nx
          lines.append([id+i,itype+nlevel+1] + [x,ylo,0.0,x,yhi,0.0])
        id += nx+1
        for i in range(ny+1):
          y = ylo + float(i)*(yhi-ylo)/ny
          lines.append([id+i,itype+nlevel+1] + [xlo,y,0.0,xhi,y,0.0])
        id += ny+1

    if self.gridflag and self.dim == 3:
      for idparent in self.idparents:
        box,subgrid = self.parents[idparent]
        levels = idparent.split('-')
        nlevel = len(levels)
        if idparent == "0": nlevel = 0

        if not box:
          if nlevel <= 1:
            idchild = int(idparent)
            idgrandparent = "0"
          else:
            idchild = int(levels[-1])
            idgrandparent = "-".join(levels[:-1])
          gbox,gsubgrid = self.parents[idgrandparent]
          # compute parent box from grandparent box and store in parents hash
          xlo = gbox[0]; xhi = gbox[1]; 
          ylo = gbox[2]; yhi = gbox[3]; 
          zlo = gbox[4]; zhi = gbox[5]; 
          nx = gsubgrid[0]; ny = gsubgrid[1]; nz = gsubgrid[2]
          ix = (idchild-1) % nx
          iy = ((idchild-1)/nx) % ny
          iz = (idchild-1) / (nx*ny)
          box = (xlo + float(ix)*(xhi-xlo)/nx, xlo + float(ix+1)*(xhi-xlo)/nx,
                 ylo + float(iy)*(yhi-ylo)/ny, ylo + float(iy+1)*(yhi-ylo)/ny,
                 zlo + float(iz)*(zhi-zlo)/nz, zlo + float(iz+1)*(zhi-zlo)/nz)
          self.parents[idparent] = [box,subgrid]
          
        xlo = box[0]; xhi = box[1]; 
        ylo = box[2]; yhi = box[3]; 
        zlo = box[4]; zhi = box[5]; 
        nx = subgrid[0]; ny = subgrid[1]; nz = subgrid[2]

        for i in range(ny+1):
          y = ylo + float(i)*(yhi-ylo)/ny
          for j in range(nz+1):
            z = zlo + float(j)*(zhi-zlo)/nz
            lines.append([id+i,itype+nlevel+1] + [xlo,y,z,xhi,y,z])
        id += (ny+1)*(nz+1)
        for i in range(nx+1):
          x = xlo + float(i)*(xhi-xlo)/nx
          for j in range(nz+1):
            z = zlo + float(j)*(zhi-zlo)/nz
            lines.append([id+i,itype+nlevel+1] + [x,ylo,z,x,yhi,z])
        id += (nx+1)*(nz+1)
        for i in range(nx+1):
          x = xlo + float(i)*(xhi-xlo)/nx
          for j in range(ny+1):
            y = ylo + float(j)*(yhi-ylo)/ny
            lines.append([id+i,itype+nlevel+1] + [x,y,zlo,x,y,zhi])
        id += (nx+1)*(ny+1)

    return 0,self.bbox(),atoms,bonds,tris,lines

  # --------------------------------------------------------------------
  # time query from other tools

  def findtime(self,n):
    if n == 0: return 0
    raise StandardError, "no step %d exists" % (n)

  # --------------------------------------------------------------------
  # return box size

  def maxbox(self):
    return self.bbox()

  # --------------------------------------------------------------------
  # return box that bounds all selected objects
  # use grid box if defined
  
  def bbox(self):
    xlo = ylo = zlo = BIG
    xhi = yhi = zhi = -BIG
    for surf in self.surfs:
      if not surf.select: continue
      box = bbox(surf.points)
      xlo = min(xlo,box[0])
      ylo = min(ylo,box[1])
      zlo = min(ylo,box[2])
      xhi = max(xhi,box[3])
      yhi = max(yhi,box[4])
      zhi = max(zhi,box[5])
      
    if self.gridflag:
      root = self.parents["0"]
      box = root[0]
      if self.dim == 2:
        xlo = box[0]; xhi = box[1]
        ylo = box[2]; yhi = box[3]
      if self.dim == 3:
        xlo = box[0]; xhi = box[1]
        ylo = box[2]; yhi = box[3]
        zlo = box[4]; zhi = box[5]

    return (xlo,ylo,zlo,xhi,yhi,zhi)

  # --------------------------------------------------------------------

  def random(self):
    k = self.seed/IQ
    self.seed = IA*(self.seed-k*IQ) - IR*k
    if self.seed < 0:
      self.seed += IM
    return AM*self.seed

# --------------------------------------------------------------------
# random number generator class

IM = 2147483647
AM = 1.0/IM
IA = 16807
IQ = 127773
IR = 2836

# --------------------------------------------------------------------
# return c = a x b

def cross(a,b):
  c = 3*[0]
  c[0] = a[1]*b[2] - a[2]*b[1]
  c[1] = a[2]*b[0] - a[0]*b[2]
  c[2] = a[0]*b[1] - a[1]*b[0]
  return c

# --------------------------------------------------------------------
# normalize vector a to unit length

def normalize(a):
  length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
  if length == 0.0: return
  a[0] /= length
  a[1] /= length
  a[2] /= length

# --------------------------------------------------------------------
# return bounding box of points

def bbox(points):
  xlo = ylo = zlo = BIG
  xhi = yhi = zhi = -BIG
  for pt in points:
    xlo = min(xlo,pt[0])
    ylo = min(ylo,pt[1])
    zlo = min(zlo,pt[2])
    xhi = max(xhi,pt[0])
    yhi = max(yhi,pt[1])
    zhi = max(zhi,pt[2])
  return (xlo,ylo,zlo,xhi,yhi,zhi)

# --------------------------------------------------------------------
# add a vertex v to vertices list unless already exists in vdict dictionary
# return index of where v is in vertices list

def vertex(v,vertices,vdict):
  if vdict.has_key(v): return vdict[v]
  n = len(vertices)
  vertices.append(v)
  vdict[v] = n
  return n

# --------------------------------------------------------------------
# compute normal for a triangle with 3 vertices

def normal(x,y,z):
 v1 = 3*[0]
 v1[0] = y[0] - x[0]
 v1[1] = y[1] - x[1]
 v1[2] = y[2] - x[2]

 v2 = 3*[0]
 v2[0] = z[0] - y[0]
 v2[1] = z[1] - y[1]
 v2[2] = z[2] - y[2]

 n = 3*[0]
 n[0] = v1[1]*v2[2] - v1[2]*v2[1]
 n[1] = v1[2]*v2[0] - v1[0]*v2[2]
 n[2] = v1[0]*v2[1] - v1[1]*v2[0]

 length = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
 n[0] /= length
 n[1] /= length
 n[2] /= length

 return n

# --------------------------------------------------------------------
# triangulate a unit box from (0,0,0) to (1,1,1) with spacings q1,q2,q3
# return lists of vertices and triangles
# insure right-hand rule for each tri points OUT of the box

def box_triangulate(q1,q2,q3):
  if q1: dx = 1.0 / q1
  if q2: dy = 1.0 / q2
  if q3: dz = 1.0 / q3
  vdict = {}
  vertices = []
  triangles = []
  for j in xrange(q2):
    for k in xrange(q3):
      v1 = (0, j*dy,     k*dz)
      v2 = (0, (j+1)*dy, k*dz)
      v3 = (0, (j+1)*dy, (k+1)*dz)
      v4 = (0, j*dy,     (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv3,iv2])
      triangles.append([iv1,iv4,iv3])
      v1 = (1, j*dy,     k*dz)
      v2 = (1, (j+1)*dy, k*dz)
      v3 = (1, (j+1)*dy, (k+1)*dz)
      v4 = (1, j*dy,     (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv2,iv3])
      triangles.append([iv1,iv3,iv4])
  for i in xrange(q1):
    for k in xrange(q3):
      v1 = (i*dx,     0, k*dz)
      v2 = ((i+1)*dx, 0, k*dz)
      v3 = ((i+1)*dx, 0, (k+1)*dz)
      v4 = (i*dx,     0, (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv2,iv3])
      triangles.append([iv1,iv3,iv4])
      v1 = (i*dx,     1, k*dz)
      v2 = ((i+1)*dx, 1, k*dz)
      v3 = ((i+1)*dx, 1, (k+1)*dz)
      v4 = (i*dx,     1, (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv3,iv2])
      triangles.append([iv1,iv4,iv3])
  for i in xrange(q1):
    for j in xrange(q2):
      v1 = (i*dx,     j*dy,     0)
      v2 = ((i+1)*dx, j*dy,     0)
      v3 = ((i+1)*dx, (j+1)*dy, 0)
      v4 = (i*dx,     (j+1)*dy, 0)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv3,iv2])
      triangles.append([iv1,iv4,iv3])
      v1 = (i*dx,     j*dy,     1)
      v2 = ((i+1)*dx, j*dy,     1)
      v3 = ((i+1)*dx, (j+1)*dy, 1)
      v4 = (i*dx,     (j+1)*dy, 1)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1,iv2,iv3])
      triangles.append([iv1,iv3,iv4])
  return vertices,triangles

# Surface class

class Surface:
  def __init__(self):
    pass
