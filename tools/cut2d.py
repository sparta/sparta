#!/usr/local/bin/python

# clip test for overlap of line segment PQ with cell
# return 1 if intersects cell, 0 if not

def cliptest(p,q,cell):
  if p[0] >= cell[1] and p[0] <= cell[2] and \
        p[1] >= cell[3] and p[1] <= cell[4]: return 1
  if q[0] >= cell[1] and q[0] <= cell[2] and \
        q[1] >= cell[3] and q[1] <= cell[4]: return 1
  c = [p[0],p[1]]
  d = [q[0],q[1]]
  if c[0] < cell[1] and d[0] < cell[1]: return 0
  if c[0] < cell[1] or d[0] < cell[1]:
    y = c[1] + (cell[1]-c[0])/(d[0]-c[0])*(d[1]-c[1])
    if c[0] < cell[1]: c[0],c[1] = cell[1],y
    else: d[0],d[1] = cell[1],y
  if c[0] > cell[2] and d[0] > cell[2]: return 0
  if c[0] > cell[2] or d[0] > cell[2]:
    y = c[1] + (cell[2]-c[0])/(d[0]-c[0])*(d[1]-c[1])
    if c[0] > cell[2]: c[0],c[1] = cell[2],y
    else: d[0],d[1] = cell[2],y
  if c[1] < cell[3] and d[1] < cell[3]: return 0
  if c[1] < cell[3] or d[1] < cell[3]:
    x = c[0] + (cell[3]-c[1])/(d[1]-c[1])*(d[0]-c[0])
    if c[1] < cell[3]: c[0],c[1] = x,cell[3]
    else: d[0],d[1] = x,cell[3]
  if c[1] > cell[4] and d[1] > cell[4]: return 0
  if c[1] > cell[4] or d[1] > cell[4]:
    x = c[0] + (cell[4]-c[1])/(d[1]-c[1])*(d[0]-c[0])
    if c[1] > cell[4]: c[0],c[1] = x,cell[4]
    else: d[0],d[1] = x,cell[4]
  return 1
