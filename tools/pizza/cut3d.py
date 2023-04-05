#!/usr/local/bin/python

# clip triangle p0,p1,p2 against cell via Sutherland-Hodgman algorithm
# clips tri path against each of 6 grid face planes
# return sequence of points in path = clipped convex polygon
# return [] if no intersection with cell
# no need to remove duplicate points,
#   caller just uses result for yes/no intersection

def clip(cell,p0,p1,p2):
  path = [p0,p1,p2]
  iface = 0
  for dim in range(3):
    iface += 1
    value = cell[iface]

    newpath = []
    s = path[-1]
    for i,e in enumerate(path):
      if e[dim] >= value:
        if s[dim] < value: newpath.append(between(s,e,dim,value))
        newpath.append(e)
      elif s[dim] >= value: newpath.append(between(e,s,dim,value))
      s = e
    path = newpath
    if not path: return []
    
    iface += 1
    value = cell[iface]

    newpath = []
    s = path[-1]
    for i,e in enumerate(path):
      if e[dim] <= value:
        if s[dim] > value: newpath.append(between(s,e,dim,value))
        newpath.append(e)
      elif s[dim] <= value: newpath.append(between(e,s,dim,value))
      s = e
    path = newpath
    if not path: return []

  return path

# return intersection pt of line segment a,b in dim with coord value
# guaranteed to intersect by caller
  
def between(a,b,dim,value):
  if dim == 0:
    y = a[1] + (value-a[dim])/(b[dim]-a[dim]) * (b[1]-a[1])
    z = a[2] + (value-a[dim])/(b[dim]-a[dim]) * (b[2]-a[2])
    pt = (value,y,z)
  elif dim == 1:
    x = a[0] + (value-a[dim])/(b[dim]-a[dim]) * (b[0]-a[0])
    z = a[2] + (value-a[dim])/(b[dim]-a[dim]) * (b[2]-a[2])
    pt = (x,value,z)
  elif dim == 2:
    x = a[0] + (value-a[dim])/(b[dim]-a[dim]) * (b[0]-a[0])
    y = a[1] + (value-a[dim])/(b[dim]-a[dim]) * (b[1]-a[1])
    pt = (x,y,value)
  return pt
