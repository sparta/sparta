# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# olog tool

oneline = "Read other log files (ChemCell, SPPARKS, SPARTA) and extract time-series data"

docstr = """
o = olog("file1")                    read in one or more log files
o = olog("log1 log2.gz")             can be gzipped
o = olog("file*")                    wildcard expands to multiple files
o = olog("log.spparks","Time")       2nd arg = start string for time section
o = olog("log.cell","",0)            3rd arg = average all runs

  incomplete and duplicate thermo entries are deleted
  if specify 2nd arg, it delimits a time section
  no 2nd arg or empty string, use default = "Step"
  if specify any 3rd arg, average all runs, assume all start at time 0
  
nvec = o.nvec                        # of vectors of thermo info
nlen = o.nlen                        length of each vectors
names = o.names                      list of vector names
a,b,... = o.get("A","B",...)         return one or more vectors of values
o.write("file.txt")	 	     write all vectors to a file
o.write("file.txt","A","B",...)      write listed vectors to a file

  get and write allow abbreviated (uniquely) vector names
"""

# History
#   1/06, Steve Plimpton (SNL): original version
#   2/09, modified to allow different firststr for different log files

# ToDo list

# Variables
#   nvec = # of vectors
#   nlen = length of each vector
#   names = list of vector names
#   ptr = dictionary, key = name, value = index into data for which column
#   data[i][j] = 2d array of floats, i = 0 to # of entries, j = 0 to nvecs-1
#   firststr = string that begins a time-series section in log file

# Imports and external programs

import sys, re, glob
from os import popen

try: tmp = PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class olog:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.nvec = 0
    self.names = []
    self.ptr = {}
    self.data = []
    self.firststr = "Step"
    self.ave = 0

    # flist = list of all log file names

    words = list[0].split()
    self.flist = []
    for word in words: self.flist += glob.glob(word)
    if len(self.flist) == 0 and len(list) == 1:
      raise StandardError,"no log file specified"

    if len(list) > 1 and len(list[1]): self.firststr = list[1]
    if len(list) == 3: self.ave = 1
    
    self.read_all()

  # --------------------------------------------------------------------
  # read all log data from all files
  
  def read_all(self):
    self.read_header(self.flist[0])
    if self.nvec == 0: raise StandardError,"log file has no values"

    # read all files

    for file in self.flist: self.read_one(file)
    print

    # if no average, sort entries by timestep, cull duplicates
    # if average, call self.average()

    if self.ave == 0:
      self.data.sort(self.compare)
      self.cull()
    else: self.average()

    self.nlen = len(self.data)
    print "read %d log entries" % self.nlen

  # --------------------------------------------------------------------

  def get(self,*keys):
    if len(keys) == 0:
      raise StandardError, "no log vectors specified"

    map = []
    for key in keys:
      if self.ptr.has_key(key):
        map.append(self.ptr[key])
      else:
        count = 0
        for i in range(self.nvec):
	  if self.names[i].find(key) == 0:
	    count += 1
	    index = i
        if count == 1:
          map.append(index)
        else:
          raise StandardError, "unique log vector %s not found" % key

    vecs = []
    for i in range(len(keys)):
      vecs.append(self.nlen * [0])
      for j in xrange(self.nlen):
        vecs[i][j] = self.data[j][map[i]]

    if len(keys) == 1: return vecs[0]
    else: return vecs

  # --------------------------------------------------------------------

  def write(self,filename,*keys):
    if len(keys):
      map = []
      for key in keys:
        if self.ptr.has_key(key):
          map.append(self.ptr[key])
        else:
          count = 0
          for i in range(self.nvec):
	    if self.names[i].find(key) == 0:
	      count += 1
	      index = i
          if count == 1:
            map.append(index)
          else:
            raise StandardError, "unique log vector %s not found" % key
    else:
      map = range(self.nvec)

    f = open(filename,"w")
    for i in xrange(self.nlen):
      for j in xrange(len(map)):
        print >>f,self.data[i][map[j]],
      print >>f
    f.close()

  # --------------------------------------------------------------------

  def compare(self,a,b):
    if a[0] < b[0]:
      return -1
    elif a[0] > b[0]:
      return 1
    else:
      return 0

  # --------------------------------------------------------------------

  def cull(self):
    i = 1
    while i < len(self.data):
      if self.data[i][0] == self.data[i-1][0]: del self.data[i]
      else: i += 1

  # --------------------------------------------------------------------

  def average(self):
    counts = []
    data = []
    nlen = 0

    i = j = 0
    while i < len(self.data):
      if self.data[i][0] == 0: j = 0
      if j >= nlen:
        counts.append(0)
        data.append(self.nvec*[0])
        nlen += 1
      counts[j] += 1
      for m in xrange(self.nvec): data[j][m] += self.data[i][m]
      j += 1
      i += 1

    for i in xrange(nlen):
      for j in xrange(self.nvec):
        data[i][j] /= counts[j]

    self.nlen = nlen
    self.data = data
  
  # --------------------------------------------------------------------

  def read_header(self,file):
    if file[-3:] == ".gz":
      txt = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r').read()
    else:
      txt = open(file).read()

    s1 = txt.find(self.firststr)
    s2 = txt.find("\n",s1)
    line = txt[s1:s2]
    words = line.split()
    for i in range(len(words)):
      self.names.append(words[i])
      self.ptr[words[i]] = i

    self.nvec = len(self.names)

  # --------------------------------------------------------------------

  def read_one(self,*list):

    # if 2nd arg exists set file ptr to that value
    # read entire (rest of) file into txt

    file = list[0]
    if file[-3:] == ".gz":
      f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'rb')
    else:
      f = open(file,'rb')

    if len(list) == 2: f.seek(list[1])
    txt = f.read()
    if file[-3:] == ".gz": eof = 0
    else: eof = f.tell()
    f.close()

    start = last = 0
    while not last:

      # chunk = contiguous set of concentration entries
      # s1 = 1st char on 1st line of chunk
      # s2 = 1st char on line after chunk
      # set last = 1 if this is last chunk in file, leave 0 otherwise
      # set start = position in file to start looking for next chunk
      # rewind eof if final entry is incomplete

      s1 = txt.find(self.firststr,start)
      s2 = txt.find("Loop time of",start+1)

      if s1 >= 0 and s2 >= 0 and s1 < s2:    # found s1,s2 with s1 before s2
        s1 = txt.find("\n",s1) + 1
      elif s1 >= 0 and s2 >= 0 and s2 < s1:  # found s1,s2 with s2 before s1
        s1 = 0
      elif s1 == -1 and s2 >= 0:             # found s2, but no s1
	last = 1
        s1 = 0
      elif s1 >= 0 and s2 == -1:             # found s1, but no s2
        last = 1
        s1 = txt.find("\n",s1) + 1
        s2 = txt.rfind("\n",s1) + 1
	eof -= len(txt) - s2
      elif s1 == -1 and s2 == -1:            # found neither
                                             # could be end-of-file section
					     # or entire read was one chunk

        if txt.find("Loop time of",start) == start:   # end of file, so exit
	  eof -= len(txt) - start                     # reset eof to "Loop"
	  break

	last = 1                                      # entire read is a chunk
        s1 = 0
        s2 = txt.rfind("\n",s1) + 1
	eof -= len(txt) - s2
	if s1 == s2: break

      chunk = txt[s1:s2-1]
      start = s2
      
      # split chunk into entries
      # parse each entry for numeric fields, append to data
  
      lines = chunk.split("\n")
      for line in lines:
        words = line.split()
        self.data.append(map(float,words))
  
      # print last timestep of chunk

      print int(self.data[len(self.data)-1][0]),
      sys.stdout.flush()

    return eof
