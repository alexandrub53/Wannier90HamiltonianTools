#!/usr/bin/python 
# ASM 7/29/13

#! /usr/bin/env python

# The purpose of this script is to lead the band structure output and
# then dump the data into the standard output with one line per k
# point and all eigenvalues sorted onto that line

import sys

try:
  fn_bands=sys.argv[1]
    
except:
  print 'Usage: ',sys.argv[0],' out.bands'
  print '       out.bands - formatted band structure data'
  sys.exit(1)

# Read band structure information from out.bands
fin = open(fn_bands,'r')
lines = fin.readlines()
fin.close()

pos_list=[] # list of positions of lines with band energy info
for i,line in enumerate(lines):
  if 'number of Kohn-Sham states' in line:
    fields=line.split('=')
    num_bnd=int(fields[1])
  if 'number of k points' in line:
    num_kpt=int(line.split('=')[1].split()[0])
  if 'band energies' in line:
    pos_list.append(i)
  if 'bands (ev)' in line:
    pos_list.append(i)

xvals=[]
for i in range(num_kpt):
  xvals.append(float(i)/float(num_kpt-1))
bnd = []
for i in range(num_bnd):
  bnd.append([])
for pos in pos_list:
  num_read=0
  j=pos+1
  while num_read<num_bnd:
    j+=1
    line=lines[j]
    fields=line.split()
    for i in range(len(fields)):
      bnd[num_read+i].append(float(fields[i]))
    num_read=num_read+len(fields)

for ik in range(num_kpt) :
   print xvals[ik] , 
   for i in range(num_bnd) :
     print bnd[i][ik] ,
   print

