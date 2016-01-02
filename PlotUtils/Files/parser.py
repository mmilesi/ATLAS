#!/usr/bin/python

out_file = open("samples2015correct.csv","w")
infile = open("samples2015test.csv")
for line in infile:
  line.strip('')
  line = line.split(',')
  idx = 0
  for token in line :
    print " line[",idx,"]: ", token
    idx = idx+1
  print "\n"
  new_line = line[0] + ',' + line[1] + ',' + line[3] + ',' + line[4] + ',' + line [5] + ',' + line[2] + ',' + line [6] + ',' + line[7]
  print " new line: ", new_line
  out_file.write(new_line)
out_file.close()
