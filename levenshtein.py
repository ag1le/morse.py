#!/usr/bin/env python
import os
import sys
import time
from optparse import OptionParser

parser = OptionParser(usage="%prog [OPTIONS] file1 file2")

parser.add_option("-v", "--verbose",
  action="store_true",
  dest="verbose",
  default=False,
  help="Prints details about errors and calls.")


def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
 
    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)
 
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
 
    return previous_row[-1]

def run_ltest(str1,str2):
	fid1 = open(str1);
	fid2 = open(str2);
	  
	s1 = file.read(fid1);
	s2 = file.read(fid2);
	ltest = levenshtein(s1, s2)
	return float(ltest)/float(len(s1))

def main(*args, **kwargs):
  
  (options, args) = parser.parse_args()

  verbosity = 0
  if options.verbose:
    verbosity = 1


  if len(args) < 2:
    print 'usage: <input file1> <input file2>' 
    exit(1)
  elif len(args) == 2:
	ltest = run_ltest(args[0],args[1])
	print ltest	
"""	
  print "L-TEST"
  print levenshtein("124","32")
  print levenshtein("1","2")
  print levenshtein("222","2")
"""  
    
if __name__ == "__main__":
  main()
