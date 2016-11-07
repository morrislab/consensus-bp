from __future__ import print_function
import sys

def main():
  for line in sys.stdin:
    guid, methods = line.split()
    methods = methods.split(',')
    if len(methods) < 5:
      print(guid, methods)
    #if len(methods) >= 5:
      #print(line.strip())

main()
