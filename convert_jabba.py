from __future__ import print_function
import sys

def main():
  outfiles = {}

  header = next(sys.stdin).strip()
  for line in sys.stdin:
    line = line.strip()
    guid = line.split('\t')[7]
    if guid not in outfiles:
      outfiles[guid] = open('%s_segments.txt' % guid, 'w')
      print(header, file=outfiles[guid])
    print(line, file=outfiles[guid])

  for outf in outfiles.values():
    outf.close()

main()
