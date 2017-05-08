import sys
import os
import re

import gat

track = 0
for line in sys.stdin:
    if line.startswith("##Iteration"):
        sys.stdout.write("track name=%i\n" % track)
        track += 1

    if not line.startswith("##Segments"):
        continue

    data = line[:-1].split("\t")[1:]
    contig = data[0]
    coords = [list(map(int, x[1:-1].split(","))) for x in data[1:]]
    # print "annotator", coords

    for start, end in coords:
        sys.stdout.write("\t".join(map(str, (contig, start, end))) + "\n")
