'''test script to check memory usage in multiprocessing.'''


import multiprocessing
from gat.SegmentList import SegmentList

nsegments = 100000000
nsegments = 10000000
ncpu = 2
nwork = 100

s = SegmentList(
    iter=[(x, x + 1) for x in range(0, nsegments, 2)], normalize=True)

print("built list")
r = input("press return")

s.share("/test")

print("shared data")
r = input("press return")


def dowork(segs):
    while 1:
        pass
    return segs.sum()

p = multiprocessing.Pool(ncpu)


print("starting mp")
r = p.map(dowork, [s for i in range(nwork)])

print(r)
