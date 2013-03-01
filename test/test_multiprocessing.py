'''test script to check memory usage in multiprocessing.'''


import multiprocessing
import GatSegmentList
import time

nsegments=100000000
nsegments=10000000
ncpu = 2
nwork = 100

s = GatSegmentList.SegmentList( iter = [(x,x+1) for x in xrange( 0,nsegments, 2)], normalize = True )

print "built list"
r = raw_input("press return")

s.share("/test")

print "shared data"
r = raw_input("press return")

def dowork( segs ):
    while 1: pass
    return segs.sum()

p = multiprocessing.Pool( ncpu )


print "starting mp"
r = p.map( dowork, [ s for i in range( nwork ) ] )

print r
