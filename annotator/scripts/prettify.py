#
# Makes ascii output of Annotator bit more readable
#

import sys
import kaks

printall = False
printenst = False
printkaks = False
sign = 1

filename = sys.argv[1]
pvalcutoff = float( sys.argv[2] )

if len(sys.argv) == 4:
    printkaks = True
    anndict = kaks.kaksdict( sys.argv[3], 1 )

print "# Output file: "+filename

lines = open(filename,'r').readlines()
data = []
dummy = "----"
skipstats = False
dataline = False

for line in lines:
    line = line[:-1]
    if (line + dummy)[:10] == "# Isochore":
        skipstats = True
    elif (line + dummy)[:3] != "#  ":
        skipstats = False

    if printall:
        printline = (line + dummy)[0] == '#' and not skipstats
    else:
        printline = (line + dummy)[:9] == "# Reading"
    if (line + dummy)[:2] == "--":
        printline = True

    if printline and not printenst:
        print line

    if dataline:
        lineelts = line.split('\t')
        if len(lineelts) == 9:
            devexp, relrep, exppval, dum1, dum2, annprop, eprop, sd, ann = lineelts
            values = (devexp, relrep, exppval, annprop, eprop, sd)
            values = map(float, values)
            data.append( (ann, values) )

    if (line + dummy)[:6] == "DevExp":
        dataline = True
        if not printenst:
            print "--"
            print "--  Using p-value cutoff: ",pvalcutoff

data2 = []
prevvalues = 0
for ann, values in data:
    (devexp, relrep, exppval, annprop, eprop, sd) = values
    if exppval <= pvalcutoff * 1.00001:
        if prevvalues == values:
            data2[-1][0] += " +"
        else:
            prevvalues = values
            data2.append( [ann, values] )

print "RelRep\tExperP\tZ-score\tExpected\tObserved\tAnnotation"
for ann, values in data2:
    (devexp, relrep, exppval, annprop, eprop, sd) = values
    if printenst:
        if relrep * sign > 0.0:
            print ann.split(' ')[1]
    else:
        if printkaks:
            if ann in anndict:
                median, quant, ann = anndict[ann]
            else:
                median, quant = -1, -1
            print "%1.2f\t%0.5f\t%0.7f\t%0.7f\t%0.3f\t%2d\t%s" % (1+relrep/100.0, exppval, eprop, annprop, median, quant, ann)
        else:
            print "%2.2f\t%0.5f\t%1.2f\t%0.7f\t%0.7f\t%s" % (relrep, exppval, devexp, eprop, annprop, ann)

    

    
    
