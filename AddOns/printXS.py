#!/usr/bin/env python

import sys
desc="Read cross-sections from results data base. All units are in pb."
usage="%s Results.db"%(sys.argv[0])

from optparse import OptionParser
parser = OptionParser(usage=usage, description=desc)
parser.add_option("-v", dest="VERBOSITY", help="VERBOSITY level", default=1, type=int)

(opts, args) = parser.parse_args()

if len(args) !=1:
    print usage
    sys.exit(1)


# Open data base and get all items
import sqlite3
c=sqlite3.connect(sys.argv[1]).cursor()
c.execute("select file from path")

raw_files=[i[0] for i in c if "XS_" in i[0] and not i[0].endswith(".bvi")]
xs_entries=[]
# We need to filter out the summary files
for i in raw_files:
    temp = i.split("/")
    if not temp[2] in temp[1]:
        xs_entries.append(i)

processes=[i.split("XS_")[-1] for i in xs_entries]
raw_groups   =list(set([i.split("XS_")[-1].split("/")[0] for i in xs_entries]))

xsgroups={}
for r in raw_groups:
    xsgroups[r] = [r]
    if "BVI" in r:
        mult, proc_core = r.split("__",1)
        proc_core=proc_core.rsplit("__",1)[0]
        # find corresponding RS xs
        a,b = mult.split("_")
        m_rs = a+"_"+str(int(b)+1)
        g_rs = [x for x in raw_groups if m_rs in x and proc_core in x and "(RS)" in x][0]
        xsgroups[r.replace("BVI", "BVI+RS")] = [r, g_rs]


def getSingleXS(cursor, name):
    # Factor to convert GeV2 to picobarn
    to_pb = 3.89379656e8
    cursor.execute("select content from path where file='%s'"%name)
    xs_raw = [l for l in cursor][0][0].split()
    return  (float(xs_raw[1]) * to_pb, float(xs_raw[3]) * to_pb)


xs_sum = {}
for k, v in xsgroups.iteritems():
    t_xs, delta_t_xs = 0 ,0
    for g in v:
        this_proc = [getSingleXS(c,p) for p in xs_entries if g in p]
        for tp in this_proc:
            t_xs += tp[0]
            delta_t_xs += tp[1]**2
    from math import sqrt
    xs_sum[k] = (t_xs, sqrt(delta_t_xs))

if opts.VERBOSITY==2:
    for k in sorted(xs_sum.keys()):
        print k, "%e"%xs_sum[k][0], "%e"%xs_sum[k][1]
elif opts.VERBOSITY==3:
    for p in sorted(xs_entries):
        print p, getSingleXS(c, p)
else:
    red = [x for x in xs_sum.keys() if not "(BVI)" in x and not "(RS)" in x]
    for k in sorted(red):
        print k, "%e"%xs_sum[k][0], "%e"%xs_sum[k][1]



sys.exit(0)
