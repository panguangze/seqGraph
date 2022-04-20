import sys
import re
inp = open(sys.argv[1])
inp2 = open(sys.argv[2])
segf = open(sys.argv[4])
outs = open(sys.argv[3],"w")
tmp = {}
segs = []
depth = 30
for line in inp:
    vs = line.split("\t")
    a = re.split(":|,|;",vs[0])
    tmp[a[0]]=a[1:]
for line in segf:
    t = re.split(r"[\+\-]",line)
    for v in t:
        segs.append(v)
for line in inp2:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        if vs[1] in segs:
            outs.write(line)
        continue
    k = vs[1]
    kc = vs[1]+"'"
    if vs[2] == '-':
        k = vs[1]+"'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3]+"'"
    if vs[4] == '-':
        v = vs[3]+"'"
        vc = vs[3]
    if vs[1] in segs and vs[3] in segs:
        outs.write(" ".join(vs)+"\n")
for item in tmp.keys():
    first = item
    fdir = "+"
    if item[-1] == "'":
        first = item[0:-1]
        fdir = "-"
    for i in tmp[item]:
        if i == "":
            continue
        second = i
        sdir = "+"
        if i[-1] == "'":
            second = i[0:-1]
            sdir = "-"
        if first in segs and second in segs:
            outs.write("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))

