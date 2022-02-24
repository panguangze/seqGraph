import sys
import re
inp = open(sys.argv[1])
inp2 = open(sys.argv[2])
outs = open(sys.argv[3],"w")
tmp = {}
for line in inp:
    vs = line.split("\t")
    a = re.split(":|,|;",vs[0])
    tmp[a[0]]=a[1:]
for line in inp2:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
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
    if (k in tmp.keys() and v in tmp[k]):
        outs.write(" ".join(vs)+"\n")
        tmp[k].remove(v)
    elif (vc in tmp.keys() and kc in tmp[vc]):
        outs.write(" ".join(vs)+"\n")
        tmp[vc].remove(kc)
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
        outs.write("JUNC {} {} {} {} 1\n".format(first,fdir,second,sdir,))
