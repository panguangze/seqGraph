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
print(tmp)
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
    if (k in tmp.keys() and v in tmp[k]) or (vc in tmp.keys() and kc in tmp[vc]):
        outs.write(" ".join(vs)+"\n")