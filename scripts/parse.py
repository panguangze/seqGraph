import sys
import re
inp = open(sys.argv[1])
inp2 = open(sys.argv[2])
outs = open(sys.argv[3],"w")
depth = sys.argv[4]
f_th = sys.argv[5]
tmp = {}
for line in inp:
    vs = line.split("\t")
    a = re.split(":|,|;",vs[0])
    tmp[a[0]]=a[1:]
all_segs = {}
write_segs = set()
write_juncs = []
for line in inp2:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        all_segs[vs[1]] = line
        # outs.write(line)
        continue
    if int(vs[-1]) < int(f_th):
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
    if k =="EDGE_5369_length_3828_cov_7.082378":
        print(tmp[k])
        print(v)
    if (k in tmp.keys() and v in tmp[k]):
        # print(v)
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        #write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        write_juncs.append(" ".join(vs)+"\n")
        tmp[k].remove(v)
    if (vc in tmp.keys() and kc in tmp[vc]):
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        #write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        write_juncs.append(" ".join(vs)+"\n")
        tmp[vc].remove(kc)

    write_segs.add(all_segs[vs[1]])
    write_segs.add(all_segs[vs[3]])
# print(tmp['EDGE_5369_length_3828_cov_7.082378'])

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
        write_juncs.append("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
        write_segs.add(all_segs[first])
        write_segs.add(all_segs[second])

for item in write_segs:
    print(item)
    outs.write(item)
for item in write_juncs:
    outs.write(item)