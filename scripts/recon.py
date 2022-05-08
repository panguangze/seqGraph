import sys
from operator import itemgetter
import re

fai = open(sys.argv[1])  # fai
graph = open(sys.argv[2])  # graph
segf = open(sys.argv[4])  # 需要从新的segs
outs = sys.argv[3]  # 输出prefix
depth = float(sys.argv[5])  # samtools 计算而来
tmp = {}
segs = []
seg_g = {}
out_juncs = []
for line in fai:
    vs = line.split("\t")
    a = re.split(":|,|;", vs[0])
    tmp[a[0]] = a[1:]
ref_name_records = {}
for idx, line in enumerate(segf):
    l = line.strip().split('\t')
    t = re.split(r"[\+\-]", l[0])
    ref_name_records[idx] = l[1]
    line_segs = []
    out_juncs.append([])
    for v in t:
        line_segs.append(v)
        seg_g[v] = idx
    segs.append(line_segs)

out_segs = [[] for i in range(0, len(out_juncs))]


def in_segs(seg1, seg2=None):
    if seg2 is None:
        for s in segs:
            if seg1 in s:
                return 1
    else:
        for i in len(segs):
            s = segs[i]
            if seg1 in s and seg2 in s:
                return i + 1
    return 0


for line in graph:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        if in_segs(vs[1]):
            out_segs[seg_g[vs[1]]].append(line)
        continue

    k = vs[1]
    kc = vs[1] + "'"
    if vs[2] == '-':
        k = vs[1] + "'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3] + "'"
    if vs[4] == '-':
        v = vs[3] + "'"
        vc = vs[3]
    if in_segs(vs[1], vs[3]):
        out_juncs[in_segs(vs[1], vs[3]) - 1].append(vs)
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
        if in_segs(first, second):
            vs = [first, fdir, second, sdir, str(int(depth) / 2)]
            out_juncs[in_segs(vs[1], vs[3]) - 1].append(vs)
            # out_juncs.append(vs)
            # out_juncs("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
# s_juncs = sorted(out_juncs, key=itemgetter(2))
i = 0
for js in out_segs:
    outs_f = open(outs + str(i) + "_" + ref_name_records[i] + ".txt", "w")
    i = i + 1
    for j in js:
        outs_f.write(" ".join(j) + "\n")
    outs_f.close()
i = 0
for js in out_juncs:
    outs_f = open(outs + str(i) + "_" + ref_name_records[i] + ".txt", "a")
    i = i + 1
    for j in sorted(js):
        outs_f.write(" ".join(j) + "\n")
    outs_f.close()
