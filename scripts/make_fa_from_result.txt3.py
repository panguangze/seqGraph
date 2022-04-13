from Bio import SeqIO
import sys
from Bio.Seq import Seq

fain = sys.argv[1]
orderin = open(sys.argv[2])
faout = open(sys.argv[3], "w")
blastin = open(sys.argv[4])
blast_ratio = float(sys.argv[5])
gene_hit = sys.argv[6]
plasscore_file = sys.argv[7]
blast_segs = set()
prev_seg = ""
prev_len = 0
for line in blastin.readlines():
    t = line.strip().split("\t")
    if prev_seg != t[0] and prev_seg != "":
        elen = prev_seg.split("_")[3]
        if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)
        prev_seg = t[0]
        prev_len = int(t[3])
    else:
        if float(t[2]) > 60:
            prev_len = prev_len + int(t[3])
        prev_seg = t[0]
elen = int(prev_seg.split("_")[3])
if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
    blast_segs.add(t[0])
# gene_res = set()
# with open(gene_hit, 'r') as gene_lst:
#     for gene_r in gene_lst:
#         gene_res.add(gene_r.strip())


plasscore = {}
with open(plasscore_file, 'r') as ps:
    for s in ps:
        item = s.strip().split('\t')
        id_ = int(item[0].split('_')[1])
        plasscore[id_] = float(item[1])

genehit = {}
for i in range(len(orderin.readlines())):
    genehit[i + 1] = False
with open(gene_hit, 'r') as gh:
    for s in gh:
        item = s.strip().split('\t')
        if len(s.strip()) == 0:
            continue
        id_ = int(item[0].split('_')[1])
        genehit[id_] = True

record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
n_seq = Seq("N" * 40)
self_tag = False
for idx, line in enumerate(orderin.readlines()):
    if line.startswith("iter") or line.startswith("self"):
        if line.startswith("self"):
            self_tag = True
        continue
    seq = ""
    tmp = line.strip().split("\t")
    if len(tmp) == 1 and self_tag:
        # if tmp[0][0:-1] not in blast_segs or plasscore[idx] < 0.7:
        #     continue
        if genehit[idx + 1] or plasscore[idx + 1] >= 0.7:
            for t in tmp:
                tmp_seq = record_dict[t[0:-1]].seq
                if t[-1] == '-':
                    tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
                seq = seq + tmp_seq
            # faout.write(">self" + "".join(tmp) + "\n" + str(seq) + "\n")
            print(">self-gene" + "".join(tmp))
        else:
            for t in tmp:
                tmp_seq = record_dict[t[0:-1]].seq
                if t[-1] == '-':
                    tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
                seq = seq + tmp_seq
                faout.write(">self" + "".join(tmp) + "\n" + str(seq) + "\n")
        continue
    flags = False
    is_gene = False
    blast_len = 0
    all_len = 0
    for t in tmp:
        if genehit[idx + 1]:
            flags = True
            is_gene = True
        all_len = all_len + elen
        if t[0:-1] in blast_segs:
            blast_len = blast_len + elen
    if (float(blast_len) / float(all_len)) > 0.1:
        flags = True
    if all_len < 1000:
        flags = False
    if not flags and plasscore[idx + 1] < 0.7:
        continue
    for t in tmp:
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    if genehit[idx + 1] and plasscore[idx + 1] >= 0.7:
        print(">gene_score" + "".join(tmp))
    else:
        if plasscore[idx + 1] >= 0.7:
            faout.write(">score" + "".join(tmp) + "\n" + str(seq) + "\n")
        else:
            faout.write(">gene" + "".join(tmp) + "\n" + str(seq) + "\n")
