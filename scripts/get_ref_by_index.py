import sys
from pyfaidx import SeqIO

fa = open(sys.argvs[1])
fai = open(sys.argvs[2])
ref_indx=open(sys.argvs[3])
out_ref = open(sys.argvs[4],"w")

fai_array = []
for line in fai.readlines():
    t = line.split("\t")
    fai_array.append(t[0])

record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

for line in ref_indx.readlines():
    t = line.split("\t")
    ref_name = fai_array[t[1]]
    out_ref.write(">"+ref_name+"\n")
    out_ref.write(str(record_dict[ref_name].seq)+"\n")