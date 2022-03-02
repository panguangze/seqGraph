from Bio import SeqIO
import sys
fain = sys.argv[1]
orderin = open(sys.argv[2])
faout= open(sys.argv[3],"w")
record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
for line in orderin.readlines():
    seq = ""
    tmp = line.strip().split("\t")
    for t in tmp:
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    faout.write(seq+"\n")
