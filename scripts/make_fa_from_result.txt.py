from Bio import SeqIO
import sys
fain = sys.argv[1]
orderin = sys.argv[2]
faout= sys.argvs[3]
record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
for line in orderin.readlines():
    seq = ""
    tmp = line.split("\t")
    for t in tmp:
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    faout.write(seq+"\n")
