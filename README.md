# seqGraph
SEQ is a conjugated assembly algorithm that constructs a conjugate map according to Barcode connection information between Contig, and finds the most likely connection relationship between Contig.
## Requirement
htslib
## Install
```
git clone https://github.com/panguangze/seqGraph
mkdir build
cd build
cmake .. && make
```
## Running
#### make input
```
overlap_end input.bam graph.txt 500,
```
The last parameters is the end length of contig for barcode statistic
```
seqGraph \
    -g  graph.txt\ // graph input
    -o  assembly_orders.txt\ //assemly contigs order and orientations
    -i  1 //iteration times, more iteration times longer scaf length but less accuracy
```
make fasta from assembly result
```
python scripts/make_fa.py assembly_orders.txt out.fasta
```
