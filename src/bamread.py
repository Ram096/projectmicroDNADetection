import pysam

# Open a BAM file for reading
bamfile = pysam.AlignmentFile("../data/SRR413984.sorted.NC_000001.10.bam", "rb")  # "rb" = read binary

# Iterate through all aligned reads
for read in bamfile:
    print(read.query_name, read.reference_name, read.reference_start, read.cigarstring)

bamfile.close()
