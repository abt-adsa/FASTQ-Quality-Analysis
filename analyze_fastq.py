"""
FASTQ analyzer with basic statistics
"""

from Bio import SeqIO
import matplotlib.pyplot as plt

# Load FASTQ file
fastq_file = "SRR562646_1.fastq"

read_lengths = []

# Parse the FASTQ file
for record in SeqIO.parse(fastq_file, "fastq"):
    read_lengths.append(len(record.seq))

# Basic statistics
print(f"Total reads: {len(read_lengths)}")
print(f"Average read length: {sum(read_lengths)/len(read_lengths):.2f}")
print(f"Min: {min(read_lengths)}, Max: {max(read_lengths)}")

# Simple histogram
plt.figure(figsize=(10, 6))
plt.hist(read_lengths, bins=20)
plt.title("Read Length Distribution")
plt.xlabel("Read Length")
plt.ylabel("Count")
plt.show()
