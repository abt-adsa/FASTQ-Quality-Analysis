"""
FASTQ Quality Analysis
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load FASTQ file
fastq_file = "SRR562646_1.fastq"

read_lengths = []
quality_scores = []

print("Processing FASTQ file...")
print("Note: This is R1 (forward reads) - should also analyze R2 later")

# Parse the FASTQ file
count = 0
for record in SeqIO.parse(fastq_file, "fastq"):
    read_lengths.append(len(record.seq))

    # Get quality scores
    qualities = record.letter_annotations["phred_quality"]
    avg_quality = sum(qualities) / len(qualities)
    quality_scores.append(avg_quality)

    count += 1
    if count % 10000 == 0:
        print(f"Processed {count} reads...")

    # Limit for testing
    if count >= 50000:
        break

# Statistics
print(f"\nTotal reads: {len(read_lengths)}")
print(f"Average read length: {np.mean(read_lengths):.2f}")
print(f"Average quality: {np.mean(quality_scores):.2f}")

# Plot both metrics
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.hist(read_lengths, bins=20)
ax1.set_title("Read Lengths")
ax1.set_xlabel("Length")

ax2.hist(quality_scores, bins=20)
ax2.set_title("Quality Scores")
ax2.set_xlabel("Average Quality")

plt.tight_layout()
plt.savefig("analysis_v3.png")
plt.show()
