"""
FASTQ Quality Analysis
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load FASTQ file
fastq_file = "SRR562646_1.fastq"

read_lengths = []
quality_scores = []

print("Analyzing FASTQ file...")
print("Analyzing R1 (forward reads) - paired-end data")

# Parse the FASTQ file with progress
count = 0
max_reads = 75000

for record in SeqIO.parse(fastq_file, "fastq"):
    read_lengths.append(len(record.seq))
    quality_scores.append(np.mean(record.letter_annotations["phred_quality"]))

    count += 1
    if count % 10000 == 0:
        print(f"Processed {count} reads...")

    if count >= max_reads:
        break

# Enhanced statistics
print(f"\nProcessed {len(read_lengths)} reads")
print(
    f"Read length - Mean: {np.mean(read_lengths):.2f}, Median: {np.median(read_lengths):.2f}"
)
print(
    f"Quality score - Mean: {np.mean(quality_scores):.2f}, Median: {np.median(quality_scores):.2f}"
)

# Quality assessment
avg_qual = np.mean(quality_scores)
if avg_qual >= 30:
    print("Quality: Excellent (Q30+)")
elif avg_qual >= 20:
    print("Quality: Good (Q20+)")
else:
    print("Quality: Poor")

# Better plots with seaborn
sns.set_style("whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Read length distribution
sns.histplot(read_lengths, bins=25, kde=True, ax=ax1)
ax1.set_title("Read Length Distribution")
ax1.set_xlabel("Read Length (bp)")
ax1.axvline(np.mean(read_lengths), color="red", linestyle="--", label="Mean")
ax1.legend()

# Quality distribution
sns.histplot(quality_scores, bins=25, kde=True, ax=ax2, color="green")
ax2.set_title("Quality Score Distribution")
ax2.set_xlabel("Average Quality Score")
ax2.axvline(np.mean(quality_scores), color="red", linestyle="--", label="Mean")
ax2.axvline(30, color="orange", linestyle=":", label="Q30")
ax2.legend()

plt.tight_layout()
plt.savefig("quality_analysis_v4.png", dpi=300)
print("Saved plot as quality_analysis_v4.png")
plt.show()
