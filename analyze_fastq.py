"""
FASTQ Quality Analysis
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


def analyze_fastq(filename):
    """Analyze a FASTQ file and return statistics"""

    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        return None

    read_lengths = []
    quality_scores = []

    print(f"Analyzing {filename}...")

    # Parse with error handling
    count = 0
    max_reads = 100000

    try:
        for record in SeqIO.parse(filename, "fastq"):
            if len(record.seq) > 0:  # Skip empty reads
                read_lengths.append(len(record.seq))
                quality_scores.append(
                    np.mean(record.letter_annotations["phred_quality"])
                )

            count += 1
            if count % 10000 == 0:
                print(f"  Processed {count:,} reads...")

            if count >= max_reads:
                print(f"  Limited to first {max_reads:,} reads")
                break

    except Exception as e:
        print(f"Warning: Error processing file: {e}")
        if len(read_lengths) == 0:
            return None

    return read_lengths, quality_scores


# Main analysis
fastq_file = "SRR562646_1.fastq"
result = analyze_fastq(fastq_file)

if result is None:
    print("Analysis failed!")
    exit(1)

read_lengths, quality_scores = result

# Statistics
print(f"\n=== RESULTS ===")
print(f"Reads processed: {len(read_lengths):,}")
print(f"Read length: {np.mean(read_lengths):.1f} ± {np.std(read_lengths):.1f} bp")
print(f"Quality score: {np.mean(quality_scores):.1f} ± {np.std(quality_scores):.1f}")

# Quality classification
avg_qual = np.mean(quality_scores)
quality_class = (
    "Excellent (Q30+)"
    if avg_qual >= 30
    else "Good (Q20+)" if avg_qual >= 20 else "Poor"
)
print(f"Overall quality: {quality_class}")

# Visualization
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Read lengths
sns.histplot(read_lengths, bins=30, kde=True, color="skyblue", ax=axes[0])
axes[0].set_title("Read Length Distribution")
axes[0].set_xlabel("Length (bp)")
axes[0].axvline(
    np.mean(read_lengths),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(read_lengths):.0f} bp",
)
axes[0].legend()

# Quality scores
sns.histplot(quality_scores, bins=30, kde=True, color="lightgreen", ax=axes[1])
axes[1].set_title("Quality Score Distribution")
axes[1].set_xlabel("Average Quality (Phred)")
axes[1].axvline(
    np.mean(quality_scores),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(quality_scores):.1f}",
)
axes[1].axvline(30, color="orange", linestyle=":", label="Q30 threshold")
axes[1].legend()

plt.tight_layout()
plt.savefig("fastq_analysis.png", dpi=300, bbox_inches="tight")
print("Saved: fastq_analysis.png")
plt.show()
