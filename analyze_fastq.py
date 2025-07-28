"""
FASTQ Quality Analysis
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


def analyze_fastq(filename, label=""):
    """Analyze a FASTQ file and return statistics"""

    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        return None

    read_lengths = []
    quality_scores = []

    print(f"Analyzing {filename} ({label})...")

    # Parse with error handling
    count = 0
    max_reads = 75000

    try:
        for record in SeqIO.parse(filename, "fastq"):
            if len(record.seq) > 0:
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
        print(f"Warning: Error processing {filename}: {e}")
        if len(read_lengths) == 0:
            return None

    return read_lengths, quality_scores


# Analyze both R1 and R2
fastq_r1 = "SRR562646_1.fastq"
fastq_r2 = "SRR562646_2.fastq"

print("=== PAIRED-END ANALYSIS ===")
result_r1 = analyze_fastq(fastq_r1, "R1 - Forward reads")
result_r2 = analyze_fastq(fastq_r2, "R2 - Reverse reads")

if result_r1 is None or result_r2 is None:
    print("Analysis failed!")
    exit(1)

r1_lengths, r1_quality = result_r1
r2_lengths, r2_quality = result_r2

# Compare statistics
print(f"\n=== COMPARISON RESULTS ===")
print(f"R1 reads processed: {len(r1_lengths):,}")
print(f"R2 reads processed: {len(r2_lengths):,}")

print(f"\nRead Lengths:")
print(f"  R1: {np.mean(r1_lengths):.1f} ± {np.std(r1_lengths):.1f} bp")
print(f"  R2: {np.mean(r2_lengths):.1f} ± {np.std(r2_lengths):.1f} bp")

print(f"\nQuality Scores:")
print(f"  R1: {np.mean(r1_quality):.1f} ± {np.std(r1_quality):.1f}")
print(f"  R2: {np.mean(r2_quality):.1f} ± {np.std(r2_quality):.1f}")

# Quality classification for both
r1_class = (
    "Excellent (Q30+)"
    if np.mean(r1_quality) >= 30
    else "Good (Q20+)" if np.mean(r1_quality) >= 20 else "Poor"
)
r2_class = (
    "Excellent (Q30+)"
    if np.mean(r2_quality) >= 30
    else "Good (Q20+)" if np.mean(r2_quality) >= 20 else "Poor"
)

print(f"\nQuality Assessment:")
print(f"  R1: {r1_class}")
print(f"  R2: {r2_class}")

# Visualization - Compare R1 vs R2
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# R1 Read lengths
sns.histplot(r1_lengths, bins=25, kde=True, color="skyblue", ax=axes[0, 0])
axes[0, 0].set_title("R1 Read Length Distribution")
axes[0, 0].set_xlabel("Length (bp)")
axes[0, 0].axvline(
    np.mean(r1_lengths),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(r1_lengths):.0f} bp",
)
axes[0, 0].legend()

# R2 Read lengths
sns.histplot(r2_lengths, bins=25, kde=True, color="lightcoral", ax=axes[0, 1])
axes[0, 1].set_title("R2 Read Length Distribution")
axes[0, 1].set_xlabel("Length (bp)")
axes[0, 1].axvline(
    np.mean(r2_lengths),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(r2_lengths):.0f} bp",
)
axes[0, 1].legend()

# R1 Quality scores
sns.histplot(r1_quality, bins=25, kde=True, color="lightgreen", ax=axes[1, 0])
axes[1, 0].set_title("R1 Quality Score Distribution")
axes[1, 0].set_xlabel("Average Quality (Phred)")
axes[1, 0].axvline(
    np.mean(r1_quality),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(r1_quality):.1f}",
)
axes[1, 0].axvline(30, color="orange", linestyle=":", label="Q30")
axes[1, 0].legend()

# R2 Quality scores
sns.histplot(r2_quality, bins=25, kde=True, color="gold", ax=axes[1, 1])
axes[1, 1].set_title("R2 Quality Score Distribution")
axes[1, 1].set_xlabel("Average Quality (Phred)")
axes[1, 1].axvline(
    np.mean(r2_quality),
    color="red",
    linestyle="--",
    label=f"Mean: {np.mean(r2_quality):.1f}",
)
axes[1, 1].axvline(30, color="orange", linestyle=":", label="Q30")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("paired_end_comparison.png", dpi=300, bbox_inches="tight")
print("\nSaved: paired_end_comparison.png")
plt.show()

print("\nNote: R1 and R2 are independent reads from paired-end sequencing")
print("Each provides complete quality information for analysis")
