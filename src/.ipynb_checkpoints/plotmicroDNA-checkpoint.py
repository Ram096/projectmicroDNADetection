import matplotlib.pyplot as plt
from detectmicro import detect_microdnas, validate_circles


def plot_microdna_scores(detected_circles, output_file="microdna_scores_distribution.png"):
    """
    Plot and save the distribution of scores for detected microDNAs.

    Parameters:
        detected_circles (list): List of detected circles.
        output_file (str): Filename to save the plot.
    """
    scores = [circle["score"] for circle in detected_circles]
    plt.figure(figsize=(8, 6))
    plt.hist(scores, bins=20, color='blue', alpha=0.7, edgecolor='black')
    plt.title("Distribution of MicroDNA Scores")
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(output_file)  # Save the plot as a file
    plt.close()  # Close the plot


def plot_evidence_count(detected_circles, output_file="evidence_count_vs_score.png"):
    """
    Plot and save the relationship between evidence count and score for detected microDNAs.

    Parameters:
        detected_circles (list): List of detected circles.
        output_file (str): Filename to save the plot.
    """
    scores = [circle["score"] for circle in detected_circles]
    evidence_counts = [circle["evidence_count"] for circle in detected_circles]

    plt.figure(figsize=(8, 6))
    plt.scatter(evidence_counts, scores, color='green', alpha=0.7)
    plt.title("Evidence Count vs. MicroDNA Score")
    plt.xlabel("Evidence Count")
    plt.ylabel("Score")
    plt.grid(True)
    plt.savefig(output_file)  # Save the plot as a file
    plt.close()  # Close the plot


def plot_positions(detected_circles, output_file="genomic_positions_vs_score.png"):
    """
    Plot and save the genomic positions of detected microDNAs.

    Parameters:
        detected_circles (list): List of detected circles.
        output_file (str): Filename to save the plot.
    """
    positions = [circle["position"][0] for circle in detected_circles]
    scores = [circle["score"] for circle in detected_circles]

    plt.figure(figsize=(10, 6))
    plt.scatter(positions, scores, alpha=0.7, color='orange', edgecolor='black')
    plt.title("Genomic Positions vs. MicroDNA Score")
    plt.xlabel("Genomic Position (Start)")
    plt.ylabel("Score")
    plt.grid(True)
    plt.savefig(output_file)  # Save the plot as a file
    plt.close()  # Close the plot


def main():
    # Path to your BAM file
    bam_file = "../data/SRR413984.sorted.NC_000001.10.bam"

    # Detect microDNAs
    print("Running microDNA detection...")
    detected_circles = detect_microdnas(bam_file)

    # Generate and save plots
    print("Saving plots...")
    plot_microdna_scores(detected_circles)
    plot_evidence_count(detected_circles)
    plot_positions(detected_circles)

    print("Plots have been saved successfully.")


if __name__ == "__main__":
    main()