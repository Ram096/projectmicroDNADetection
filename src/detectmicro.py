import pysam
from collections import defaultdict

def detect_microdnas(bam_file, soft_clip_threshold=10, circle_score_threshold=5):
    """
    Detect and quantify microDNAs from a BAM file.

    Parameters:
        bam_file (str): Path to the input BAM file.
        soft_clip_threshold (int): Minimum length of soft clips to consider evidence for circles.
        circle_score_threshold (int): Minimum score to report a detected circle.

    Returns:
        list: A list of detected circles, each as a dictionary with position and score.
    """
    # Store potential circles
    circles = defaultdict(list)

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Check if the read has soft clips in the CIGAR string
            if any(cigar_op[0] == 4 and cigar_op[1] >= soft_clip_threshold for cigar_op in read.cigartuples):
                # Record the evidence for a circle
                circles[read.reference_start].append({
                    "query_name": read.query_name,
                    "reference_end": read.reference_end,
                    "cigarstring": read.cigarstring,
                    "mapping_quality": read.mapping_quality,
                })

    # Aggregate and score circles
    detected_circles = []
    for start, evidence in circles.items():
        score = sum(e["mapping_quality"] for e in evidence) / len(evidence)
        if score >= circle_score_threshold:
            detected_circles.append({
                "position": (start, evidence[0]["reference_end"]),
                "score": score,
                "evidence_count": len(evidence),
            })

    return detected_circles


def validate_circles(bam_file, detected_circles, sample_size=5):
    """
    Validate a subset of detected circles by reanalyzing their evidence.

    Parameters:
        bam_file (str): Path to the input BAM file.
        detected_circles (list): List of detected circles.
        sample_size (int): Number of circles to validate.

    Returns:
        list: Validation results for the sampled circles.
    """
    import random
    sample = random.sample(detected_circles, min(sample_size, len(detected_circles)))
    validation_results = []

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for circle in sample:
            start, end = circle["position"]
            evidence_count = 0
            bam.reset()
            for read in bam.fetch():
                if read.reference_start == start and read.reference_end == end:
                    evidence_count += 1

            validation_results.append({
                "position": circle["position"],
                "original_score": circle["score"],
                "original_evidence_count": circle["evidence_count"],
                "validation_evidence_count": evidence_count,
            })

    return validation_results


if __name__ == "__main__":
    # Path to your BAM file
    bam_file = "../data/SRR413984.sorted.NC_000001.10.bam"
    output_file = "microdna_summary.txt"

    # Detect microDNAs
    print("Detecting microDNAs...")
    detected_circles = detect_microdnas(bam_file)

    # Print and save summary
    num_detected = len(detected_circles)
    summary_lines = [f"Identified total of {num_detected} candidate microDNAs.\n"]

    for circle in detected_circles:
        line = f"Circle at {circle['position']} with score {circle['score']:.2f} and evidence count {circle['evidence_count']}"
        print(line)
        summary_lines.append(line)

    # Write to file
    with open(output_file, "w") as f:
        f.write("\n".join(summary_lines))

    print(f"\nSummary saved to {output_file}")

    # Validate a subset of circles
    print("\nValidating detected circles...")
    validation_results = validate_circles(bam_file, detected_circles)
    for result in validation_results:
        print(f"Circle {result['position']} - Original Evidence: {result['original_evidence_count']}, "
              f"Validation Evidence: {result['validation_evidence_count']}")
