import re

def internal_rbs_checker(dna):
    """
    Finds potential ribosome binding sites in a DNA sequence.
    Searches for potential Shine-Dalgarno sequences with proper spacing to a start codon (AUG/ATG)
    Spacing ranges and Shine-Dalgarno sequences taken from "Determination of the optimal aligned 
        spacing between the Shine-Dalgarno sequence and the translation initiation codon of 
        Escherichia coli mRNAs." (Chen et al, 1994, in Nucleic Acids Research)

    Parameters:
        dna (str): Input DNA sequence.

    Returns:
        tuple: (True, None) if no sites found, (False, potential_sites) if sites are found.
    """
    # Convert T to U for RNA
    rna_sequence = dna.upper().replace('T', 'U')

    # Define Shine-Dalgarno sequences
    sd_sequences = [
        "UAAGGAGGU", "AAGGG", "GGA", "AGGA", "AAGG", "AAGGA", "AAGGGU", "AAGGAGGU"
    ]
    # Spacing ranges
    spacing_ranges = {
        "UAAGGAGGU": (3, 11),
        "AAGGG": (2, 15),
        "GGA": (5, 8),
        "AGGA": (7, 17),
        "AAGG": (6, 23),
        "AAGGA": (4, 12),
        "AAGGGU": (7, 13),
        "AAGGAGGU": (3, 11)
    }

    # Start codon
    start_codon = "AUG"

    # Find all matches
    potential_sites = []
    for sd in sd_sequences:
        spacing_min, spacing_max = spacing_ranges.get(sd, (0, 0))
        # Look for Shine-Dalgarno sequence matches
        for match in re.finditer(sd, rna_sequence):
            sd_start = match.start()
            sd_end = match.end()
            # Look for start codon within the spacing range
            start_search_start = sd_end + spacing_min
            start_search_end = sd_end + spacing_max + len(start_codon)
            if start_search_end > len(rna_sequence):
                start_search_end = len(rna_sequence)
            search_region = rna_sequence[start_search_start:start_search_end]
            start_codon_match = search_region.find(start_codon)
            while(start_codon_match) != -1:
                start_codon_position = start_search_start + start_codon_match
                potential_sites.append({
                    "Shine-Dalgarno Sequence": sd,
                    "SD Start": sd_start,
                    "SD End": sd_end,
                    "Start Codon Position": start_codon_position,
                    "Spacing": start_codon_position - sd_end
                })
                start_codon_match = search_region.find(start_codon, start_codon_match + 3)

    # Return based on the presence of potential sites
    if potential_sites:
        return False, potential_sites
    else:
        return True, None


# Example usage
if __name__ == "__main__":
    dna_seq = "ATGGAATGGTGUAAGGAGGTTTACCATGAGGATAGGGAAGGTTATGAAGGGTT"
    status, results = internal_rbs_checker(dna_seq)

    if status:
        print("No potential ribosome binding sites found.")
    else:
        print("Potential ribosome binding sites found:")
        for site in results:
            print(site)
