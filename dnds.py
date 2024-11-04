import argparse
import math

# Amino acid array to map codon numbers to amino acids
amino_acid_array = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L',
    'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P',
    'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q',
    'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', '*', 'W', 'R', 'R', 'R', 'R',
    'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G'
]

# Function to convert a codon to a unique number (0-63)
def codon_to_number(codon):
    base_number = {'T': 0, 'C': 1, 'A': 2, 'G': 3}
    return base_number[codon[1]] * 16 + base_number[codon[0]] * 4 + base_number[codon[2]]

# Function to calculate potential synonymous sites for a codon
def syn_sites(codon):
    syn_site_array = [
        1, 1, 2, 2, 3, 3, 4, 4, 2, 2, 2, 0, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 3, 3, 4, 4, 1, 1, 2, 2, 3, 3, 3, 3
    ]
    return syn_site_array[codon_to_number(codon)]

# Function to process pairwise comparison of sequences and write to TSV
def process_pairwise_comparisons(names, sequences, output_file):
    num_sequences = len(sequences)
    seq_length = len(sequences[0])

    with open(output_file, 'w') as out:
        # Write header for TSV
        out.write("Seq1\tSeq2\tSynonymous Changes\tNon-synonymous Changes\tdS\tdN\tdS/dN Ratio\n")

        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                A = sequences[i]
                B = sequences[j]
                syn_codons = 0.0
                nonsyn_codons = 0.0
                SA_Nei = 0.0
                SB_Nei = 0.0
                count_insertions = 0
                count_Ns = 0
                total_codons = 0

                for k in range(0, seq_length, 3):
                    codonA = A[k:k + 3]
                    codonB = B[k:k + 3]
                    if '-' in codonA + codonB:
                        count_insertions += 1
                        continue
                    if 'N' in codonA + codonB:
                        count_Ns += 1
                        continue

                    total_codons += 1
                    SA_Nei += syn_sites(codonA)
                    SB_Nei += syn_sites(codonB)

                    aaA = amino_acid_array[codon_to_number(codonA)]
                    aaB = amino_acid_array[codon_to_number(codonB)]

                    if aaA == aaB:
                        if codonA != codonB:
                            syn_codons += 1
                    else:
                        nonsyn_codons += 1

                # Calculations for Nei's method
                compared_codons = total_codons - (count_insertions + count_Ns)
                potential_syn = (SA_Nei / 3 + SB_Nei / 3) / 2
                potential_nonsyn = 3 * compared_codons - potential_syn

                ps = syn_codons / potential_syn if potential_syn else 0
                pn = nonsyn_codons / potential_nonsyn if potential_nonsyn else 0

                ds = -3.0 / 4.0 * math.log(1 - (4.0 * ps / 3.0)) if ps < 0.75 else 'NA'
                dn = -3.0 / 4.0 * math.log(1 - (4.0 * pn / 3.0)) if pn < 0.75 else 'NA'

                ds_dn_ratio = (ds / dn) if ds != 'NA' and dn != 'NA' and dn != 0 else 'NA'

                # Write results to TSV
                out.write(f"{names[i]}\t{names[j]}\t{syn_codons:.2f}\t{nonsyn_codons:.2f}\t{ds}\t{dn}\t{ds_dn_ratio}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate dS/dN for sequences and output to TSV.")
    parser.add_argument("input_file", help="Input FASTA-like formatted file")
    parser.add_argument("output_file", help="Output TSV file")
    args = parser.parse_args()

    # Read sequences from file
    names = []
    sequences = []
    with open(args.input_file, 'r') as file:
        for line in file:
            if line.strip():
                name, seq = line.strip().split()
                names.append(name)
                sequences.append(seq)

    process_pairwise_comparisons(names, sequences, args.output_file)

if __name__ == "__main__":
    main()
