import argparse
import logging
import math
import itertools
from multiprocessing import Pool, Lock, Manager

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Dictionary for base to number conversion
base_to_number = {'T': 0, 'C': 1, 'A': 2, 'G': 3}

# Array of amino acids corresponding to codon indices
aa_array = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L',
    'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P',
    'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', 'Z', 'Z', 'H', 'H', 'Q', 'Q',
    'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', 'Z', 'W', 'R', 'R', 'R', 'R',
    'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G'
]

# Array of synonymous sites
syn_site_array = [
    1, 1, 2, 2, 3, 3, 4, 4, 2, 2, 2, 0,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 2, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 0, 3, 3, 4, 4, 1, 1, 2, 2,
    3, 3, 3, 3
]

lock = Lock()

def codon_conversion(codon):
    """Converts a codon into a numeric index (0-63)."""
    return base_to_number[codon[1]] * 16 + base_to_number[codon[0]] * 4 + base_to_number[codon[2]]

def syn_site(codon):
    """Returns the number of potential synonymous sites for a codon."""
    return syn_site_array[codon_conversion(codon)]

def calculate_syn_nonsyn(seq1, seq2):
    count_codons = 0
    syn_codons = 0.0
    nonsyn_codons = 0.0
    SA_Nei = 0.0
    SB_Nei = 0.0

    for i in range(0, len(seq1) - 2, 3):
        codonA = seq1[i:i+3]
        codonB = seq2[i:i+3]
        if '-' in codonA or '-' in codonB or 'N' in codonA or 'N' in codonB:
            continue

        count_codons += 1
        SA_Nei += syn_site(codonA)
        SB_Nei += syn_site(codonB)

        codon_index_A = codon_conversion(codonA)
        codon_index_B = codon_conversion(codonB)

        if codonA != codonB:
            if aa_array[codon_index_A] == aa_array[codon_index_B]:
                syn_codons += 1.0
            else:
                nonsyn_codons += 1.0

    potential_syn = (SA_Nei / 3 + SB_Nei / 3) / 2
    potential_nonsyn = 3 * count_codons - potential_syn

    ps = syn_codons / potential_syn if potential_syn > 0 else 0
    pn = nonsyn_codons / potential_nonsyn if potential_nonsyn > 0 else 0

    ds = -3.0 / 4.0 * math.log(1 - 4.0 * ps / 3.0) if ps < 0.75 else 'NA'
    dn = -3.0 / 4.0 * math.log(1 - 4.0 * pn / 3.0) if pn < 0.75 else 'NA'
    ratio = ds / dn if ds != 'NA' and dn != 'NA' and dn != 0 else 'NA'

    return syn_codons, nonsyn_codons, potential_syn, potential_nonsyn, ps, pn, ds, dn, ratio

def compare_and_write(seq1, seq2, name1, name2, output_file):
    """Compares two sequences, calculates dN/dS, and writes the result."""
    result = calculate_syn_nonsyn(seq1, seq2)
    line = f"{name1}\t{name2}\t" + "\t".join(map(str, result)) + "\n"

    # Ensure that only one process writes to the file at a time
    with lock:
        with open(output_file, 'a') as f:
            f.write(line)
            f.flush()  # Ensures the line is written to the disk immediately

    return name1

def process_pairwise_comparisons(names, sequences, output_file, max_workers=4):
    # Initialize output file with headers
    with open(output_file, 'w') as summary:
        summary.write("Seq1\tSeq2\tSd\tSn\tS\tN\tps\tpn\tds\tdn\tds/dn\n")

    processed_names = set()
    with Pool(processes=max_workers) as pool:
        results = []
        for i, j in itertools.combinations(range(len(sequences)), 2):
            name1, name2 = names[i], names[j]
            seq1, seq2 = sequences[i], sequences[j]

            if name1 not in processed_names:
                logging.info(f"Processing comparisons for sequence: {name1}")
                processed_names.add(name1)

            result = pool.apply_async(compare_and_write, (seq1, seq2, name1, name2, output_file))
            results.append(result)

        # Wait for all results to complete
        for result in results:
            result.wait()

def main():
    parser = argparse.ArgumentParser(description='Calculates dN/dS for sequences.')
    parser.add_argument('input_file', help='Input file with sequences.')
    parser.add_argument('--output', default='output.tsv', help='Output file (TSV).')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel processes.')
    args = parser.parse_args()

    # Read sequences
    names, sequences = [], []
    with open(args.input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                names.append(parts[0])
                sequences.append(parts[1])
            else:
                logging.warning(f"Malformed line: {line.strip()}")

    if len(names) != len(sequences):
        raise ValueError("The number of names does not match the number of sequences.")

    process_pairwise_comparisons(names, sequences, args.output, args.workers)
    logging.info("Process completed. Results saved in " + args.output)

if __name__ == '__main__':
    main()
