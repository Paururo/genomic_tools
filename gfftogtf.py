#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to convert GFF files to GTF format, ensuring the inclusion of necessary attributes,
adding exon entries for each CDS in prokaryotic annotations, and removing lines with unknown strand orientation ("?").
Enhanced to add attributes such as gene_name, product, Dbxref, protein_id, and locus_tag when available.
Handles cases where ID is 'nan' by generating unique IDs based on functional categories.
Additionally, for exon features, the end coordinate is reduced by 3.
"""

import argparse
import logging
from collections import defaultdict
import re

# Configure logging
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger("gff_to_gtf_enhanced")

def parse_attributes(attributes_str):
    """Parse the attributes column of GFF/GTF and return a dictionary."""
    attrs = {}
    for attr in attributes_str.split(";"):
        if attr.strip():
            key_value = attr.strip().split("=")
            if len(key_value) == 2:
                key, value = key_value
                attrs[key] = value
    return attrs

def sanitize_string(s):
    """Sanitize a string to be used in IDs by replacing non-alphanumeric characters with underscores."""
    return re.sub(r'\W+', '_', s)

def gff_to_gtf_with_transcripts(input_gff, output_gtf):
    """
    Convert a GFF file to GTF format, ensuring the presence of 'gene_id' and 'transcript_id',
    adding exon entries to match each CDS entry for prokaryotic GTF compatibility.
    Generates unique IDs when 'ID' is 'nan'.
    Additionally, reduces the end coordinate by 3 for exon features.
    
    Args:
        input_gff (str): Path to the input GFF file.
        output_gtf (str): Path to the output GTF file.
    """
    # Data structure to hold transcript information
    transcripts = defaultdict(lambda: {
        "gene_id": None,
        "transcript_id": None,
        "exons": [],
        "cds": []
    })

    # Counter for unique ID generation
    unique_id_counter = defaultdict(int)

    with open(input_gff, 'r') as infile:
        for line_number, line in enumerate(infile, 1):
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            if len(columns) != 9:
                logger.warning(f"Line {line_number}: Skipping malformed line (expected 9 columns).")
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = columns

            # Filter out lines with unknown strand
            if strand == "?":
                logger.info(f"Line {line_number}: Skipping line with unknown strand '?'.")
                continue

            attr_dict = parse_attributes(attributes)

            # Check for 'ID' attribute
            if "ID" not in attr_dict:
                logger.warning(f"Line {line_number}: Missing 'ID' attribute. Line skipped.")
                continue

            original_id = attr_dict["ID"]

            # Generate unique ID if ID is 'nan'
            if original_id.lower() == "nan":
                # Attempt to use 'product' or 'Name' for functional category
                functional_category = attr_dict.get("product") or attr_dict.get("Name") or "unknown_function"
                sanitized_category = sanitize_string(functional_category)
                unique_id_counter[sanitized_category] += 1
                unique_id = f"{sanitized_category}_{unique_id_counter[sanitized_category]}"
                logger.info(f"Line {line_number}: Generated unique ID '{unique_id}' based on functional category.")
            else:
                unique_id = original_id

            gene_id = unique_id
            transcript_id = unique_id  # Assuming one transcript per gene; modify if multiple transcripts per gene

            if feature == "transcript":
                transcripts[transcript_id]["gene_id"] = gene_id
                transcripts[transcript_id]["transcript_id"] = transcript_id
            elif feature in {"CDS", "exon", "ncRNA-region"}:
                parent = attr_dict.get("Parent", gene_id)
                # If 'Parent' is 'nan', use the generated unique_id
                if parent.lower() == "nan":
                    parent = unique_id
                    logger.info(f"Line {line_number}: 'Parent' attribute was 'nan', set to '{parent}'.")

                transcripts[parent]["gene_id"] = gene_id
                transcripts[parent]["transcript_id"] = transcript_id

                feature_info = {
                    "seqname": seqname,
                    "source": source,
                    "feature": feature,
                    "start": start,
                    "end": end,
                    "score": score,
                    "strand": strand,
                    "frame": frame,
                    "attributes": attributes
                }

                if feature == "CDS":
                    transcripts[parent]["cds"].append(feature_info)
                    # Add a corresponding exon entry for each CDS to ensure compatibility
                    exon_info = feature_info.copy()
                    exon_info["feature"] = "exon"
                    exon_info["frame"] = "."  # Exons typically don't have frame information
                    transcripts[parent]["exons"].append(exon_info)
                elif feature == "exon" or feature == "ncRNA-region":
                    exon_info = feature_info.copy()
                    exon_info["frame"] = "."  # Exons typically don't have frame information
                    transcripts[parent]["exons"].append(exon_info)
            else:
                pass

    with open(output_gtf, 'w') as outfile:
        for transcript_id, info in transcripts.items():
            if not info["gene_id"] or not info["transcript_id"]:
                logger.warning(f"Transcript {transcript_id} is missing 'gene_id' or 'transcript_id'. Skipped.")
                continue

            # Determine transcript boundaries
            all_starts = [int(entry["start"]) for entry in info["exons"] + info["cds"]]
            all_ends = [int(entry["end"]) for entry in info["exons"] + info["cds"]]
            transcript_start = min(all_starts)
            transcript_end = max(all_ends)

            # Prepare transcript attributes
            transcript_attrs = f'gene_id "{info["gene_id"]}"; transcript_id "{info["transcript_id"]}";'

            # Write transcript line
            transcript_line = (
                f'{info["exons"][0]["seqname"]}\t'
                f'{info["exons"][0]["source"]}\t'
                f'transcript\t{transcript_start}\t{transcript_end}\t.\t'
                f'{info["exons"][0]["strand"]}\t.\t{transcript_attrs}\n'
            )
            outfile.write(transcript_line)

            # Write exon lines with end coordinate reduced by 3
            for exon in sorted(info["exons"], key=lambda x: int(x["start"])):
                # Calcular el nuevo end restando 3
                original_end = int(exon["end"])
                adjusted_end = original_end - 3
                # Asegurar que el end no sea menor que start
                if adjusted_end < int(exon["start"]):
                    logger.warning(
                        f"Transcript {transcript_id}: Exon end ({adjusted_end}) is less than start ({exon['start']}). "
                        f"Setting end to start."
                    )
                    adjusted_end = int(exon["start"])
                # Convert back a string
                adjusted_end_str = str(adjusted_end)

                exon_attrs = f'gene_id "{info["gene_id"]}"; transcript_id "{info["transcript_id"]}";'
                exon_attr_dict = parse_attributes(exon["attributes"])
                if "Name" in exon_attr_dict and exon_attr_dict["Name"].lower() != "nan":
                    exon_attrs += f' gene_name "{exon_attr_dict["Name"]}";'
                elif "product" in exon_attr_dict and exon_attr_dict["product"].lower() != "nan":
                    exon_attrs += f' gene_name "{exon_attr_dict["product"]}";'
                if "product" in exon_attr_dict and exon_attr_dict["product"].lower() != "nan":
                    exon_attrs += f' product "{exon_attr_dict["product"]}";'
                if "Dbxref" in exon_attr_dict and exon_attr_dict["Dbxref"].lower() != "nan":
                    exon_attrs += f' Dbxref "{exon_attr_dict["Dbxref"]}";'
                if "protein_id" in exon_attr_dict and exon_attr_dict["protein_id"].lower() != "nan":
                    exon_attrs += f' protein_id "{exon_attr_dict["protein_id"]}";'
                if "locus_tag" in exon_attr_dict and exon_attr_dict["locus_tag"].lower() != "nan":
                    exon_attrs += f' locus_tag "{exon_attr_dict["locus_tag"]}";'

                # Exons típicamente no tienen frame, se establece como '.'
                exon_line = (
                    f'{exon["seqname"]}\t'
                    f'{exon["source"]}\t'
                    f'exon\t{exon["start"]}\t{adjusted_end_str}\t{exon["score"]}\t'
                    f'{exon["strand"]}\t.\t{exon_attrs}\n'
                )
                outfile.write(exon_line)

            # Write CDS lines
            for cds in sorted(info["cds"], key=lambda x: int(x["start"])):
                cds_attrs = f'gene_id "{info["gene_id"]}"; transcript_id "{info["transcript_id"]}";'
                cds_attr_dict = parse_attributes(cds["attributes"])
                if "Name" in cds_attr_dict and cds_attr_dict["Name"].lower() != "nan":
                    cds_attrs += f' gene_name "{cds_attr_dict["Name"]}";'
                elif "product" in cds_attr_dict and cds_attr_dict["product"].lower() != "nan":
                    cds_attrs += f' gene_name "{cds_attr_dict["product"]}";'
                if "product" in cds_attr_dict and cds_attr_dict["product"].lower() != "nan":
                    cds_attrs += f' product "{cds_attr_dict["product"]}";'
                if "Dbxref" in cds_attr_dict and cds_attr_dict["Dbxref"].lower() != "nan":
                    cds_attrs += f' Dbxref "{cds_attr_dict["Dbxref"]}";'
                if "protein_id" in cds_attr_dict and cds_attr_dict["protein_id"].lower() != "nan":
                    cds_attrs += f' protein_id "{cds_attr_dict["protein_id"]}";'
                if "locus_tag" in cds_attr_dict and cds_attr_dict["locus_tag"].lower() != "nan":
                    cds_attrs += f' locus_tag "{cds_attr_dict["locus_tag"]}";'

                cds_line = (
                    f'{cds["seqname"]}\t'
                    f'{cds["source"]}\t'
                    f'CDS\t{cds["start"]}\t{cds["end"]}\t{cds["score"]}\t'
                    f'{cds["strand"]}\t{cds["frame"]}\t{cds_attrs}\n'
                )
                outfile.write(cds_line)

    logger.info(f"Conversion completed: {input_gff} -> {output_gtf}")

def main():
    parser = argparse.ArgumentParser(description="Convert GFF to GTF with transcript, exon, and CDS annotations.")
    parser.add_argument("-i", "--input", required=True, help="Input GFF file.")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file.")

    args = parser.parse_args()
    gff_to_gtf_with_transcripts(args.input, args.output)

if __name__ == "__main__":
    main()
