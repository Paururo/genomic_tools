import os
import json
import shutil
import logging
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Root directory containing the folders
root_dir = '/scr/ruizro/foldseek/process_proteins'
root_dir_out = '/scr/ruizro/foldseek/process_proteins/proteins'
# Regex pattern to match the portion to be replaced
pattern = r"_scores_rank_001_alphafold2_ptm_model_[1-5]_seed_000"

# Iterate through the directories in root_dir
for folder_name in os.listdir(root_dir):
    folder_path = os.path.join(root_dir, folder_name)
    
    # Check if it's a directory that starts with "prot_"
    if os.path.isdir(folder_path) and folder_name.startswith("prot_"):
        json_file = None
        pdb_file = None

        # Search for .json and .pdb files in the directory
        for file_name in os.listdir(folder_path):
            if "scores_rank_001" in file_name:
                full_path = os.path.join(folder_path, file_name)
                
                # Identify JSON and PDB file paths
                json_file = full_path
                logging.info(f"JSON file found: {json_file}")
                
                pdb_file = full_path.replace(".json", ".pdb").replace("scores", "unrelaxed")
                logging.info(f"PDB file found: {pdb_file}")

                # Calculate mean pLDDT from JSON file
                with open(json_file) as f:
                    data = json.load(f)
                    plddt_values = data["plddt"]
                    mean_plddt = sum(plddt_values) / len(plddt_values)
                    logging.info(f"Mean pLDDT calculated: {mean_plddt}")

                # Only proceed if mean pLDDT is 70 or higher
                if mean_plddt >= 70:
                    # Build the new name for the PDB file by removing the dynamic portion using regex
                    base_name = re.sub(pattern, "", file_name).replace("-AMAP.panaroo", "").replace(".json", "")
                    new_pdb_name = f"{folder_name.replace('prot_', '')}-{base_name}-{mean_plddt:.2f}.pdb"
                    new_pdb_path = os.path.join(root_dir_out, new_pdb_name)
                    
                    # Copy the renamed PDB file to the desired location
                    shutil.copy(pdb_file, new_pdb_path)
                    logging.info(f"File {pdb_file} copied and renamed to {new_pdb_path}")
                else:
                    logging.info(f"File {pdb_file} not copied as mean pLDDT {mean_plddt:.2f} is below 70")
