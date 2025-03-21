import os
from Bio import SeqIO
from os import listdir, path
from shutil import move

def rename_region_gbk(folder_path: str): 
    """Renames .gbk files in a given antiSMASH output folder."""
    assembly = folder_path[20:38]  # Extract the assembly name - change according to name (mine is ./antismash_results/xxx.barcodeXX_assembly)

    for file in listdir(folder_path):
        if ".region" in file:  # Process only region files
            print(f"\nProcessing file: {file}")  # Debugging output
            file_path = path.join(folder_path, file)  # Original file path

            if not path.exists(file_path):  
                print(f"⚠️ File not found (skipping): {file_path}")
                continue  # Avoid FileNotFoundError

            try:
                with open(file_path, "r") as gbk:
                    records = list(SeqIO.parse(gbk, "genbank"))  # Convert iterator to list
                    
                    if not records:  # If parsing failed or file is empty
                        print(f"⚠️ Empty or invalid GenBank file (skipping): {file}")
                        continue

                    for record in records:
                        for feature in record.features:
                            if feature.type == "protocluster":
                                # Get 'product' qualifier, default to 'N/A' if not found
                                product = feature.qualifiers.get("product", ["N/A"])[0]

                                # Sanitize product name for safe filenames
                                product = product.replace(" ", "_").replace("/", "_")

                                # Ensure filename format is correct
                                filename_parts = file.split(".")
                                if len(filename_parts) < 4:
                                    print(f"⚠️ Unexpected filename format (skipping): {file}")
                                    continue  # Avoid IndexError

                                contig = filename_parts[0]
                                region = filename_parts[1]

                                # Construct new filename
                                newfilename = f"{product}_{assembly}_{contig}_{region}.gbk"
                                new_file_path = path.join(folder_path, newfilename)

                                # Check if the target file already exists
                                if path.exists(new_file_path):
                                    print(f"⚠️ Target filename already exists (skipping): {newfilename}")
                                    continue

                                # Rename file
                                move(file_path, new_file_path)
                                print(f"✅ Renamed: {file} -> {newfilename}")

            except Exception as e:
                print(f"❌ Error processing {file}: {e}")

def process_all_folders(base_path: str):
    """Loops through all subdirectories in 'antismash_results' and processes them."""
    for folder in listdir(base_path):
        if folder == ".DS_Store":  # Skip system files (MacOS)
            continue

        folder_path = path.join(base_path, folder)
        if path.isdir(folder_path):  # Only process directories
            print(f"\n--- Processing Folder: {folder} ---")
            rename_region_gbk(folder_path)

# Example Usage
# process_all_folders("path/to/antismash_results")
