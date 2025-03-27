from os import listdir, rename, path
import pandas as pd
from Bio import SeqIO

antismash_results = "../antismash_results"

""" Clean up folder names from antiSMASH raw output: remove '_assembly' """
for assembly in listdir(antismash_results):
    if assembly == ".DS_Store": # Skip system files (MacOS thing)
        continue
    if assembly.endswith("_assembly"):
        new_name = assembly.replace("_assembly", "")
        current_path = path.join(antismash_results, assembly)
        new_path = path.join(antismash_results, new_name)

        if path.exists(new_path):
            continue
        else:
            rename(current_path, new_path)

""" Modify genbank region filenames to avoid duplicate names when running BiG-SCAPE """
for assembly in listdir(antismash_results):
    assembly_path = f"{antismash_results}/{assembly}"
    if assembly in [".DS_Store", "all_bgc_annotations.csv"]:
        continue

    for file in listdir(assembly_path):
        if ".region" in file:
            new_name = f"{assembly}_{file}"
            current_path = path.join(assembly_path, file)
            new_path = path.join(assembly_path, new_name)
            
            if path.exists(new_path):
                continue
            else:
                rename(current_path, new_path)

""" Extract annotations from each BGC predicted by antiSMASH """
def extract_bgc_annos(assembly_path: str, assembly: str):
    all_annos = []  # Master list to store annotations from all GenBank files

    for file in listdir(assembly_path):
        if ".region" in file:
            file_path = f"{assembly_path}/{file}"
            
            with open(file_path, "r") as gbk:
                records = SeqIO.parse(gbk, "genbank")

                for record in records:
                    for feature in record.features:
                        if feature.type == "protocluster":
                            tmp_dict = {}
                            tmp_dict["Position"] = file.split(".gbk")[0]

                            number = feature.qualifiers.get("protocluster_number").pop()
                            category = feature.qualifiers.get("category").pop()
                            product = feature.qualifiers.get("product").pop()

                            tmp_dict["Protocluster"] = number
                            tmp_dict["Category"] = category
                            tmp_dict["Product"] = product
                            
                            references = []
                            for q in feature.qualifiers:
                                if q.endswith("product_classes"):
                                    references = feature.qualifiers.get(q, [])  # Get reference(s) if found
                            
                            if not references:
                                references = [product]  # Default to product if no reference found

                            for ref in references:
                                new_entry = tmp_dict.copy()
                                new_entry["Reference"] = ref
                                new_entry["Assembly"] = assembly
                                all_annos.append(new_entry)

    if all_annos:
        df = pd.DataFrame(all_annos)

        column_order = ["Assembly", "Position", "Protocluster", "Category", "Product", "Reference"]
        df = df[column_order]

        # Save individual assembly BGC annotations
        output_path = f"{assembly_path}/{assembly}_bgc_annotations.csv"
        df.to_csv(output_path, index=False)

        return df
    
    return None

def process_all_folders(antismash_results_path: str):
    """Loops through all subdirectories in 'antismash_results' and processes them."""
    all_dfs = []  # List to collect all individual DataFrames

    for assembly in listdir(antismash_results_path):
        if assembly == ".DS_Store":  
            continue
        
        assembly_path = f"{antismash_results_path}/{assembly}"
        df = extract_bgc_annos(assembly_path=assembly_path, assembly=assembly)

        if df is not None:
            all_dfs.append(df)

    # Merge all DataFrames into one
    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
        merged_df.to_csv(f"{antismash_results_path}/all_bgc_annotations.csv", index=False)
        return merged_df

    return "No BGC annotations found in any assembly!"

# Run the processing
merged_annotations = process_all_folders(antismash_results)
