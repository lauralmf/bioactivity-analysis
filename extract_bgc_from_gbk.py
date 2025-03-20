from os import listdir
import pandas as pd
from Bio import SeqIO

# For a single antiSMASH output folder
def extract_bgc_from_gbk(path: str, dictionary: dict):
    Category = []
    Product = []
    Assembly = path[20:38]
    for file in listdir(path):
        if ".region" in file:
            file_path = f"{path}/{file}"
            with open(file_path, "r") as gbk:
                records = SeqIO.parse(gbk, "genbank")
                for record in records:
                    for feature in record.features:
                        if feature.type == "protocluster":
                            category = feature.qualifiers.get("category", ["N/A"])[0]
                            Category.append(category)

                            product = feature.qualifiers.get("product", ["N/A"])[0]
                            Product.append(product)
                            
    dictionary["Assembly"].extend([Assembly] * len(Product))
    dictionary["Category"].extend(Category)
    dictionary["Product"].extend(Product)
    
    return dictionary

# For many antiSMASH output folders
def process_all_folders(base_path: str, dictionary: dict):
    """Loop through all subdirectories in 'antismash_results' and process them."""
    for folder in listdir(base_path):
        if folder == ".DS_Store":
            continue

        folder_path = f"{base_path}/{folder}"
        results = extract_bgc_from_gbk(folder_path, dictionary)
    return results

d = {"Assembly": [], "Category": [], "Product": []} # Should always be used!

antismash_folder_path = "path/to/antismash_results_folder"
bgc_output_csv = "path/to/output_file/output_file.csv"

# Example use
antismash_results = process_all_folders(antismash_folder_path, d)
antismash_results = pd.DataFrame.from_dict(antismash_results)

pd.DataFrame.to_csv(antismash_results, bgc_output_csv, index=False)