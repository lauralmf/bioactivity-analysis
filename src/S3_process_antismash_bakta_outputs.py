from os import listdir, rename, path
import pandas as pd
from Bio import SeqIO
from bs4 import BeautifulSoup
import csv
from collections import defaultdict
import json

""" antiSMASH """
# antiSMASH results path
antismash_results = "../antismash_results"

# Trim folder names from antiSMASH raw output
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

# Rename antiSMASH gbk regions to include sequencing ID and BGC class
for assembly in listdir(antismash_results):
    assembly_path = f"{antismash_results}/{assembly}"
    if assembly == ".DS_Store": # Skip system files (MacOS thing)
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

# Get all BGC annotations from the antiSMASH JSON output file
def get_bgcs_from_json(assembly, assembly_path):
    df = {
        "SEQID": [],
        "Contig": [],
        "Region": [],
        "Category": [],
        "Number of hits": [],
        "Most similar known cluster": [],
        "Similarity": []
    }

    for file in listdir(assembly_path):
        if file.endswith("assembly.json"):
            file_path = f"{assembly_path}/{file}"
            with open(file_path) as f:
                data = json.load(f)
                records = data.get("records")

                for i in range(len(records)):
                    record = records[i]
                    contig = record["id"]
                    module = record["modules"]

                    for contig_region in record["features"]:
                        if contig_region["type"] == "region":
                            qualifier = contig_region["qualifiers"]
                            region = qualifier["region_number"]
                            p = qualifier["product"]
                            df["SEQID"].append(assembly)    # Get SEQID
                            df["Contig"].append(contig)     # Get BGC contig number
                            df["Region"].append(region[0])  # Get BGC region number
                            df["Category"].append(p[0])     # Get BGC category

                    for modkey in module.keys():

                        # Getting ClusterBlast hits (all BGC-like regions):
                        if modkey == "antismash.modules.clusterblast":
                            line = module[modkey]

                            c = line["knowncluster"]
                            c_results = c["results"]

                            for region_result in c_results:
                                ranks = region_result["ranking"]
                                n_hits = region_result["total_hits"]
                                df["Number of hits"].append(n_hits)

                                bgc_list = []
                                similarity_list = []

                                for l in ranks:
                                    for hits_dict in l:
                                        if "description" in hits_dict:
                                            bgc = hits_dict["description"]
                                            bgc_list.append(bgc)

                                        if "similarity" in hits_dict:
                                            similarity = hits_dict["similarity"]/100
                                            similarity_list.append(similarity)

                                refs_dict = {bgc_list[i]: similarity_list[i] for i in range(len(bgc_list))}
                                if len(refs_dict):
                                    best_bgc = max(refs_dict, key=refs_dict.get)
                                    best_bgc_similarity = refs_dict[best_bgc]

                                    df["Most similar known cluster"].append(best_bgc)
                                    df["Similarity"].append(best_bgc_similarity)
                                
                                if not len(refs_dict):
                                    df["Most similar known cluster"].append(None)
                                    df["Similarity"].append(0)

    return pd.DataFrame(df)

def get_all_json_bgcs(antismash_results_path):
    master_df = []

    for assembly in listdir(antismash_results_path):
        if assembly != ".DS_Store": # MacOS thing
            assembly_path = f"{antismash_results_path}/{assembly}"
            df = get_bgcs_from_json(assembly, assembly_path)
            master_df.append(df)

    master_df = pd.concat(master_df, ignore_index=True)

    # Replace NAs for missing most similar known clusters with category
    master_df["Most similar known cluster"] = master_df["Most similar known cluster"].fillna(master_df["Category"])

    bgc_annos_path = "../antismash_results/all_bgc_annotations_from_json.csv" # Desired output path
    return pd.DataFrame.to_csv(master_df, bgc_annos_path, index=False)        # Write to csv file to desired output path

""" bakta """
# Convert a bakta whole-genome annotation in .gbk or .embl format to a csv with gene/product info per contig and totals:
def bakta_to_csv(input_file, output_file, file_format="embl"):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        headers = ["SEQID", "Gene", "Product"]
        writer.writerow(headers)

        global_totals = defaultdict(int)

        for record in SeqIO.parse(infile, file_format):
            seqid = record.id

            for feature in record.features:
                if feature.type == "CDS":
                    gene = feature.qualifiers.get("gene", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]

                    writer.writerow([seqid, gene, product])

                    # Optional: keep number of total genes if needed
                    global_totals["total_genes"] += 1
                    if gene:
                        global_totals[gene] += 1

# Apply bakta_to_csv to whole-genome annotations for all assemblies:
for assembly in listdir("../bakta_results/"):
    embl_file = f"../bakta_results/{assembly}/{assembly}_assembly.embl"
    if embl_file == "../bakta_results/.DS_Store/.DS_Store.embl": # MacOS thing
        continue
    bakta_to_csv(embl_file, f"../bakta_results/{assembly}/{assembly}_bakta.csv", file_format="embl")

# Merge all bakta annotations into a single csv file:
master = []
for assembly in listdir("../bakta_results"):
    assembly_path = f"../bakta_results/{assembly}"
    if assembly == ".DS_Store":
        continue
    
    for file in listdir(assembly_path):
        if file.endswith("bakta.csv"):
            bakta_csv_path = f"{assembly_path}/{file}"
            bakta_df = pd.read_csv(bakta_csv_path)
            bakta_df["SEQID"] = assembly
            master.extend(bakta_df[["SEQID", "Gene", "Product"]].to_dict(orient="records")) # Append results for a single assembly to master list

master_df = pd.DataFrame(master) # Convert all results to pandas dataframe

bakta_annos_path = "../bakta_results/all_bakta_annos.csv"     # Desired output path
pd.DataFrame.to_csv(master_df, bakta_annos_path, index=False) # Write results to csv to desired path
