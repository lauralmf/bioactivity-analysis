rule run_antismash: # Before running, remember to install antiSMASH correctly.
    input:
        "./assembly_files"
    
    output:
        "./antismash_results"

    conda:
        "./envs/antismash.yml"

    shell:
        """
        mkdir -p antismash_results

        for file in {input}/*.fna; do
            base=$(basename "$file" .fasta)  # Extract filename without extension
            antismash "$file" --output-dir {output}/"$base" \
                -t bacteria \
                --genefinding-tool prodigal-m --asf \
                --cc-mibig --cb-knownclusters --pfam2go --rre --tfbs --cpu 16
        done
        """

rule run_gtdbtk: # Before running, remember to install GTDB-Tk correctly.
    input:
        "./assembly_files"
    
    output:
        "./gtdbtk_results"

    conda:
        "./envs/gtdbtk.yml"
    
    shell:
        """
        mkdir -p gtdbtk_results

        gtdbtk classify_wf --genome_dir {input} -x fasta --out_dir {output} --mash_db path/to/mash_db
        """

rule run_bakta: # Before running, remember to install bakta correctly.
    input:
        "./assembly_files"
    
    output:
        "./bakta_results"

    conda:
        "./envs/bakta.yml"

    shell:
        """
        mkdir -p bakta_results

        for file in {input}/*.fasta; do bakta --db db "$file" --output {output} --threads 8; done
        """

rule run_busco: # Before running, remember to install BUSCO correctly.
    input:
        "./assembly_files"
    
    output:
        "./busco_results"

    conda:
        "./envs/busco.yml"

    shell:
        """
        mkdir -p busco_results

        for file in {input}/*.fasta; do busco -i "$file" -m genome --lineage_dataset bacteria_odb12 -c 8 -o {output}; done
        """
    
rule run_busco_phylo:
    input:
        busco_results = "./busco_results"
        assemblies = "./assembly_files"
    
    output:
        phylo = "./busco_phylo_results",
        aln = "./busco_phylo_results/supermatrix/alignments"
        aln_trimmed = "./busco_phylo_results/supermatrix/alignments_trimmed"
    conda:
        pass

    shell:
        """
        mkdir -p busco_phylo_results

        BUSCO_phylogenetics.py -i {input.busco_results} -o {output.phylo} -t 8

        for file in {output.phylo}/supermatrix/sequences/*; do muscle -in muscle "$file" -out {output.aln}/"${file%.faa}"_alignment.fasta; done
        for file in {output.aln}; do trimal -in "$file" -automated1 -out {output.aln_trimmed}"${file%.fasta}"_trimmed.fasta; done

        AMAS.py concat -i {input.assemblies}/*.fasta -f fasta -d aa -p {output.phylo}/supermatrix/concatenated_partition.txt -t {output.phylo/supermatrix/concatenated_alignment.fasta

        fasttree {output.phylo/supermatrix/concatenated_alignment.fasta > {output.phylo}/supermatrix/busco_tree.tree
        """