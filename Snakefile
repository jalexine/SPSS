import yaml

config = yaml.safe_load(open("config.yaml"))
FASTA_FILES = [f.split("/")[-1].replace(".fasta", "") for f in config["fasta_files"]]
K_VALUES = config["k"]
T_VALUES = config["t"]
MODES = ["simplitig", "unitig"]

rule all:
    input:
        expand("benchmark/{dataset}/fmi_{dataset}_{mode}_k{k}_t{t}.dump",
               dataset=FASTA_FILES, mode=MODES, k=K_VALUES, t=T_VALUES),
        expand("benchmark/stats/stats_{dataset}_{mode}.csv",
               dataset=FASTA_FILES, mode=MODES),
        expand("benchmark/res_query/{dataset}/query_{mode}_k{k}_t{t}.tsv",
               dataset=FASTA_FILES, mode=MODES, k=K_VALUES, t=T_VALUES),
        "benchmark/plots/spss_reads_05x.png",
        "benchmark/plots/threshold_reads_80x.png"



rule generate_dump:
    input:
        fasta=lambda wildcards: f"test-data/{wildcards.dataset}.fasta"
    output:
        "benchmark/{dataset}/fmi_{dataset}_{mode}_k{k}_t{t}.dump"
    shell:
        """
        mkdir -p $(dirname {output})
        python src/sequences_to_indexed_spss.py -i {input.fasta} -k {wildcards.k} -t {wildcards.t} -m {wildcards.mode} \
            -q {output}
        """
rule generate_stats:
    input:
        fasta=lambda wildcards: f"test-data/{wildcards.dataset}.fasta"
    output:
        "benchmark/stats/stats_{dataset}_{mode}.csv"
    shell:
        """
        mkdir -p $(dirname {output})
        for k in {K_VALUES}; do
            for t in {T_VALUES}; do
                python src/sequences_to_indexed_spss.py -i {input.fasta} -k $k -t $t -m {wildcards.mode} -stats {output}
            done
        done
        """


rule query_fmi:
    input:
        index="benchmark/{dataset}/fmi_{dataset}_{mode}_k{k}_t{t}.dump",
        query="test-data/queries.fasta"
    output:
        "benchmark/res_query/{dataset}/query_{mode}_k{k}_t{t}.tsv"
    log:
        "benchmark/res_query/logs/query_{dataset}_{mode}_k{k}_t{t}.log"
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        echo "Running query: dataset={wildcards.dataset}, mode={wildcards.mode}, k={wildcards.k}, t={wildcards.t}" > {log}
        echo "Resources Used:" >> {log}
        
        # Utilise gtime au lieu de time
        (
            gtime -v python src/query_indexed_spss.py \\
                -q {input.query} \\
                -i {input.index} \\
                -k {wildcards.k} \\
                -o {output}
        ) 2>&1 | tee -a {log}
        
        printf "Disk usage:\t" >> {log}
        du -sh "{output}" >> {log}
        
        printf "Hostname:\t" >> {log}
        hostname >> {log}
        """


rule plot_spss:
    input:
        unitig="benchmark/stats/stats_{dataset}_unitig.csv",
        simplitig="benchmark/stats/stats_{dataset}_simplitig.csv"
    output:
        "benchmark/plots/spss_{dataset}.png"
    shell:
        """
        mkdir -p benchmark/plots
        python benchmark/scripts/plot_spss.py {input.unitig} {input.simplitig} {output}
        """


rule plot_threshold:
    input:
        simplitig="benchmark/stats/stats_{dataset}_simplitig.csv"
    output:
        "benchmark/plots/threshold_{dataset}.png"
    shell:
        """
        mkdir -p benchmark/plots
        python benchmark/scripts/plot_threshold.py {input.simplitig} {output}
        """
