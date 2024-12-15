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
        expand("benchmark/plots/spss_{mode}.png", mode=MODES),
        expand("benchmark/plots/threshold_05x_{mode}.png", mode=MODES),
        expand("benchmark/plots/threshold_30xs_{mode}.png", mode=MODES)






rule generate_dump:
    input:
        fasta=lambda wildcards: f"test-data/{wildcards.dataset}.fasta"
    output:
        "benchmark/{dataset}/fmi_{dataset}_{mode}_k{k}_t{t}.dump"
    shell:
        """
        mkdir -p $(dirname {output})
        python src/sequences_to_indexed_spss.py -i {input.fasta} -k {wildcards.k} -t {wildcards.t} -m {wildcards.mode} \
            -o {output}
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
        small ="benchmark/stats/stats_reads_05x_{mode}.csv",
        big ="benchmark/stats/stats_reads_80x_{mode}.csv"
    output:
        "benchmark/plots/spss_{mode}.png"
    shell:
        """
        mkdir -p benchmark/plots
        python benchmark/scripts/plot_spss.py {input.small} {input.big} {output}
        """


rule plot_threshold_05x:
    input:
        data ="benchmark/stats/stats_reads_05x_{mode}.csv",
    output:
        "benchmark/plots/threshold_05x_{mode}.png"
    shell:
        """
        mkdir -p benchmark/plots
        python benchmark/scripts/plot_threshold_05x.py {input.data} {output}
        """

rule plot_threshold_30x:
    input:
        data ="benchmark/stats/stats_reads_30x_{mode}.csv",
        datapoly ="benchmark/stats/stats_reads_30x_poly_{mode}.csv"
    output:
        "benchmark/plots/threshold_30xs_{mode}.png"
    shell:
        """
        mkdir -p benchmark/plots
        python benchmark/scripts/plot_threshold_30x.py {input.data} {input.datapoly} {output}
        """
