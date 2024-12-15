import os
import re
import pandas as pd

# Chemin vers le dossier contenant les logs
log_dir = "mydirhehe"

# Initialisation des données
data = []

# Parser chaque fichier de log
for log_file in os.listdir(log_dir):
    if log_file.endswith(".log"):
        file_path = os.path.join(log_dir, log_file)
        with open(file_path, "r") as f:
            content = f.read()

            # Extraire les informations principales depuis le contenu du fichier
            dataset_match = re.search(r"Running query: dataset=(.+?), mode=(.+?), k=(\d+), t=(\d+)", content)
            if dataset_match:
                dataset = dataset_match.group(1)
                mode = dataset_match.group(2)
                k = int(dataset_match.group(3))
                t = int(dataset_match.group(4))

                # Extraire le temps d'exécution (Elapsed time, m:ss)
                time_match = re.search(r"Elapsed \(wall clock\) time.*?:\s+(\d+):(\d+\.\d+)", content)
                if time_match:
                    minutes = int(time_match.group(1))
                    seconds = float(time_match.group(2))
                    elapsed_time_formatted = f"{minutes:02}:{seconds:.2f}"
                else:
                    elapsed_time_formatted = "00:00.00"

                # Extraire la mémoire maximale utilisée (Peak Memory)
                mem_match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", content)
                if mem_match:
                    memory = int(mem_match.group(1)) / 1024  # Convertir en MB
                else:
                    memory = None

                # Ajouter les données collectées
                data.append({
                    "Dataset": dataset,
                    "k": k,
                    "t": t,
                    "Mode": mode,
                    "Elapsed Time (mm:ss)": elapsed_time_formatted,
                    "Peak Mem. (MB)": round(memory, 2) if memory else None
                })

# Convertir les données en DataFrame
df = pd.DataFrame(data)

# Appliquer les filtres : garder uniquement k=21 ou k=31, et t=2 ou t=3
df_filtered = df[(df["k"].isin([21, 31])) & (df["t"].isin([2, 3]))]

# Trier les données par Dataset, puis k, t, et enfin Mode
df_filtered = df_filtered.sort_values(by=["Dataset", "k", "t", "Mode"])

# Sauvegarder les données filtrées dans un CSV propre
df_filtered.to_csv("logs_table.csv", index=False)

# Afficher un aperçu des données
print(df_filtered)
