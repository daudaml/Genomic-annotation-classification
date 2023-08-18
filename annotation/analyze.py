from glob import glob

import numpy as np
import pandas as pd
from scipy.io import savemat

column_names = ["Protein accession", "Sequence", "length", "Analysis",
                "Signature accession", "signature description", "Start location",
                "Stop location", "Score", "Status", "Date", "InterPro accession",
                "InterPro description"]
data_dict = {"all_proteins": [], "scan": [], "class": [], "bacteria_name": []}
directories = ["interproscan_results/biocontrol", "interproscan_results/pathogen"]
data_class = [0, 1]


def tsv_to_dict(data):
    data_dict = {}
    for ind in data.index:
        key = data["InterPro accession"][ind]
        if key in data_dict:
            data_dict[key] += 1
        else:
            data_dict[key] = 1
    return data_dict


def top_ten_proteins(data):
    keys = list(data.keys())
    values = list(data.values())
    sorted_value_index = np.argsort(values)[::-1]
    for k, i in enumerate(sorted_value_index):
        print(keys[i], values[i])
        if k == 10:
            break


print("TSVs to dict representation")
count = 0
for idx in range(len(directories)):
    directory = directories[idx]
    for file in glob("{}/*.tsv".format(directory)):
        scan_file = pd.read_csv(file, sep='\t', names=column_names)
        protein_dict = tsv_to_dict(scan_file)
        data_dict["all_proteins"].extend(key for key in protein_dict.keys() if key not in data_dict["all_proteins"])
        data_dict["scan"].append(protein_dict)
        data_dict["class"].append(idx)
        data_dict["bacteria_name"].append(file)
        # print(len(data_dict["all_proteins"]), count)
        # count += 1

print("dict to standard data representation")
data_rep = {"X": [], "Y": [], "all_proteins": [], "bacteria_name": []}
n_samples = len(data_dict["scan"])
n_features = len(data_dict["all_proteins"])
x_data = np.zeros((n_samples, n_features), dtype=int)
y_data = np.zeros(n_samples, dtype=int)
print(n_samples)
for i in range(n_samples):
    for j in range(n_features):
        protein_label = data_dict["all_proteins"][j]
        if protein_label in data_dict["scan"][i]:
            x_data[i, j] = data_dict["scan"][i][protein_label]
    y_data[i] = data_dict["class"][i]

data_rep["X"] = x_data
data_rep["Y"] = y_data
data_rep["all_proteins"] = data_dict["all_proteins"]
data_rep["bacteria_name"] = data_dict["bacteria_name"]
savemat("pseudomonas_data.mat", data_rep)
