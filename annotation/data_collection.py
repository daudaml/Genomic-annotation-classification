import os
import shutil
import urllib.request

import pandas as pd
from sh import gunzip
from tqdm import tqdm

data = pd.read_csv('pseudomonas_report_all.csv', sep='\t')
plant_pathogen = data.loc[data["category"] == 3]
bio_control = data.loc[data["category"] == 1]

plant_pathogen.to_csv("plant_pathogen.csv")
bio_control.to_csv("bio_control.csv")


def links_to_download(links_list, save_dir):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    else:
        shutil.rmtree(save_dir)
        os.mkdir(save_dir)
    for i in tqdm(range(len(links_list)), desc="Writing to {}".format(save_dir)):
        try:
            link = links_list[i]
            file_name = link.split("/")[-1] + '_genomic.fna.gz'
            file_link = link + '/' + file_name
            file_path = os.path.join(save_dir, file_name)
            urllib.request.urlretrieve(file_link, file_path)
            gunzip(file_path)
        except:
            print("Something went wrong with this file {}".format(file_name))


links_to_download(list(plant_pathogen["GenBank FTP"]), 'pathogen')
links_to_download(list(bio_control["GenBank FTP"]), 'biocontrol')
