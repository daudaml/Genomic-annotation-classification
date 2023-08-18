# Genomic Annotation and Classification Pipeline

![Network Image](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/39/4/10.1093_bioinformatics_btad135/1/btad135f1.jpeg?Expires=1695346411&Signature=T7QYCBbUHK0l0vOiQwjo5WPP0Ja9RU0~NSgFnD9tYIXyUfabAXC43cXtY1nYNKf4pinMRI-ZEUB9DebKpqmFb~Zig4TqJUxqScQlgmujvXMb5DSpkxFg00uYTZQFrBW-3DkTSU7t3nnPZE2c0cLmeZ8VijzJFM3Nh-f9dzzRCR7dK0ixxudeIRvB4z37-xpouCsDogbMyKhOpGbEIAs7omDM~KX0rrp5PnZ3UVbL36XWXNYkTosu2cAgL3SJnqxKSaOLXzmPB2yTWEh~tQyhitZP2~mntTdjN1cgUTg8rFuzfIBIpMN52RTexCGYWo92TaAncjDY95he~8rQduAw8g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

This repository contains a Genomic Annotation and Classification Pipeline that utilizes the GRACES (GRAph Convolutional nEtwork feature Selector) algorithm for feature selection. The GRACES algorithm is detailed in the paper: ["Graph convolutional network-based feature selection for high-dimensional and low-sample size data"](https://academic.oup.com/bioinformatics/article/39/4/btad135/7135826), authored by Can Chen, Scott T Weiss, and Yang-Yu Liu.

## Overview

The primary goal of this pipeline is to perform genomic annotation and classification on high-dimensional and low-sample size (HDLSS) data. The pipeline follows these main steps:

1. **Genomic Annotation**: The pipeline starts by annotating genomic data using tools like Prodigal and InterProScan. Prodigal is used for predicting protein-coding genes, and InterProScan is used to assign functional annotations to these genes.

2. **Feature Extraction**: Based on the annotation results, a feature dataset is constructed. The features include information about predicted proteins and their functional annotations.

3. **Feature Selection with GRACES**: The GRACES algorithm, a deep learning-based method, is employed to perform feature selection. GRACES leverages graph convolutional networks and overfitting-reducing techniques to identify the most relevant features for classification.

4. **Classification**: The selected features are then used to train a classification model. In this project, we've performed classification on the pathogenicity of different samples within the Pseudomonas genus.

5. **Accuracy Evaluation and Graph Plotting**: The pipeline calculates the accuracy of the classification model and plots a graph that illustrates the relationship between the number of selected features and the achieved accuracy.

## Tools Used

### Prodigal

Prodigal is a tool for predicting protein-coding genes in DNA sequences. It identifies the location and orientation of genes within the input genomic sequences.

Learn more about Prodigal: [Prodigal GitHub Repository](https://github.com/hyattpd/Prodigal)

### InterProScan

InterProScan is a sequence analysis tool that searches for functional domains, sites, and motifs in protein sequences. It provides insights into the functional characteristics of predicted proteins.

Learn more about InterProScan: [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/)

## Results and Graph

In this project, we ran the genomic annotation and classification pipeline on the Pseudomonas genus to classify pathogenicity in different samples. The achieved accuracy was around 96%. To further understand the relationship between the number of selected features and accuracy, we reran the pipeline for different feature counts and plotted a graph.

![Accuracy Graph](https://github.com/daudaml/Genomic-annotation-classification/blob/main/images/1.jpg?raw=true)

## Acknowledgments

Special thanks to the collaborators from the USDA Agricultural Research Service (ARS) for their support during the development of this pipeline.

## Draft Manuscript

Please note that the manuscript related to this work is currently in the draft writing stage. The relevant citation for the GRACES algorithm is:
```bibtex
@article{chen2023graph,
  title={Graph convolutional network-based feature selection for high-dimensional and low-sample size data},
  author={Chen, Can and Weiss, Scott T and Liu, Yang-Yu},
  journal={Bioinformatics},
  volume={39},
  number={4},
  pages={btad135},
  year={2023},
  publisher={Oxford University Press}
}
```

If you find this Pipeline useful in your research or work, we encourage you to cite this repository.
If you have any questions or suggestions, please open an issue.