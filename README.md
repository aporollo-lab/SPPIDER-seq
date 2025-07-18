# SPPIDER-seq
## Sequence-based partner-aware prediction of PPI interaction sites

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/sppider_seq.ipynb)

This project predicts protein-peptide interaction sites using pretrained models based on ESM-2 embeddings of the query protein and its interaction partner(s).

## 📂 Repository Structure

- `notebooks/` – main analysis notebook and utility scripts
- `models/` – pretrained models
- `data/` – optional input/output examples

## 🚀 Run the Model

Click the Colab badge above to launch the notebook in your browser. The notebook will automatically download the pretrained models from this repository.

Run each cell sequentially. When you reach the **Input Form**, you can either:

- Provide your own protein of interest along with one or more potential interaction partners, or  
- Start with the example inputs available in the `data/` folder.

After specifying the sequences, click the **Run PPI predictions** button to initiate the interface prediction process.

Once the prediction is complete, the results can be reviewed directly within the notebook using interactive visualizations and structured tables. Additionally, all output files, including raw probability scores and metadata, can be downloaded locally for downstream analysis or record-keeping.

Optionally, users may choose to assess the statistical significance of predicted interaction sites by computing a null-background distribution of scores based on scrambled partner sequences. This step helps identify regions with significantly elevated interface probabilities beyond random expectations. However, this analysis is computationally intensive and may be time-consuming, especially for longer proteins or large sets of interaction partners.


## 🧬 FASTA Sequence Utilities

SPPIDER-seq includes a suite of utilities for FASTA sequence manipulation to aid in hypothesis testing and sensitivity analyses:

### 1. Sequence Scrambling
Randomly shuffles each protein sequence while preserving amino acid composition. Useful for generating negative controls to assess prediction specificity.

[![Launch in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/scramble_fasta.ipynb)

### 2. Sequence Alanine Scanning
Introduces windows of alanine residues in a sliding manner across the query sequence to assess the impact of local disruptions on predicted interface probabilities.

[![Launch in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/alanine_scanning.ipynb)

### 3. Sequence Signal Amplification
Appends selected subsequences from the query to its C-terminus to identify regions that may amplify interaction signals in downstream prediction.

[![Launch in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/amplify_sequence_signal.ipynb)

## 📦 Requirements

Listed in the main notebook.

## 📣 Citation

While the manuscript is under preparation, please cite the GitHub version if using this tool in your work:

A. Porollo, O. Jadhav, A. Alvarez, J. Chen
SPPIDER-seq: Sequence-based partner-aware predictor of protein-protein interaction sites.
https://github.com/aporollo-lab/SPPIDER-seq/
