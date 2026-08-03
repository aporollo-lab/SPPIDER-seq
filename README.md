# SPPIDER-seq

## Sequence-based partner-aware prediction of protein–protein interaction sites

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/sppider_seq.ipynb)

**SPPIDER-seq** is a deep learning framework for **partner-aware prediction of protein–protein interaction (PPI) sites directly from amino acid sequences**. Unlike conventional interface predictors that estimate a residue's general propensity to participate in protein interactions, SPPIDER-seq predicts **which residues are likely to bind a specific interaction partner**, enabling context-specific interface prediction without requiring experimental or predicted protein structures.

The method combines **ESM-2 protein language model embeddings** with a **cross-attention neural network** that jointly models the query protein and its interaction partner. By explicitly learning partner-dependent sequence representations, SPPIDER-seq is particularly well suited for identifying interaction interfaces within **intrinsically disordered regions (IDRs)** and **short linear motifs (SLiMs)**, where partner context is often essential for determining binding specificity.

<p align="center">
  <img src="docs/figures/workflow.png" width="900">
  <br>
  <em>Overview of the SPPIDER-seq workflow.</em>
</p>


## 📂 Repository Structure

- `notebooks/` – Google Colab notebooks providing the main SPPIDER-seq workflow and sequence analysis utilities
- `cli/` – command-line interface and supporting scripts for local and batch execution
- `docs/` – documentation assets and figures
- `examples/` – optional input/output examples

The current SPPIDER-seq resources are maintained on Hugging Face:

- **Pretrained models:** https://huggingface.co/aporollo-lab/SPPIDER-seq
- **Training and benchmark datasets:** https://huggingface.co/datasets/aporollo-lab/SPPIDER-seq-datasets

The Google Colab notebook and the command-line version of SPPIDER-seq automatically download the current pretrained receptor-centric and peptide-centric models from Hugging Face. No manual model download or configuration is required for standard use. Custom model paths only need to be specified when using user-trained models or models retained from previous SPPIDER-seq releases.


## 🚀 Running SPPIDER-seq

### Option 1: Google Colab (recommended for interactive use)

Click the Colab badge above to launch the notebook in your browser.
Change runtime type to T4 GPU (or any other GPU if available).

Run each cell sequentially. When you reach the **Input Form**, you can either:

- Provide your own protein of interest along with one or more interaction partners, or  
- Use example inputs from the `examples/` folder  

After specifying sequences, click **Run PPI predictions** to initiate inference.

Results are displayed as:
- PPI site probability plots  
- structured tables  
- downloadable output files (zipped png and tsv files)

Optionally, users may enable statistical significance estimation using a null-background distribution based on scrambled partner sequences. This helps identify residues with statistically significant prediction probabilities. Note that this step can be computationally intensive, particularly when executed on a CPU runtime.

---

### Option 2: Command-Line Interface (CLI)

SPPIDER-seq also provides CLI scripts for batch processing and reproducibility.

To view all available options, run:
```
python cli/sppider_seq.py -h
```

**Required arguments**

At minimum, the following inputs must be provided:
```
--query path/to/query.fasta \
--partners path/to/partners.fasta \
--output-dir path/to/outputs
```

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


## ⚙️ Environment and Requirements

SPPIDER-seq was developed and tested using Python 3.12.* with PyTorch and HuggingFace Transformers.

### Option 1: Google Colab

No manual installation of libraries is required. The notebook installs dependencies automatically and runs either on GPU- or CPU-enabled runtime.

### Option 2: Conda environment

Create an environment using:

```
conda env create -f environment.yml
conda activate sppider-seq
```

### CUDA / GPU notes
- The model benefits significantly from GPU acceleration
- Tested on NVIDIA GPUs (e.g., T4, A100) with CUDA 11.x
- CPU execution is supported but slower


## 📣 Citation

If you use **SPPIDER-seq** in your research, please cite:

> Porollo A, Jadhav O, Alvarez A, Chen J. **SPPIDER-seq: sequence-based partner-aware predictor of protein-protein interaction sites.** *Bioinformatics*. 2026;42(7):btag388. https://doi.org/10.1093/bioinformatics/btag388

