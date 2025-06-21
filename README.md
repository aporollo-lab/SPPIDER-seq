# SPPIDER-seq
## Sequence-based partner-aware prediction of PPI interaction sites

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aporollo-lab/SPPIDER-seq/blob/main/notebooks/sppider_seq.ipynb)

This project predicts protein-peptide interaction sites using pretrained models based on ESM-2 embeddings of the query protein and its interaction partner(s).

## 📂 Repository Structure

- `notebooks/` – main analysis notebook
- `models/` – pretrained models
- `data/` – optional input/output examples

## 🚀 Run the Model

Click the Colab badge above to launch the notebook in your browser. The notebook will automatically download the pretrained models from this repository.

Run each cell sequentially. When you reach the **Input Form**, you can either:

- Provide your own protein of interest along with one or more potential interaction partners, or

- Start with the example inputs available in the `data/` folder.

After specifying the sequences, click the `Run PPI predictions` button to initiate the interface prediction process.

## 📦 Requirements

Listed in the notebook

## 📣 Citation

While the manuscript is under preparation, please cite the GitHub version (https://github.com/aporollo-lab/SPPIDER-seq/) if using this tool in your work.
