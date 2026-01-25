# Dataset File Description

This directory contains curated protein–protein interaction (PPI) datasets used for training, validation, and benchmarking of partner-aware PPI site prediction models.

---

## 1. Dataset prefixes and model usage

Files with different prefixes correspond to different model variants:

* Files starting with the prefix **`rec_`** were used for training, validation, and benchmarking of the **receptor-centric model**.
* Files starting with the prefix **`pep_`** were used for training, validation, and benchmarking of the **peptide-centric model**.

---

## 2. Pair definition files (`.tsv`)

Files with the `.tsv` extension define protein–peptide pairs used for different stages of model development:

* `*_train.tsv` — training dataset
* `*_val.tsv` — validation dataset used during training
* `*_blind.tsv` — blind benchmark dataset not seen during training or validation

Each `.tsv` file contains **exactly two columns**:

1. Receptor protein sequence
2. Peptide / ligand protein sequence

Each row corresponds to a single receptor–ligand interaction pair.

---

## 3. Sequence naming convention

Protein sequences in the `.tsv` files follow the naming convention:

```
PDBcode_ChainID|UniProtID|ChainPair
```

where:

* **PDBcode**
  Indicates the PDB entry from which the interacting protein pair was derived.

* **ChainID**
  Specifies the chain identifier corresponding to the sequence listed.

* **UniProtID**
  Identifies the UniProt accession from which the full-length reference sequence was obtained.

* **ChainPair**
  Describes the interaction context within the PDB complex.
  The first chain ID always corresponds to the **receptor** chain, and the second chain ID corresponds to the **ligand (peptide)** chain.

If multiple chain pairs within the same PDB entry are identical in sequence and interaction context, they are listed as multiple `ChainPair` entries separated by the `=` character.

### Example protein–protein pairs

```
1ca9_A|Q12933|AG=EG    1ca9_G|P20333|AG=EG
1cg9_A|P01889|AC       1cg9_C|P03204|AC
1d8d_A|Q04631|AP       1d8d_P|P01116|AP
1d8d_B|Q02293|BP       1d8d_P|P01116|BP
1ddv_A|Q9Z214|AB       1ddv_B|P31424|AB
1ds5_A|P28523|AE=BE    1ds5_E|P67870|AE=BE
```

---

## 4. Ground-truth PPI site annotation files (`.txt`)

Files with the `.txt` extension contain residue-level ground-truth PPI site annotations derived from PDB structures.

Each protein entry consists of **four consecutive lines**:

1. **Header line**
   Starts with the `>` character and contains the protein sequence name.

2. **Full-length sequence**
   The complete protein sequence derived from the corresponding UniProt ID.

3. **PDB-derived sequence**
   The protein sequence extracted from the PDB entry and specified chain, aligned to the full-length UniProt reference sequence.
   Missing or unresolved residues are represented by `-` characters.

4. **Ground-truth PPI annotation**
   A position-wise binary annotation indicating the PPI site state for each residue in the given interaction context (chain pair).
   Interacting residues are marked with `1`, non-interacting residues with `0`.

### Example PPI site annotation entries

```
>1be3_K|P07552|JK
MLTRFLGPRYRQLARNWVPTASLWGAVGAVGLVWATDWRLILDWVPYINGKFKKDD
--------------RNWVPTAQLWGAVGAVGLVSAT--------------------
00000000000000010011001101100110010000000000000000000000

>1be3_J|P00130|JK
MVAPTLTARLYSLLFRRTSTFALTIVVGALFFERAFDQGADAIYEHINEGKLWKHIKHKYENKE
-VAPTLTARLYSLLFRRTSTFALTIVVGALFFERAFDQGADAIYEHINEGKLWKHIKHKYENK-
0000000000000000111101100110010000000000000000000000000000000000

>1vf5_G|P83797|BG
MVEPLLDGLVLGLVFATLGGLFYAAYQQYKRPNELGG
--------LVLGLVFATLGGLFYAAYQQYKR------
0000000000000110010011010011001000000

>1vf5_B|P83792|BG
MATLKKPDLSDPKLRAKLAKGMGHNYYGEPAWPNDLLYVFPVVIMGTFACIVALSVLDPAMVGEPADPFATPLEILPEWYLYPVFQILRSVPNKLLGVLLMASVPLGLILVPFIENVNKFQNPFRRPVATTIFLFGTLVTIWLGIGATFPLDKTLTLGLF
-----------------LAKGMGHNYYGEPAWPNDLLYVFPVVIMGTFACIVALSVLDPAMVGEPANPFATPLEILPEWYLYPVFQILRSLPNKLLGVLLMASVPLGLILVPFIENVNKFQNPFRRPVATTIFLFGTLVTIWLGIGAALPLDKTL-----
0000000000000000000000000000000100000010011001001100110010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
```

