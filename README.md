# 🧬 3D Protein Visualization and Analysis Tool

A web app built with **Streamlit** that lets you search, load, visualize, and analyze protein and DNA structures from the **RCSB Protein Data Bank (PDB)**. You can also upload your own PDB file, compare sequences to detect mutations, and download detailed analysis reports.

---

## ✨ Features

- **Multiple input methods**
  - Search by protein name (uses the RCSB Search API)
  - Enter a PDB ID directly (e.g. `1CRN`)
  - Upload your own `.pdb` file
- **Interactive 3D visualization** (powered by py3Dmol)
  - Styles: cartoon, stick, sphere, surface
  - Color schemes: spectrum, greenCarbon, chain
  - Ligands highlighted in red
- **Structure information**
  - PDB ID, title, experimental method
  - Number of atoms, chains, and chain IDs
  - Ligand detection (HETATM residues)
  - Secondary structure count (alpha helices and beta sheets)
- **Sequence analysis**
  - Automatic detection of Protein vs DNA structure
  - Sequence extraction from SEQRES records
  - 3-letter to 1-letter code conversion
  - DNA sequence extraction
- **Protein statistics**
  - Sequence length and estimated molecular weight
  - Hydrophobic / polar / charged residue composition (with pie chart)
  - Amino acid composition (counts + bar graph)
- **Mutation analysis**
  - Paste a query sequence and compare it with the reference
  - Lists all point mutations in standard notation (e.g. `A12G`)
  - Download the mutation report as a `.txt` file
- **Downloads**
  - Download the PDB file
  - Download the complete analysis report
  - Download the mutation report

---

## 🛠️ Tech Stack

| Tool | Purpose |
|------|---------|
| Python 3 | Core language |
| Streamlit | Web interface |
| py3Dmol | 3D molecular visualization |
| Requests | Fetching data from the RCSB PDB API |
| Matplotlib | Pie and bar charts |

---

## 📦 Installation

1. **Clone the repository**
   ```bash
   git clone <your-repo-url>
   cd <your-project-folder>
   ```

2. **(Optional) Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate      # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install streamlit requests py3Dmol matplotlib
   ```

---

## ▶️ Usage

Run the app (replace `app.py` with your file name):

```bash
streamlit run app.py
```

Then open the link shown in the terminal (usually `http://localhost:8501`).

### Steps
1. Choose an input method in the sidebar: search a protein name, enter a PDB ID, or upload a `.pdb` file.
2. Select the visualization style and color scheme.
3. Click **Load Structure**.
4. Explore the 3D model, structure details, statistics, and charts.
5. (Optional) Paste a sequence in the **Mutation Analysis** section and click **Compare Sequence**.
6. Download the reports you need.

> **Note:** If a protein name is entered, it takes priority over the PDB ID field. If a file is uploaded, it takes priority over both.

---

## 🧪 Example

- PDB ID: `1CRN` (Crambin) – default example
- Try `1BNA` for a DNA structure
- Try `1LYZ` (Lysozyme) for a protein with ligands

---

## 📁 Project Structure

```
├── app.py          # Main Streamlit application
└── README.md       # Project documentation
```

---

## ⚠️ Limitations

- Mutation analysis only works when the query sequence has the **same length** as the reference sequence (point mutations only, no insertions/deletions).
- Molecular weight is a simple sum of residue weights (water loss from peptide bonds is not subtracted), so it is an estimate.
- Sequence information is read from `SEQRES` records, so files without them will show an empty sequence.
- Protein-name search returns the top match from RCSB full-text search, which may not always be the exact protein you want.
- Requires an internet connection to fetch structures from RCSB.

---

## 🚀 Future Improvements

- Sequence alignment to support insertions and deletions
- Accurate molecular weight and isoelectric point (pI) calculation
- Support for mmCIF file format
- Ramachandran plot and B-factor visualization
- Option to choose from multiple search results

---

## 🙏 Acknowledgements

- [RCSB Protein Data Bank](https://www.rcsb.org/) for structure data and search API
- [py3Dmol](https://pypi.org/project/py3Dmol/) and [3Dmol.js](https://3dmol.csb.pitt.edu/) for visualization
- [Streamlit](https://streamlit.io/) for the web framework

---

## 👤 Author

Made by Priyanshi Gupta
