import streamlit as st
import requests
import py3Dmol
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
def get_pdb(pdb_id):
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        return response.text
    else:
        return None
    
#Search protien Manualy    
def search_protein(name):

    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": name
            }
        },
        "return_type": "entry"
    }

    response = requests.post(url, json=query)

    if response.status_code == 200:
        return response.json()

    return None



# Sidebar
st.sidebar.title("Controls")
st.sidebar.markdown("### Input Method")
protein_name = st.sidebar.text_input(
    "Search Protein Name"
)
pdb_id = st.sidebar.text_input(
    "Enter PDB ID",
    "1CRN"
)

uploaded_file = st.sidebar.file_uploader(
    "Upload PDB File",
    type=["pdb"]
)
style = st.sidebar.selectbox(
    "Visualization Style",
    ["cartoon", "stick", "sphere", "surface"]
)

color = st.sidebar.selectbox(
    "Color Scheme",
    ["spectrum", "greenCarbon", "chain"]
)


def show_protein(pdb_text, style, color):

    viewer = py3Dmol.view(
        width=800,
        height=600
    )

    viewer.addModel(pdb_text, "pdb")


    if style == "cartoon":
        viewer.setStyle({ "cartoon": {"color": color} })

    elif style == "stick":
        viewer.setStyle({ "stick": { "colorscheme": color }})

    elif style == "sphere":
        viewer.setStyle({ "sphere": { "colorscheme": color }})

    elif style == "surface":
        viewer.addSurface(
            py3Dmol.VDW,
            {
                "opacity": 0.8,
                "color": "lightblue"
            } 
        )

    #legand color
    viewer.setStyle(
    {'hetflag': True},
    {'stick': {'color': 'red'}}
    )

    viewer.zoomTo()

    return viewer._make_html()


#Secondary Structure count function
def count_secondary_structure(pdb_text):

    helix_count = 0
    sheet_count = 0

    for line in pdb_text.split("\n"):

        if line.startswith("HELIX"):
            helix_count += 1

        elif line.startswith("SHEET"):
            sheet_count += 1

    return helix_count, sheet_count



#Protien Statistics function
def protein_statistics(sequence):

    # Amino acid molecular weights (Da)
    aa_weight = {
        "ALA": 89.09,
        "ARG": 174.20,
        "ASN": 132.12,
        "ASP": 133.10,
        "CYS": 121.16,
        "GLN": 146.15,
        "GLU": 147.13,
        "GLY": 75.07,
        "HIS": 155.16,
        "ILE": 131.17,
        "LEU": 131.17,
        "LYS": 146.19,
        "MET": 149.21,
        "PHE": 165.19,
        "PRO": 115.13,
        "SER": 105.09,
        "THR": 119.12,
        "TRP": 204.23,
        "TYR": 181.19,
        "VAL": 117.15
    }
    amino_acids = sequence.split()

    length = len(amino_acids)


    # Molecular weight
    molecular_weight = 0

    for aa in amino_acids:
        if aa in aa_weight:
            molecular_weight += aa_weight[aa]


    # Groups
    hydrophobic = [
        "ALA","VAL","ILE","LEU",
        "MET","PHE","TRP","PRO"
    ]

    polar = [
        "SER","THR","ASN",
        "GLN","TYR","CYS"
    ]

    charged = [
        "ASP","GLU",
        "LYS","ARG","HIS"
    ]


    hydrophobic_count = 0
    polar_count = 0
    charged_count = 0


    for aa in amino_acids:

        if aa in hydrophobic:
            hydrophobic_count += 1

        elif aa in polar:
            polar_count += 1

        elif aa in charged:
            charged_count += 1

    return (
        length,
        molecular_weight,
        hydrophobic_count,
        polar_count,
        charged_count
    )


if "loaded" not in st.session_state:
    st.session_state.loaded = False
if st.button("Load Structure"):
    st.session_state.loaded = True
if st.session_state.loaded:

    if protein_name:
        result = search_protein(protein_name)

        if result and result["result_set"]:
            pdb_id = result["result_set"][0]["identifier"]
            st.info(f"Found PDB ID: {pdb_id}")
        else:
            st.error("Protein not found!")
            st.stop()


    if uploaded_file is not None:
        pdb_text = uploaded_file.read().decode("utf-8")
    else:
        pdb_text = get_pdb(pdb_id)


    if pdb_text:

        st.success("Structure loaded successfully!")
        report_button_placeholder = st.container()

        html = show_protein(pdb_text, style,color)
        #html = None
        components.html(
            html,
            height=600
        )
        st.download_button(
            label="Download PDB File",
            data=pdb_text,
            file_name=f"{pdb_id}.pdb",
            mime="text/plain"
        )

        # Metadata Extraction
        title = ""
        method = ""

        for line in pdb_text.split("\n"):
            if line.startswith("TITLE"):
                title += line[10:].strip() +" "
            elif line.startswith("EXPDTA"):
                method = line[10:].strip()

        #Atom count
        atom_count = 0
        for line in pdb_text.split("\n"):

            if line.startswith("ATOM"):
                atom_count += 1


        # Chain Count
        chains = set()
        for line in pdb_text.split("\n"):

            if line.startswith("ATOM"):
                chains.add(line[21])
    
        


        # Ligand Detection
        ligands = set()

        for line in pdb_text.split("\n"):
            if line.startswith("HETATM"):
                ligand = line[17:20].strip()

                if ligand:
                    ligands.add(ligand)


        #protien sequence extraction
        sequence = ""

        for line in pdb_text.split("\n"):

            if line.startswith("SEQRES"):
                sequence += line[19:].strip() + " "

        aa_map = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
            "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
            "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
            "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
            "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
        }
        reference_sequence = ""

        for aa in sequence.split():

            if aa in aa_map:
                reference_sequence += aa_map[aa]
      

        # Sequence Length
        sequence_length = len(sequence.split())


        # detect sequence DNA or Protien
        dna_codes = ["DA", "DT", "DG", "DC"]

        is_dna = False

        for code in dna_codes:
            if code in sequence:
                is_dna = True
                break
        if is_dna:
            st.warning("This is a DNA structure 🧬")

            dna_map = {
                "DA": "A",
                "DT": "T",
                "DG": "G",
                "DC": "C"
            }

            dna_sequence = ""

            for code in sequence.split():
                if code in dna_map:
                    dna_sequence += dna_map[code]

            st.subheader("DNA Sequence")
            st.code(dna_sequence)

        else:
            st.success("This is a Protein structure 🧪")


        # Protien Information
        structure_type = "DNA" if is_dna else "Protein"
        st.subheader(f"{structure_type} Information")
        st.write("**PDB ID** :", pdb_id)
        st.write("**Title**:", title)
        st.write("**Experimental Method** :", method)
        st.write("**Number of characters**:", len(pdb_text))
        st.write("**Number of Atoms** :", atom_count)
        st.write("**Number of Chains** :", len(chains))
        st.write("**Chains** :", ", ".join(sorted(chains)))
        st.write("**Sequence Length** :", sequence_length)

        # Ligands
        st.subheader("Ligands Found")

        if ligands:
            for ligand in sorted(ligands):
                st.write("•", ligand)
        else:
            st.write("No ligands found")

        # secondary structure count
        helix_count, sheet_count = count_secondary_structure(pdb_text)

        st.subheader("Secondary Structure")

        st.write("Alpha Helices:", helix_count)
        st.write("Beta Sheets:", sheet_count)

        if is_dna:
            st.subheader("Nucleotide Sequence")
        else:
            st.subheader("Protein Sequence")
        st.text_area(
            "Sequence",
            sequence,
            height=150
        )


        # Mutation Analysis
        st.subheader("Reference Sequence (1-letter)")
        st.code(reference_sequence) 
        if not is_dna:
            st.subheader("Mutation Analysis")

            query_sequence = st.text_area(
                "Paste your Protein Sequence (1-letter code)",
                height=150,
                placeholder="Example: TTCCPSIVARSNFNV..."
            )
            compare = st.button("Compare Sequence")

            if compare:

                # Remove spaces and convert to uppercase
                query_sequence = query_sequence.replace(" ", "").upper()

                if query_sequence == "":
                    st.warning("Please paste a protein sequence.")

                elif len(query_sequence) != len(reference_sequence):
                    st.error("Sequence length does not match the reference sequence.")

                else:
                    mutations = []

                    for i in range(len(reference_sequence)):

                        if reference_sequence[i] != query_sequence[i]:
                            mutations.append(
                                (i + 1, reference_sequence[i], query_sequence[i]))

                    if len(mutations) == 0:
                        st.success("✅ No Mutation Detected")

                    else:
                        st.error(f"⚠ Mutation Detected ({len(mutations)} mutation(s))")
                        st.subheader("Mutation Details")
                        for pos, ref, mut in mutations:
                            st.write(
                            f"Position {pos}: {ref} → {mut} ({ref}{pos}{mut})"
                        )
                        # Report text
                        report = "MUTATION ANALYSIS REPORT\n"
                        report += "========================\n\n"

                        report += f"PDB ID: {pdb_id}\n"
                        report += f"Reference Sequence Length: {len(reference_sequence)}\n"
                        report += f"Query Sequence Length: {len(query_sequence)}\n"
                        report += f"Total Mutations: {len(mutations)}\n\n"

                        report += "Mutation Details\n"
                        report += "----------------\n"

                        for pos, ref, mut in mutations:

                            mutation_name = f"{ref}{pos}{mut}"

                            st.write(
                                f"Position {pos}: {ref} → {mut} ({mutation_name})"
                            )

                            report += (
                                f"Position: {pos} | "
                                f"Reference: {ref} | "
                                f"Mutated: {mut} | "
                                f"Mutation: {mutation_name}\n"
                            )

                        # Download report
                        st.download_button(
                            label="📄 Download Mutation Report",
                            data=report,
                            file_name=f"{pdb_id}_mutation_report.txt",
                            mime="text/plain"
                        )

        #Amino count
        if not is_dna:

            # Protein Statistics
            length, molecular_weight, hydrophobic, polar, charged = protein_statistics(sequence)

            st.subheader("Protein Statistics")

            st.write("Sequence Length:", length, "amino acids")
            st.write("Molecular Weight:", round(molecular_weight, 2), "Da")

            # Residue Percentage

            total = length

            hydrophobic_percent = (hydrophobic / total) * 100
            polar_percent = (polar / total) * 100
            charged_percent = (charged / total) * 100

            st.subheader("Residue Composition")

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric(
                    "Hydrophobic",
                    f"{hydrophobic}",
                    f"{hydrophobic_percent:.1f}%"
                )
            with col2:
                st.metric(
                    "Polar",
                    f"{polar}",
                    f"{polar_percent:.1f}%"
                )
            with col3:
                st.metric(
                    "Charged",
                    f"{charged}",
                    f"{charged_percent:.1f}%"
            )
           
 
            labels = [
                "Hydrophobic", "Polar", "Charged"
            ]

            values = [
                hydrophobic_percent, polar_percent, charged_percent
            ]
            
            fig, ax = plt.subplots()

            ax.pie(
                values, labels=labels, autopct="%1.1f%%"
            )
            ax.set_title("Amino Acid Property Distribution")

            st.pyplot(fig)



            # Amino Acid Composition
            amino_counts = {}

            for aa in sequence.split():
                if aa in amino_counts:
                    amino_counts[aa] += 1
                else:
                    amino_counts[aa] = 1

            # Two columns: composition on left, graph on right
            col1, col2 = st.columns([1, 1.5])

            #left columns
            with col1:
                st.subheader("Amino Acid Composition")

                for aa, count in sorted(amino_counts.items()):
                    st.write(f"**{aa}:** {count}")

            #Right columns
            with col2:
                st.subheader("Amino Acid Composition Graph")
                aa_names = list(amino_counts.keys())
                aa_counts = list(amino_counts.values())

                fig, ax = plt.subplots(figsize=(6, 4))
                ax.bar(aa_names, aa_counts)
                ax.set_xlabel("Amino Acids")
                ax.set_ylabel("Count")
                ax.set_title("Amino Acid Composition")

                plt.xticks(rotation=45)

                st.pyplot(fig)

        # Download Complete Analysis Report
        report = ""

        report += "3D PROTEIN VISUALIZATION AND ANALYSIS REPORT\n"
        report += "=" * 50 + "\n\n"

        # 1. Structure Information
        report += "1. STRUCTURE INFORMATION\n"
        report += "-" * 30 + "\n"

        report += f"PDB ID: {pdb_id}\n"
        report += f"Structure Type: {'DNA' if is_dna else 'Protein'}\n"
        report += f"Title: {title.strip() or 'Not available'}\n"
        report += f"Experimental Method: {method or 'Not available'}\n"
        report += f"Number of Atoms: {atom_count}\n"
        report += f"Number of Chains: {len(chains)}\n"
        report += f"Chains: {', '.join(sorted(chains)) or 'Not available'}\n"
        report += f"Sequence Length: {sequence_length}\n\n"

        # 2. Ligand Information
        report += "2. LIGAND INFORMATION\n"
        report += "-" * 30 + "\n"

        if ligands:
            report += "Ligands found: " + ", ".join(sorted(ligands)) + "\n"
        else:
            report += "No HETATM residue names detected.\n"

        report += "\n"

        # 3. Secondary Structure
        report += "3. SECONDARY STRUCTURE\n"
        report += "-" * 30 + "\n"
        report += f"Alpha Helices: {helix_count}\n"
        report += f"Beta Sheets: {sheet_count}\n\n"

        # 4. Sequence Information
        report += "4. SEQUENCE INFORMATION\n"
        report += "-" * 30 + "\n"

        if is_dna:
            report += "DNA Sequence:\n"
            report += dna_sequence + "\n\n"
        else:
            report += "Protein Sequence (one-letter code):\n"
            report += reference_sequence + "\n\n"

        # 5. Protein Statistics
        if not is_dna:

            length, molecular_weight, hydrophobic, polar, charged = (
                protein_statistics(sequence)
            )

            report += "5. PROTEIN STATISTICS\n"
            report += "-" * 30 + "\n"

            report += f"Sequence Length: {length} amino acids\n"
            report += (
                f"Estimated Molecular Weight: "
                f"{molecular_weight:.2f} Da\n"
            )

            total = length

            if total > 0:
                report += (
                    f"Hydrophobic Residues: {hydrophobic} "
                    f"({hydrophobic / total * 100:.2f}%)\n"
                )
                report += (
                    f"Polar Residues: {polar} "
                    f"({polar / total * 100:.2f}%)\n"
                )
                report += (
                    f"Charged Residues: {charged} "
                    f"({charged / total * 100:.2f}%)\n"
                )

            # Amino Acid Composition
            amino_counts = {}

            for aa in sequence.split():
                amino_counts[aa] = amino_counts.get(aa, 0) + 1

            report += "\n6. AMINO ACID COMPOSITION\n"
            report += "-" * 30 + "\n"

            for aa, count in sorted(amino_counts.items()):
                report += f"{aa}: {count}\n"

        with report_button_placeholder:
            st.subheader("Download Analysis Report")
            st.download_button(
                label="📄 Download Complete Analysis Report",
                data=report,
                file_name=f"{pdb_id}_analysis_report.txt",
                mime="text/plain",
                key="download_complete_analysis_report"
            )

    else:
        st.error("Invalid PDB ID")