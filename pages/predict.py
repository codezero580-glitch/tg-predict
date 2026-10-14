import streamlit as st
from streamlit_ketcher import st_ketcher

import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdDepictor

import navigation

from utils.rdkit_structure_draw import draw_rdkit_structure
from utils.describe_atom import describe_atom
from utils.make_pmiles import make_psmiles


def show_predict():
    col1, col2, col3, col4 = st.columns([1, 3, 3, 1])

    with col1:
        st.button("Back", on_click=navigation.buka_landing,
                  use_container_width=True,
                  key="about_kembali")

    DEFAULT_SMILES = "Oc1ccccc1C"

    BOND_TYPES = {
        "Single (-)": Chem.BondType.SINGLE,
        "Double (=)": Chem.BondType.DOUBLE,
        "Triple (#)": Chem.BondType.TRIPLE,
    }

    if "base_smiles" not in st.session_state:
        st.session_state.base_smiles = DEFAULT_SMILES

    if "generated_psmiles" not in st.session_state:
        st.session_state.generated_psmiles = ""

    saved_mol = Chem.MolFromSmiles(
        st.session_state.base_smiles or ""
    )

    if saved_mol is None:
        st.session_state.base_smiles = DEFAULT_SMILES
        st.session_state.generated_psmiles = ""

    st.header("Draw a Chemical Structure")

    base_mol = Chem.MolFromSmiles(st.session_state.base_smiles)

    if base_mol is None:
        st.toast("The initial structure is invalid.", icon="❌")
        st.stop()

    rdDepictor.Compute2DCoords(base_mol)
    base_molfile = Chem.MolToMolBlock(base_mol)

    ALLOWED_ELEMENTS = {
        "C", "H", "O", "N", "S", "P",
        "F", "Cl", "Br", "Si"
    }

    # Mencegah struktur awal diproses ulang saat aplikasi pertama kali dibuka
    if "_last_processed_molfile" not in st.session_state:
        st.session_state["_last_processed_molfile"] = base_molfile

    drawn_molfile = st_ketcher(
        value=base_molfile,
        molecule_format="MOLFILE",
        key="polymer_ketcher",
        height=450,
    )

    if (
        isinstance(drawn_molfile, str)
        and not drawn_molfile.strip()
        and st.session_state["_last_processed_molfile"]
    ):
        st.session_state["_last_processed_molfile"] = ""
        st.session_state.base_smiles = ""
        st.session_state.generated_psmiles = ""

        for key in ("atom_1", "atom_2", "bond_1", "bond_2"):
            st.session_state.pop(key, None)

        st.info("Structure cleared.")
        st.stop()

    # Apply: editor mengirim struktur baru
    elif (
        isinstance(drawn_molfile, str)
        and drawn_molfile.strip()
        and drawn_molfile != st.session_state["_last_processed_molfile"]
    ):
        st.session_state["_last_processed_molfile"] = drawn_molfile

        drawn_mol = Chem.MolFromMolBlock(
            drawn_molfile,
            removeHs=False,
        )

        if drawn_mol is None:
            st.toast("RDKit could not parse the structure.", icon="❌")

        elif any(
            atom.GetAtomicNum() == 0
            for atom in drawn_mol.GetAtoms()
        ):
            st.toast(
                "For Mode 1, draw the molecule without wildcard atoms.",
                icon="❌",
            )

        else:
            unsupported = sorted({
                atom.GetSymbol()
                for atom in drawn_mol.GetAtoms()
                if atom.GetAtomicNum() != 0
                and atom.GetSymbol() not in ALLOWED_ELEMENTS
            })

            if unsupported:
                st.toast(
                    f"Unsupported elements: {', '.join(unsupported)}",
                    icon="❌",
                )

            else:
                st.session_state.base_smiles = Chem.MolToSmiles(drawn_mol)
                st.session_state.generated_psmiles = ""

                st.toast("Structure saved successfully.", icon="✅")

    drawn_mol = Chem.MolFromMolBlock(
                drawn_molfile,
                removeHs=False,
                )

    if drawn_mol is None:
        st.toast("""RDKit could not parse the
                 structure from the editor.""", icon="❌")
        st.stop()

    smiles = st.session_state.base_smiles
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        st.toast("The saved structure is invalid.", icon="❌")
        st.stop()

    candidates = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() not in (0, 1)
    ]

    st.header("Select Connection Atoms")

    table_placeholder = st.container()

    if len(candidates) < 2:
        st.warning("The structure must contain at least two candidate atoms.")
        st.stop()

    labels = {
        describe_atom(mol, idx): idx
        for idx in candidates
    }

    st.subheader("Select Two Connection Atoms")

    col5, col6 = st.columns(2)

    with col5:
        label1 = st.selectbox(
            "First Connection Atom",
            list(labels.keys()),
            key="atom_1",
        )

        bond1_label = st.selectbox(
            "First Bond Type",
            list(BOND_TYPES.keys()),
            index=0,
            key="bond_1",
        )

    with col6:
        other_labels = [
            label for label in labels if label != label1
        ]

        label2 = st.selectbox(
            "Second Connection Atom",
            other_labels,
            key="atom_2",
        )

        bond2_label = st.selectbox(
            "Second Bond Type",
            list(BOND_TYPES.keys()),
            index=0,
            key="bond_2",
        )

    atom_idx_1 = labels[label1]
    atom_idx_2 = labels[label2]

    atom_rows = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()

        neighbors = [
            f"{neighbor.GetSymbol()} ({neighbor.GetIdx()})"
            for neighbor in atom.GetNeighbors()
        ]

        atom_rows.append({
            "Index": idx,
            "Element": atom.GetSymbol(),
            "Atomic Number": atom.GetAtomicNum(),
            "Total Hydrogens": atom.GetTotalNumHs(),
            "Aromatic": atom.GetIsAromatic(),
            "Formal Charge": atom.GetFormalCharge(),
            "Neighbors": ", ".join(neighbors),
            "Connection Point": (
                "Connection 1" if idx == atom_idx_1
                else "Connection 2" if idx == atom_idx_2
                else "-"
            ),
        })

    atom_df = pd.DataFrame(atom_rows)

    with table_placeholder:
        if drawn_mol is not None and drawn_mol.GetNumAtoms() > 0:
            st.subheader("Atom Information Table")
            st.dataframe(
                atom_df,
                use_container_width=True,
                hide_index=True,
            )

    col11, col12, col13 = st.columns([1, 5, 1])

    with col13:
        if st.button("Generate Structure"):
            try:
                psmiles = make_psmiles(
                    smiles=smiles,
                    atom_idx_1=atom_idx_1,
                    atom_idx_2=atom_idx_2,
                    bond1=BOND_TYPES[bond1_label],
                    bond2=BOND_TYPES[bond2_label],
                )
    
                st.session_state.generated_psmiles = psmiles
    
            except Exception as exc:
                st.session_state.generated_psmiles = ""
                st.toast(str(exc))

    if st.session_state.generated_psmiles:
        st.header("Candidate Structure for Prediction")

        psmiles = st.session_state.generated_psmiles

        psmiles_mol = Chem.MolFromSmiles(psmiles)

        if psmiles_mol is not None:
            st.image(
                draw_rdkit_structure(psmiles_mol),
                caption="""Please inspect the structure
                above to ensure that it matches your
                intended design""",
                use_container_width=True,
            )
        else:
            st.warning("The structure could not be visualized.")

        col7, col8, col9, col10 = st.columns([1, 3, 3, 1])
        with col7:
            st.button(
                "Back",
                on_click=navigation.buka_landing,
                use_container_width=True,
                key="predict_kembali"
            )

        with col10:
            if psmiles_mol is not None:
                st.button(
                    "Predict",
                    on_click=navigation.prediksi_output,
                    use_container_width=True,
                    key="predict_output"
                )
