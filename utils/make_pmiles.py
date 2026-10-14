import re
import streamlit as st

from rdkit import Chem


def make_psmiles(
        smiles,
        atom_idx_1,
        atom_idx_2,
        bond1,
        bond2,):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        st.toast("The input SMILES is invalid.", icon="❌")

    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        st.toast(
            "The input structure already contains wildcard atoms.",
            icon="❌"
        )

    if atom_idx_1 == atom_idx_2:
        st.toast(
            "Please select two different connection atoms.",
            icon="❌"
        )

    atom_count = mol.GetNumAtoms()

    for idx in (atom_idx_1, atom_idx_2):
        if not 0 <= idx < atom_count:
            st.toast(f"Atom index {idx} is invalid.", icon="❌")

    rw = Chem.RWMol(mol)

    wildcard1 = rw.AddAtom(Chem.Atom(0))

    rw.AddBond(atom_idx_1, wildcard1, bond1)

    wildcard2 = rw.AddAtom(Chem.Atom(0))

    rw.AddBond(atom_idx_2, wildcard2, bond2)

    result = rw.GetMol()

    try:
        Chem.SanitizeMol(result)

    except Exception as exc:
        st.toast(
            "Adding wildcard atoms resulted in a structure "
            "that is invalid according to RDKit. Check the "
            "selected atoms and bond types. Details: " + str(exc),
            icon="❌"
        )

    wildcard_count = sum(
        atom.GetAtomicNum() == 0
        for atom in result.GetAtoms()
    )

    if wildcard_count != 2:
        st.toast(
            f"Expected two wildcard atoms, but found {wildcard_count}.",
            icon="❌"
        )

    if len(Chem.GetMolFrags(result)) != 1:
        st.toast("The resulting structure is disconnected.", icon="❌")

    output = Chem.MolToSmiles(result, canonical=True)

    output = re.sub(r"(?<!\[)\*(?!\])", "[*]", output)

    if len(re.findall(r"\[\*\]", output)) != 2:
        st.toast("The pSMILES format is invalid.", icon="❌")

    check = Chem.MolFromSmiles(
        output.replace("[*]", "*")
    )

    if check is None:
        st.toast("pSMILES revalidation failed.", icon="❌")

    return output
