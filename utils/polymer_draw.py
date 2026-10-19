
from rdkit import Chem
from rdkit.Chem import AllChem


ALLOWED_ELEMENTS = {
    "C", "H", "O", "N", "S", "P",
    "F", "Cl", "Br",
}


def build_polymer(psmiles: str, repeat_count: int = 5):
    """
    Build a polymer chain from a pSMILES repeat unit.

    Assumptions:
    - The repeat unit contains exactly two dummy atoms [*].
    - Each dummy atom has exactly one neighboring atom.
    - Inter-unit connections use single bonds.
    - Only elements listed in ALLOWED_ELEMENTS are supported.

    Returns:
        A sanitized RDKit molecule with explicit hydrogen atoms
        and 3D coordinates, along with the name of the force field
        used for geometry optimization.
    """

    if not isinstance(psmiles, str) or not psmiles.strip():
        raise ValueError("pSMILES must not be empty.")

    if (
        not isinstance(repeat_count, int)
        or isinstance(repeat_count, bool)
        or repeat_count < 1
    ):
        raise ValueError(
            "The repeat count must be a positive integer."
        )

    # 1. Parse the repeat unit
    unit = Chem.MolFromSmiles(psmiles)

    if unit is None:
        raise ValueError(
            "Invalid pSMILES or unsupported SMILES syntax."
        )

    # 2. Validate supported elements before building the chain
    unsupported = sorted({
        atom.GetSymbol()
        for atom in unit.GetAtoms()
        if atom.GetAtomicNum() != 0
        and atom.GetSymbol() not in ALLOWED_ELEMENTS
    })

    if unsupported:
        elements = ", ".join(unsupported)
        raise ValueError(
            f"The following elements are not supported: {elements}. "
            "Please use elements supported by the application."
        )

    # 3. Identify the dummy atoms [*]
    dummy_atoms = [
        atom for atom in unit.GetAtoms()
        if atom.GetAtomicNum() == 0
    ]

    if len(dummy_atoms) != 2:
        raise ValueError(
            f"Exactly two dummy atoms [*] are required, "
            f"but {len(dummy_atoms)} were found."
        )

    for atom in dummy_atoms:
        if atom.GetDegree() != 1:
            raise ValueError(
                "Each dummy atom [*] must have exactly one neighbor."
            )

        # For this linear-chain implementation, each dummy atom
        # must be connected to its neighboring atom by a single bond.
        bond = unit.GetBondBetweenAtoms(
            atom.GetIdx(),
            atom.GetNeighbors()[0].GetIdx()
        )

        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError(
                "Each dummy atom [*] must be connected to its "
                "neighboring atom by a single bond in the "
                "polymer representation supported by this application."
            )

    # 4. Store the dummy atom indices and connection endpoints
    dummy_ids = {atom.GetIdx() for atom in dummy_atoms}

    left_dummy, right_dummy = sorted(
        dummy_atoms,
        key=lambda atom: atom.GetIdx()
    )

    left_endpoint = left_dummy.GetNeighbors()[0].GetIdx()
    right_endpoint = right_dummy.GetNeighbors()[0].GetIdx()

    if left_endpoint == right_endpoint:
        raise ValueError(
            "The two connection points must not refer "
            "to the same atom."
        )

    # 5. Create copies of the repeat unit
    rw = Chem.RWMol()
    endpoint_pairs = []
    all_dummy_ids = []

    for repeat_idx in range(repeat_count):
        offset = rw.GetNumAtoms()

        for atom in unit.GetAtoms():
            rw.AddAtom(Chem.Atom(atom))

        for bond in unit.GetBonds():
            rw.AddBond(
                offset + bond.GetBeginAtomIdx(),
                offset + bond.GetEndAtomIdx(),
                bond.GetBondType(),
            )

        endpoint_pairs.append((
            offset + left_endpoint,
            offset + right_endpoint,
        ))

        all_dummy_ids.extend(
            offset + idx for idx in dummy_ids
        )

    # 6. Connect adjacent repeat units
    for i in range(repeat_count - 1):
        right_idx = endpoint_pairs[i][1]
        next_left_idx = endpoint_pairs[i + 1][0]

        rw.AddBond(
            right_idx,
            next_left_idx,
            Chem.BondType.SINGLE,
        )

    # 7. Remove the dummy atoms
    for idx in sorted(all_dummy_ids, reverse=True):
        rw.RemoveAtom(idx)

    mol = rw.GetMol()

    # 8. Validate the resulting polymer structure
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ValueError(
            "The resulting polymer structure is invalid. "
            "Check the valences, bond types, and pSMILES representation."
        ) from exc

    # 9. Add hydrogen atoms
    mol = Chem.AddHs(mol)

    # 10. Generate 3D coordinates
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxIterations = 1000

    status = AllChem.EmbedMolecule(mol, params)

    if status != 0:
        raise RuntimeError(
            "Failed to generate 3D coordinates using ETKDG."
        )

    # 11. Optimize the molecular geometry
    if AllChem.MMFFHasAllMoleculeParams(mol):
        result = AllChem.MMFFOptimizeMolecule(
            mol, maxIters=1000
        )
        force_field = "MMFF94"

    elif AllChem.UFFHasAllMoleculeParams(mol):
        result = AllChem.UFFOptimizeMolecule(
            mol, maxIters=1000
        )
        force_field = "UFF"

    else:
        raise RuntimeError(
            "Neither MMFF94 nor UFF supports all atoms "
            "in this molecular structure."
        )

    if result < 0:
        raise RuntimeError(
            "Geometry optimization failed."
        )

    return mol, force_field
