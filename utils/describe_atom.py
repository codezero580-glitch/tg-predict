def describe_atom(mol, idx):
    atom = mol.GetAtomWithIdx(idx)

    neighbors = [
        f"{n.GetSymbol()}{n.GetIdx()}"
        for n in atom.GetNeighbors()
    ]

    aromaticity = "aromatic" if atom.GetIsAromatic() else "non-aromatic"

    return (
        f"Atom {idx} | {atom.GetSymbol()} | {aromaticity} | "
        f"H={atom.GetTotalNumHs()} | "
        f"Neighbors={', '.join(neighbors) or '-'}"
    )
