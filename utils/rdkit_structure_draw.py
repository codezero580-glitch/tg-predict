from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def draw_rdkit_structure(mol, width=750, height=450):
    mol_copy = Chem.Mol(mol)

    if mol_copy.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(mol_copy)

    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.addAtomIndices = False
    options.bondLineWidth = 2

    drawer.DrawMolecule(mol_copy)
    drawer.FinishDrawing()

    return drawer.GetDrawingText()
