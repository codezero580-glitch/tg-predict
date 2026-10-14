
from rdkit import Chem
from rdkit.Chem import AllChem


ALLOWED_ELEMENTS = {
    "C", "H", "O", "N", "S", "P",
    "F", "Cl", "Br",
}


def build_polymer(psmiles: str, repeat_count: int = 5):
    """
    Membentuk rantai polimer dari repeat unit pSMILES.

    Asumsi:
    - Terdapat tepat dua atom dummy [*].
    - Setiap atom dummy memiliki satu tetangga.
    - Sambungan antarunit menggunakan ikatan tunggal.
    - Hanya unsur dalam ALLOWED_ELEMENTS yang didukung.

    Returns:
        Molekul RDKit dengan atom hidrogen dan koordinat 3D,
        serta nama force field yang digunakan.
    """

    if not isinstance(psmiles, str) or not psmiles.strip():
        raise ValueError("pSMILES tidak boleh kosong.")

    if not isinstance(repeat_count, int) or isinstance(
        repeat_count, bool
    ) or repeat_count < 1:
        raise ValueError(
            "Jumlah unit pengulangan harus berupa integer minimal 1."
        )

    # 1. Parsing repeat unit
    unit = Chem.MolFromSmiles(psmiles)

    if unit is None:
        raise ValueError(
            "pSMILES tidak valid atau tidak dapat dibaca RDKit."
        )

    # 2. Validasi unsur SEBELUM membangun rantai
    unsupported = sorted({
        atom.GetSymbol()
        for atom in unit.GetAtoms()
        if atom.GetAtomicNum() != 0
        and atom.GetSymbol() not in ALLOWED_ELEMENTS
    })

    if unsupported:
        elements = ", ".join(unsupported)
        raise ValueError(
            f"Unsur berikut belum didukung: {elements}. "
            "Silakan gunakan unsur yang didukung aplikasi."
        )

    # 3. Identifikasi atom dummy [*]
    dummy_atoms = [
        atom for atom in unit.GetAtoms()
        if atom.GetAtomicNum() == 0
    ]

    if len(dummy_atoms) != 2:
        raise ValueError(
            f"Diperlukan tepat 2 atom [*], "
            f"ditemukan {len(dummy_atoms)}."
        )

    for atom in dummy_atoms:
        if atom.GetDegree() != 1:
            raise ValueError(
                "Setiap atom [*] harus memiliki tepat satu tetangga."
            )

        # Untuk implementasi rantai linear ini, dummy harus
        # memiliki ikatan tunggal ke atom tetangganya.
        bond = unit.GetBondBetweenAtoms(
            atom.GetIdx(),
            atom.GetNeighbors()[0].GetIdx()
        )

        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError(
                "Ikatan dari [*] ke atom tetangga harus tunggal "
                "pada format yang didukung aplikasi ini."
            )

    # 4. Simpan indeks atom dummy dan titik ujung
    dummy_ids = {atom.GetIdx() for atom in dummy_atoms}

    left_dummy, right_dummy = sorted(
        dummy_atoms,
        key=lambda atom: atom.GetIdx()
    )

    left_endpoint = left_dummy.GetNeighbors()[0].GetIdx()
    right_endpoint = right_dummy.GetNeighbors()[0].GetIdx()

    if left_endpoint == right_endpoint:
        raise ValueError(
            "Kedua titik sambungan tidak boleh menunjuk "
            "ke atom yang sama."
        )

    # 5. Bangun salinan repeat unit
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

    # 6. Hubungkan ujung antarrepeat unit
    for i in range(repeat_count - 1):
        right_idx = endpoint_pairs[i][1]
        next_left_idx = endpoint_pairs[i + 1][0]

        rw.AddBond(
            right_idx,
            next_left_idx,
            Chem.BondType.SINGLE,
        )

    # 7. Hapus atom dummy
    for idx in sorted(all_dummy_ids, reverse=True):
        rw.RemoveAtom(idx)

    mol = rw.GetMol()

    # 8. Validasi struktur hasil penyambungan
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ValueError(
            "Struktur hasil penyambungan tidak valid. "
            "Periksa valensi, jenis ikatan, dan pSMILES."
        ) from exc

    # 9. Tambahkan hidrogen
    mol = Chem.AddHs(mol)

    # 10. Buat koordinat 3D
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxIterations = 1000

    status = AllChem.EmbedMolecule(mol, params)

    if status != 0:
        raise RuntimeError(
            "Gagal membuat koordinat 3D dengan ETKDG."
        )

    # 11. Optimasi geometri
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
            "MMFF94 dan UFF tidak mendukung seluruh atom "
            "dalam struktur ini."
        )

    if result < 0:
        raise RuntimeError(
            "Terjadi kegagalan saat optimasi geometri."
        )

    return mol, force_field
