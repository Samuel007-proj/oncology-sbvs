# %%
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, SDMolSupplier, SDWriter

PREP = Path("PDBs/prepared")
COCRYSTAL = Path("Ligands/co-crystal ligand")
OUT = Path("Ligands/validation")
OUT.mkdir(parents=True, exist_ok=True)

tags = ["9TC5", "7EWA", "9ZL3", "8Q68", "6UYB", "7ZJQ", "5OAQ", "6HAZ"]

# %%
# 1. Clean RCSB crystal ligands (keep original coordinates)
for tag in tags:
    src = COCRYSTAL / f"{tag}_LIG.sdf"
    dst_dir = OUT / tag
    dst_dir.mkdir(exist_ok=True)

    mol = next((m for m in SDMolSupplier(str(src), removeHs=False, sanitize=False) if m), None)
    if mol is None:
        print(f"{tag}: could not read")
        continue

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    clean_sdf = dst_dir / f"{tag}_crystal.sdf"
    w = SDWriter(str(clean_sdf))
    w.write(mol)
    w.close()
    print(tag, "-> cleaned")

# %%
# 2. Convert crystal ligands to PDBQT (OpenBabel)
for tag in tags:
    clean_sdf = OUT / tag / f"{tag}_crystal.sdf"
    pdbqt = OUT / tag / f"{tag}_crystal.pdbqt"

    if not clean_sdf.exists():
        print(f"{tag}: missing cleaned SDF")
        continue

    cmd = ["obabel", "-isdf", str(clean_sdf), "-opdbqt", "-O", str(pdbqt)]
    print(f"OpenBabel {tag} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("  FAILED:", result.stderr[:300])
    else:
        print("  OK")

# %%
# 3. Local minimization of the crystal pose
for tag in tags:
    receptor = PREP / f"{tag}_prep.pdb"
    ligand = OUT / tag / f"{tag}_crystal.pdbqt"
    out_sdf = OUT / tag / f"{tag}_min.sdf"
    log_file = OUT / tag / f"{tag}_min.log"

    if not ligand.exists():
        print(f"{tag}: missing PDBQT - skip")
        continue

    cmd = [
        "smina",
        "-r", str(receptor),
        "-l", str(ligand),
        "--autobox_ligand", str(OUT / tag / f"{tag}_crystal.sdf"),
        "--autobox_add", "2",
        "--minimize",
        "--num_modes", "1",
        "-o", str(out_sdf),
        "--log", str(log_file),
    ]
    print(f"=== Minimizing {tag} ===")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("  FAILED:", result.stderr[:300])
    else:
        for line in result.stdout.splitlines():
            if "Affinity" in line:
                print(" ", line.strip())

# %%
# 4. Heavy-atom RMSD
def calc_rmsd_heavy(docked_sdf, ref_sdf):
    docked = next((m for m in SDMolSupplier(str(docked_sdf), removeHs=True, sanitize=False) if m), None)
    ref = next((m for m in SDMolSupplier(str(ref_sdf), removeHs=True, sanitize=False) if m), None)

    if docked is None or ref is None:
        return "read failed"

    print(f"  heavy atoms: docked={docked.GetNumAtoms()}  ref={ref.GetNumAtoms()}")

    if docked.GetNumAtoms() != ref.GetNumAtoms():
        return "atom count mismatch"

    try:
        Chem.SanitizeMol(docked)
        Chem.SanitizeMol(ref)
        rmsd = rdMolAlign.GetBestRMS(docked, ref)
        return f"{rmsd:.2f} A"
    except Exception as e:
        return f"failed: {e}"

print("\n=== RMSD results ===")
for tag in tags:
    docked = OUT / tag / f"{tag}_min.sdf"
    ref = OUT / tag / f"{tag}_crystal.sdf"
    if docked.exists() and ref.exists():
        result = calc_rmsd_heavy(docked, ref)
        print(f"{tag}: {result}")
    else:
        print(f"{tag}: files missing")
