# %%
import time
import requests
import pandas as pd
from pathlib import Path
from urllib.parse import quote
from rdkit import Chem
from chembl_webresource_client.new_client import new_client
from rdkit.Chem import AllChem, DataStructs

SEEDS = Path("Ligands/cohort_1_seeds.csv")
OUT = Path("Ligands/cohorts")
OUT.mkdir(parents=True, exist_ok=True)

CHEMOTYPES = ["K1", "P1", "P2", "P3", "T1", "S1"]
N_PER_COHORT = 500
SIM_PCT = 40
SLEEP = 0.3

similarity_api = new_client.similarity
substructure_api = new_client.substructure
molecule_api = new_client.molecule

# %%
df = pd.read_csv(SEEDS)
print("columns:", list(df.columns))
print(df)

seeds = {}
for _, row in df.iterrows():
    ct = str(row.chemotype).strip() if "chemotype" in df.columns else str(row.Mol_name)[:2]
    mol = Chem.MolFromSmiles(row.SMILES)
    if mol is None:
        print("BAD SMILES:", row.Mol_name, row.SMILES[:80])
        continue
    seeds.setdefault(ct, []).append({"name": row.Mol_name, "smiles": Chem.MolToSmiles(mol)})
    print(f"  [{ct}] {row.Mol_name}")
print("loaded:", {k: len(v) for k, v in seeds.items()})

# %%
def chembl_similarity(smiles, pct=SIM_PCT, limit=500):
    hits = []
    try:
        res = similarity_api.filter(smiles=smiles, similarity=pct).only(
            ["molecule_chembl_id", "similarity", "molecule_structures"]
        )
        for i, m in enumerate(res):
            if i >= limit:
                break
            st = m.get("molecule_structures") or {}
            smi = st.get("canonical_smiles")
            if not smi:
                try:
                    mol = molecule_api.get(m["molecule_chembl_id"])
                    st = mol.get("molecule_structures") or {}
                    smi = st.get("canonical_smiles")
                except Exception:
                    continue
            if not smi:
                continue
            hits.append({
                "name": m.get("molecule_chembl_id", smi),
                "smiles": smi,
                "score": float(m.get("similarity", pct)) / 100.0,
                "method": "chembl_sim",
                "source": "ChEMBL",
            })
    except Exception as e:
        print("    chembl sim error:", type(e).__name__, e)
    print(f"    chembl sim -> {len(hits)}")
    time.sleep(SLEEP)
    return hits


def chembl_substructure(smiles, limit=300):
    hits = []
    try:
        res = substructure_api.filter(smiles=smiles).only(
            ["molecule_chembl_id", "molecule_structures"]
        )
        for i, m in enumerate(res):
            if i >= limit:
                break
            st = m.get("molecule_structures") or {}
            smi = st.get("canonical_smiles")
            if not smi:
                continue
            hits.append({
                "name": m.get("molecule_chembl_id", smi),
                "smiles": smi,
                "score": 0.99,
                "method": "chembl_sub",
                "source": "ChEMBL",
            })
    except Exception as e:
        print("    chembl sub error:", type(e).__name__, e)
    print(f"    chembl sub -> {len(hits)}")
    time.sleep(SLEEP)
    return hits

# %%
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def pubchem_similarity(smiles, threshold=SIM_PCT, max_cids=500):
    hits = []
    try:
        r = requests.post(
            f"{PUG}/compound/fastsimilarity_2d/smiles/cids/JSON",
            data={"smiles": smiles},
            params={"Threshold": threshold, "MaxRecords": max_cids},
            timeout=120,
        )
        if r.status_code != 200:
            print(f"    pubchem sim HTTP {r.status_code}: {r.text[:200]}")
            return hits
        cids = r.json().get("IdentifierList", {}).get("CID", [])[:max_cids]
        print(f"    pubchem cids -> {len(cids)}")
    except Exception as e:
        print("    pubchem sim error:", e)
        return hits

    time.sleep(SLEEP)
    for i in range(0, len(cids), 50):
        batch = ",".join(map(str, cids[i:i+50]))
        try:
            r = requests.get(
                f"{PUG}/compound/cid/{batch}/property/ConnectivitySMILES,CanonicalSMILES/JSON",
                timeout=60,
            )
            if r.status_code != 200:
                print(f"    pubchem props HTTP {r.status_code}: {r.text[:120]}")
                continue
            props = r.json().get("PropertyTable", {}).get("Properties", [])
            for prop in props:
                smi = prop.get("CanonicalSMILES") or prop.get("ConnectivitySMILES")
                if not smi:
                    continue
                hits.append({
                    "name": f"CID{prop['CID']}",
                    "smiles": smi,
                    "score": threshold / 100.0,
                    "method": "pubchem_sim",
                    "source": "PubChem",
                })
        except Exception as e:
            print("    pubchem props error:", e)
        time.sleep(SLEEP)
    print(f"    pubchem sim -> {len(hits)}")
    return hits

# %%
def canonicalize(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else None


def recruit_api(seed_list, n_target=N_PER_COHORT):
    pool = {}
    rank = {"seed": 0, "chembl_sub": 1, "chembl_sim": 2, "pubchem_sim": 3}

    for s in seed_list:
        can = canonicalize(s["smiles"])
        if can:
            pool[can] = {"name": s["name"], "smiles": can, "score": 1.0,
                         "method": "seed", "source": "seed"}

    for s in seed_list:
        print(f"  query {s['name']}")
        print(f"    SMILES len={len(s['smiles'])}")
        for hit in chembl_similarity(s["smiles"]) + chembl_substructure(s["smiles"]) + pubchem_similarity(s["smiles"]):
            can = canonicalize(hit["smiles"])
            if not can:
                continue
            prev = pool.get(can)
            if prev is None or rank.get(hit["method"], 9) < rank.get(prev["method"], 9):
                hit["smiles"] = can
                pool[can] = hit
            elif hit["score"] > prev["score"]:
                hit["smiles"] = can
                pool[can] = hit

    ranked = sorted(pool.values(), key=lambda x: (
        rank.get(x["method"], 9), -x["score"], x["name"]
    ))
    return ranked[:n_target]

# %%
cohorts = {}
for ct in CHEMOTYPES:
    if ct not in seeds:
        print(f"{ct}: no seed - skip")
        continue
    print(f"\n{'='*40}\n=== {ct} ===")
    hits = recruit_api(seeds[ct], N_PER_COHORT)
    cohorts[ct] = hits
    methods, sources = {}, {}
    for h in hits:
        methods[h["method"]] = methods.get(h["method"], 0) + 1
        sources[h["source"]] = sources.get(h["source"], 0) + 1
    print(f"{ct}: {len(hits)} unique SMILES | methods={methods} | sources={sources}")

# %%
def tanimoto(smiles_a, smiles_b):
    ma, mb = Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b)
    if ma is None or mb is None:
        return 0.0
    fa = AllChem.GetMorganFingerprintAsBitVect(ma, 2, 2048)
    fb = AllChem.GetMorganFingerprintAsBitVect(mb, 2, 2048)
    return DataStructs.TanimotoSimilarity(fa, fb)

# %%
for ct in cohorts:
    out_dir = OUT / ct
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = pd.DataFrame(cohorts[ct])
    seed_smi = seeds[ct][0]["smiles"]
    pdf["score"] = [
        1.0 if m == "seed" else tanimoto(smi, seed_smi)
        for smi, m in zip(pdf["smiles"], pdf["method"])
    ]
    pdf = pdf.sort_values("score", ascending=False)
    pdf.to_csv(out_dir / f"{ct}_cohort.csv", index=False)
    print(ct, "score range:", pdf["score"].min(), "-", pdf["score"].max())
