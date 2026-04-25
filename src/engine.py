import functools
import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

PARETO_COLORS = {1:"#FF4136", 2:"#FF851B", 3:"#FFDC00", 4:"#2ECC40", 5:"#0074D9", 0:"#BBBBBB"}

def compute_pareto_ranks(df, x_col, y_col, x_dir, y_dir, n_ranks):
    vals = df[[x_col, y_col]].copy().astype(float)
    if x_dir == "maximize": vals[x_col] = -vals[x_col]
    if y_dir == "maximize": vals[y_col] = -vals[y_col]
    ranks = pd.Series(0, index=df.index)
    remaining = vals.copy()
    rank = 1
    while len(remaining) > 0 and rank <= n_ranks:
        arr = remaining.values
        front_idx = []
        for i in range(len(arr)):
            dominated = any(
                arr[j,0] <= arr[i,0] and arr[j,1] <= arr[i,1] and
                (arr[j,0] < arr[i,0] or arr[j,1] < arr[i,1])
                for j in range(len(arr)) if j != i
            )
            if not dominated: front_idx.append(remaining.index[i])
        ranks[front_idx] = rank
        remaining = remaining.drop(front_idx)
        rank += 1
    return ranks

@functools.lru_cache(maxsize=20000)
def mol_to_fp(smiles):
    if not RDKIT_OK or not smiles: return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
    except: 
        return None

def compute_tanimoto(query_smiles, smiles_list):
    qfp = mol_to_fp(query_smiles)
    if qfp is None: return None
    
    fps = []
    valid_indices = []
    
    for i, s in enumerate(smiles_list):
        fp = mol_to_fp(s)
        if fp is not None:
            fps.append(fp)
            valid_indices.append(i)
            
    results = [-1.0] * len(smiles_list)
    
    if fps:
        sims = DataStructs.BulkTanimotoSimilarity(qfp, fps)
        for idx, sim in zip(valid_indices, sims):
            results[idx] = sim
            
    return results

def normalize_smiles(smiles):
    if not RDKIT_OK or not smiles: return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else smiles
    except: return smiles

def canonicalize_smiles(smiles):
    """回傳 (canonical_smiles, is_valid)"""
    if not RDKIT_OK or not smiles: return None, False
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol:
            return Chem.MolToSmiles(mol), True
        return None, False
    except:
        return None, False
