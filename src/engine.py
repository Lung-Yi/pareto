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

def compute_pareto_ranks(df, x_col, y_col, x_dir, y_dir, n_ranks=None):
    """
    計算所有分子的 Pareto Rank。
    如果 n_ranks 為 None，則計算所有層級。
    """
    vals = df[[x_col, y_col]].copy().astype(float)
    if x_dir == "maximize": vals[x_col] = -vals[x_col]
    if y_dir == "maximize": vals[y_col] = -vals[y_col]
    
    # 轉為 numpy 加快速度
    data = vals.values
    n = len(data)
    ranks = np.zeros(n, dtype=int)
    
    # indices 追蹤還沒被分配 rank 的點
    remaining_indices = np.arange(n)
    current_rank = 1
    
    while len(remaining_indices) > 0:
        current_data = data[remaining_indices]
        is_efficient = np.ones(len(current_data), dtype=bool)
        
        # 尋找目前的 Pareto Front
        for i, c in enumerate(current_data):
            if is_efficient[i]:
                # 檢查是否有其他點支配點 i
                # 支配條件：x_j <= x_i AND y_j <= y_i 且 (x_j < x_i OR y_j < y_i)
                dominated = np.any(
                    (current_data[:, 0] <= c[0]) & (current_data[:, 1] <= c[1]) &
                    ((current_data[:, 0] < c[0]) | (current_data[:, 1] < c[1]))
                )
                if dominated:
                    is_efficient[i] = False
        
        # 分配 Rank
        front_indices = remaining_indices[is_efficient]
        ranks[front_indices] = current_rank
        
        # 移除已分配的點
        remaining_indices = remaining_indices[~is_efficient]
        
        # 如果使用者有指定 n_ranks 且已達到，則跳出（雖然我們現在傾向算完）
        if n_ranks is not None and current_rank >= n_ranks:
            break
        current_rank += 1
        
    return pd.Series(ranks, index=df.index)

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
