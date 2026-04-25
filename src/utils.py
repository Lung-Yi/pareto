import base64
import io
import pandas as pd
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

def smiles_to_image_b64(smiles, size=(120, 120)):
    if not RDKIT_OK or not smiles: return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO(); img.save(buf, format="PNG", optimize=True)
        return base64.b64encode(buf.getvalue()).decode()
    except: return None

def parse_csv(contents, filename):
    try:
        _, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string)
        if not filename.lower().endswith(".csv"):
            return None, "請上傳 .csv 格式的檔案。"
        df = pd.read_csv(io.StringIO(decoded.decode("utf-8")))
    except Exception as e:
        return None, f"讀取失敗：{e}"
    if "SMILES" in df.columns and "smiles" not in df.columns:
        df = df.rename(columns={"SMILES":"smiles"})
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) < 2:
        return None, "CSV 中需要至少 2 個數值欄位。"
    return df, None
