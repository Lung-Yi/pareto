# Pareto Front Scatter Explorer 🧬

這是一個專為化學分子數據分析設計的互動式探索工具。透過多目標優化（Pareto Front）演算法，協助研究人員從大量分子數據中篩選出性質最優的候選結構，並整合了 RDKit 進行分子結構的視覺化與標準化處理。

## 🎯 使用目的

在藥物設計或材料科學中，通常需要同時優化多個性質（例如：降低毒性同時提高活性）。本工具旨在：
- **視覺化權衡 (Trade-offs)**：透過 2D 散布圖直觀觀察不同性質間的關係。
- **Pareto 篩選**：自動計算分子的 Pareto Rank，識別出不被其他分子「支配」的最佳候選者。
- **結構關聯分析**：結合 Tanimoto 相似度，尋找與目標分子結構相似且性質優異的化合物。
- **高效管理**：透過分頁圖庫瀏覽成千上萬的分子結構，並能跨頁勾選、整合成新的清單匯出。

## ✨ 核心功能

1. **數據導入與預處理**：
   - 支援上傳包含 SMILES 的 `.csv` 檔案。
   - **自動標準化 (Canonicalization)**：上傳時自動校正 SMILES 格式並過濾無效列。
2. **互動式 Pareto 散布圖**：
   - 自定義 X/Y 軸性質及優化方向（最大化或最小化）。
   - 懸停 (Hover) 即時預覽分子結構及所有 CSV 登載之性質。
3. **高效分子圖庫 (Gallery)**：
   - **分頁顯示**：支援大數據集流暢預覽。
   - **全量 Rank 計算**：顯示所有分子的真實 Pareto Rank 順序。
   - **選取連動**：在圖庫中勾選分子，散布圖會同步高亮標示位置。
4. **相似度搜尋**：
   - 輸入目標 SMILES，計算全資料集的 Tanimoto 相似度。
5. **數據匯出**：
   - 支援將選取的分子完整性質匯出為 `pareto_selected.csv`。

## 🚀 安裝與執行

### 1. 環境需求
建議使用 Python 3.9+。

### 2. 安裝套件
請在終端機執行以下指令安裝所需依賴：
```bash
pip install -r requirements.txt
```

### 3. 啟動程式
執行主程式後，系統會自動開啟瀏覽器視窗：
```bash
python main.py
```

## 📦 打包成 .exe 執行檔

本專案已經過優化，支援使用 `PyInstaller` 打包為獨立執行檔。

### 打包步驟：
1. **安裝 PyInstaller**：
   ```bash
   pip install pyinstaller
   ```
2. **執行打包指令**：
   在專案根目錄下執行以下指令：
   ```bash
   pyinstaller --onefile --noconsole --name "ParetoExplorer" --add-data "src:src" main.py
   ```
   - `--onefile`: 打包成單一執行檔。
   - `--noconsole`: 執行時不顯示黑色的控制台視窗。
   - `--add-data "src:src"`: 將原始碼目錄作為資源打包進去（Linux/Mac 使用 `:` 分隔，Windows 建議使用 `;`）。

3. **取得檔案**：
   打包完成後，可在生成的 `dist/` 資料夾中找到 `ParetoExplorer.exe`。

## 📂 專案結構

- `main.py`: 程式進入點，處理路徑相容性與瀏覽器啟動。
- `src/`:
  - `app.py`: UI 介面佈局定義。
  - `callbacks.py`: 互動邏輯與分頁處理。
  - `engine.py`: Pareto 演算法與 RDKit 計算核心。
  - `utils.py`: CSV 解析與高效圖片轉換。
- `data/`: 存放測試資料夾。
- `GEMINI.md`: 專案開發進度與技術備忘錄。

---
**開發備註**：本工具使用 Dash (Plotly) 框架開發，計算核心依賴 RDKit。對於超過 5000 筆的數據，圖庫分頁功能可顯著節省客戶端記憶體消耗。
