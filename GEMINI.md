# Pareto Front Scatter Explorer - 專案導覽

## 🎯 專案目標
開發一個互動式的化學分子 Pareto Front 散布圖工具。
- 支援上傳 CSV 進行多目標優化 (Pareto) 分析。
- 整合 RDKIT 進行 SMILES 結構預覽與 Tanimoto 相似度計算。
- 提供直觀的 UI/UX，並最終支援打包為獨立的 `.exe` 執行檔。

## 🏗️ 專案架構 (Modular Structure)
- `main.py`: 程式進入點，處理資源路徑與自動開啟瀏覽器。
- `src/app.py`: UI 佈局 (Layout) 定義，包含分子圖庫與分頁功能。
- `src/callbacks.py`: 互動邏輯與數據流處理，包含高亮連動、分頁控制與 CSV 匯出。
- `src/engine.py`: 核心計算邏輯 (Pareto Ranks, RDKit Functions, Canonicalization)。
- `src/utils.py`: 輔助工具 (CSV Parsing, Base64 Image Conversion with optimization)。
- `data/`: 存放測試用的資料集 (例如 `data_1000.csv`)。

## 📈 目前進度
- [x] **專案初始化**: 完成單一檔案重構為模組化架構。
- [x] **核心功能**: 實作 Pareto 分層計算與 Tanimoto 相似度比對。
- [x] **效能優化 (新)**:
    - **分頁顯示**: 分子圖庫採用分頁機制 (每頁 48 筆)，大幅降低瀏覽器負擔。
    - **解析度調整**: 降低圖庫分子圖片解析度，減少 Base64 數據傳輸量。
    - **Lazy Loading**: 圖片支援延遲載入。
- [x] **SMILES 處理**: 
    - 上傳時自動進行 SMILES Canonicalization。
    - 自動過濾無法解析的無效 SMILES 列。
    - 匯出時保留 `orig_smiles` 與標準化後的 `smiles`。
- [x] **UI 優化**: 完成分佈圖與分子結構懸停 (Hover) 顯示，動態呈現性質欄位。
- [x] **分子圖庫與導出**: 
    - 顯示所有分子結構，並支援分頁切換與統計資訊。
    - 每個分子卡片均顯示對應的 **X 軸與 Y 軸性質數值**。
    - 支援**跨頁勾選**分子，並在 2D 散布圖中即時高亮顯示。
    - 支援匯出勾選分子的完整性質為 CSV。

## 🛠️ 下一步計畫
- [ ] **功能增強**: 加入過濾器 (Filter) 功能，篩選特定數值範圍的分子。
- [ ] **打包測試**: 實際使用 PyInstaller 進行 Windows 執行檔打包驗證。

## 💡 開發備忘錄
- 修改計算邏輯請優先查看 `src/engine.py`。
- 修改 UI 佈局請調整 `src/app.py`。
- 圖庫分頁大小可在 `src/callbacks.py` 中的 `PAGE_SIZE` 調整。
- 匯出時 `orig_smiles` 為原始輸入，`smiles` 為標準化結果。
- 執行方式: `python main.py`
