import dash
from dash import dcc, html

app = dash.Dash(__name__, title="Pareto Scatter Explorer")
app.config.suppress_callback_exceptions = True

SMILES_PANEL = html.Div([
    # 左側：搜尋與控制區
    html.Div([
        html.Label("🔎 查詢 SMILES", style={"fontWeight":"bold","fontSize":"14px"}),
        html.Div([
            dcc.Input(
                id="query-smiles",
                type="text",
                placeholder="輸入 SMILES，例如：c1ccccc1",
                style={
                    "width":"340px", "padding":"7px 10px",
                    "fontSize":"13px", "border":"1px solid #aac",
                    "borderRadius":"6px", "fontFamily":"monospace",
                    "marginRight":"8px",
                },
            ),
            html.Button("✔ 檢查是否存在", id="btn-check", n_clicks=0,
                style={"padding":"7px 14px","fontSize":"13px","background":"#0074D9",
                       "color":"white","border":"none","cursor":"pointer",
                       "borderRadius":"6px","marginRight":"6px"}),
            html.Button("⬡ 計算 Tanimoto", id="btn-tanimoto", n_clicks=0,
                style={"padding":"7px 14px","fontSize":"13px","background":"#2ECC40",
                       "color":"white","border":"none","cursor":"pointer",
                       "borderRadius":"6px","marginRight":"6px"}),
            html.Button("✕ 清除", id="btn-clear-sim", n_clicks=0,
                style={"padding":"7px 12px","fontSize":"13px","background":"#999",
                       "color":"white","border":"none","cursor":"pointer",
                       "borderRadius":"6px"}),
        ], style={"display":"flex","alignItems":"center","marginTop":"6px","flexWrap":"wrap","gap":"4px"}),
        html.Div(id="query-status",
                 style={"marginTop":"6px","fontSize":"13px","color":"#555","minHeight":"20px"}),
    ], style={"flex":"1"}),

    # 右側：固定顯示當前查詢分子
    html.Div([
        html.Div("當前查詢分子", style={"fontSize":"12px", "color":"#666", "textAlign":"center", "marginBottom":"4px", "fontWeight":"bold"}),
        html.Img(id="query-mol-img", 
                 style={"maxWidth":"120px", "maxHeight":"120px", "display":"none", 
                        "borderRadius":"6px", "border":"1px solid #ccc", "background":"white"}),
    ], style={"width":"140px", "display":"flex", "flexDirection":"column", "alignItems":"center", 
              "justifyContent":"center", "marginLeft":"20px", "paddingLeft":"20px", "borderLeft":"1px solid #b0d4f1"})

], style={
    "display": "flex", "alignItems": "stretch", # 使用 Flex 讓左右並排
    "background":"#f0f8ff","border":"1px solid #b0d4f1",
    "borderRadius":"8px","padding":"12px 16px","marginBottom":"12px",
})

app.layout = html.Div([
    dcc.Store(id="stored-data"),
    dcc.Store(id="similarity-data"),
    dcc.Store(id="selected-indices", data=[]),
    dcc.Store(id="gallery-page", data=0), # 新增：記錄目前圖庫頁碼
    dcc.Download(id="download-csv"),

    html.H2("🧬 Pareto Front Scatter Explorer",
            style={"textAlign":"center","fontFamily":"monospace",
                   "color":"#222","margin":"20px 0 12px"}),

    # 上傳區（永遠顯示）
    html.Div([
        dcc.Upload(
            id="upload-data",
            children=html.Div([
                html.Span("📂 拖曳或點擊上傳 CSV 檔案",
                          style={"fontSize":"16px","color":"#555"}),
                html.Br(),
                html.Span("支援欄位：smiles/SMILES、Name（選填）、數值性質欄位",
                          style={"fontSize":"12px","color":"#aaa"}),
            ]),
            style={
                "width":"100%","height":"80px",
                "borderWidth":"2px","borderStyle":"dashed",
                "borderRadius":"10px","borderColor":"#ccc",
                "display":"flex","alignItems":"center","justifyContent":"center",
                "cursor":"pointer","background":"#fafafa","boxSizing":"border-box",
            },
            multiple=False,
        ),
        html.Div(id="upload-status",
                 style={"marginTop":"6px","fontSize":"13px",
                        "textAlign":"center","minHeight":"20px"}),
    ], style={"marginBottom":"12px"}),

    # Pareto 控制列（上傳後顯示）
    html.Div(id="pareto-controls", style={"display":"none"}, children=[
        html.Div([
            html.Div([
                html.Label("X 軸性質", style={"fontWeight":"bold"}),
                dcc.Dropdown(id="x-col", clearable=False),
                dcc.RadioItems([
                    {"label":" Minimize ↓","value":"minimize"},
                    {"label":" Maximize ↑","value":"maximize"},
                ], value="minimize", id="x-dir", inline=True,
                   style={"marginTop":"4px"}),
            ], style={"flex":"1","minWidth":"180px"}),
            html.Div([
                html.Label("Y 軸性質", style={"fontWeight":"bold"}),
                dcc.Dropdown(id="y-col", clearable=False),
                dcc.RadioItems([
                    {"label":" Minimize ↓","value":"minimize"},
                    {"label":" Maximize ↑","value":"maximize"},
                ], value="minimize", id="y-dir", inline=True,
                   style={"marginTop":"4px"}),
            ], style={"flex":"1","minWidth":"180px"}),
            html.Div([
                html.Label("Pareto Front 層數", style={"fontWeight":"bold"}),
                dcc.Slider(1, 5, 1, value=3, id="n-ranks",
                           marks={i:str(i) for i in range(1,6)}),
                html.Div(id="rank-legend", style={"marginTop":"6px","fontSize":"13px"}),
            ], style={"flex":"2","minWidth":"260px"}),
            html.Div([
                html.Button("🖼 開啟分子圖庫", id="btn-open-gallery", n_clicks=0,
                            style={"width":"100%","height":"50px","background":"#6f42c1",
                                   "color":"white","border":"none","borderRadius":"8px",
                                   "cursor":"pointer","fontSize":"16px","fontWeight":"bold"}),
            ], style={"flex":"1","minWidth":"150px","display":"flex","alignItems":"center"}),
        ], style={"display":"flex","flexWrap":"wrap","gap":"24px",
                  "padding":"12px 16px","background":"#f7f7f7",
                  "borderRadius":"8px","marginBottom":"12px"}),
    ]),

    # SMILES 查詢列（上傳後顯示）
    html.Div(id="smiles-panel", style={"display":"none"}, children=[SMILES_PANEL]),

    # 圖 + 分子圖（上傳後顯示）
    html.Div(id="chart-area", style={"display":"none"}, children=[
        html.Div([
            dcc.Graph(id="scatter-plot",
                      style={"flex":"3","minWidth":"400px","height":"600px"},
                      config={"displayModeBar":True,"scrollZoom":True}),
            html.Div([
                html.P("🔍 Hover 某個點以顯示分子結構", id="mol-hint",
                       style={"color":"#888","fontSize":"13px",
                              "textAlign":"center","marginTop":"40px"}),
                html.Img(id="mol-img",
                         style={"maxWidth":"100%","display":"none",
                                "borderRadius":"8px","border":"1px solid #ddd"}),
                html.Div(id="mol-info",
                         style={"fontSize":"12px","color":"#555",
                                "marginTop":"8px","padding":"0 8px",
                                "wordBreak":"break-all"}),
            ], style={"flex":"1","minWidth":"240px","maxWidth":"320px",
                      "display":"flex","flexDirection":"column","alignItems":"center",
                      "background":"#fafafa","borderRadius":"8px",
                      "border":"1px solid #e0e0e0","padding":"12px"}),
        ], style={"display":"flex","gap":"16px","alignItems":"flex-start",
                  "paddingBottom":"20px"}),
    ]),

    # 分子圖庫彈出視窗 (Modal-like)
    html.Div(id="gallery-modal", style={"display":"none"}, children=[
        html.Div([
            html.Div([
                html.Div([
                    html.H3("🖼 分子圖庫 (依 Rank 排序)", style={"margin":"0"}),
                    html.Span(id="gallery-stats", style={"marginLeft":"15px", "color":"#666", "fontSize":"14px"}),
                ], style={"display":"flex", "alignItems":"center"}),
                html.Div([
                    html.Button("◀ 上一頁", id="btn-prev-page", n_clicks=0,
                                style={"padding":"8px 12px","background":"#eee","border":"1px solid #ccc","borderRadius":"4px","cursor":"pointer","marginRight":"8px"}),
                    html.Span(id="page-display", style={"marginRight":"8px", "fontWeight":"bold"}),
                    html.Button("下一頁 ▶", id="btn-next-page", n_clicks=0,
                                style={"padding":"8px 12px","background":"#eee","border":"1px solid #ccc","borderRadius":"4px","cursor":"pointer","marginRight":"20px"}),
                    html.Button("🗑 清除勾選", id="btn-clear-selection", n_clicks=0,
                                style={"padding":"8px 16px","background":"#6c757d","color":"white",
                                       "border":"none","borderRadius":"4px","cursor":"pointer","marginRight":"12px"}),
                    html.Button("📥 匯出勾選 CSV", id="btn-export", n_clicks=0,
                                style={"padding":"8px 16px","background":"#28a745","color":"white",
                                       "border":"none","borderRadius":"4px","cursor":"pointer","marginRight":"12px"}),
                    html.Button("✕ 關閉", id="btn-close-gallery", n_clicks=0,
                                style={"padding":"8px 16px","background":"#dc3545","color":"white",
                                       "border":"none","borderRadius":"4px","cursor":"pointer"}),
                ], style={"display":"flex", "alignItems":"center"}),
            ], style={"display":"flex","justifyContent":"space-between","alignItems":"center",
                      "paddingBottom":"15px","borderBottom":"1px solid #ddd","marginBottom":"15px"}),
            
            html.Div(id="gallery-content", style={
                "display":"grid","gridTemplateColumns":"repeat(auto-fill, minmax(200px, 1fr))",
                "gap":"15px","maxHeight":"70vh","overflowY":"auto","padding":"10px"
            }),
        ], style={
            "position":"fixed","top":"50%","left":"50%","transform":"translate(-50%, -50%)",
            "width":"90%","maxWidth":"1100px","background":"white","padding":"25px",
            "borderRadius":"12px","boxShadow":"0 0 20px rgba(0,0,0,0.3)","zIndex":"1000"
        }),
        html.Div(style={
            "position":"fixed","top":"0","left":"0","width":"100%","height":"100%",
            "background":"rgba(0,0,0,0.5)","zIndex":"999"
        }, id="gallery-overlay")
    ]),

], style={"maxWidth":"1200px","margin":"0 auto",
          "padding":"0 20px","fontFamily":"sans-serif"})
