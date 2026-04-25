import io
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, html, dcc, ALL
from .app import app
from .engine import (
    compute_pareto_ranks, compute_tanimoto, normalize_smiles, 
    canonicalize_smiles, PARETO_COLORS, RDKIT_OK
)
from .utils import parse_csv, smiles_to_image_b64

# ── Callback 1：解析 CSV ──────────────────────────────────────────────────────
@app.callback(
    Output("stored-data","data"),
    Output("upload-status","children"),
    Output("upload-status","style"),
    Output("pareto-controls","style"),
    Output("smiles-panel","style"),
    Output("chart-area","style"),
    Output("x-col","options"), Output("x-col","value"),
    Output("y-col","options"), Output("y-col","value"),
    Input("upload-data","contents"),
    State("upload-data","filename"),
    prevent_initial_call=True,
)
def on_upload(contents, filename):
    hidden  = {"display":"none"}
    visible = {"display":"block"}
    ok_style  = {"marginTop":"6px","fontSize":"13px","color":"#2ECC40","textAlign":"center"}
    err_style = {"marginTop":"6px","fontSize":"13px","color":"#FF4136","textAlign":"center"}

    if not contents:
        return None,"",{},hidden,hidden,hidden,[],None,[],None

    df, err = parse_csv(contents, filename)
    if err:
        return None,f"❌ {err}",err_style,hidden,hidden,hidden,[],None,[],None

    # SMILES 規範化與過濾
    if "smiles" in df.columns:
        df = df.rename(columns={"smiles": "orig_smiles"})
        res = df["orig_smiles"].apply(canonicalize_smiles)
        df["smiles"] = res.apply(lambda x: x[0])
        df["_valid"] = res.apply(lambda x: x[1])
        dropped = len(df) - df["_valid"].sum()
        df = df[df["_valid"]].drop(columns=["_valid"])
    else:
        dropped = 0

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) < 2:
        return None, "❌ CSV 需要至少兩個數值欄位", err_style, hidden, hidden, hidden, [], None, [], None

    options = [{"label":c,"value":c} for c in numeric_cols]
    x_val, y_val = numeric_cols[0], numeric_cols[1] if len(numeric_cols)>1 else numeric_cols[0]

    note = f"✅ 已載入 {filename}｜{len(df)} 筆"
    if dropped > 0: note += f" (已過濾 {dropped} 筆無效 SMILES)"
    note += f"｜{len(numeric_cols)} 個數值欄位"
    
    return (df.to_json(date_format="iso", orient="split"),
            note, ok_style,
            visible, visible, visible,
            options, x_val, options, y_val)


# ── Callback 2：處理 SMILES 查詢按鈕 ─────────────────────────────────────────
@app.callback(
    Output("similarity-data","data"),
    Output("query-status","children"),
    Output("query-mol-img", "src"),
    Output("query-mol-img", "style"),
    Input("btn-check","n_clicks"),
    Input("btn-tanimoto","n_clicks"),
    Input("btn-clear-sim","n_clicks"),
    State("query-smiles","value"),
    State("stored-data","data"),
    prevent_initial_call=True,
)
def handle_query(n_check, n_tanimoto, n_clear, query_smi, stored):
    triggered = ctx.triggered_id
    red   = {"color":"#FF4136"}
    green = {"color":"#2ECC40"}
    gray  = {"color":"#888"}
    
    img_hidden = {"maxWidth":"120px", "maxHeight":"120px", "display":"none"}
    img_visible = {"maxWidth":"120px", "maxHeight":"120px", "display":"block", "borderRadius":"6px", "border":"1px solid #ccc", "background":"white"}

    if triggered == "btn-clear-sim":
        return None, html.Span("已清除 Similarity 模式", style=gray), "", img_hidden
    if not stored:
        return None, html.Span("❌ 請先上傳 CSV", style=red), "", img_hidden
    if not query_smi or not query_smi.strip():
        return None, html.Span("❌ 請先輸入 SMILES", style=red), "", img_hidden
    if not RDKIT_OK:
        return None, html.Span("❌ 未安裝 rdkit，無法計算", style=red), "", img_hidden

    q = normalize_smiles(query_smi.strip())
    from rdkit import Chem 
    if not q or Chem.MolFromSmiles(q) is None:
        return None, html.Span("❌ SMILES 格式無效", style=red), "", img_hidden

    df = pd.read_json(io.StringIO(stored), orient="split")
    if "smiles" not in df.columns:
        return None, html.Span("❌ 資料中無 SMILES 欄位", style=red), "", img_hidden

    b64_query = smiles_to_image_b64(q, size=(120, 120))
    img_src = f"data:image/png;base64,{b64_query}" if b64_query else ""
    img_style = img_visible if b64_query else img_hidden

    smiles_list = df["smiles"].fillna("").tolist()
    matched = [i for i,s in enumerate(smiles_list) if normalize_smiles(s)==q]

    if triggered == "btn-check":
        if matched:
            names = [str(df.iloc[i].get("Name","")) for i in matched
                     if "Name" in df.columns and df.iloc[i].get("Name","")]
            label = "、".join(names) if names else f"第 {[m+1 for m in matched]} 筆"
            msg = html.Span([html.B(f"✅ 找到 {len(matched)} 筆相符："), f" {label}"], style=green)
        else:
            msg = html.Span("🔍 資料中不存在此分子", style={"color":"#FF851B"})
        return {"mode":"check","query_smiles":q,"matched_indices":matched,"similarities":None}, msg, img_src, img_style

    if triggered == "btn-tanimoto":
        sims = compute_tanimoto(q, smiles_list)
        if sims is None:
            return None, html.Span("❌ Fingerprint 計算失敗", style=red), img_src, img_style
        valid = [s for s in sims if s >= 0]
        msg = html.Span([
            html.B("✅ Tanimoto 計算完成"),
            f"  最高: {max(valid):.3f}  平均: {np.mean(valid):.3f}" if valid else "",
            html.Span("  ★ = 完全相符" if matched else "", style={"color":"#FF4136","marginLeft":"8px"}),
        ], style=green)
        return {"mode":"tanimoto","query_smiles":q,"matched_indices":matched,"similarities":sims}, msg, img_src, img_style

    return None, "", "", img_hidden


# ── Callback 3：更新散布圖 ────────────────────────────────────────────────────
@app.callback(
    Output("scatter-plot","figure"),
    Output("rank-legend","children"),
    Input("x-col","value"), Input("y-col","value"),
    Input("x-dir","value"), Input("y-dir","value"),
    Input("n-ranks","value"),
    Input("similarity-data","data"),
    Input("selected-indices", "data"),
    State("stored-data","data"),
    prevent_initial_call=True,
)
def update_scatter(x_col, y_col, x_dir, y_dir, n_ranks, sim_data, selected_indices, stored):
    if not stored or not x_col or not y_col:
        return go.Figure(), []

    df = pd.read_json(io.StringIO(stored), orient="split")
    has_name   = "Name"   in df.columns
    has_smiles = "smiles" in df.columns

    sub     = df[[x_col,y_col]].dropna()
    work_df = df.loc[sub.index].copy()
    work_df["_rank"] = compute_pareto_ranks(work_df, x_col, y_col, x_dir, y_dir, n_ranks)

    orig_indices = list(df.index)

    def get_labels(d):
        return d["Name"].fillna("") if has_name else pd.Series("", index=d.index)

    def get_cdata(d, smap=None):
        return d.apply(lambda row: {
            "all_data": row.to_dict(),
            "smiles": row.get("smiles", "") if has_smiles else "",
            "name": row.get("Name", "") if has_name else "",
            "rank": row.get("_rank", 0),
            "sim": smap.get(row.name, -1) if smap else -1,
            "index": row.name
        }, axis=1).tolist()

    def draw_pareto(fig):
        m0 = work_df["_rank"]==0
        if m0.any():
            s0 = work_df[m0]
            fig.add_trace(go.Scatter(
                x=s0[x_col], y=s0[y_col], mode="markers+text",
                marker=dict(color=PARETO_COLORS[0],size=7,opacity=0.5,
                            line=dict(width=0.5,color="white")),
                text=get_labels(s0), textposition="top center",
                textfont=dict(size=9,color="#999"),
                name="其他", customdata=get_cdata(s0),
                hovertemplate=(f"<b>{x_col}</b>: %{{x:.3f}}<br>"
                               f"<b>{y_col}</b>: %{{y:.3f}}<br>"
                               "<b>Name</b>: %{customdata.name}<extra></extra>"),
            ))
        for rank in range(1, n_ranks+1):
            mr = work_df["_rank"]==rank
            if not mr.any(): continue
            sr = work_df[mr]; color = PARETO_COLORS.get(rank,"#888")
            fig.add_trace(go.Scatter(
                x=sr[x_col], y=sr[y_col], mode="markers+text",
                marker=dict(color=color,size=10,opacity=0.9,
                            line=dict(width=1,color="white")),
                text=get_labels(sr), textposition="top center",
                textfont=dict(size=9,color=color),
                name=f"Rank {rank}", customdata=get_cdata(sr),
                hovertemplate=(f"<b>Rank {rank}</b><br>"
                               f"<b>{x_col}</b>: %{{x:.3f}}<br>"
                               f"<b>{y_col}</b>: %{{y:.3f}}<br>"
                               "<b>Name</b>: %{customdata.name}<extra></extra>"),
            ))

    def add_star(fig, matched_orig):
        sm = work_df[work_df.index.isin(matched_orig)]
        if sm.empty: return
        fig.add_trace(go.Scatter(
            x=sm[x_col], y=sm[y_col], mode="markers",
            marker=dict(symbol="star",size=20,color="#FF4136",
                        line=dict(width=1.5,color="#8B0000")),
            name="★ 完全相符",
            customdata=get_cdata(sm),
            hovertemplate=(f"<b>★ 完全相符</b><br>"
                           f"<b>{x_col}</b>: %{{x:.3f}}<br>"
                           f"<b>{y_col}</b>: %{{y:.3f}}<extra></extra>"),
        ))

    fig = go.Figure()
    legend_items = []

    if sim_data and sim_data.get("mode")=="tanimoto" and sim_data.get("similarities"):
        sims_list   = sim_data["similarities"]
        matched_idx = set(sim_data.get("matched_indices",[]))
        smap = {orig_indices[i]: sims_list[i]
                for i in range(min(len(orig_indices),len(sims_list)))}
        work_df["_sim"] = work_df.index.map(lambda i: smap.get(i,-1))

        invalid = work_df[work_df["_sim"] < 0]
        if not invalid.empty:
            fig.add_trace(go.Scatter(
                x=invalid[x_col], y=invalid[y_col], mode="markers",
                marker=dict(color="#ddd",size=6,opacity=0.4),
                name="無法計算", showlegend=True,
                customdata=get_cdata(invalid, smap),
                hovertemplate="無 SMILES 資料<extra></extra>",
            ))

        valid = work_df[work_df["_sim"] >= 0]
        if not valid.empty:
            fig.add_trace(go.Scatter(
                x=valid[x_col], y=valid[y_col], mode="markers+text",
                marker=dict(
                    color=valid["_sim"].values, colorscale="Viridis",
                    cmin=0, cmax=1, size=9, opacity=0.9,
                    line=dict(width=0.5,color="white"),
                    colorbar=dict(
                        title=dict(text="Tanimoto<br>Similarity",side="right"),
                        thickness=16, len=0.75, x=1.02, tickformat=".2f",
                    ),
                    showscale=True,
                ),
                text=get_labels(valid), textposition="top center",
                textfont=dict(size=9,color="#444"),
                name="Similarity", customdata=get_cdata(valid, smap),
                hovertemplate=(f"<b>{x_col}</b>: %{{x:.3f}}<br>"
                               f"<b>{y_col}</b>: %{{y:.3f}}<br>"
                               "<b>Name</b>: %{customdata.name}<br>"
                               "<b>Tanimoto</b>: %{customdata.sim:.3f}"
                               "<extra></extra>"),
            ))

        matched_orig = {orig_indices[m] for m in matched_idx if m < len(orig_indices)}
        add_star(fig, matched_orig)
        legend_items = [html.Span("● Tanimoto Similarity 模式",
                                   style={"color":"#555","fontStyle":"italic"})]

    elif sim_data and sim_data.get("mode")=="check":
        draw_pareto(fig)
        matched_orig = {orig_indices[m] for m in sim_data.get("matched_indices",[])
                        if m < len(orig_indices)}
        add_star(fig, matched_orig)
        legend_items = ([html.Span(f"■ Rank {r}  ",
                          style={"color":PARETO_COLORS.get(r,"#888"),
                                 "fontWeight":"bold","marginRight":"8px"})
                         for r in range(1, n_ranks+1)]
                        + [html.Span("■ 其他",style={"color":PARETO_COLORS[0]})])

    else:
        draw_pareto(fig)
        legend_items = ([html.Span(f"■ Rank {r}  ",
                          style={"color":PARETO_COLORS.get(r,"#888"),
                                 "fontWeight":"bold","marginRight":"8px"})
                         for r in range(1, n_ranks+1)]
                        + [html.Span("■ 其他",style={"color":PARETO_COLORS[0]})])

    if selected_indices:
        sel_df = work_df[work_df.index.isin(selected_indices)]
        if not sel_df.empty:
            fig.add_trace(go.Scatter(
                x=sel_df[x_col], y=sel_df[y_col], mode="markers",
                marker=dict(symbol="circle-open", size=18, color="#FF4136", line=dict(width=3)),
                name="已勾選", showlegend=True, hoverinfo="skip"
            ))

    fig.update_layout(
        xaxis_title=f"{x_col}  ({'→ maximize' if x_dir=='maximize' else '← minimize'})",
        yaxis_title=f"{y_col}  ({'↑ maximize' if y_dir=='maximize' else '↓ minimize'})",
        legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="right",x=1),
        margin=dict(l=60,r=80,t=50,b=50),
        plot_bgcolor="white", paper_bgcolor="white",
        hovermode="closest", uirevision="constant",
    )
    fig.update_xaxes(showgrid=True, gridcolor="#eee", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="#eee", zeroline=False)
    return fig, legend_items


# ── Callback 4：Hover 顯示分子圖 ──────────────────────────────────────────────
@app.callback(
    Output("mol-img","src"), Output("mol-img","style"),
    Output("mol-info","children"), Output("mol-hint","style"),
    Input("scatter-plot","hoverData"),
    State("similarity-data","data"),
    prevent_initial_call=True,
)
def show_molecule(hover_data, sim_data):
    hidden  = {"maxWidth":"100%","display":"none","borderRadius":"8px","border":"1px solid #ddd"}
    visible = {"maxWidth":"100%","display":"block","borderRadius":"8px","border":"1px solid #ddd"}
    hint_s  = {"color":"#888","fontSize":"13px","textAlign":"center","marginTop":"40px"}
    hint_h  = {"display":"none"}

    if not hover_data: return "","",hint_s ,""

    try:
        cd = hover_data["points"][0].get("customdata", {})
        if not cd: return "", hidden, "", hint_s
        
        all_data = cd.get("all_data", {})
        smiles = cd.get("smiles", "")
        name = cd.get("name", "")
        rank = cd.get("rank", 0)
        sim = cd.get("sim", -1)
    except: return "",hidden,"",hint_s

    if not smiles: return "",hidden,"（無 SMILES 資料）",hint_h
    if not RDKIT_OK: return "",hidden,f"未安裝 rdkit",hint_h

    b64 = smiles_to_image_b64(smiles)
    if b64 is None: return "",hidden,f"無法解析 SMILES",hint_h

    props = []
    for k, v in all_data.items():
        if k in ["smiles", "Name", "_rank", "_sim", "orig_smiles"]: continue
        val_str = f"{v:.4f}" if isinstance(v, (float, np.float64)) else str(v)
        props.append(html.Div([
            html.Span(f"{k}: ", style={"color":"#888", "fontSize":"11px"}),
            html.Span(val_str, style={"fontWeight":"bold", "color":"#333"})
        ], style={"borderBottom":"1px solid #f0f0f0", "padding":"2px 0"}))

    info = [
        html.Div(name if name else "（無名稱）", style={"fontWeight":"bold", "fontSize":"14px", "marginBottom":"4px"}),
        html.Div(f"Rank {rank}" if rank>0 else "非 Pareto Front",
                  style={"color":PARETO_COLORS.get(rank,"#888"),"fontWeight":"bold", "fontSize":"13px"}),
    ]
    if sim_data and sim_data.get("mode")=="tanimoto" and sim>=0:
        info += [html.Div(f"Tanimoto: {sim:.4f}",
                           style={"color":"#0074D9","fontWeight":"bold","fontSize":"13px"})]
    info += [
        html.Div(props, style={"marginTop":"8px", "maxHeight":"160px", "overflowY":"auto"}),
        html.Div(f"SMILES: {smiles}", style={"fontSize":"10px","color":"#bbb", "marginTop":"8px", "wordBreak":"break-all"})
    ]

    return f"data:image/png;base64,{b64}", visible, info, hint_h


# ── Callback 5：圖庫彈窗切換 ───────────────────────────────────────────────
@app.callback(
    Output("gallery-modal", "style"),
    Output("gallery-page", "data", allow_duplicate=True),
    Input("btn-open-gallery", "n_clicks"),
    Input("btn-close-gallery", "n_clicks"),
    Input("gallery-overlay", "n_clicks"),
    prevent_initial_call=True,
)
def toggle_gallery(n_open, n_close, n_overlay):
    triggered = ctx.triggered_id
    if triggered == "btn-open-gallery":
        return {"display":"block"}, 0
    return {"display":"none"}, 0


# ── Callback 6：處理分頁切換 ─────────────────────────────────────────────────
@app.callback(
    Output("gallery-page", "data"),
    Input("btn-prev-page", "n_clicks"),
    Input("btn-next-page", "n_clicks"),
    State("gallery-page", "data"),
    prevent_initial_call=True,
)
def change_page(n_prev, n_next, current_page):
    triggered = ctx.triggered_id
    if triggered == "btn-prev-page":
        return max(0, current_page - 1)
    return current_page + 1


# ── Callback 7：生成圖庫內容與分頁顯示 ─────────────────────────────────────────
@app.callback(
    Output("gallery-content", "children"),
    Output("page-display", "children"),
    Output("gallery-stats", "children"),
    Input("gallery-page", "data"),
    Input("gallery-modal", "style"),
    Input("btn-clear-selection", "n_clicks"), 
    State("stored-data", "data"),
    State("x-col", "value"), State("y-col", "value"),
    State("x-dir", "value"), State("y-dir", "value"),
    State("n-ranks", "value"),
    State("selected-indices", "data"),
    prevent_initial_call=True,
)
def render_gallery(page, modal_style, n_clear, stored, x_col, y_col, x_dir, y_dir, n_ranks, selected):
    if not modal_style or modal_style.get("display") != "block":
        return "", "", ""
        
    if not stored:
        return "請先上傳資料", "", ""
    
    triggered = ctx.triggered_id
    effective_selected = selected if triggered != "btn-clear-selection" else []
    
    try:
        df = pd.read_json(io.StringIO(stored), orient="split")
        work_df = df.copy()
        work_df["_rank"] = compute_pareto_ranks(work_df, x_col, y_col, x_dir, y_dir, n_ranks)
        
        work_df["_sort_weight"] = work_df["_rank"].apply(lambda r: r if r > 0 else 999)
        all_df = work_df.sort_values("_sort_weight")
        
        PAGE_SIZE = 50
        total_mols = len(all_df)
        total_pages = max(1, (total_mols + PAGE_SIZE - 1) // PAGE_SIZE)
        
        current_page = min(max(0, page if page is not None else 0), total_pages - 1)
        
        start_idx = current_page * PAGE_SIZE
        end_idx = start_idx + PAGE_SIZE
        gallery_df = all_df.iloc[start_idx:end_idx]
        
        cards = []
        for idx, row in gallery_df.iterrows():
            smiles = row.get("smiles", "")
            name = row.get("Name", "")
            rank = row["_rank"]
            x_val = row[x_col]
            y_val = row[y_col]
            b64 = smiles_to_image_b64(smiles, size=(120, 120))
            
            x_str = f"{x_val:.3f}" if isinstance(x_val, (float, np.float64)) else str(x_val)
            y_str = f"{y_val:.3f}" if isinstance(y_val, (float, np.float64)) else str(y_val)

            cards.append(html.Div([
                html.Img(src=f"data:image/png;base64,{b64}" if b64 else "", 
                         style={"width":"100%", "borderRadius":"4px"}),
                html.Div([
                    dcc.Checklist(
                        options=[{"label": "", "value": idx}],
                        value=[idx] if idx in effective_selected else [],
                        id={"type": "gallery-check", "index": idx},
                        style={"display":"inline-block", "marginRight":"8px"}
                    ),
                    html.Span(name if name else f"ID: {idx}", 
                              style={"fontSize":"11px", "fontWeight":"bold", "overflow":"hidden", "textOverflow":"ellipsis", "whiteSpace":"nowrap"}),
                ], style={"display":"flex", "alignItems":"center", "marginTop":"5px"}),
                html.Div(f"Rank {rank}" if rank > 0 else "非 Pareto", 
                         style={"fontSize":"10px", "color":PARETO_COLORS.get(rank, "#888"), "fontWeight":"bold"}),
                html.Div([
                    html.Div(f"{x_col}: {x_str}", style={"fontSize":"9px", "color":"#555"}),
                    html.Div(f"{y_col}: {y_str}", style={"fontSize":"9px", "color":"#555"}),
                ], style={"marginTop":"4px", "borderTop":"1px solid #f9f9f9", "paddingTop":"4px"})
            ], style={
                "border":"1px solid #eee", "padding":"8px", "borderRadius":"8px", 
                "background":"#fff", "boxShadow":"0 2px 5px rgba(0,0,0,0.05)"
            }))
        
        page_info = f"第 {current_page + 1} / {total_pages} 頁"
        stats_info = f"共 {total_mols} 個分子 (已選取 {len(effective_selected)} 個)"
            
        return cards, page_info, stats_info
    except Exception as e:
        return f"渲染失敗: {str(e)}", "錯誤", ""


# ── Callback 8：同步勾選狀態 (跨頁勾選支援) ──────────────────────────────────────────
@app.callback(
    Output("selected-indices", "data"),
    Input({"type": "gallery-check", "index": ALL}, "value"),
    Input("btn-clear-selection", "n_clicks"),
    State("selected-indices", "data"),
    prevent_initial_call=True,
)
def update_selection(values, n_clear, current_selected):
    triggered = ctx.triggered_id
    if triggered == "btn-clear-selection":
        return []
    if values:
        triggered_ids = [item["id"]["index"] for item in ctx.inputs_list[0]]
        page_selected = [v[0] for v in values if v]
        new_selected = set(current_selected) if current_selected else set()
        for tid in triggered_ids:
            if tid in new_selected:
                new_selected.remove(tid)
        for psid in page_selected:
            new_selected.add(psid)
        return list(new_selected)
    return current_selected


# ── Callback 9：匯出 CSV ─────────────────────────────────────────────────────
@app.callback(
    Output("download-csv", "data"),
    Input("btn-export", "n_clicks"),
    State("selected-indices", "data"),
    State("stored-data", "data"),
    prevent_initial_call=True,
)
def export_selected(n_clicks, selected, stored):
    if not selected or not stored: return None
    df = pd.read_json(io.StringIO(stored), orient="split")
    export_df = df.loc[selected]
    return dcc.send_data_frame(export_df.to_csv, "pareto_selected.csv", index=False)
