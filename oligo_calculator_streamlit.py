"""
Oligonucleotide Impurity MW Calculator - Streamlit Version
运行方式: streamlit run oligo_calculator_streamlit.py
"""

import streamlit as st
import pandas as pd

# ─────────────────────────────────────────────
# 数据字典
# ─────────────────────────────────────────────

def get_oligo_dict():
    return {
        # MOE nucleotides
        "MOE G": 419.06646, "MOE G P=O": 403.0893, "MOE G nucleoside": 341.13353,
        "MOE A": 403.0715, "MOE A P=O": 387.09438, "MOE A nucleoside": 325.13861,
        "MOE MeU": 394.05997, "MOE MeU P=O": 378.08282, "MOE MeU nucleoside": 316.12705,
        "MOE MeC": 393.07596, "MOE MeC P=O": 377.0988, "MOE MeC nucleoside": 315.14304,
        # DNA nucleotides
        "dG": 345.02968, "dG P=O": 239.05252, "dG nucleoside": 267.09675,
        "dA": 329.03476, "dA P=O": 313.05761, "dA nucleoside": 251.10184,
        "dT": 320.02319, "dT P=O": 304.04604, "dT nucleoside": 242.09027,
        "dMeC": 319.03918, "dMeC P=O": 303.06202, "dMeC nucleoside": 241.10626,
        "dC P=O": 289.04637, "dC nucleoside": 227.09061,
        # NMA nucleotides
        "NMA G": 432.0617, "NMA A": 416.06679, "NMA G P=O": 416.08455, "NMA G nucleoside": 354.12878,
        "NMA MeU": 407.05522, "NMA A P=O": 400.08963, "NMA A nucleoside": 338.13387,
        "NMA MeC": 406.07121, "NMA MeU P=O": 391.07807, "NMA MeU nucleoside": 329.1223,
        "NMA MeC P=O": 390.09405, "NMA MeC nucleoside": 328.13828,
        # m-nucleotides
        "mG": 375.05024, "mA": 359.04533, "mG P=O": 359.06308, "mG nucleoside": 297.10732,
        "mU": 336.01811, "mA P=O": 343.06817, "mA nucleoside": 281.1124,
        "mC": 335.03409, "mU P=O": 320.04095, "mU nucleoside": 258.08519,
        "mC P=O": 319.05694, "mC nucleoside": 257.10117,
        # cEt nucleotides
        "cEt G": 387.04024, "cEt A": 371.04533, "cEt G P=O": 371.06308, "cEt G nucleoside": 309.10732,
        "cEt MeU": 362.03376, "cEt A P=O": 355.06817, "cEt A nucleoside": 293.1124,
        "cEt MeC": 361.04974, "cEt MeU P=O": 346.0566, "cEt MeU nucleoside": 284.10084,
        "cEt MeC P=O": 345.07259, "cEt MeC nucleoside": 283.11682,
        # f-nucleotides
        "fG": 363.02025, "fA": 347.02534, "fG P=O": 347.0431, "fG nucleoside": 285.08733,
        "fU": 323.99812, "fA P=O": 331.04818, "fA nucleoside": 269.09242,
        "fC": 323.01411, "fU P=O": 308.02097, "fU nucleoside": 246.0652,
        "fC P=O": 307.03695, "fC nucleoside": 245.08118,
        # Other
        "AH": 195.04829, "AH P=O": 179.07113,
        "GalNAc": 1518.782075,
    }

# ─────────────────────────────────────────────
# 计算逻辑（与原版完全一致）
# ─────────────────────────────────────────────

def calculate_molecular_weight(sequence):
    oligo_dict = get_oligo_dict()
    components = [c.strip() for c in sequence.split("-")]
    total_weight = 0
    component_weights = []
    missing_components = []
    for component in components:
        if component in oligo_dict:
            weight = oligo_dict[component]
            total_weight += weight
            component_weights.append((component, weight))
        else:
            missing_components.append(component)
    found_all = len(missing_components) == 0
    return total_weight, component_weights, found_all, missing_components

def check_po_impurity(components):
    po_eligible = {
        "MOE G","MOE A","MOE MeU","MOE MeC",
        "dG","dA","dT","dMeC",
        "NMA G","NMA A","NMA MeU","NMA MeC",
        "mG","mA","mU","mC",
        "cEt G","cEt A","cEt MeU","cEt MeC",
        "fG","fA","fU","fC",
    }
    return any(c in po_eligible for c in components)

def calculate_impurities(sequence, average_mw):
    oligo_dict = get_oligo_dict()
    components = [c.strip() for c in sequence.split("-")]
    if not all(c in oligo_dict for c in components):
        return [], [], []

    unique_components = []
    for c in components:
        if c not in unique_components and c in oligo_dict:
            unique_components.append(c)

    n_minus_1 = [(c, average_mw - oligo_dict[c]) for c in unique_components]
    n_plus_1  = [(c, average_mw + oligo_dict[c]) for c in unique_components]

    if check_po_impurity(components):
        n_plus_1.append(("P=O", average_mw - 15.977156))

    moe_set = {"MOE G","MOE A","MOE MeU","MOE MeC",
               "MOE G P=O","MOE A P=O","MOE MeU P=O","MOE MeC P=O",
               "MOE G nucleoside","MOE A nucleoside","MOE MeU nucleoside","MOE MeC nucleoside"}
    if any(c in moe_set for c in components):
        n_plus_1.append(("2-O-CH3", average_mw - 44.026215))

    a_nuc = {"MOE A","dA","NMA A","mA","cEt A","fA",
             "MOE A P=O","dA P=O","NMA A P=O","mA P=O","cEt A P=O","fA P=O",
             "MOE A nucleoside","dA nucleoside","NMA A nucleoside","mA nucleoside","cEt A nucleoside","fA nucleoside"}
    c_nuc = {"MOE MeC","dMeC","NMA MeC","cEt MeC",
             "MOE MeC P=O","dMeC P=O","NMA MeC P=O","cEt MeC P=O",
             "MOE MeC nucleoside","dMeC nucleoside","NMA MeC nucleoside","cEt MeC nucleoside"}
    g_nuc = {"MOE G","dG","NMA G","mG","cEt G","fG",
             "MOE G P=O","dG P=O","NMA G P=O","mG P=O","cEt G P=O","fG P=O",
             "MOE G nucleoside","dG nucleoside","NMA G nucleoside","mG nucleoside","cEt G nucleoside","fG nucleoside"}
    t_nuc = {"dT","dT P=O","dT nucleoside"}
    u_nuc = {"MOE MeU","NMA MeU","cEt MeU",
             "MOE MeU P=O","NMA MeU P=O","cEt MeU P=O",
             "MOE MeU nucleoside","NMA MeU nucleoside","cEt MeU nucleoside"}

    has_a = any(c in a_nuc for c in components)
    has_c = any(c in c_nuc for c in components)
    has_g = any(c in g_nuc for c in components)
    has_t = any(c in t_nuc for c in components)
    has_u = any(c in u_nuc for c in components)
    has_galnac = "GalNAc" in components

    if has_a:
        n_plus_1 += [("Debase-A", average_mw - 117.051755), ("AMPA", average_mw + 98.0368)]
    if has_c:
        n_plus_1 += [("Debase-C", average_mw - 107.056172), ("OPC", average_mw + 80.0262),
                     ("Deamination", average_mw + 0.984),
                     ("DMT-C-phosphonates (Thiolation incomplete)", average_mw + 270.1586),
                     ("DMT-C-phosphonates (Oxygen generation incomplete)", average_mw + 286.1358)]
    if has_g:
        n_plus_1 += [("Debase-G", average_mw - 133.04667), ("ADP", average_mw + 41.026549),
                     ("n+iBu", average_mw + 70.0419), ("IDP", average_mw + 69.0579)]
    if has_t:
        n_plus_1 += [("Debase-T", average_mw - 108.040188), ("CNET", average_mw + 53.034374)]
    if has_u:
        n_plus_1.append(("Debase-U", average_mw - 108.040188))
    if has_galnac:
        n_plus_1 += [
            ("Penta GalNAc",       average_mw + 923.4951),
            ("Tetra GalNAc",       average_mw + 549.2898),
            ("THA Branch Lost",    average_mw + 347.2053),
            ("HexGal",             average_mw - 303.1682),
            ("GalNAc+H2O",         average_mw - 203.0794),
            ("n-TrisGalNAc",       average_mw - 1226.6632),
            ("AHI+Ac",             average_mw - 1297.7003),
            ("n+p(AH)",            average_mw - 179.0712),
            ("n-Ade+TEA",          average_mw - 33.9341),
            ("GalNAc-Ac",          average_mw - 42.0105),
            ("Aminohexyl Phosphate", average_mw - 1339.7109),
        ]

    n_plus_1 += [("Disulfate/Sulfate", average_mw + 15.977156), ("Uny-TP", average_mw + 276.9811)]

    mobile_phase = [
        ("Mobile phase false peak 1", average_mw + 219.2),
        ("Mobile phase false peak 2", average_mw + 234.8),
        ("Mobile phase false peak 3", average_mw + 255.2),
        ("Mobile phase false peak 4", average_mw + 467.2),
        ("TBuA",    average_mw + 185.2143),
        ("TBuA+O",  average_mw + 201.2092),
    ]
    return n_minus_1, n_plus_1, mobile_phase

def get_all_results(sequence, average_mw, charge):
    SPECIAL_ORDER = [
        "P=O","2-O-CH3","CNET","Debase-A","Debase-C","Debase-G","Debase-T","Debase-U",
        "Disulfate/Sulfate","ADP","AMPA","OPC","dGalNAc","Deamination",
        "DMT-C-phosphonates (Thiolation incomplete)","DMT-C-phosphonates (Oxygen generation incomplete)",
        "n+iBu","IDP","Uny-TP",
        "Penta GalNAc","Tetra GalNAc","THA Branch Lost","HexGal","GalNAc+H2O",
        "n-TrisGalNAc","AHI+Ac","n+p(AH)","n-Ade+TEA","GalNAc-Ac","Aminohexyl Phosphate",
    ]
    total_weight, component_weights, found_all, missing = calculate_molecular_weight(sequence)
    if not found_all:
        return {"success": False, "missing": missing}

    n_minus_1, n_plus_1, mobile_phase = calculate_impurities(sequence, average_mw)

    special_set = set(SPECIAL_ORDER)
    all_impurities = []
    special_added = set()
    for t in SPECIAL_ORDER:
        for name, w in n_plus_1:
            if name == t and name not in special_added:
                all_impurities.append((name, w))
                special_added.add(name)
    seq_added = set()
    for name, w in n_plus_1:
        if name not in special_set and name not in seq_added:
            all_impurities.append((f"N+{name}", w))
            seq_added.add(name)
    for name, w in n_minus_1:
        all_impurities.append((f"N-{name}", w))

    # charge state calculations
    p, na, k, fe, n = 1.007825, 22.989769, 38.963707, 55.9349, charge
    def cs(mw, adduct=0): return (mw + adduct - n * p) / n

    std  = [("Average MW", cs(average_mw))]  + [(nm, cs(mw))        for nm, mw in all_impurities]
    naA  = [("Average MW", cs(average_mw, na))] + [(nm, cs(mw, na)) for nm, mw in all_impurities]
    kA   = [("Average MW", cs(average_mw, k))]  + [(nm, cs(mw, k))  for nm, mw in all_impurities]
    feA  = [("Average MW", cs(average_mw, fe))] + [(nm, cs(mw, fe)) for nm, mw in all_impurities]
    mob  = [(nm, cs(mw)) for nm, mw in mobile_phase]

    return {
        "success": True,
        "sequence": sequence,
        "calculated_mw": total_weight,
        "average_mw": average_mw,
        "charge": charge,
        "component_weights": component_weights,
        "all_impurities": all_impurities,
        "mobile_phase": mobile_phase,
        "std": std, "naA": naA, "kA": kA, "feA": feA, "mob": mob,
    }

# ─────────────────────────────────────────────
# Streamlit 界面
# ─────────────────────────────────────────────

def to_df(rows, prefix=""):
    return pd.DataFrame(
        [(f"{prefix}{name}", f"{mw:.5f}") for name, mw in rows],
        columns=["杂质名称", "m/z (Da)"]
    )

def main():
    st.set_page_config(
        page_title="寡核苷酸杂质MW计算器",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 寡核苷酸杂质分子量计算器")
    st.caption("Oligonucleotide Impurity Molecular Weight Calculator · v2.8 Simple Addition Mode")

    oligo_dict = get_oligo_dict()

    # ── 初始化 session_state ──
    if "seq_input" not in st.session_state:
        st.session_state["seq_input"] = ""

    # 处理按钮动作（在 text_input 渲染之前修改其 session_state）
    if "_action" in st.session_state:
        action = st.session_state.pop("_action")
        current = st.session_state.get("seq_input", "").strip()
        if action == "clear":
            st.session_state["seq_input"] = ""
        elif action == "undo":
            parts = [p.strip() for p in current.split("-") if p.strip()]
            if parts:
                parts.pop()
            st.session_state["seq_input"] = "-".join(parts)
        else:
            # action 是要追加的组分名称
            if current:
                st.session_state["seq_input"] = current + "-" + action
            else:
                st.session_state["seq_input"] = action

    # ── 侧边栏：字典参考 ──
    with st.sidebar:
        st.header("📖 核苷酸字典")
        groups = {
            "MOE":       [k for k in oligo_dict if k.startswith("MOE")],
            "DNA (d)":   [k for k in oligo_dict if k.startswith("d")],
            "NMA":       [k for k in oligo_dict if k.startswith("NMA")],
            "m-nuc":     [k for k in oligo_dict if k.startswith("m") and not k.startswith("MOE")],
            "cEt":       [k for k in oligo_dict if k.startswith("cEt")],
            "f-nuc":     [k for k in oligo_dict if k.startswith("f")],
            "Other":     ["AH", "AH P=O", "GalNAc"],
        }
        for group_name, keys in groups.items():
            with st.expander(group_name):
                df = pd.DataFrame(
                    [(k, f"{oligo_dict[k]:.5f}") for k in sorted(keys)],
                    columns=["名称", "MW (Da)"]
                )
                st.dataframe(df, hide_index=True, use_container_width=True)

    # ── 主区域：输入 ──
    st.subheader("输入参数")
    col1, col2, col3 = st.columns([3, 1.5, 1])

    with col1:
        sequence = st.text_input(
            "序列（用连字符 - 分隔各组分）",
            placeholder="例如: MOE G-dA-dT-MOE MeC-dG",
            help="组分名称必须与字典中完全一致。可手动输入，也可用下方按钮快速构建。",
            key="seq_input",
        )
    with col2:
        average_mw = st.number_input(
            "平均分子量 Average MW (Da)",
            min_value=0.0, value=0.0, step=0.00001, format="%.5f",
        )
    with col3:
        charge = st.number_input(
            "电荷态 Charge (z)",
            min_value=1, value=5, step=1,
        )

    # 快速添加组分
    with st.expander("🔧 快速构建序列（点击添加组分）"):
        # 撤销 & 清空按钮
        undo_col, clear_col, preview_col = st.columns([1, 1, 4])
        with undo_col:
            if st.button("↩️ 撤销上一个", use_container_width=True):
                st.session_state["_action"] = "undo"
                st.rerun()
        with clear_col:
            if st.button("🗑️ 清空序列", use_container_width=True):
                st.session_state["_action"] = "clear"
                st.rerun()

        current_seq = st.session_state.get("seq_input", "").strip()
        if current_seq:
            parts = [p.strip() for p in current_seq.split("-") if p.strip()]
            st.success(f"当前序列 ({len(parts)} 个组分): `{current_seq}`")

        groups = {
            "MOE":     [k for k in oligo_dict if k.startswith("MOE")],
            "DNA (d)": [k for k in oligo_dict if k.startswith("d")],
            "NMA":     [k for k in oligo_dict if k.startswith("NMA")],
            "m-nuc":   [k for k in oligo_dict if k.startswith("m") and not k.startswith("MOE")],
            "cEt":     [k for k in oligo_dict if k.startswith("cEt")],
            "f-nuc":   [k for k in oligo_dict if k.startswith("f")],
            "Other":   ["AH", "AH P=O", "GalNAc"],
        }
        for gname, keys in groups.items():
            st.caption(gname)
            cols = st.columns(6)
            for i, k in enumerate(sorted(keys)):
                if cols[i % 6].button(k, key=f"btn_{k}", use_container_width=True):
                    st.session_state["_action"] = k
                    st.rerun()

    st.divider()

    # ── 计算按钮 ──
    calc_btn = st.button("🔬 开始计算", type="primary", use_container_width=False)

    if calc_btn:
        if not sequence.strip():
            st.error("请输入序列。")
        elif average_mw <= 0:
            st.error("请输入有效的平均分子量。")
        else:
            results = get_all_results(sequence.strip(), average_mw, int(charge))

            if not results["success"]:
                missing = results["missing"]
                suggestions = []
                for m in missing:
                    sugg = [k for k in oligo_dict if m.lower() in k.lower()]
                    suggestions.append(f"**{m}** → 建议: {', '.join(sugg[:4]) if sugg else '未找到匹配项'}")
                st.error("以下组分在字典中未找到：\n\n" + "\n\n".join(suggestions))
            else:
                # ── 摘要卡片 ──
                st.subheader("计算结果")
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("计算分子量 (Da)", f"{results['calculated_mw']:.5f}")
                m2.metric("输入平均MW (Da)", f"{results['average_mw']:.5f}")
                m3.metric("电荷态 z", results["charge"])
                m4.metric("序列组分数", len(results["component_weights"]))

                # ── 组分明细 ──
                with st.expander("组分分子量明细"):
                    cw_df = pd.DataFrame(
                        results["component_weights"],
                        columns=["组分", "MW (Da)"]
                    )
                    cw_df["MW (Da)"] = cw_df["MW (Da)"].apply(lambda x: f"{x:.5f}")
                    st.dataframe(cw_df, hide_index=True, use_container_width=True)

                st.divider()

                # ── 分Tab展示各离子系列 ──
                z = results["charge"]
                tab_std, tab_na, tab_k, tab_fe, tab_mob, tab_ref = st.tabs([
                    f"Standard (z={z})",
                    "Na-addition",
                    "K-addition",
                    "Fe-addition",
                    "Mobile Phase",
                    "Na/K/Fe Reference",
                ])

                with tab_std:
                    st.caption(f"标准去质子化离子：(MW - {z}×1.007825) / {z}")
                    df = to_df(results["std"], prefix=f"{z}-")
                    st.dataframe(df, hide_index=True, use_container_width=True)

                with tab_na:
                    st.caption("Na加成峰（仅显示 Average MW 和 P=O）")
                    filtered = [(n, mw) for n, mw in results["naA"] if n in ("Average MW", "P=O")]
                    st.dataframe(to_df(filtered, "Na-"), hide_index=True, use_container_width=True)

                with tab_k:
                    st.caption("K加成峰（仅显示 Average MW 和 P=O）")
                    filtered = [(n, mw) for n, mw in results["kA"] if n in ("Average MW", "P=O")]
                    st.dataframe(to_df(filtered, "K-"), hide_index=True, use_container_width=True)

                with tab_fe:
                    st.caption("Fe加成峰（仅显示 Average MW 和 P=O）")
                    filtered = [(n, mw) for n, mw in results["feA"] if n in ("Average MW", "P=O")]
                    st.dataframe(to_df(filtered, "Fe-"), hide_index=True, use_container_width=True)

                with tab_mob:
                    st.caption("流动相相关杂质（仅标准电荷态）")
                    st.dataframe(to_df(results["mob"], f"{z}-"), hide_index=True, use_container_width=True)

                with tab_ref:
                    st.caption("Na/K/Fe 全杂质列表（供详细分析参考，已排除 Average MW 和 P=O）")
                    excl = {"Average MW", "P=O"}
                    for label, data, pfx in [
                        ("Na-addition", results["naA"], "Na-"),
                        ("K-addition",  results["kA"],  "K-"),
                        ("Fe-addition", results["feA"], "Fe-"),
                    ]:
                        st.markdown(f"**{label}**")
                        filtered = [(n, mw) for n, mw in data if n not in excl]
                        st.dataframe(to_df(filtered, pfx), hide_index=True, use_container_width=True)

                # ── 导出 CSV ──
                st.divider()
                st.subheader("导出结果")
                all_export = []
                for name, mw in results["std"]:
                    all_export.append({"Section": f"Standard (z={z})", "Name": f"{z}-{name}", "MW (Da)": f"{mw:.5f}"})
                for name, mw in results["mob"]:
                    all_export.append({"Section": "Mobile Phase", "Name": f"{z}-{name}", "MW (Da)": f"{mw:.5f}"})
                for name, mw in results["naA"]:
                    all_export.append({"Section": "Na-addition", "Name": f"Na-{name}", "MW (Da)": f"{mw:.5f}"})
                for name, mw in results["kA"]:
                    all_export.append({"Section": "K-addition", "Name": f"K-{name}", "MW (Da)": f"{mw:.5f}"})
                for name, mw in results["feA"]:
                    all_export.append({"Section": "Fe-addition", "Name": f"Fe-{name}", "MW (Da)": f"{mw:.5f}"})

                export_df = pd.DataFrame(all_export)
                csv = export_df.to_csv(index=False, encoding="utf-8-sig")
                st.download_button(
                    label="⬇️ 下载完整结果 CSV",
                    data=csv,
                    file_name=f"oligo_results_{sequence[:20].replace('-','_')}.csv",
                    mime="text/csv",
                )


if __name__ == "__main__":
    main()
