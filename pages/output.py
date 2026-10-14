import streamlit as st
import py3Dmol
import navigation
from rdkit import Chem
from streamlit.components.v1 import html

from utils.polymer_draw import build_polymer
from rdkit.Chem import Descriptors
from utils.tg_predict import predict_tg

import pickle
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

token2idx_path = BASE_DIR / "model" / "token2idx.pkl"


def show_output():

    repeat_count = None
    mol = None
    force_field = None

    col10, col11 = st.columns([1, 1])

    with col10:
        st.subheader("Prediction Results")

    with col11:
        st.markdown("""
            <style>
            /* Container */
            div[data-testid="stVerticalBlockBorderWrapper"] {
                min-height: 100px;
            }

            /* Metric value */
            [data-testid="stMetricValue"] {
                font-size: 20px;
            }

            /* Metric label */
            [data-testid="stMetricLabel"] {
                font-size: 8px;
            }
            </style>
            """, unsafe_allow_html=True)

    col1, col2 = st.columns([4, 1])

    with col1:
        with st.container(height=400, border=True):
            with st.container(height=250, border=True):
                viewer_placeholder = st.empty()

            with st.form("polymer_form"):
                col3, col4 = st.columns([3, 1])

                with col3:
                    repeat_count = st.number_input(
                        "Number of Repeating Units",
                        min_value=1,
                        max_value=20,
                        value=5,
                        step=1,
                    )

                with col4:
                    st.markdown(
                        "<div style='height: 28px'></div>",
                        unsafe_allow_html=True,
                    )
                    submitted = st.form_submit_button(
                        "Generate 3D Structure",
                        use_container_width=True,
                    )

            if submitted or st.session_state.get("active_tab") == "output":
                try:
                    with st.spinner(
                        "Generating and optimizing the 3D structure..."
                    ):
                        mol, force_field = build_polymer(
                            psmiles=st.session_state.get("psmiles", ""),
                            repeat_count=int(repeat_count),
                        )

                    st.toast("3D structure generated successfully.", icon="✅")

                    # Convert the RDKit molecule to MOL format
                    mol_block = Chem.MolToMolBlock(mol)

                    viewer = py3Dmol.view(width=900, height=210)
                    viewer.addModel(mol_block, "mol")

                    # Ball-and-stick representation with element-based colors
                    viewer.setStyle({
                        "stick": {
                            "radius": 0.16,
                            "colorscheme": "Jmol",
                        },
                        "sphere": {
                            "scale": 0.28,
                            "colorscheme": "Jmol",
                        },
                    })

                    viewer.setBackgroundColor("#0E1117")
                    viewer.zoomTo()

                    # Display the 3D structure in Streamlit
                    with viewer_placeholder:
                        html(
                            viewer._make_html(),
                            height=210,
                            scrolling=False,
                        )

                except (ValueError, RuntimeError) as exc:
                    st.toast(str(exc), icon="❌")

                except Exception as exc:
                    st.toast(
                        f"An error occurred during visualization: {exc}",
                        icon="❌"
                    )

    with col2:
        psmiles = st.session_state.get("psmiles", "")

        with open(token2idx_path, "rb") as f:
            token2idx = pickle.load(f)

        with st.container(height=88, border=True):
            tg_pred = predict_tg(
                psmiles,
                token2idx,
                maxLen=120
            )

            st.metric(
                "Tg",
                f"{tg_pred:.2f} K" if repeat_count is not None else "-"
            )

        with st.container(height=88, border=True):
            st.metric(
                "Number of Atoms",
                mol.GetNumAtoms() if mol is not None else "-"
            )

        with st.container(height=88, border=True):
            mw = Descriptors.MolWt(mol) if mol is not None else None

            st.metric(
                "Molecular Weight",
                f"{mw:.2f} g/mol" if mw is not None else "-"
            )

        with st.container(height=88, border=True):
            logp = Descriptors.MolLogP(mol) if mol is not None else None

            st.metric(
                "Hydrophobicity (LogP)",
                f"{logp:.2f}" if logp is not None else "-"
            )

    col6, col7, col8, col9 = st.columns([1, 2, 2, 1])

    with col6:
        st.button(
            "Back",
            on_click=navigation.buka_predict,
            use_container_width=True,
            key="about_kembali"
        )

    with col9:
        st.button(
            "Finish",
            on_click=navigation.buka_landing,
            use_container_width=True,
            key="about_mulai"
        )
