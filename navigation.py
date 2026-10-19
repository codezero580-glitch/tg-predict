
import streamlit as st
from utils.polymer_draw import build_polymer


def buka_landing():
    st.session_state.active_tab = "landing"


def buka_about():
    st.session_state.active_tab = "about"


def buka_predict():
    st.session_state.active_tab = "predict"


def buka_output():
    st.session_state.active_tab = "output"


def prediksi_output():
    psmiles = st.session_state.get(
        "generated_psmiles", ""
    )

    repeat_count = int(
        st.session_state.get("repeat_count", 5)
    )

    try:
        mol, force_field = build_polymer(
            psmiles=psmiles,
            repeat_count=repeat_count,
        )

        st.session_state["psmiles"] = psmiles
        st.session_state["polymer_mol"] = mol
        st.session_state["force_field"] = force_field
        st.session_state["polymer_error"] = None

    except Exception as exc:
        st.session_state["polymer_mol"] = None
        st.session_state["force_field"] = None
        st.session_state["polymer_error"] = str(exc)

    # Pindah halaman setelah proses selesai
    st.session_state.active_tab = "output"


def finish_process():
    st.session_state.clear()

    st.session_state.active_tab = "landing"
