
import streamlit as st

from pages import landing, about, predict, output


st.set_page_config(
    page_title="Polymer Tg Predictor",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

if "active_tab" not in st.session_state:
    st.session_state.active_tab = "landing"

if st.session_state.active_tab == "landing":
    landing.show_landing()

elif st.session_state.active_tab == "about":
    about.show_about()

elif st.session_state.active_tab == "predict":
    predict.show_predict()

elif st.session_state.active_tab == "output":
    output.show_output()
