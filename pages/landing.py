import streamlit as st
import navigation


def show_landing():
    st.write("""
        <h1 style="text-align: center;
        font-size: 48px; color: deepskyblue;">
        Polymer Tg Predictor</h1>
    """, unsafe_allow_html=True)

    st.write("""
        <p style="text-align: center;">
            An application for predicting the glass transition
            temperature (Tg) of polymers based on their chemical structures.
        </p>
    """, unsafe_allow_html=True)

    col1, col2, col3, col4 = st.columns([4, 3, 2, 4])

    with col2:
        st.button(
            "About the Application",
            on_click=navigation.buka_about,
            use_container_width=True,
            key="landing_tentang"
        )

    with col3:
        st.button(
            "Start Prediction",
            on_click=navigation.buka_predict,
            use_container_width=True,
            key="landing_mulai"
        )
