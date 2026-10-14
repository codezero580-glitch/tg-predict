import streamlit as st
import navigation


def show_about():
    col1, col2, col3, col4 = st.columns([1, 3, 3, 1])

    with col1:
        st.button(
            "Back",
            on_click=navigation.buka_landing,
            use_container_width=True,
            key="about_kembali"
        )

    with col4:
        st.button(
            "Start Prediction",
            on_click=navigation.buka_predict,
            use_container_width=True,
            key="about_mulai"
        )

    st.markdown("# About the Application")

    st.write("""
        <p style="text-align: center;">
            Polymer Tg Predictor is a machine learning-based application
            designed to predict the glass transition temperature (Tg) of
            polymers from their chemical structures. It explores the
            integration of chemistry, data science, and computational
            technology for polymer property prediction.
        </p>
    """, unsafe_allow_html=True)

    st.write("""
        <p style="text-align: center;">
            This application is still under development and has limitations
            related to data availability, model performance, and prediction
            accuracy. The predicted values should therefore be interpreted
            as estimates rather than substitutes for experimental measurements.
        </p>
    """, unsafe_allow_html=True)

    st.write("""
        <p style="text-align: center;">
            The developer is continuously working to improve the application,
            enhance its predictive capabilities, and address its existing
            limitations as development progresses.
        </p>
    """, unsafe_allow_html=True)

    st.write("""
        <p style="text-align: center;">
            Thank you for visiting and exploring this application.
        </p>
    """, unsafe_allow_html=True)

    st.write("""
        <footer style="text-align: center;">
            — CodeZero —
        </footer>
    """, unsafe_allow_html=True)
