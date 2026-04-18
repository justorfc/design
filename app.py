import streamlit as st

st.set_page_config(
    page_title="Diseños Experimentales",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

inicio = st.Page("pages/00_inicio.py", title="Inicio", icon="🏠", default=True)
dca = st.Page("pages/01_dca.py", title="DCA", icon="🎲")
dbca = st.Page("pages/02_dbca.py", title="DBCA", icon="🧱")
dcl = st.Page("pages/03_dcl.py", title="DCL", icon="🔲")
dpd = st.Page("pages/04_dpd.py", title="DPD", icon="📐")
factoriales = st.Page("pages/05_factoriales.py", title="Factoriales", icon="🧮")
superficie = st.Page("pages/06_superficie_respuesta.py", title="Superficie de Respuesta", icon="📈")

pg = st.navigation(
    {
        "General": [inicio],
        "Diseños clásicos": [dca, dbca, dcl, dpd],
        "Diseños avanzados": [factoriales, superficie],
    },
    position="sidebar",
    expanded=True,
)

with st.sidebar:
    st.markdown("---")
    st.caption("App docente para el estudio paso a paso de diseños experimentales.")
    st.caption("Universidad de Sucre · Ingeniería · Estadística Aplicada")

pg.run()
