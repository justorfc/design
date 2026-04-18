import streamlit as st

st.title("📐 DPD — Diseño en Parcelas Divididas")
st.write(
    "Página base para desarrollar el módulo del DPD. Aquí podrá organizar factores de parcela principal, "
    "subparcela, estructura de error, análisis manual y validación en R."
)

st.subheader("Estructura sugerida del futuro módulo")
st.markdown(
    """
    - Contexto experimental.
    - Factor de parcela principal.
    - Factor de subparcela.
    - Bloques o repeticiones.
    - Tabla de observaciones.
    - Modelo lineal.
    - ANOVA con errores diferenciados.
    - Interpretación de efectos principales e interacción.
    - Script de validación en R.
    """
)
