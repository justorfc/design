import streamlit as st

st.title("🧱 DBCA — Diseño en Bloques Completos al Azar")
st.write(
    "Página base para desarrollar el módulo del DBCA. Aquí podrá incorporar el análisis con tratamientos y bloques, "
    "explicación del bloqueo, procedimiento manual y validación computacional."
)

st.subheader("Estructura sugerida del futuro módulo")
st.markdown(
    """
    - Planteamiento del problema.
    - Tratamientos y bloques.
    - Matriz de observaciones.
    - Totales por tratamiento y por bloque.
    - Sumas de cuadrados.
    - Tabla ANOVA.
    - F calculado y F crítico.
    - Conclusión técnica.
    - Script de validación en R.
    """
)
