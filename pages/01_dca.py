import streamlit as st

st.title("🎲 DCA — Diseño Completamente Aleatorizado")
st.write(
    "Página base para desarrollar el módulo del DCA. Aquí podrá incorporar más adelante: "
    "formulario del problema, ingreso de tratamientos, réplicas, cálculos manuales, ANOVA y validación en R."
)

st.subheader("Estructura sugerida del futuro módulo")
st.markdown(
    """
    - Contexto experimental.
    - Número de tratamientos.
    - Número de repeticiones.
    - Tabla de datos.
    - Cálculo de totales y medias.
    - Sumas de cuadrados.
    - Tabla ANOVA.
    - Prueba F y comparación con F crítico.
    - Script de validación en R.
    """
)
