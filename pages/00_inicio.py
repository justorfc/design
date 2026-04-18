import streamlit as st

st.title("🧪 App multipágina de Diseños Experimentales")
st.markdown(
    """
    Esta aplicación organiza en páginas separadas los principales diseños experimentales:

    - **DCA**: Diseño Completamente Aleatorizado.
    - **DBCA**: Diseño en Bloques Completos al Azar.
    - **DCL**: Diseño en Cuadro Latino.
    - **DPD**: Diseño en Parcelas Divididas.
    - **Diseños Factoriales**.
    - **Metodología de Superficie de Respuesta**.

    ### Propósito pedagógico
    La idea es que cada página pueda evolucionar hacia una herramienta de apoyo docente que:

    1. Reciba el planteamiento del problema.
    2. Organice los datos experimentales.
    3. Muestre el desarrollo manual paso a paso.
    4. Construya la tabla ANOVA.
    5. Muestre la validación en **R** y/o **Python**.

    En esta versión ya se incluye una página funcional base para **DCL** y se deja la estructura lista
    para ampliar las demás páginas de forma ordenada.
    """
)

st.info(
    "Use el menú lateral para entrar a cada diseño. "
    "La arquitectura quedó lista para crecer de forma modular."
)
