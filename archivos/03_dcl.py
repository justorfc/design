import math
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import streamlit as st
from scipy.stats import f


# ============================================================
# CONFIGURACIÓN GENERAL DE LA PÁGINA
# ============================================================
st.title("Diseño en Cuadro Latino (DCL)")
st.markdown("""
Esta página guía al estudiante en el desarrollo **manual paso a paso** de un
**Diseño en Cuadro Latino (DCL)**, mostrando:

- El planteamiento del problema.
- La estructura del cuadro latino.
- Los cálculos manuales.
- La construcción de la tabla ANOVA.
- El cálculo de los valores de **F**.
- El valor de **F crítico**.
- Un script en **R** para validación computacional.
""")


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================
def generar_etiquetas(prefijo: str, k: int) -> List[str]:
    """Genera etiquetas tipo F1, F2, ..., o C1, C2, ..."""
    return [f"{prefijo}{i}" for i in range(1, k + 1)]


def generar_tratamientos(k: int) -> List[str]:
    """
    Genera etiquetas de tratamientos:
    A, B, C, ..., Z, T27, T28, ...
    """
    etiquetas = []
    for i in range(k):
        if i < 26:
            etiquetas.append(chr(65 + i))
        else:
            etiquetas.append(f"T{i + 1}")
    return etiquetas


def construir_latin_square_base(tratamientos: List[str]) -> pd.DataFrame:
    """
    Construye un cuadro latino base por desplazamiento cíclico.
    Cada fila es una rotación de la anterior.
    """
    k = len(tratamientos)
    datos = []
    for i in range(k):
        fila = tratamientos[i:] + tratamientos[:i]
        datos.append(fila)
    return pd.DataFrame(datos, columns=generar_etiquetas("C", k), index=generar_etiquetas("F", k))


def validar_cuadro_latino(df_trat: pd.DataFrame, tratamientos: List[str]) -> Tuple[bool, List[str]]:
    """
    Verifica que cada tratamiento aparezca exactamente una vez por fila y una vez por columna.
    """
    errores = []
    k = len(tratamientos)
    conjunto_esperado = sorted(tratamientos)

    # Validación por filas
    for idx in df_trat.index:
        fila = sorted(df_trat.loc[idx].astype(str).tolist())
        if fila != conjunto_esperado:
            errores.append(
                f"La fila '{idx}' no contiene exactamente una vez cada tratamiento."
            )

    # Validación por columnas
    for col in df_trat.columns:
        columna = sorted(df_trat[col].astype(str).tolist())
        if columna != conjunto_esperado:
            errores.append(
                f"La columna '{col}' no contiene exactamente una vez cada tratamiento."
            )

    # Validación de tamaño
    if df_trat.shape != (k, k):
        errores.append("La matriz de tratamientos no tiene dimensión k × k.")

    return len(errores) == 0, errores


def convertir_respuestas_a_numerico(df_resp: pd.DataFrame) -> pd.DataFrame:
    """Convierte la tabla de respuestas a formato numérico."""
    return df_resp.apply(pd.to_numeric, errors="coerce")


def calcular_totales(
    df_resp: pd.DataFrame,
    df_trat: pd.DataFrame,
    tratamientos: List[str]
) -> Dict[str, object]:
    """
    Calcula:
    - Gran total
    - Totales por fila
    - Totales por columna
    - Totales por tratamiento
    """
    gran_total = float(df_resp.to_numpy().sum())
    totales_fila = df_resp.sum(axis=1)
    totales_columna = df_resp.sum(axis=0)

    totales_trat = {}
    for t in tratamientos:
        mascara = (df_trat == t)
        total_t = df_resp.where(mascara).sum().sum()
        totales_trat[t] = float(total_t)

    return {
        "gran_total": gran_total,
        "totales_fila": totales_fila,
        "totales_columna": totales_columna,
        "totales_trat": totales_trat,
    }


def expresion_suma_cuadrados(valores: List[float]) -> str:
    """Devuelve una expresión tipo a² + b² + c²."""
    return " + ".join([f"{v:.4f}²" if isinstance(v, float) and not float(v).is_integer() else f"{int(v)}²" for v in valores])


def calcular_anova_dcl(
    df_resp: pd.DataFrame,
    df_trat: pd.DataFrame,
    tratamientos: List[str],
    alpha: float = 0.05
) -> Dict[str, object]:
    """
    Realiza todos los cálculos del DCL:
    - TC
    - SCT
    - SCF
    - SCC
    - SCTr
    - SCE
    - gl
    - CM
    - F
    - F crítico
    """
    k = len(tratamientos)
    n = k * k

    y = df_resp.to_numpy(dtype=float)

    gran_total = y.sum()
    suma_cuadrados_observaciones = np.sum(y ** 2)

    totales_fila = df_resp.sum(axis=1).to_numpy(dtype=float)
    totales_columna = df_resp.sum(axis=0).to_numpy(dtype=float)

    totales_trat = []
    for t in tratamientos:
        mascara = (df_trat == t).to_numpy()
        total_t = df_resp.where(df_trat == t).sum().sum()
        totales_trat.append(float(total_t))
    totales_trat = np.array(totales_trat, dtype=float)

    tc = (gran_total ** 2) / n
    sct = suma_cuadrados_observaciones - tc
    scf = np.sum(totales_fila ** 2) / k - tc
    scc = np.sum(totales_columna ** 2) / k - tc
    sctr = np.sum(totales_trat ** 2) / k - tc
    sce = sct - scf - scc - sctr

    gl_filas = k - 1
    gl_columnas = k - 1
    gl_trat = k - 1
    gl_error = (k - 1) * (k - 2)
    gl_total = n - 1

    cm_filas = scf / gl_filas
    cm_columnas = scc / gl_columnas
    cm_trat = sctr / gl_trat
    cm_error = sce / gl_error

    f_filas = cm_filas / cm_error
    f_columnas = cm_columnas / cm_error
    f_trat = cm_trat / cm_error

    f_crit = f.ppf(1 - alpha, gl_trat, gl_error)

    tabla_anova = pd.DataFrame(
        {
            "Fuente de variación": ["Filas", "Columnas", "Tratamientos", "Error", "Total"],
            "SC": [scf, scc, sctr, sce, sct],
            "gl": [gl_filas, gl_columnas, gl_trat, gl_error, gl_total],
            "CM": [cm_filas, cm_columnas, cm_trat, cm_error, np.nan],
            "F calculado": [f_filas, f_columnas, f_trat, np.nan, np.nan],
        }
    )

    return {
        "k": k,
        "n": n,
        "gran_total": gran_total,
        "suma_cuadrados_observaciones": suma_cuadrados_observaciones,
        "totales_fila": totales_fila,
        "totales_columna": totales_columna,
        "totales_trat": totales_trat,
        "tc": tc,
        "sct": sct,
        "scf": scf,
        "scc": scc,
        "sctr": sctr,
        "sce": sce,
        "gl_filas": gl_filas,
        "gl_columnas": gl_columnas,
        "gl_trat": gl_trat,
        "gl_error": gl_error,
        "gl_total": gl_total,
        "cm_filas": cm_filas,
        "cm_columnas": cm_columnas,
        "cm_trat": cm_trat,
        "cm_error": cm_error,
        "f_filas": f_filas,
        "f_columnas": f_columnas,
        "f_trat": f_trat,
        "f_crit": f_crit,
        "tabla_anova": tabla_anova,
    }


def formatear_num(x: float, dec: int = 4) -> str:
    """Formatea números de forma amigable."""
    if pd.isna(x):
        return ""
    if float(x).is_integer():
        return str(int(x))
    return f"{x:.{dec}f}"


def generar_script_r(
    nombre_respuesta: str,
    filas: List[str],
    columnas: List[str],
    tratamientos: List[str],
    df_trat: pd.DataFrame,
    df_resp: pd.DataFrame,
) -> str:
    """
    Genera script en R con easyanova y con aov().
    """
    trat_vector = []
    fila_vector = []
    col_vector = []
    resp_vector = []

    for i, fila in enumerate(df_resp.index):
        for j, col in enumerate(df_resp.columns):
            trat_vector.append(str(df_trat.loc[fila, col]))
            fila_vector.append(str(fila))
            col_vector.append(str(col))
            resp_vector.append(float(df_resp.loc[fila, col]))

    trat_r = ", ".join([f'"{x}"' for x in trat_vector])
    fila_r = ", ".join([f'"{x}"' for x in fila_vector])
    col_r = ", ".join([f'"{x}"' for x in col_vector])
    resp_r = ", ".join([formatear_num(x, 4) for x in resp_vector])

    nombre_respuesta_r = nombre_respuesta.strip().replace(" ", "_").lower()
    if not nombre_respuesta_r:
        nombre_respuesta_r = "respuesta"

    script = f'''# ==========================================================
# VALIDACIÓN COMPUTACIONAL EN R
# Diseño en Cuadro Latino (DCL)
# ==========================================================

# Instale el paquete si es necesario:
# install.packages("easyanova")

library(easyanova)

tratamiento <- as.factor(c({trat_r}))
fila <- as.factor(c({fila_r}))
columna <- as.factor(c({col_r}))
{nombre_respuesta_r} <- c({resp_r})

datos <- data.frame(tratamiento, fila, columna, {nombre_respuesta_r})

cat("\\n--- Datos del experimento ---\\n")
print(datos)

# ----------------------------------------------------------
# Opción 1: easyanova
# ----------------------------------------------------------
cat("\\n--- ANOVA con easyanova ---\\n")
anova_dcl <- eal(datos, design = 3)
print(anova_dcl$`Analysis of variance`)

# ----------------------------------------------------------
# Opción 2: aov() base R
# ----------------------------------------------------------
cat("\\n--- ANOVA con aov() ---\\n")
modelo <- aov({nombre_respuesta_r} ~ fila + columna + tratamiento, data = datos)
summary(modelo)
'''
    return script


def cargar_ejemplo_1() -> Tuple[str, str, str, str, int, List[str], List[str], List[str], pd.DataFrame, pd.DataFrame]:
    """
    Ejemplo 1 del material del usuario: Secador Solar de Convección Natural.
    """
    titulo = "Secador Solar de Convección Natural"
    variable = "Eficiencia térmica (%)"
    filas_desc = "Filas = Horarios de medición"
    columnas_desc = "Columnas = Posición en la cámara"
    k = 4

    filas = ["F1", "F2", "F3", "F4"]
    columnas = ["C1", "C2", "C3", "C4"]
    tratamientos = ["A", "B", "C", "D"]

    df_trat = pd.DataFrame(
        [
            ["A", "B", "C", "D"],
            ["B", "C", "D", "A"],
            ["C", "D", "A", "B"],
            ["D", "A", "B", "C"],
        ],
        index=filas,
        columns=columnas,
    )

    df_resp = pd.DataFrame(
        [
            [25, 37, 42, 33],
            [35, 46, 35, 36],
            [48, 37, 39, 48],
            [36, 41, 47, 55],
        ],
        index=filas,
        columns=columnas,
    )

    return titulo, variable, filas_desc, columnas_desc, k, filas, columnas, tratamientos, df_trat, df_resp


def cargar_ejemplo_2() -> Tuple[str, str, str, str, int, List[str], List[str], List[str], pd.DataFrame, pd.DataFrame]:
    """
    Ejemplo 2 del material del usuario: Secador de Túnel de Convección Forzada.
    """
    titulo = "Secador de Túnel de Convección Forzada"
    variable = "Pérdida de humedad (%)"
    filas_desc = "Filas = Densidad de carga"
    columnas_desc = "Columnas = Lote de cosecha"
    k = 4

    filas = ["D1", "D2", "D3", "D4"]
    columnas = ["L1", "L2", "L3", "L4"]
    tratamientos = ["V1", "V2", "V3", "V4"]

    df_trat = pd.DataFrame(
        [
            ["V1", "V2", "V3", "V4"],
            ["V2", "V3", "V4", "V1"],
            ["V3", "V4", "V1", "V2"],
            ["V4", "V1", "V2", "V3"],
        ],
        index=filas,
        columns=columnas,
    )

    df_resp = pd.DataFrame(
        [
            [44, 51, 62, 63],
            [45, 56, 59, 48],
            [42, 53, 40, 48],
            [49, 36, 47, 52],
        ],
        index=filas,
        columns=columnas,
    )

    return titulo, variable, filas_desc, columnas_desc, k, filas, columnas, tratamientos, df_trat, df_resp


# ============================================================
# PANEL DE CONFIGURACIÓN
# ============================================================
st.subheader("1. Configuración del problema")

modo = st.radio(
    "Seleccione el modo de trabajo:",
    options=[
        "Ejemplo 1: Secador Solar de Convección Natural",
        "Ejemplo 2: Secador de Túnel de Convección Forzada",
        "Problema personalizado",
    ],
    horizontal=False,
)

alpha = st.selectbox(
    "Nivel de significancia (α):",
    options=[0.10, 0.05, 0.01],
    index=1,
)

if modo == "Ejemplo 1: Secador Solar de Convección Natural":
    (
        titulo_problema,
        nombre_respuesta,
        desc_filas,
        desc_columnas,
        k,
        filas,
        columnas,
        tratamientos,
        df_trat_inicial,
        df_resp_inicial,
    ) = cargar_ejemplo_1()

elif modo == "Ejemplo 2: Secador de Túnel de Convección Forzada":
    (
        titulo_problema,
        nombre_respuesta,
        desc_filas,
        desc_columnas,
        k,
        filas,
        columnas,
        tratamientos,
        df_trat_inicial,
        df_resp_inicial,
    ) = cargar_ejemplo_2()

else:
    col1, col2 = st.columns(2)

    with col1:
        titulo_problema = st.text_input(
            "Título o planteamiento breve del problema:",
            value="Problema personalizado en DCL",
        )
        nombre_respuesta = st.text_input(
            "Nombre de la variable respuesta:",
            value="Respuesta",
        )
        desc_filas = st.text_input(
            "¿Qué representan las filas?",
            value="Filas = Bloque 1",
        )

    with col2:
        desc_columnas = st.text_input(
            "¿Qué representan las columnas?",
            value="Columnas = Bloque 2",
        )
        k = st.number_input(
            "Orden del Cuadro Latino (k):",
            min_value=3,
            max_value=8,
            value=4,
            step=1,
        )

    filas = generar_etiquetas("F", k)
    columnas = generar_etiquetas("C", k)
    tratamientos = generar_tratamientos(k)

    df_trat_inicial = construir_latin_square_base(tratamientos)
    df_resp_inicial = pd.DataFrame(
        np.zeros((k, k), dtype=float),
        index=filas,
        columns=columnas,
    )

st.markdown("---")

# ============================================================
# INFORMACIÓN DEL PLANTEAMIENTO
# ============================================================
st.subheader("2. Planteamiento del problema")

st.markdown(f"**Problema:** {titulo_problema}")
st.markdown(f"**Variable respuesta:** {nombre_respuesta}")
st.markdown(f"**Bloqueo por filas:** {desc_filas}")
st.markdown(f"**Bloqueo por columnas:** {desc_columnas}")
st.markdown(f"**Tratamientos:** {', '.join(tratamientos)}")
st.markdown(f"**Orden del diseño:** k = {k}")

st.info(
    "En un DCL de orden k, cada tratamiento debe aparecer exactamente una vez en cada fila "
    "y exactamente una vez en cada columna."
)

# ============================================================
# INGRESO DE ESTRUCTURA Y DATOS
# ============================================================
st.subheader("3. Estructura del Cuadro Latino")

st.markdown("### 3.1. Matriz de tratamientos")
st.caption("Puede editar la disposición de los tratamientos, pero debe conservar la estructura de cuadro latino.")

df_trat = st.data_editor(
    df_trat_inicial,
    use_container_width=True,
    num_rows="fixed",
    key="tratamientos_editor",
)

st.markdown("### 3.2. Matriz de respuestas")
st.caption("Ingrese aquí los valores observados de la variable respuesta.")

df_resp = st.data_editor(
    df_resp_inicial,
    use_container_width=True,
    num_rows="fixed",
    key="respuesta_editor",
)

df_trat = pd.DataFrame(df_trat, index=filas, columns=columnas).astype(str)
df_resp = pd.DataFrame(df_resp, index=filas, columns=columnas)
df_resp = convertir_respuestas_a_numerico(df_resp)

# ============================================================
# VALIDACIONES
# ============================================================
st.subheader("4. Validación de la estructura")

valido, errores = validar_cuadro_latino(df_trat, tratamientos)

if not valido:
    st.error("La estructura ingresada no corresponde a un Cuadro Latino válido.")
    for err in errores:
        st.write(f"- {err}")
else:
    st.success("La matriz de tratamientos cumple la estructura de un Cuadro Latino.")

if df_resp.isna().any().any():
    st.warning("Hay celdas de respuesta vacías o no numéricas. Complete todos los datos para continuar.")
    st.stop()

if not valido:
    st.stop()

# ============================================================
# BOTÓN DE CÁLCULO
# ============================================================
calcular = st.button("Calcular desarrollo manual y tabla ANOVA", type="primary")

if calcular:
    resultados = calcular_anova_dcl(df_resp, df_trat, tratamientos, alpha=alpha)

    # --------------------------------------------------------
    # MOSTRAR TABLA COMBINADA DEL EXPERIMENTO
    # --------------------------------------------------------
    st.markdown("---")
    st.subheader("5. Tabla del experimento")

    tabla_mixta = df_resp.copy().astype(object)
    for i in filas:
        for j in columnas:
            tabla_mixta.loc[i, j] = f"{formatear_num(float(df_resp.loc[i, j]), 4)} ({df_trat.loc[i, j]})"

    st.dataframe(tabla_mixta, use_container_width=True)

    # --------------------------------------------------------
    # CÁLCULO DE TOTALES
    # --------------------------------------------------------
    st.subheader("6. Paso 1: Cálculo de totales")

    totales = calcular_totales(df_resp, df_trat, tratamientos)

    st.markdown(f"### 6.1. Gran total $Y_{{..}}$")
    valores_obs = df_resp.to_numpy().flatten().tolist()
    expr_obs = " + ".join([formatear_num(v, 4) for v in valores_obs])
    st.latex(rf"Y_{{..}} = {expr_obs} = {formatear_num(totales['gran_total'], 4)}")

    st.markdown("### 6.2. Totales por fila $Y_{i.}$")
    for idx, total in zip(df_resp.index, totales["totales_fila"]):
        elementos = df_resp.loc[idx].tolist()
        expr = " + ".join([formatear_num(v, 4) for v in elementos])
        st.latex(rf"{idx}: \quad {expr} = {formatear_num(float(total), 4)}")

    st.markdown("### 6.3. Totales por columna $Y_{.j}$")
    for col, total in zip(df_resp.columns, totales["totales_columna"]):
        elementos = df_resp[col].tolist()
        expr = " + ".join([formatear_num(v, 4) for v in elementos])
        st.latex(rf"{col}: \quad {expr} = {formatear_num(float(total), 4)}")

    st.markdown("### 6.4. Totales por tratamiento $Y_{..t}$")
    for t in tratamientos:
        valores_t = []
        for i in filas:
            for j in columnas:
                if df_trat.loc[i, j] == t:
                    valores_t.append(float(df_resp.loc[i, j]))
        expr = " + ".join([formatear_num(v, 4) for v in valores_t])
        st.latex(rf"{t}: \quad {expr} = {formatear_num(float(totales['totales_trat'][t]), 4)}")

    # --------------------------------------------------------
    # SUMAS DE CUADRADOS
    # --------------------------------------------------------
    st.subheader("7. Paso 2: Cálculo de sumas de cuadrados")

    k_val = resultados["k"]
    n_val = resultados["n"]

    st.markdown("### 7.1. Fórmulas generales del DCL")
    st.latex(r"TC = \frac{(Y_{..})^2}{k^2}")
    st.latex(r"SCT = \sum y_{ij(t)}^2 - TC")
    st.latex(r"SCF = \frac{\sum Y_{i.}^2}{k} - TC")
    st.latex(r"SCC = \frac{\sum Y_{.j}^2}{k} - TC")
    st.latex(r"SCTr = \frac{\sum Y_{..t}^2}{k} - TC")
    st.latex(r"SCE = SCT - SCF - SCC - SCTr")

    st.markdown("### 7.2. Término de corrección (TC)")
    st.latex(
        rf"TC = \frac{{({formatear_num(resultados['gran_total'], 4)})^2}}{{{n_val}}} = {formatear_num(resultados['tc'], 4)}"
    )

    st.markdown("### 7.3. Suma de Cuadrados Total (SCT)")
    expr_sq_obs = expresion_suma_cuadrados(valores_obs)
    st.latex(
        rf"SCT = ({expr_sq_obs}) - {formatear_num(resultados['tc'], 4)} = {formatear_num(resultados['suma_cuadrados_observaciones'], 4)} - {formatear_num(resultados['tc'], 4)} = {formatear_num(resultados['sct'], 4)}"
    )

    st.markdown("### 7.4. Suma de Cuadrados de Filas (SCF)")
    expr_sq_filas = expresion_suma_cuadrados(resultados["totales_fila"].tolist())
    suma_sq_filas = np.sum(resultados["totales_fila"] ** 2)
    st.latex(
        rf"SCF = \frac{{{expr_sq_filas}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = \frac{{{formatear_num(suma_sq_filas, 4)}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = {formatear_num(resultados['scf'], 4)}"
    )

    st.markdown("### 7.5. Suma de Cuadrados de Columnas (SCC)")
    expr_sq_cols = expresion_suma_cuadrados(resultados["totales_columna"].tolist())
    suma_sq_cols = np.sum(resultados["totales_columna"] ** 2)
    st.latex(
        rf"SCC = \frac{{{expr_sq_cols}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = \frac{{{formatear_num(suma_sq_cols, 4)}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = {formatear_num(resultados['scc'], 4)}"
    )

    st.markdown("### 7.6. Suma de Cuadrados de Tratamientos (SCTr)")
    expr_sq_trat = expresion_suma_cuadrados(resultados["totales_trat"].tolist())
    suma_sq_trat = np.sum(resultados["totales_trat"] ** 2)
    st.latex(
        rf"SCTr = \frac{{{expr_sq_trat}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = \frac{{{formatear_num(suma_sq_trat, 4)}}}{{{k_val}}} - {formatear_num(resultados['tc'], 4)} = {formatear_num(resultados['sctr'], 4)}"
    )

    st.markdown("### 7.7. Suma de Cuadrados del Error (SCE)")
    st.latex(
        rf"SCE = {formatear_num(resultados['sct'], 4)} - {formatear_num(resultados['scf'], 4)} - {formatear_num(resultados['scc'], 4)} - {formatear_num(resultados['sctr'], 4)} = {formatear_num(resultados['sce'], 4)}"
    )

    # --------------------------------------------------------
    # TABLA ANOVA
    # --------------------------------------------------------
    st.subheader("8. Paso 3: Construcción de la tabla ANOVA")

    st.markdown("### 8.1. Grados de libertad")
    st.latex(r"gl_{Filas} = k - 1")
    st.latex(r"gl_{Columnas} = k - 1")
    st.latex(r"gl_{Tratamientos} = k - 1")
    st.latex(r"gl_{Error} = (k - 1)(k - 2)")
    st.latex(r"gl_{Total} = k^2 - 1")

    st.latex(rf"gl_{{Filas}} = {k_val} - 1 = {resultados['gl_filas']}")
    st.latex(rf"gl_{{Columnas}} = {k_val} - 1 = {resultados['gl_columnas']}")
    st.latex(rf"gl_{{Tratamientos}} = {k_val} - 1 = {resultados['gl_trat']}")
    st.latex(rf"gl_{{Error}} = ({k_val} - 1)({k_val} - 2) = {resultados['gl_error']}")
    st.latex(rf"gl_{{Total}} = {k_val}^2 - 1 = {resultados['gl_total']}")

    st.markdown("### 8.2. Cuadrados medios")
    st.latex(r"CM_{Filas} = \frac{SCF}{gl_{Filas}}")
    st.latex(r"CM_{Columnas} = \frac{SCC}{gl_{Columnas}}")
    st.latex(r"CM_{Tratamientos} = \frac{SCTr}{gl_{Tratamientos}}")
    st.latex(r"CM_{Error} = \frac{SCE}{gl_{Error}}")

    st.latex(
        rf"CM_{{Filas}} = \frac{{{formatear_num(resultados['scf'], 4)}}}{{{resultados['gl_filas']}}} = {formatear_num(resultados['cm_filas'], 4)}"
    )
    st.latex(
        rf"CM_{{Columnas}} = \frac{{{formatear_num(resultados['scc'], 4)}}}{{{resultados['gl_columnas']}}} = {formatear_num(resultados['cm_columnas'], 4)}"
    )
    st.latex(
        rf"CM_{{Tratamientos}} = \frac{{{formatear_num(resultados['sctr'], 4)}}}{{{resultados['gl_trat']}}} = {formatear_num(resultados['cm_trat'], 4)}"
    )
    st.latex(
        rf"CM_{{Error}} = \frac{{{formatear_num(resultados['sce'], 4)}}}{{{resultados['gl_error']}}} = {formatear_num(resultados['cm_error'], 4)}"
    )

    st.markdown("### 8.3. Estadísticos F calculados")
    st.latex(r"F_{Filas} = \frac{CM_{Filas}}{CM_{Error}}")
    st.latex(r"F_{Columnas} = \frac{CM_{Columnas}}{CM_{Error}}")
    st.latex(r"F_{Tratamientos} = \frac{CM_{Tratamientos}}{CM_{Error}}")

    st.latex(
        rf"F_{{Filas}} = \frac{{{formatear_num(resultados['cm_filas'], 4)}}}{{{formatear_num(resultados['cm_error'], 4)}}} = {formatear_num(resultados['f_filas'], 4)}"
    )
    st.latex(
        rf"F_{{Columnas}} = \frac{{{formatear_num(resultados['cm_columnas'], 4)}}}{{{formatear_num(resultados['cm_error'], 4)}}} = {formatear_num(resultados['f_columnas'], 4)}"
    )
    st.latex(
        rf"F_{{Tratamientos}} = \frac{{{formatear_num(resultados['cm_trat'], 4)}}}{{{formatear_num(resultados['cm_error'], 4)}}} = {formatear_num(resultados['f_trat'], 4)}"
    )

    st.markdown("### 8.4. Valor crítico de F")
    st.latex(
        rf"F_{{crítico}} = F_{{1-\alpha}}({resultados['gl_trat']}, {resultados['gl_error']}) = {formatear_num(resultados['f_crit'], 4)}"
    )
    st.markdown(
        f"Para **α = {alpha}**, con **gl1 = {resultados['gl_trat']}** y **gl2 = {resultados['gl_error']}**, "
        f"el valor crítico es **{formatear_num(resultados['f_crit'], 4)}**."
    )

    tabla_anova = resultados["tabla_anova"].copy()
    for col in ["SC", "CM", "F calculado"]:
        tabla_anova[col] = tabla_anova[col].apply(lambda x: np.nan if pd.isna(x) else round(float(x), 4))

    st.markdown("### 8.5. Tabla ANOVA final")
    st.dataframe(tabla_anova, use_container_width=True)

    # --------------------------------------------------------
    # DECISIONES
    # --------------------------------------------------------
    st.subheader("9. Interpretación de decisiones con F calculado y F crítico")

    def decision(f_calc: float, f_crit: float, nombre: str) -> str:
        if f_calc > f_crit:
            return (
                f"Como **F calculado de {nombre} = {formatear_num(f_calc, 4)}** es mayor que "
                f"**F crítico = {formatear_num(f_crit, 4)}**, se **rechaza H₀** para {nombre}."
            )
        return (
            f"Como **F calculado de {nombre} = {formatear_num(f_calc, 4)}** no supera "
            f"**F crítico = {formatear_num(f_crit, 4)}**, no se rechaza H₀ para {nombre}."
        )

    st.markdown(f"- {decision(resultados['f_filas'], resultados['f_crit'], 'Filas')}")
    st.markdown(f"- {decision(resultados['f_columnas'], resultados['f_crit'], 'Columnas')}")
    st.markdown(f"- {decision(resultados['f_trat'], resultados['f_crit'], 'Tratamientos')}")

    # Mejor tratamiento según promedio observado
    promedios_trat = {}
    for t in tratamientos:
        vals = []
        for i in filas:
            for j in columnas:
                if df_trat.loc[i, j] == t:
                    vals.append(float(df_resp.loc[i, j]))
        promedios_trat[t] = np.mean(vals)

    mejor_trat = max(promedios_trat, key=promedios_trat.get)
    mejor_prom = promedios_trat[mejor_trat]

    st.markdown("### 9.1. Promedio observado por tratamiento")
    df_prom = pd.DataFrame(
        {
            "Tratamiento": list(promedios_trat.keys()),
            "Promedio observado": [round(v, 4) for v in promedios_trat.values()],
        }
    )
    st.dataframe(df_prom, use_container_width=True)

    st.markdown(
        f"El tratamiento con mayor promedio observado es **{mejor_trat}**, "
        f"con un valor promedio de **{formatear_num(mejor_prom, 4)}**."
    )

    # --------------------------------------------------------
    # SCRIPT EN R
    # --------------------------------------------------------
    st.subheader("10. Script de validación en R")

    script_r = generar_script_r(
        nombre_respuesta=nombre_respuesta,
        filas=filas,
        columnas=columnas,
        tratamientos=tratamientos,
        df_trat=df_trat,
        df_resp=df_resp,
    )

    st.code(script_r, language="r")

    # --------------------------------------------------------
    # DESCARGAS
    # --------------------------------------------------------
    st.subheader("11. Exportación")

    csv_anova = tabla_anova.to_csv(index=False).encode("utf-8")
    csv_datos = df_resp.to_csv().encode("utf-8")

    st.download_button(
        label="Descargar tabla ANOVA en CSV",
        data=csv_anova,
        file_name="tabla_anova_dcl.csv",
        mime="text/csv",
    )

    st.download_button(
        label="Descargar matriz de respuestas en CSV",
        data=csv_datos,
        file_name="datos_dcl_respuestas.csv",
        mime="text/csv",
    )

    st.download_button(
        label="Descargar script R",
        data=script_r.encode("utf-8"),
        file_name="validacion_dcl.R",
        mime="text/plain",
    )