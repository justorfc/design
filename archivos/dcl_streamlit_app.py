import math
from typing import List, Tuple

import pandas as pd
import streamlit as st
from scipy.stats import f as f_dist


st.set_page_config(page_title="DCL paso a paso", layout="wide")


def parse_labels(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def cyclic_latin_square(labels: List[str]) -> List[List[str]]:
    k = len(labels)
    return [[labels[(i + j) % k] for j in range(k)] for i in range(k)]


def validate_latin_square(square: List[List[str]], treatments: List[str]) -> Tuple[bool, List[str]]:
    errors = []
    k = len(treatments)
    expected = set(treatments)

    if len(square) != k or any(len(row) != k for row in square):
        errors.append("La matriz de tratamientos no tiene dimensiones k × k.")
        return False, errors

    for i, row in enumerate(square, start=1):
        if set(row) != expected:
            errors.append(f"La fila {i} no contiene exactamente una vez cada tratamiento.")

    for j in range(k):
        col = [square[i][j] for i in range(k)]
        if set(col) != expected:
            errors.append(f"La columna {j + 1} no contiene exactamente una vez cada tratamiento.")

    return len(errors) == 0, errors


def flatten_by_rows(matrix: List[List[float]]) -> List[float]:
    return [value for row in matrix for value in row]


def format_sum(values: List[float], decimals: int = 4) -> str:
    return " + ".join(f"{v:.{decimals}f}" if isinstance(v, float) and not float(v).is_integer() else f"{int(v)}" if float(v).is_integer() else str(v) for v in values)


def format_square_sum(values: List[float], decimals: int = 4) -> str:
    parts = []
    for v in values:
        if float(v).is_integer():
            parts.append(f"{int(v)}²")
        else:
            parts.append(f"({v:.{decimals}f})²")
    return " + ".join(parts)


def compute_dcl(values: List[List[float]], treatment_square: List[List[str]], row_labels: List[str], col_labels: List[str], treatment_labels: List[str], alpha: float):
    k = len(values)
    flat_values = flatten_by_rows(values)
    grand_total = sum(flat_values)
    n = k * k
    tc = (grand_total ** 2) / n
    sct = sum(v ** 2 for v in flat_values) - tc

    row_totals = [sum(row) for row in values]
    col_totals = [sum(values[i][j] for i in range(k)) for j in range(k)]

    treatment_totals = []
    for tr in treatment_labels:
        total_tr = sum(values[i][j] for i in range(k) for j in range(k) if treatment_square[i][j] == tr)
        treatment_totals.append(total_tr)

    scf = sum(rt ** 2 for rt in row_totals) / k - tc
    scc = sum(ct ** 2 for ct in col_totals) / k - tc
    sctr = sum(tt ** 2 for tt in treatment_totals) / k - tc
    sce = sct - scf - scc - sctr

    gl_filas = k - 1
    gl_columnas = k - 1
    gl_trat = k - 1
    gl_error = (k - 1) * (k - 2)
    gl_total = k * k - 1

    cm_filas = scf / gl_filas
    cm_columnas = scc / gl_columnas
    cm_trat = sctr / gl_trat
    cm_error = sce / gl_error

    f_filas = cm_filas / cm_error
    f_columnas = cm_columnas / cm_error
    f_trat = cm_trat / cm_error

    f_crit = f_dist.ppf(1 - alpha, gl_trat, gl_error)

    return {
        "k": k,
        "n": n,
        "grand_total": grand_total,
        "tc": tc,
        "sct": sct,
        "row_totals": row_totals,
        "col_totals": col_totals,
        "treatment_totals": treatment_totals,
        "scf": scf,
        "scc": scc,
        "sctr": sctr,
        "sce": sce,
        "gl": {
            "Filas": gl_filas,
            "Columnas": gl_columnas,
            "Tratamientos": gl_trat,
            "Error": gl_error,
            "Total": gl_total,
        },
        "cm": {
            "Filas": cm_filas,
            "Columnas": cm_columnas,
            "Tratamientos": cm_trat,
            "Error": cm_error,
        },
        "f": {
            "Filas": f_filas,
            "Columnas": f_columnas,
            "Tratamientos": f_trat,
        },
        "f_crit": f_crit,
    }


def build_r_script(problem_name: str, response_name: str, row_name: str, col_name: str, treatment_name: str,
                   row_labels: List[str], col_labels: List[str], treatment_labels: List[str],
                   values: List[List[float]], treatment_square: List[List[str]], alpha: float) -> str:
    treatment_vector = [treatment_square[i][j] for i in range(len(row_labels)) for j in range(len(col_labels))]
    row_vector = [row_labels[i] for i in range(len(row_labels)) for _ in range(len(col_labels))]
    col_vector = [col_labels[j] for _ in range(len(row_labels)) for j in range(len(col_labels))]
    y_vector = [values[i][j] for i in range(len(row_labels)) for j in range(len(col_labels))]

    tr_str = ", ".join(f'"{x}"' for x in treatment_vector)
    row_str = ", ".join(f'"{x}"' for x in row_vector)
    col_str = ", ".join(f'"{x}"' for x in col_vector)
    y_str = ", ".join(str(int(x)) if float(x).is_integer() else str(x) for x in y_vector)

    return f'''# ==========================================
# Validación en R del DCL: {problem_name}
# ==========================================

# install.packages("easyanova")
library(easyanova)

tratamiento <- as.factor(c({tr_str}))
fila <- as.factor(c({row_str}))
columna <- as.factor(c({col_str}))
respuesta <- c({y_str})

datos <- data.frame(tratamiento, fila, columna, respuesta)
colnames(datos) <- c("{treatment_name}", "{row_name}", "{col_name}", "{response_name}")

cat("\\n--- Datos del problema: {problem_name} ---\\n")
print(datos)

cat("\\n--- ANOVA con easyanova (design = 3, Cuadro Latino) ---\\n")
resultado_easyanova <- eal(datos, design = 3)
print(resultado_easyanova$`Analysis of variance`)

cat("\\n--- Validación adicional con aov() ---\\n")
modelo <- aov({response_name} ~ {treatment_name} + {row_name} + {col_name}, data = datos)
print(summary(modelo))

cat("\\n--- F crítico al nivel alpha = {alpha} ---\\n")
gl1 <- length(levels(datos${treatment_name})) - 1
gl2 <- (length(levels(datos${treatment_name})) - 1) * (length(levels(datos${treatment_name})) - 2)
F_critico <- qf(1 - {alpha}, df1 = gl1, df2 = gl2)
print(F_critico)
'''


def default_example_1():
    row_labels = ["F1", "F2", "F3", "F4"]
    col_labels = ["C1", "C2", "C3", "C4"]
    treatment_labels = ["A", "B", "C", "D"]
    treatment_square = [
        ["A", "B", "C", "D"],
        ["B", "C", "D", "A"],
        ["C", "D", "A", "B"],
        ["D", "A", "B", "C"],
    ]
    values = [
        [25, 37, 42, 33],
        [35, 46, 35, 36],
        [48, 37, 39, 48],
        [36, 41, 47, 55],
    ]
    return {
        "problem_name": "Secador Solar de Convección Natural",
        "response_name": "eficiencia",
        "row_name": "horario",
        "col_name": "posicion",
        "treatment_name": "material",
        "row_labels": row_labels,
        "col_labels": col_labels,
        "treatment_labels": treatment_labels,
        "treatment_square": treatment_square,
        "values": values,
    }


def default_example_2():
    row_labels = ["D1", "D2", "D3", "D4"]
    col_labels = ["L1", "L2", "L3", "L4"]
    treatment_labels = ["V1", "V2", "V3", "V4"]
    treatment_square = [
        ["V1", "V2", "V3", "V4"],
        ["V2", "V3", "V4", "V1"],
        ["V3", "V4", "V1", "V2"],
        ["V4", "V1", "V2", "V3"],
    ]
    values = [
        [44, 51, 62, 63],
        [45, 56, 59, 48],
        [42, 53, 40, 48],
        [49, 36, 47, 52],
    ]
    return {
        "problem_name": "Secador de Túnel de Convección Forzada",
        "response_name": "humedad",
        "row_name": "densidad",
        "col_name": "lote",
        "treatment_name": "velocidad",
        "row_labels": row_labels,
        "col_labels": col_labels,
        "treatment_labels": treatment_labels,
        "treatment_square": treatment_square,
        "values": values,
    }


def load_example(example_name: str):
    data = default_example_1() if example_name == "Ejercicio 1" else default_example_2()
    st.session_state.problem_name = data["problem_name"]
    st.session_state.response_name = data["response_name"]
    st.session_state.row_name = data["row_name"]
    st.session_state.col_name = data["col_name"]
    st.session_state.treatment_name = data["treatment_name"]
    st.session_state.row_labels_text = ", ".join(data["row_labels"])
    st.session_state.col_labels_text = ", ".join(data["col_labels"])
    st.session_state.treatment_labels_text = ", ".join(data["treatment_labels"])
    st.session_state.alpha = 0.05
    st.session_state.values_df = pd.DataFrame(data["values"], index=data["row_labels"], columns=data["col_labels"])
    st.session_state.treatment_df = pd.DataFrame(data["treatment_square"], index=data["row_labels"], columns=data["col_labels"])


if "initialized" not in st.session_state:
    load_example("Ejercicio 1")
    st.session_state.initialized = True

st.title("Diseño en Cuadro Latino (DCL): solución manual guiada y validación en R")
st.markdown(
    "Esta aplicación guía al estudiante desde el planteamiento del problema hasta la tabla ANOVA, "
    "mostrando fórmulas, operaciones aritméticas y la validación computacional en R."
)

with st.sidebar:
    st.header("Configuración")
    ejemplo = st.selectbox("Cargar ejemplo base", ["Ejercicio 1", "Ejercicio 2", "Personalizado"])
    if st.button("Cargar selección"):
        if ejemplo in ["Ejercicio 1", "Ejercicio 2"]:
            load_example(ejemplo)
    alpha = st.number_input("Nivel de significancia α", min_value=0.001, max_value=0.2, value=float(st.session_state.get("alpha", 0.05)), step=0.001, format="%.3f")
    st.session_state.alpha = alpha
    st.info("Para un DCL válido, cada tratamiento debe aparecer exactamente una vez en cada fila y en cada columna.")

with st.form("problem_form"):
    st.subheader("1) Planteamiento del problema")
    col1, col2 = st.columns(2)
    with col1:
        problem_name = st.text_input("Nombre del problema", value=st.session_state.get("problem_name", ""))
        response_name = st.text_input("Variable respuesta", value=st.session_state.get("response_name", "respuesta"))
        row_name = st.text_input("Nombre de las filas (bloqueo 1)", value=st.session_state.get("row_name", "filas"))
        row_labels_text = st.text_input("Etiquetas de filas separadas por coma", value=st.session_state.get("row_labels_text", "F1, F2, F3, F4"))
    with col2:
        treatment_name = st.text_input("Nombre de los tratamientos", value=st.session_state.get("treatment_name", "tratamientos"))
        col_name = st.text_input("Nombre de las columnas (bloqueo 2)", value=st.session_state.get("col_name", "columnas"))
        col_labels_text = st.text_input("Etiquetas de columnas separadas por coma", value=st.session_state.get("col_labels_text", "C1, C2, C3, C4"))
        treatment_labels_text = st.text_input("Etiquetas de tratamientos separadas por coma", value=st.session_state.get("treatment_labels_text", "A, B, C, D"))

    submitted = st.form_submit_button("Actualizar estructura")

if submitted:
    st.session_state.problem_name = problem_name
    st.session_state.response_name = response_name
    st.session_state.row_name = row_name
    st.session_state.col_name = col_name
    st.session_state.treatment_name = treatment_name
    st.session_state.row_labels_text = row_labels_text
    st.session_state.col_labels_text = col_labels_text
    st.session_state.treatment_labels_text = treatment_labels_text

row_labels = parse_labels(st.session_state.get("row_labels_text", "F1, F2, F3, F4"))
col_labels = parse_labels(st.session_state.get("col_labels_text", "C1, C2, C3, C4"))
treatment_labels = parse_labels(st.session_state.get("treatment_labels_text", "A, B, C, D"))

sizes = {len(row_labels), len(col_labels), len(treatment_labels)}
if len(sizes) != 1 or list(sizes)[0] < 2:
    st.error("Las filas, columnas y tratamientos deben tener el mismo número de niveles (k) y k debe ser al menos 2.")
    st.stop()

k = len(row_labels)

if "values_df" not in st.session_state or list(st.session_state.values_df.index) != row_labels or list(st.session_state.values_df.columns) != col_labels:
    st.session_state.values_df = pd.DataFrame([[0.0] * k for _ in range(k)], index=row_labels, columns=col_labels)

if "treatment_df" not in st.session_state or list(st.session_state.treatment_df.index) != row_labels or list(st.session_state.treatment_df.columns) != col_labels:
    st.session_state.treatment_df = pd.DataFrame(cyclic_latin_square(treatment_labels), index=row_labels, columns=col_labels)

st.subheader("2) Ingreso del cuadro latino")
col_a, col_b = st.columns(2)
with col_a:
    st.markdown("**Matriz de respuestas numéricas**")
    values_df = st.data_editor(
        st.session_state.values_df,
        num_rows="fixed",
        use_container_width=True,
        key="values_editor"
    )
with col_b:
    st.markdown("**Matriz de tratamientos**")
    treatment_df = st.data_editor(
        st.session_state.treatment_df,
        num_rows="fixed",
        use_container_width=True,
        key="treatment_editor"
    )

st.session_state.values_df = values_df
st.session_state.treatment_df = treatment_df

values = [[float(values_df.loc[r, c]) for c in col_labels] for r in row_labels]
treatment_square = [[str(treatment_df.loc[r, c]).strip() for c in col_labels] for r in row_labels]

valid, validation_errors = validate_latin_square(treatment_square, treatment_labels)
for i in range(k):
    for j in range(k):
        if treatment_square[i][j] not in treatment_labels:
            validation_errors.append(
                f"La celda ({row_labels[i]}, {col_labels[j]}) contiene '{treatment_square[i][j]}', que no está en la lista de tratamientos."
            )
            valid = False

if not valid:
    st.error("La disposición de tratamientos no corresponde a un Cuadro Latino válido.")
    for err in validation_errors:
        st.write(f"- {err}")
    st.stop()

results = compute_dcl(
    values=values,
    treatment_square=treatment_square,
    row_labels=row_labels,
    col_labels=col_labels,
    treatment_labels=treatment_labels,
    alpha=st.session_state.alpha,
)

st.subheader("3) Fórmulas generales del DCL")
st.latex(r"TC = \frac{(Y_{..})^2}{k^2}")
st.latex(r"SCT = \sum y_{ij(t)}^2 - TC")
st.latex(r"SCF = \frac{\sum Y_{i.}^2}{k} - TC")
st.latex(r"SCC = \frac{\sum Y_{.j}^2}{k} - TC")
st.latex(r"SCTr = \frac{\sum Y_{..t}^2}{k} - TC")
st.latex(r"SCE = SCT - SCF - SCC - SCTr")
st.latex(r"gl_{filas}=k-1,\; gl_{columnas}=k-1,\; gl_{trat}=k-1,\; gl_{error}=(k-1)(k-2)")
st.latex(r"CM = \frac{SC}{gl},\qquad F = \frac{CM_{fuente}}{CM_{error}}")

st.subheader("4) Paso a paso de la solución manual")
flat_values = flatten_by_rows(values)
row_totals = results["row_totals"]
col_totals = results["col_totals"]
tr_totals = results["treatment_totals"]

with st.expander("Paso 1. Gran total y totales marginales", expanded=True):
    st.markdown(f"**Gran total $Y_{{..}}$** = {format_sum(flat_values)} = **{results['grand_total']:.4f}**")

    st.markdown(f"**Totales por {st.session_state.row_name}**")
    for label, row in zip(row_labels, values):
        st.write(f"{label}: {format_sum(row)} = {sum(row):.4f}")

    st.markdown(f"**Totales por {st.session_state.col_name}**")
    for j, label in enumerate(col_labels):
        col_vals = [values[i][j] for i in range(k)]
        st.write(f"{label}: {format_sum(col_vals)} = {sum(col_vals):.4f}")

    st.markdown(f"**Totales por {st.session_state.treatment_name}**")
    for tr, tot in zip(treatment_labels, tr_totals):
        vals = [values[i][j] for i in range(k) for j in range(k) if treatment_square[i][j] == tr]
        st.write(f"{tr}: {format_sum(vals)} = {tot:.4f}")

with st.expander("Paso 2. Término de corrección", expanded=True):
    st.latex(rf"TC = \frac{{({results['grand_total']:.4f})^2}}{{{results['n']}}} = {results['tc']:.4f}")

with st.expander("Paso 3. Suma de cuadrados total", expanded=True):
    sum_sq = sum(v ** 2 for v in flat_values)
    st.markdown(f"$\\sum y^2$ = {format_square_sum(flat_values)} = **{sum_sq:.4f}**")
    st.latex(rf"SCT = {sum_sq:.4f} - {results['tc']:.4f} = {results['sct']:.4f}")

with st.expander(f"Paso 4. Suma de cuadrados de {st.session_state.row_name}", expanded=True):
    sq_rows = sum(rt ** 2 for rt in row_totals)
    st.markdown(f"$\\sum Y_{{i.}}^2$ = {format_square_sum(row_totals)} = **{sq_rows:.4f}**")
    st.latex(rf"SCF = \frac{{{sq_rows:.4f}}}{{{k}}} - {results['tc']:.4f} = {results['scf']:.4f}")

with st.expander(f"Paso 5. Suma de cuadrados de {st.session_state.col_name}", expanded=True):
    sq_cols = sum(ct ** 2 for ct in col_totals)
    st.markdown(f"$\\sum Y_{{.j}}^2$ = {format_square_sum(col_totals)} = **{sq_cols:.4f}**")
    st.latex(rf"SCC = \frac{{{sq_cols:.4f}}}{{{k}}} - {results['tc']:.4f} = {results['scc']:.4f}")

with st.expander(f"Paso 6. Suma de cuadrados de {st.session_state.treatment_name}", expanded=True):
    sq_trs = sum(tt ** 2 for tt in tr_totals)
    st.markdown(f"$\\sum Y_{{..t}}^2$ = {format_square_sum(tr_totals)} = **{sq_trs:.4f}**")
    st.latex(rf"SCTr = \frac{{{sq_trs:.4f}}}{{{k}}} - {results['tc']:.4f} = {results['sctr']:.4f}")

with st.expander("Paso 7. Suma de cuadrados del error", expanded=True):
    st.latex(rf"SCE = {results['sct']:.4f} - {results['scf']:.4f} - {results['scc']:.4f} - {results['sctr']:.4f} = {results['sce']:.4f}")

with st.expander("Paso 8. Grados de libertad, cuadrados medios y F calculado", expanded=True):
    st.write(f"gl Filas = {results['gl']['Filas']}")
    st.write(f"gl Columnas = {results['gl']['Columnas']}")
    st.write(f"gl Tratamientos = {results['gl']['Tratamientos']}")
    st.write(f"gl Error = {results['gl']['Error']}")
    st.write(f"gl Total = {results['gl']['Total']}")

    st.latex(rf"CM_{{Filas}} = \frac{{{results['scf']:.4f}}}{{{results['gl']['Filas']}}} = {results['cm']['Filas']:.4f}")
    st.latex(rf"CM_{{Columnas}} = \frac{{{results['scc']:.4f}}}{{{results['gl']['Columnas']}}} = {results['cm']['Columnas']:.4f}")
    st.latex(rf"CM_{{Tratamientos}} = \frac{{{results['sctr']:.4f}}}{{{results['gl']['Tratamientos']}}} = {results['cm']['Tratamientos']:.4f}")
    st.latex(rf"CM_{{Error}} = \frac{{{results['sce']:.4f}}}{{{results['gl']['Error']}}} = {results['cm']['Error']:.4f}")

    st.latex(rf"F_{{Filas}} = \frac{{{results['cm']['Filas']:.4f}}}{{{results['cm']['Error']:.4f}}} = {results['f']['Filas']:.4f}")
    st.latex(rf"F_{{Columnas}} = \frac{{{results['cm']['Columnas']:.4f}}}{{{results['cm']['Error']:.4f}}} = {results['f']['Columnas']:.4f}")
    st.latex(rf"F_{{Tratamientos}} = \frac{{{results['cm']['Tratamientos']:.4f}}}{{{results['cm']['Error']:.4f}}} = {results['f']['Tratamientos']:.4f}")
    st.latex(rf"F_{{crítico}} = F_{{1-\alpha;\, gl_1={results['gl']['Tratamientos']},\, gl_2={results['gl']['Error']}}} = {results['f_crit']:.4f}")

anova_df = pd.DataFrame(
    {
        "Fuente": ["Filas", "Columnas", "Tratamientos", "Error", "Total"],
        "SC": [results["scf"], results["scc"], results["sctr"], results["sce"], results["sct"]],
        "gl": [results["gl"]["Filas"], results["gl"]["Columnas"], results["gl"]["Tratamientos"], results["gl"]["Error"], results["gl"]["Total"]],
        "CM": [results["cm"]["Filas"], results["cm"]["Columnas"], results["cm"]["Tratamientos"], results["cm"]["Error"], None],
        "F calculado": [results["f"]["Filas"], results["f"]["Columnas"], results["f"]["Tratamientos"], None, None],
        "F crítico": [results["f_crit"], results["f_crit"], results["f_crit"], None, None],
    }
)

st.subheader("5) Tabla ANOVA")
st.dataframe(anova_df, use_container_width=True)

st.subheader("6) Decisión estadística")
for source in ["Filas", "Columnas", "Tratamientos"]:
    f_calc = results["f"][source]
    decision = "Rechazar H₀" if f_calc > results["f_crit"] else "No rechazar H₀"
    st.write(f"**{source}:** F calculado = {f_calc:.4f}, F crítico = {results['f_crit']:.4f}. Decisión: **{decision}**.")

r_script = build_r_script(
    problem_name=st.session_state.problem_name,
    response_name=st.session_state.response_name,
    row_name=st.session_state.row_name,
    col_name=st.session_state.col_name,
    treatment_name=st.session_state.treatment_name,
    row_labels=row_labels,
    col_labels=col_labels,
    treatment_labels=treatment_labels,
    values=values,
    treatment_square=treatment_square,
    alpha=st.session_state.alpha,
)

st.subheader("7) Script en R para validación")
st.code(r_script, language="r")
st.download_button(
    "Descargar script de R",
    data=r_script,
    file_name="validacion_dcl.R",
    mime="text/x-r-source"
)

python_requirements = """streamlit>=1.40
pandas>=2.2
scipy>=1.13
"""

st.subheader("8) Ejecución")
st.code("pip install -r requirements.txt\nstreamlit run dcl_streamlit_app.py", language="bash")
st.code(python_requirements, language="text")

st.download_button(
    "Descargar requirements.txt",
    data=python_requirements,
    file_name="requirements.txt",
    mime="text/plain"
)
