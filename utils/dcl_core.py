from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import f


@dataclass
class DCLResults:
    k: int
    grand_total: float
    correction_term: float
    sst: float
    ss_rows: float
    ss_cols: float
    ss_trt: float
    ss_error: float
    df_rows: int
    df_cols: int
    df_trt: int
    df_error: int
    df_total: int
    ms_rows: float
    ms_cols: float
    ms_trt: float
    ms_error: float
    f_rows: float
    f_cols: float
    f_trt: float
    f_crit: float
    row_totals: pd.Series
    col_totals: pd.Series
    trt_totals: pd.Series
    data_long: pd.DataFrame


def validate_latin_square(treatments_df: pd.DataFrame) -> Tuple[bool, List[str]]:
    messages: List[str] = []
    k = treatments_df.shape[0]

    if treatments_df.shape[0] != treatments_df.shape[1]:
        return False, ["La matriz de tratamientos debe ser cuadrada (k × k)."]

    unique_treatments = sorted(pd.unique(treatments_df.values.ravel()))
    if len(unique_treatments) != k:
        messages.append(
            f"Debe haber exactamente {k} tratamientos distintos y se encontraron {len(unique_treatments)}."
        )

    for i in range(k):
        row_values = list(treatments_df.iloc[i, :])
        if len(set(row_values)) != k:
            messages.append(f"La fila {i+1} no contiene los {k} tratamientos una sola vez.")

    for j in range(k):
        col_values = list(treatments_df.iloc[:, j])
        if len(set(col_values)) != k:
            messages.append(f"La columna {j+1} no contiene los {k} tratamientos una sola vez.")

    return len(messages) == 0, messages


def build_long_dataframe(
    response_df: pd.DataFrame,
    treatments_df: pd.DataFrame,
    row_labels: List[str],
    col_labels: List[str],
    row_factor_name: str,
    col_factor_name: str,
    treatment_name: str,
    response_name: str,
) -> pd.DataFrame:
    records = []
    for i, row in enumerate(row_labels):
        for j, col in enumerate(col_labels):
            records.append(
                {
                    row_factor_name: row,
                    col_factor_name: col,
                    treatment_name: str(treatments_df.iloc[i, j]),
                    response_name: float(response_df.iloc[i, j]),
                }
            )
    return pd.DataFrame(records)


def compute_dcl_anova(
    response_df: pd.DataFrame,
    treatments_df: pd.DataFrame,
    row_labels: List[str],
    col_labels: List[str],
    row_factor_name: str,
    col_factor_name: str,
    treatment_name: str,
    response_name: str,
    alpha: float = 0.05,
) -> DCLResults:
    k = response_df.shape[0]
    y = response_df.astype(float)

    data_long = build_long_dataframe(
        response_df,
        treatments_df,
        row_labels,
        col_labels,
        row_factor_name,
        col_factor_name,
        treatment_name,
        response_name,
    )

    grand_total = float(y.to_numpy().sum())
    correction_term = grand_total**2 / (k**2)
    sst = float((y.to_numpy() ** 2).sum() - correction_term)

    row_totals = y.sum(axis=1)
    col_totals = y.sum(axis=0)
    trt_totals = data_long.groupby(treatment_name)[response_name].sum().sort_index()

    ss_rows = float(((row_totals**2).sum() / k) - correction_term)
    ss_cols = float(((col_totals**2).sum() / k) - correction_term)
    ss_trt = float(((trt_totals**2).sum() / k) - correction_term)
    ss_error = float(sst - ss_rows - ss_cols - ss_trt)

    df_rows = k - 1
    df_cols = k - 1
    df_trt = k - 1
    df_error = (k - 1) * (k - 2)
    df_total = k**2 - 1

    ms_rows = ss_rows / df_rows
    ms_cols = ss_cols / df_cols
    ms_trt = ss_trt / df_trt
    ms_error = ss_error / df_error

    f_rows = ms_rows / ms_error
    f_cols = ms_cols / ms_error
    f_trt = ms_trt / ms_error
    f_crit = float(f.ppf(1 - alpha, df_trt, df_error))

    return DCLResults(
        k=k,
        grand_total=grand_total,
        correction_term=correction_term,
        sst=sst,
        ss_rows=ss_rows,
        ss_cols=ss_cols,
        ss_trt=ss_trt,
        ss_error=ss_error,
        df_rows=df_rows,
        df_cols=df_cols,
        df_trt=df_trt,
        df_error=df_error,
        df_total=df_total,
        ms_rows=ms_rows,
        ms_cols=ms_cols,
        ms_trt=ms_trt,
        ms_error=ms_error,
        f_rows=f_rows,
        f_cols=f_cols,
        f_trt=f_trt,
        f_crit=f_crit,
        row_totals=row_totals,
        col_totals=col_totals,
        trt_totals=trt_totals,
        data_long=data_long,
    )


def anova_table_df(res: DCLResults, row_factor_name: str, col_factor_name: str, treatment_name: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Fuente": [row_factor_name, col_factor_name, treatment_name, "Error", "Total"],
            "SC": [res.ss_rows, res.ss_cols, res.ss_trt, res.ss_error, res.sst],
            "gl": [res.df_rows, res.df_cols, res.df_trt, res.df_error, res.df_total],
            "CM": [res.ms_rows, res.ms_cols, res.ms_trt, res.ms_error, np.nan],
            "F calculado": [res.f_rows, res.f_cols, res.f_trt, np.nan, np.nan],
        }
    )


def build_r_script(
    data_long: pd.DataFrame,
    row_factor_name: str,
    col_factor_name: str,
    treatment_name: str,
    response_name: str,
) -> str:
    row_values = ", ".join([f'"{v}"' for v in data_long[row_factor_name]])
    col_values = ", ".join([f'"{v}"' for v in data_long[col_factor_name]])
    trt_values = ", ".join([f'"{v}"' for v in data_long[treatment_name]])
    resp_values = ", ".join([str(v) for v in data_long[response_name]])

    return f'''# ==========================================
# VALIDACIÓN EN R DEL DISEÑO EN CUADRO LATINO
# ==========================================

# install.packages("easyanova")
library(easyanova)

{treatment_name} <- as.factor(c({trt_values}))
{row_factor_name} <- as.factor(c({row_values}))
{col_factor_name} <- as.factor(c({col_values}))
{response_name} <- c({resp_values})

datos <- data.frame({treatment_name}, {row_factor_name}, {col_factor_name}, {response_name})

cat("\\n--- ANOVA con easyanova ---\\n")
resultado_easyanova <- eal(datos, design = 3)
print(resultado_easyanova$`Analysis of variance`)

cat("\\n--- ANOVA con aov() ---\\n")
modelo <- aov({response_name} ~ {row_factor_name} + {col_factor_name} + {treatment_name}, data = datos)
summary(modelo)

cat("\\n--- F crítico para tratamientos ---\\n")
alpha <- 0.05
k <- length(levels({treatment_name}))
gl1 <- k - 1
gl2 <- (k - 1) * (k - 2)
qf(1 - alpha, df1 = gl1, df2 = gl2)
'''
