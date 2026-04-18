# App multipágina de Diseños Experimentales

Aplicación desarrollada con **Streamlit** para apoyar la enseñanza y el aprendizaje de los principales diseños experimentales, con énfasis en el desarrollo **manual paso a paso**, la **tabla ANOVA** y la **validación computacional**.

## Propósito

Esta app busca sistematizar progresivamente los siguientes diseños:

- DCA: Diseño Completamente Aleatorizado
- DBCA: Diseño en Bloques Completos al Azar
- DCL: Diseño en Cuadro Latino
- DPD: Diseño de Parcelas Divididas
- Diseños Factoriales
- Metodología de Superficie de Respuesta

## Enfoque pedagógico

Cada página está pensada para:

1. Recibir el planteamiento del problema.
2. Permitir el ingreso estructurado de datos.
3. Guiar el desarrollo manual paso a paso.
4. Construir la tabla ANOVA.
5. Calcular e interpretar los valores de F.
6. Mostrar una validación computacional en R.

## Estado actual del proyecto

Actualmente se encuentra funcional el módulo:

- **DCL (Diseño en Cuadro Latino)**

Los demás módulos están preparados dentro de la arquitectura multipágina y se desarrollarán progresivamente.

## Estructura del proyecto

```text
Design/
│
├── app.py
├── requirements.txt
├── README.md
├── .gitignore
├── pages/
│   ├── 00_inicio.py
│   ├── 01_dca.py
│   ├── 02_dbca.py
│   ├── 03_dcl.py
│   ├── 04_dpd.py
│   ├── 05_factoriales.py
│   └── 06_superficie_respuesta.py
└── utils/
    └── dcl_core.py