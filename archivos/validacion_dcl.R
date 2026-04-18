# ==========================================================
# VALIDACIÓN COMPUTACIONAL EN R
# Diseño en Cuadro Latino (DCL)
# ==========================================================

# Instale el paquete si es necesario:
# install.packages("easyanova")

library(easyanova)

tratamiento <- as.factor(c("A", "B", "C", "D", "B", "C", "D", "A", "C", "D", "A", "B", "D", "A", "B", "C"))
fila <- as.factor(c("F1", "F1", "F1", "F1", "F2", "F2", "F2", "F2", "F3", "F3", "F3", "F3", "F4", "F4", "F4", "F4"))
columna <- as.factor(c("C1", "C2", "C3", "C4", "C1", "C2", "C3", "C4", "C1", "C2", "C3", "C4", "C1", "C2", "C3", "C4"))
eficiencia_térmica_(%) <- c(25, 37, 42, 33, 35, 46, 35, 36, 48, 37, 39, 48, 36, 41, 47, 55)

datos <- data.frame(tratamiento, fila, columna, eficiencia_térmica_(%))

cat("\n--- Datos del experimento ---\n")
print(datos)

# ----------------------------------------------------------
# Opción 1: easyanova
# ----------------------------------------------------------
cat("\n--- ANOVA con easyanova ---\n")
anova_dcl <- eal(datos, design = 3)
print(anova_dcl$`Analysis of variance`)

# ----------------------------------------------------------
# Opción 2: aov() base R
# ----------------------------------------------------------
cat("\n--- ANOVA con aov() ---\n")
modelo <- aov(eficiencia_térmica_(%) ~ fila + columna + tratamiento, data = datos)
summary(modelo)
