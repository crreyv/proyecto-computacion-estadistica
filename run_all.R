# Ejecuta todo el proyecto de principio a fin (cada script fija su propia semilla: set.seed(2026))
# Uso: abrir R en la carpeta del proyecto y ejecutar  source("run_all.R")
# Archivos de funciones (se cargan con source() desde los scripts, no se ejecutan solos):
#   funciones_A.R (generadores, Fase A) y funciones_EM.R (pasos E y M del EM, Fase C)
scripts <- c("01_datos.R",            # datos, ceros y variable de trabajo
             "02_mle.R",              # Fase A: máxima verosimilitud y AIC
             "03_generadores.R",      # Fase A: comprobación de generadores
             "04_validacion.R",       # Fase A: QQ, KS a mano, cuantiles
             "05_B1.R",               # Fase B1: Monte Carlo y muestreo por importancia
             "06_B2.R",               # Fase B2: bootstrap percentil
             "07_B2_completo.R",      # Fase B2: jackknife, BCa e intervalos asintóticos
             "08_EM.R",               # Fase C: algoritmo EM y convergencia
             "09_EM_inicios.R",       # Fase C: inicios aleatorios, tolerancias y degenerados
             "10_EM_interpretacion.R",# Fase C: BIC, regímenes y cola implicada
             "11_FaseD.R")            # Fase D: VaR por modelo y autocorrelación
for (s in scripts) {
  cat("\n######## Ejecutando", s, "########\n")
  source(s, echo = FALSE)
}
cat("\nListo. Ahora compile el informe (Knit en RStudio o rmarkdown::render()).\n")
