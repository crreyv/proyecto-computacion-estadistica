# Ejecuta todo el proyecto de principio a fin (cada script fija su propia semilla: set.seed(2026))
# Uso: abrir R en la carpeta del proyecto y ejecutar  source("run_all.R")
scripts <- c("01_datos.R", "02_mle.R", "03_generadores.R", "04_validacion.R",
             "05_B1.R", "06_B2.R")
for (s in scripts) {
  cat("\n######## Ejecutando", s, "########\n")
  source(s, echo = FALSE)
}
cat("\nListo. Ahora compile 'informe_avance.Rmd' (Knit en RStudio o rmarkdown::render()).\n")
