# TAREA 1 - Carga de datos, ceros y variable de trabajo (GBP/USD)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2
# Ejecutar desde la carpeta raíz del proyecto, con esta estructura:
#   proyecto/
#     data/retornos_diarios_12_series.csv
#     R/01_datos.R
# =====================================================================

set.seed(2026)   # semilla fija (aquí aún no se simula; se deja desde el inicio)

# ---- 1. Lectura con ruta RELATIVA -----------------------------------
archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)   # ambas son rutas relativas
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) {
  stop("No encuentro '", archivo, "'. Carpeta de trabajo actual: ", getwd(),
       "\nPonga el CSV en 'data/' o en esa carpeta, o cambie la carpeta de trabajo ",
       "(RStudio: Session > Set Working Directory > To Source File Location).")
}
cat("Leyendo:", ruta, "\n")
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")

cat("Columnas:", names(datos), "\n")
cat("Series disponibles en el archivo:\n")
print(table(datos$serie))

# ---- 2. Filtro: SOLO GBP/USD, ordenada en el tiempo -----------------
gbp <- datos[datos$serie == "GBP/USD", ]
gbp <- gbp[order(gbp$t), ]
r   <- gbp$retorno_pct          # retornos originales (CON ceros): se usan en B2

# ---- 3. Verificaciones contra la guía --------------------------------
n <- length(r)
cat("\n== Verificaciones ==\n")
cat("n observaciones (guía: 1866):", n, "\n")
cat("t va de", min(gbp$t), "a", max(gbp$t), "\n")
cat("t consecutivo sin saltos:", all(diff(gbp$t) == 1), "\n")
cat("valores faltantes (NA):", sum(is.na(r)), "\n")

# ---- 4. Ceros exactos -------------------------------------------------
es_cero  <- (r == 0)            # TRUE si el retorno es EXACTAMENTE 0
n_ceros  <- sum(es_cero)
pct_ceros <- 100 * n_ceros / n
cat("\n== Ceros exactos ==\n")
cat("cantidad de ceros:", n_ceros, "\n")
cat(sprintf("porcentaje de ceros: %.3f %%\n", pct_ceros))
cat("menor |r| entre los no nulos:", min(abs(r[!es_cero])), "\n")
cat("pares de ceros en días consecutivos:",
    sum(es_cero[-1] & es_cero[-n]), "\n")

# ---- 5. Variable de trabajo de la Fase A: |r| sin ceros -------------
x <- abs(r[!es_cero])
cat("\n== Variable de trabajo x = |r| (sin ceros) ==\n")
cat("tamaño m:", length(x), "\n")
cat("media:", mean(x), " desv.est.:", sd(x), "\n")
cat("mínimo:", min(x), " máximo:", max(x), "\n")
cat("cuantiles 50, 90, 99:", quantile(x, c(0.50, 0.90, 0.99), type = 7), "\n")
cat("retornos > 0:", sum(r > 0), "  retornos < 0:", sum(r < 0), "\n")

# ---- 6. Estadísticos descriptivos de los retornos originales --------
cat("\n== Retornos originales (con ceros) ==\n")
cat("media:", mean(r), " desv.est.:", sd(r), "\n")
m2 <- mean((r - mean(r))^2)
cat("asimetría:", mean((r - mean(r))^3) / m2^1.5, "\n")
cat("curtosis en exceso:", mean((r - mean(r))^4) / m2^2 - 3, "\n")

# ---- 7. Gráfico exploratorio ----------------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/01_exploracion.png", width = 1400, height = 500, res = 130)
par(mfrow = c(1, 2))
plot(gbp$t, r, type = "h", xlab = "t (día)", ylab = "retorno (%)",
     main = "GBP/USD: retornos diarios")
hist(x, breaks = 60, freq = FALSE, xlab = "|r| (%)",
     main = "Variable de trabajo |r| sin ceros")
dev.off()
cat("\nFigura guardada en figuras/01_exploracion.png\n")
