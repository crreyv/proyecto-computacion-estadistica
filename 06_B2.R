# =====================================================================
# B2 (avance) - Bootstrap percentil sobre los retornos ORIGINALES (con ceros)
# Estadísticos: desviación estándar y VaR 99% diario (cuantil 99% de la pérdida -r)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# =====================================================================
set.seed(2026)

archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct
n <- length(r)
cat("n =", n, "| ceros incluidos:", sum(r == 0), "\n")

# ---- Estadísticos ---------------------------------------------------------------------
est_sd  <- function(v) sqrt(sum((v - mean(v))^2) / (length(v) - 1))
est_var <- function(v) quantile(-v, 0.99, type = 7, names = FALSE)   # VaR 99% = cuantil 99% de -r
theta_sd  <- est_sd(r)
theta_var <- est_var(r)
cat("Estimaciones puntuales: sd =", theta_sd, " VaR99 =", theta_var, "\n")

# ---- Bootstrap no paramétrico (remuestreo con reemplazo), programado a mano -----------------
B <- 5000
rep_sd <- numeric(B); rep_var <- numeric(B)
for (b in 1:B) {
  idx <- floor(runif(n) * n) + 1          # n índices uniformes en {1,...,n}
  rb  <- r[idx]                           # muestra bootstrap
  rep_sd[b]  <- est_sd(rb)
  rep_var[b] <- est_var(rb)
}

resumen <- function(rep, theta, nombre) {
  li <- quantile(rep, 0.025, type = 7, names = FALSE)
  ls <- quantile(rep, 0.975, type = 7, names = FALSE)
  data.frame(estadistico = nombre, estimacion = theta,
             sesgo_boot = mean(rep) - theta,
             EE_boot = sqrt(sum((rep - mean(rep))^2) / (length(rep) - 1)),
             IC95_percentil_inf = li, IC95_percentil_sup = ls)
}
tabla <- rbind(resumen(rep_sd, theta_sd, "Desv. estándar"),
               resumen(rep_var, theta_var, "VaR 99%"))
cat("\n== Bootstrap percentil, B =", B, "==\n")
print(tabla, digits = 5, row.names = FALSE)

# ---- Gráfico de la distribución bootstrap ---------------------------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/04_bootstrap.png", width = 1400, height = 560, res = 130)
par(mfrow = c(1, 2))
for (j in 1:2) {
  rep <- if (j == 1) rep_sd else rep_var
  th  <- if (j == 1) theta_sd else theta_var
  tt  <- if (j == 1) "Desviación estándar" else "VaR 99% diario"
  hist(rep, breaks = 40, main = paste("Bootstrap:", tt), xlab = tt, col = "grey85")
  abline(v = th, lwd = 2)
  abline(v = quantile(rep, c(0.025, 0.975), type = 7), lty = 2, col = "red", lwd = 2)
}
dev.off()
dir.create("resultados", showWarnings = FALSE)
saveRDS(list(tabla = tabla, rep_sd = rep_sd, rep_var = rep_var, B = B), "resultados/B2.rds")
cat("\nFigura: figuras/04_bootstrap.png\n")