# B2 (completo) - Jackknife, intervalos BCa e intervalos asintóticos
# Estadísticos: desviación estándar y VaR 99% diario, sobre los retornos ORIGINALES (con ceros)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: resultados/B2.rds (de 06_B2.R)
set.seed(2026)

# ---- 0. Datos y réplicas bootstrap de 06_B2.R -----------------------------------------
archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct; n <- length(r)
b2 <- readRDS("resultados/B2.rds")
rep_sd <- b2$rep_sd; rep_var <- b2$rep_var; B <- b2$B

est_sd  <- function(v) sqrt(sum((v - mean(v))^2) / (length(v) - 1))
est_var <- function(v) quantile(-v, 0.99, type = 7, names = FALSE)
theta_sd <- est_sd(r); theta_var <- est_var(r)

# ---- 1. Normal estándar propia: FDA (vía erf por serie de Taylor) e inversa (bisección) ----
erf_propia <- function(z) {
  k <- 0:60; fact <- c(1, cumprod(1:60))
  sapply(z, function(zz) 2 / sqrt(pi) * sum((-1)^k * zz^(2 * k + 1) / (fact * (2 * k + 1))))
}
pnorm_propia <- function(z) 0.5 * (1 + erf_propia(z / sqrt(2)))
qnorm_propia <- function(p) {
  sapply(p, function(pp) {
    lo <- -6; hi <- 6
    for (i in 1:80) { mid <- (lo + hi) / 2; if (pnorm_propia(mid) < pp) lo <- mid else hi <- mid }
    (lo + hi) / 2
  })
}
z975 <- qnorm_propia(0.975)
cat("z_0.975 propio:", z975, "| verificación con qnorm:", qnorm(0.975), "\n")

# ---- 2. Jackknife (dejar uno fuera) ------------------------------------------------------
jackknife <- function(r, est) {
  n  <- length(r)
  th <- est(r)
  jk <- sapply(1:n, function(i) est(r[-i]))          # estadístico sin la observación i
  jm <- mean(jk)
  list(theta = th, jk = jk,
       sesgo = (n - 1) * (jm - th),                   # sesgo jackknife
       ee    = sqrt((n - 1) / n * sum((jk - jm)^2)),  # error estándar jackknife
       a     = sum((jm - jk)^3) / (6 * sum((jm - jk)^2)^1.5),  # aceleración para BCa
       distintos = length(unique(round(jk, 10))))     # valores distintos de las réplicas
}
jk_sd  <- jackknife(r, est_sd)
jk_var <- jackknife(r, est_var)

# ---- 3. BCa (la aceleración sale de las réplicas jackknife de la muestra ORIGINAL) --------
bca <- function(rep, theta, a, alfa = 0.05) {
  B  <- length(rep)
  p0 <- (sum(rep < theta) + 0.5 * sum(rep == theta)) / B   # proporción de réplicas por debajo (empates: mitad)
  z0 <- qnorm_propia(p0)                                    # corrección de sesgo
  za <- qnorm_propia(c(alfa / 2, 1 - alfa / 2))
  niveles <- pnorm_propia(z0 + (z0 + za) / (1 - a * (z0 + za)))   # niveles ajustados
  list(z0 = z0, niveles = niveles, ic = quantile(rep, niveles, type = 7, names = FALSE))
}
bca_sd  <- bca(rep_sd,  theta_sd,  jk_sd$a)
bca_var <- bca(rep_var, theta_var, jk_var$a)

# ---- 4. Intervalos asintóticos --------------------------------------------------------------
# (a) Desviación estándar, método delta: Var(s) ~ (m4 - s^4) / (4 s^2 n)
m4 <- mean((r - mean(r))^4); s2 <- theta_sd^2
ee_sd_as <- sqrt((m4 - s2^2) / (4 * s2 * n))
ic_sd_as <- theta_sd + c(-1, 1) * z975 * ee_sd_as

# (b) VaR 99%, cuantil muestral: Var(q) ~ p(1-p) / (n f(q)^2), con f estimada por núcleo gaussiano
L <- -r
dens_nucleo <- function(q, L, h) sum(exp(-0.5 * ((q - L) / h)^2)) / (length(L) * h * sqrt(2 * pi))
iqr <- quantile(L, 0.75, names = FALSE) - quantile(L, 0.25, names = FALSE)
h0  <- 0.9 * min(est_sd(L), iqr / 1.34) * n^(-1 / 5)         # regla de Silverman
ee_var_as <- function(h) sqrt(0.99 * 0.01 / n) / dens_nucleo(theta_var, L, h)
sens <- data.frame(multiplo_h = c(0.5, 1, 2), h = h0 * c(0.5, 1, 2))
sens$densidad <- sapply(sens$h, function(h) dens_nucleo(theta_var, L, h))
sens$EE       <- sapply(sens$h, ee_var_as)
cat("\n== Sensibilidad del EE asintótico del VaR al ancho de banda ==\n")
print(sens, digits = 5, row.names = FALSE)
ee_var_as0 <- ee_var_as(h0)
ic_var_as  <- theta_var + c(-1, 1) * z975 * ee_var_as0

# ---- 5. Resultados -----------------------------------------------------------------------------
cat("\n== Jackknife ==\n")
tab_jk <- data.frame(
  estadistico = c("Desv. estándar", "VaR 99%"),
  estimacion  = c(theta_sd, theta_var),
  sesgo_jack  = c(jk_sd$sesgo, jk_var$sesgo),
  EE_jack     = c(jk_sd$ee, jk_var$ee),
  corregida   = c(theta_sd - jk_sd$sesgo, theta_var - jk_var$sesgo),
  valores_distintos = c(jk_sd$distintos, jk_var$distintos))
print(tab_jk, digits = 5, row.names = FALSE)
cat("Aceleración a: sd =", jk_sd$a, " VaR =", jk_var$a, "\n")
cat("BCa: z0 sd =", bca_sd$z0, " VaR =", bca_var$z0, "\n")
cat("Niveles BCa ajustados: sd =", bca_sd$niveles, " VaR =", bca_var$niveles, "\n")

perc <- function(rep) quantile(rep, c(0.025, 0.975), type = 7, names = FALSE)
cat("\n== Tabla comparativa de intervalos al 95% ==\n")
tab_ic <- data.frame(
  estadistico = rep(c("Desv. estándar", "VaR 99%"), each = 3),
  metodo = rep(c("Percentil", "BCa", "Asintótico"), 2),
  inferior = c(perc(rep_sd)[1], bca_sd$ic[1], ic_sd_as[1],
               perc(rep_var)[1], bca_var$ic[1], ic_var_as[1]),
  superior = c(perc(rep_sd)[2], bca_sd$ic[2], ic_sd_as[2],
               perc(rep_var)[2], bca_var$ic[2], ic_var_as[2]))
tab_ic$amplitud <- tab_ic$superior - tab_ic$inferior
print(tab_ic, digits = 5, row.names = FALSE)

cat("\n== Errores estándar por método ==\n")
tab_ee <- data.frame(
  estadistico = c("Desv. estándar", "VaR 99%"),
  bootstrap = c(sqrt(sum((rep_sd - mean(rep_sd))^2) / (B - 1)),
                sqrt(sum((rep_var - mean(rep_var))^2) / (B - 1))),
  jackknife = c(jk_sd$ee, jk_var$ee),
  asintotico = c(ee_sd_as, ee_var_as0))
print(tab_ee, digits = 5, row.names = FALSE)

# ---- 6. Gráfico de los intervalos ----------------------------------------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/05_intervalos.png", width = 1400, height = 560, res = 130)
par(mfrow = c(1, 2), mar = c(4, 7, 3, 1))
for (j in 1:2) {
  sub <- tab_ic[tab_ic$estadistico == c("Desv. estándar", "VaR 99%")[j], ]
  th  <- c(theta_sd, theta_var)[j]
  plot(NA, xlim = range(sub$inferior, sub$superior), ylim = c(0.5, 3.5), yaxt = "n",
       xlab = c("Desviación estándar", "VaR 99% diario")[j], ylab = "",
       main = paste("Intervalos 95%:", c("desv. estándar", "VaR 99%")[j]))
  axis(2, at = 1:3, labels = sub$metodo, las = 1)
  segments(sub$inferior, 1:3, sub$superior, 1:3, lwd = 3)
  abline(v = th, lty = 2)
}
dev.off()

saveRDS(list(tab_jk = tab_jk, tab_ic = tab_ic, tab_ee = tab_ee, sens = sens,
             a = c(sd = jk_sd$a, var = jk_var$a)), "resultados/B2_completo.rds")
cat("\nGuardado: resultados/B2_completo.rds y figuras/05_intervalos.png\n")

