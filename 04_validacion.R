# =====================================================================
# TAREA 4 - Validación: gráfico de ajuste (eje log), QQ, KS a mano, cuantiles
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_A.R y resultados/params_A.rds
# =====================================================================
set.seed(2026)
source("funciones_A.R")
p <- readRDS("resultados/params_A.rds")
k <- p$k_w; lam <- p$lam_w; mu <- p$mu_ln; sg <- p$sigma_ln

# ---- 0. Datos ------------------------------------------------------------
archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct
x <- abs(r[r != 0]); m <- length(x)
xs <- sort(x)                                   # observados ordenados

# ---- 1. Muestras sintéticas (20 veces el tamaño de la real) ----------------
N  <- 20 * m
yw <- sort(rweibull_ti(N, k, lam))              # Weibull (familia elegida por AIC)
yl <- sort(rlnorm_bm(N, mu, sg))                # lognormal (se valida también, para comparar)

# ---- 2. Densidades teóricas (fórmulas propias) --------------------------------
dweibull_propia <- function(v, k, lam) (k / lam) * (v / lam)^(k - 1) * exp(-(v / lam)^k)
dlnorm_propia   <- function(v, mu, s) exp(-(log(v) - mu)^2 / (2 * s^2)) / (v * s * sqrt(2 * pi))

# ---- 3. Gráfico de ajuste con eje vertical logarítmico ---------------------------
dir.create("figuras", showWarnings = FALSE)
h <- hist(x, breaks = seq(0, ceiling(max(x) * 4) / 4, by = 0.25), plot = FALSE)
ok <- h$density > 0                              # el log no admite densidad 0
grid_x <- seq(0.01, max(x), length.out = 400)
png("figuras/02_ajuste_log.png", width = 1100, height = 750, res = 130)
plot(h$mids[ok], h$density[ok], log = "y", pch = 16,
     xlab = "|r| (%)", ylab = "densidad (escala log)",
     main = "GBP/USD: ajuste de densidades, eje vertical logarítmico")
lines(grid_x, dweibull_propia(grid_x, k, lam), lwd = 2, col = "blue")
lines(grid_x, dlnorm_propia(grid_x, mu, sg),  lwd = 2, col = "red")
legend("topright", c("Datos (histograma)", "Weibull", "Lognormal"),
       pch = c(16, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 2, 2),
       col = c("black", "blue", "red"), bty = "n")
dev.off()

# ---- 4. QQ-plot (cuantiles observados vs sintéticos) -------------------------------
pi_ <- (1:m - 0.5) / m                           # probabilidades de los m observados
cuant_sint <- function(ysorted, pp) ysorted[ceiling(pp * length(ysorted))]
qw <- cuant_sint(yw, pi_)
ql <- cuant_sint(yl, pi_)
png("figuras/03_qq.png", width = 1300, height = 1100, res = 130)
par(mfrow = c(2, 2))
plot(qw, xs, pch = 16, cex = 0.5, xlab = "cuantil Weibull simulado", ylab = "cuantil observado",
     main = "QQ Weibull"); abline(0, 1, col = "red")
plot(ql, xs, pch = 16, cex = 0.5, xlab = "cuantil lognormal simulado", ylab = "cuantil observado",
     main = "QQ lognormal"); abline(0, 1, col = "red")
plot(qw, xs, pch = 16, cex = 0.5, log = "xy", xlab = "Weibull (log)", ylab = "observado (log)",
     main = "QQ Weibull, ejes log"); abline(0, 1, col = "red")
plot(ql, xs, pch = 16, cex = 0.5, log = "xy", xlab = "lognormal (log)", ylab = "observado (log)",
     main = "QQ lognormal, ejes log"); abline(0, 1, col = "red")
dev.off()

# ---- 5. Distancia de Kolmogorov-Smirnov a mano ---------------------------------------
# (a) Dos muestras: D = max |F_m(t) - G_N(t)| sobre todos los puntos t de ambas muestras
ks_dos_muestras <- function(a, b) {
  a <- sort(a); b <- sort(b)
  t <- sort(c(a, b))
  Fa <- findInterval(t, a) / length(a)           # fracción de a <= t
  Gb <- findInterval(t, b) / length(b)           # fracción de b <= t
  max(abs(Fa - Gb))
}
D_w <- ks_dos_muestras(xs, yw)
D_l <- ks_dos_muestras(xs, yl)

# (b) Una muestra contra la FDA teórica (extra, determinista).
#     La FDA normal se construye con la función de error erf (permitida), programada por serie de Taylor
erf_propia <- function(z) {
  n <- 0:60; fact <- c(1, cumprod(1:60))
  sapply(z, function(zz) 2 / sqrt(pi) * sum((-1)^n * zz^(2 * n + 1) / (fact * (2 * n + 1))))
}
pnorm_propia  <- function(z) 0.5 * (1 + erf_propia(z / sqrt(2)))
plnorm_propia <- function(v, mu, s) pnorm_propia((log(v) - mu) / s)
ks_una_muestra <- function(xsorted, Fteo) {
  m <- length(xsorted); i <- 1:m
  max(max(i / m - Fteo), max(Fteo - (i - 1) / m))
}
D1_w <- ks_una_muestra(xs, pweibull_propia(xs, k, lam))
D1_l <- ks_una_muestra(xs, plnorm_propia(xs, mu, sg))
# Valor crítico asintótico aproximado al 5% (referencia; ver advertencia en el texto)
crit <- 1.358 * sqrt((m + N) / (m * N))
crit1 <- 1.358 / sqrt(m)

# ---- 6. Cuantiles 50, 90 y 99 ------------------------------------------------------------
pp <- c(0.50, 0.90, 0.99)
tabla_q <- data.frame(
  p = pp,
  observado      = quantile(x, pp, type = 7),
  Weibull_sint   = cuant_sint(yw, pp),
  Weibull_teo    = lam * (-log(1 - pp))^(1 / k),
  lognormal_sint = cuant_sint(yl, pp)
)
rownames(tabla_q) <- NULL

cat("\n== Distancia KS (dos muestras: observada vs sintética de tamaño", N, ") ==\n")
cat("Weibull  :", D_w, "\nLognormal:", D_l, "\n")
cat("Valor crítico aproximado 5% (dos muestras):", crit, "\n")
cat("\n== Distancia KS (una muestra, contra la FDA teórica) ==\n")
cat("Weibull  :", D1_w, "\nLognormal:", D1_l, "\n")
cat("Valor crítico aproximado 5% (una muestra):", crit1, "\n")
cat("ADVERTENCIA: los parámetros se estimaron con los mismos datos, así que el valor crítico\n",
    "es solo una referencia; D se interpreta como medida descriptiva de distancia.\n")
cat("\n== Cuantiles 50, 90, 99 ==\n")
print(tabla_q, digits = 5)

saveRDS(list(D_w = D_w, D_l = D_l, D1_w = D1_w, D1_l = D1_l, tabla_q = tabla_q),
        "resultados/validacion_A.rds")

# ---- 7. Chequeo de robustez de los generadores (muestra muy grande) ----------------------
N2 <- 2e6
z2 <- rnorm_bm(N2)
cat("\n== Robustez Box-Muller, N =", N2, "==\n")
cat("media (SE ~", 1 / sqrt(N2), "):", mean(z2), "\n")
cat("varianza (teórica 1):", mean((z2 - mean(z2))^2), "\n")
w2 <- sort(rweibull_ti(N2, k, lam))
q99_teo <- lam * (-log(0.01))^(1 / k)
se_q99  <- sqrt(0.99 * 0.01 / N2) / dweibull_propia(q99_teo, k, lam)
cat("Weibull q99 simulado / teórico / SE:", cuant_sint(w2, 0.99), "/", q99_teo, "/", se_q99, "\n")

