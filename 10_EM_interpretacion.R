# =====================================================================
# FASE C (parte 3) - BIC frente a una normal, interpretación de los regímenes y cola implicada
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_EM.R y resultados/params_A.rds (de 02_mle.R)
# =====================================================================
set.seed(2026)
source("funciones_EM.R")
pA <- readRDS("resultados/params_A.rds")

archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct; n <- length(r); s <- sd(r)

# ---- 1. Modelo final: EM con tolerancia 1e-10; componente 1 = calma (menor desviación) -----------
ini_base <- list(pi = 0.7, mu = c(mean(r), mean(r)), sigma = c(0.5 * s, 2 * s))
fit <- em_mezcla(r, ini_base, tol = 1e-10)
th <- fit$theta
if (th$sigma[1] > th$sigma[2])
  th <- list(pi = 1 - th$pi, mu = rev(th$mu), sigma = rev(th$sigma))
ll_mez <- loglik(r, th)
cat("EM final: iteraciones =", fit$iter, "| log-verosimilitud =", ll_mez, "\n")

# ---- 2. Comparación con UNA sola normal (BIC y AIC) ---------------------------------------------------
mu0 <- mean(r); sg0 <- sqrt(mean((r - mu0)^2))               # MLE de la normal única
ll_nor <- sum(log(dnorm_propia(r, mu0, sg0)))
tab_bic <- data.frame(
  modelo = c("Normal única", "Mezcla de 2 normales"),
  parametros = c(2, 5),
  logLik = c(ll_nor, ll_mez))
tab_bic$AIC <- 2 * tab_bic$parametros - 2 * tab_bic$logLik
tab_bic$BIC <- tab_bic$parametros * log(n) - 2 * tab_bic$logLik
cat("\n== Normal única frente a mezcla ==\n")
print(tab_bic, digits = 7, row.names = FALSE)
cat("Delta BIC (normal - mezcla):", tab_bic$BIC[1] - tab_bic$BIC[2], "\n")
cat("Estadístico de razón de verosimilitud 2*(ll_mezcla - ll_normal):", 2 * (ll_mez - ll_nor),
    "(NO se compara con una chi-cuadrado: ver informe)\n")

# ---- 3. Interpretación de los regímenes ---------------------------------------------------------------------
gamma  <- paso_E(r, th)           # P(calma | dato)
p_turb <- 1 - gamma               # P(turbulencia | dato)
tab_reg <- data.frame(
  regimen = c("Calma", "Turbulencia"),
  proporcion = c(th$pi, 1 - th$pi),
  media = th$mu, desviacion = th$sigma)
cat("\n== Regímenes ==\n")
print(tab_reg, digits = 6, row.names = FALSE)
cat("Probabilidad media de turbulencia:", mean(p_turb), "(debe ser 1 - pi)\n")
cat("Días con P(turbulencia) > 0.5:", sum(p_turb > 0.5), "(", mean(p_turb > 0.5), ")\n")
cat("Días con P(turbulencia) > 0.9:", sum(p_turb > 0.9), "(", mean(p_turb > 0.9), ")\n")
cat("Días con P(turbulencia) > 0.99:", sum(p_turb > 0.99), "(", mean(p_turb > 0.99), ")\n")
cat("P(turbulencia) media en los días de retorno exactamente cero:", mean(p_turb[r == 0]), "\n")
cambios <- sum((p_turb[-1] > 0.5) != (p_turb[-n] > 0.5))
cat("Cambios de régimen día a día (corte 0.5):", cambios, "| racha media (días):", n / (cambios + 1), "\n")
orden <- order(-p_turb)
cat("15 días (t) con mayor P(turbulencia):", sort(orden[1:15]), "\n")

# ---- 4. Gráfico de la probabilidad de turbulencia en el tiempo ------------------------------------------------
dir.create("figuras", showWarnings = FALSE)
w <- 60
cs <- cumsum(c(0, p_turb))
media_movil <- (cs[(w + 1):(n + 1)] - cs[1:(n - w + 1)]) / w     # media móvil de w días
png("figuras/08_prob_turbulencia.png", width = 1400, height = 560, res = 130)
plot(gbp$t, p_turb, pch = 16, cex = 0.35, col = "grey60", ylim = c(0, 1),
     xlab = "t (día)", ylab = "P(turbulencia | retorno)",
     main = "Probabilidad diaria de pertenecer al régimen turbulento")
lines((1:(n - w + 1)) + (w - 1) / 2, media_movil, col = "red", lwd = 2)
abline(h = 0.5, lty = 2)
legend("topleft", c("Probabilidad diaria", paste("Media móvil de", w, "días")),
       pch = c(16, NA), lty = c(NA, 1), col = c("grey60", "red"), lwd = c(NA, 2), bty = "n")
dev.off()

# ---- 5. Cola implicada: P(|r| > percentil 99 de los datos) ------------------------------------------------------
erf_propia <- function(z) {                       # serie de Taylor (101 términos): precisa para |z| moderado
  k <- 0:100; fact <- c(1, cumprod(1:100))
  sapply(z, function(zz) 2 / sqrt(pi) * sum((-1)^k * zz^(2 * k + 1) / (fact * (2 * k + 1))))
}
# 1 - Phi(z): serie de erf para z < 3 y fracción continua (cociente de Mills) para z >= 3,
# porque 1 - erf pierde precisión en las colas lejanas por cancelación
cola_normal <- function(z) {
  sapply(z, function(zz) {
    if (zz < 3) return(0.5 * (1 - erf_propia(zz / sqrt(2))))
    t <- zz
    for (j in 60:1) t <- zz + j / t
    exp(-0.5 * zz^2) / sqrt(2 * pi) / t
  })
}
# P(|X| > q) si X ~ N(mu, sigma^2):  P(X > q) + P(X < -q)
p_abs <- function(q, mu, sigma) cola_normal((q - mu) / sigma) + cola_normal((q + mu) / sigma)

q99 <- quantile(abs(r), 0.99, type = 7, names = FALSE)          # percentil 99 de |r| (retornos originales)
obs <- mean(abs(r) > q99)
p_mezcla <- th$pi * p_abs(q99, th$mu[1], th$sigma[1]) + (1 - th$pi) * p_abs(q99, th$mu[2], th$sigma[2])
p_normal <- p_abs(q99, mu0, sg0)
p_weib   <- exp(-(q99 / pA$lam_w)^pA$k_w)                       # modelo de la Fase A (sin ceros)
p0       <- mean(r == 0)
tab_cola <- data.frame(
  fuente = c("Frecuencia observada", "Mezcla de 2 normales", "Weibull (Fase A)",
             "Weibull (Fase A) x (1 - proporción de ceros)", "Normal única"),
  probabilidad = c(obs, p_mezcla, p_weib, p_weib * (1 - p0), p_normal))
tab_cola$razon_vs_observada <- tab_cola$probabilidad / obs
cat("\n== Cola implicada: P(|r| >", q99, ") ==\n")
cat("Datos que superan el percentil 99:", sum(abs(r) > q99), "de", n, "\n")
print(tab_cola, digits = 5, row.names = FALSE)
# Verificación con la función de la librería (NO oficial)
pn <- function(z) pnorm(z)
v_mez <- th$pi * ((1 - pn((q99 - th$mu[1]) / th$sigma[1])) + pn((-q99 - th$mu[1]) / th$sigma[1])) +
         (1 - th$pi) * ((1 - pn((q99 - th$mu[2]) / th$sigma[2])) + pn((-q99 - th$mu[2]) / th$sigma[2]))
cat("Verificación con pnorm (solo comprobación) - mezcla:", v_mez, " propio:", p_mezcla, "\n")

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(theta = th, iter = fit$iter, tab_bic = tab_bic, tab_reg = tab_reg,
             tab_cola = tab_cola, p_turb = p_turb, q99 = q99), "resultados/EM_final.rds")
cat("\nGuardado: resultados/EM_final.rds y figuras/08_prob_turbulencia.png\n")
