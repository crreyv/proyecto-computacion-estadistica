# =====================================================================
# FASE D (parte 1) - VaR por modelo y autocorrelación (dependencia temporal)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_EM.R, resultados/params_A.rds, resultados/B2.rds, resultados/EM_final.rds
# =====================================================================
set.seed(2026)
source("funciones_EM.R")
pA <- readRDS("resultados/params_A.rds")
b2 <- readRDS("resultados/B2.rds")
em <- readRDS("resultados/EM_final.rds")
th <- em$theta; p_turb <- em$p_turb

archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct; n <- length(r); L <- -r          # L = pérdida diaria

# ---- 0. Normal estándar propia (cola y cuantil) -------------------------------------------------------
erf_propia <- function(z) {
  k <- 0:100; fact <- c(1, cumprod(1:100))
  sapply(z, function(zz) 2 / sqrt(pi) * sum((-1)^k * zz^(2 * k + 1) / (fact * (2 * k + 1))))
}
cola_normal <- function(z) {                    # 1 - Phi(z)
  sapply(z, function(zz) {
    if (zz < 3) return(0.5 * (1 - erf_propia(zz / sqrt(2))))
    t <- zz
    for (j in 60:1) t <- zz + j / t             # fracción continua del cociente de Mills
    exp(-0.5 * zz^2) / sqrt(2 * pi) / t
  })
}
bisec <- function(f, lo, hi) {                  # raíz de f (decreciente) por bisección
  for (i in 1:100) { mid <- (lo + hi) / 2; if (f(mid) > 0) lo <- mid else hi <- mid }
  (lo + hi) / 2
}
qnorm_propia <- function(p) bisec(function(z) cola_normal(z) - (1 - p), -6, 6)   # Q(z) decrece: Q(z) > 1-p -> subir z
cat("z_0.99 propio:", qnorm_propia(0.99), "| verificación con qnorm:", qnorm(0.99), "\n")

# ---- 1. VaR por modelo ----------------------------------------------------------------------------------
niveles <- c(0.95, 0.99, 0.995, 0.999)
mu0 <- mean(r); sg0 <- sqrt(mean((r - mu0)^2)); p0 <- mean(r == 0)
var_emp <- sapply(niveles, function(p) quantile(L, p, type = 7, names = FALSE))
var_nor <- sapply(niveles, function(p) -mu0 + qnorm_propia(p) * sg0)         # L ~ N(-mu0, sg0^2)
cola_perdida_mezcla <- function(q)                                           # P(L > q) bajo la mezcla
  th$pi * cola_normal((q + th$mu[1]) / th$sigma[1]) + (1 - th$pi) * cola_normal((q + th$mu[2]) / th$sigma[2])
var_mez <- sapply(niveles, function(p) bisec(function(q) cola_perdida_mezcla(q) - (1 - p), 0, 10))
var_wei <- pA$lam_w * (-log(2 * (1 - niveles)))^(1 / pA$k_w)                 # Weibull con signos +-1/2 (sin ceros)
var_wei0 <- pA$lam_w * (-log((1 - niveles) / (0.5 * (1 - p0))))^(1 / pA$k_w) # idem, con la masa de ceros
tab_var <- data.frame(nivel = niveles, empirico = var_emp, normal = var_nor, mezcla = var_mez,
                      weibull_signos = var_wei, weibull_signos_ceros = var_wei0)
cat("\n== VaR diario (cuantil de la pérdida -r) por modelo ==\n")
print(tab_var, digits = 5, row.names = FALSE)
dif <- data.frame(nivel = niveles,
                  normal = 100 * (var_nor / var_emp - 1), mezcla = 100 * (var_mez / var_emp - 1),
                  weibull_signos = 100 * (var_wei / var_emp - 1))
cat("\nDiferencia porcentual frente al VaR empírico:\n"); print(dif, digits = 3, row.names = FALSE)
cat("\nIC 95% percentil del VaR 99% empírico (B2):", unlist(b2$tabla[2, c("IC95_percentil_inf", "IC95_percentil_sup")]), "\n")
cat("Pérdidas observadas por encima del nivel (n*(1-p)):", n * (1 - niveles), "\n")

# ---- 2. Autocorrelación -----------------------------------------------------------------------------------
acf_propia <- function(v, h) {
  v <- v - mean(v)
  sum(v[1:(length(v) - h)] * v[(1 + h):length(v)]) / sum(v^2)
}
lags <- 1:30
series <- list("r" = r, "|r|" = abs(r), "r^2" = r^2, "P(turbulencia)" = p_turb)
tab_acf <- data.frame(lag = lags,
                      sapply(series, function(v) sapply(lags, function(h) acf_propia(v, h))),
                      check.names = FALSE)
banda <- 1.96 / sqrt(n)
cat("\n== Autocorrelación muestral (banda de ruido blanco +-", banda, ") ==\n")
print(tab_acf[c(1, 2, 3, 5, 10, 20, 30), ], digits = 3, row.names = FALSE)
cat("Rezagos (de 30) fuera de la banda:\n")
print(sapply(tab_acf[, -1], function(a) sum(abs(a) > banda)))

# Efecto de diseño aproximado (ventana de Bartlett, K = 20): n_eff = n / (1 + 2 sum (1 - h/(K+1)) rho_h)
efecto <- function(v, K = 20)
  1 + 2 * sum((1 - (1:K) / (K + 1)) * sapply(1:K, function(h) acf_propia(v, h)))
f_abs <- efecto(abs(r)); f_r2 <- efecto(r^2)
cat("\nEfecto de diseño aproximado: |r| =", f_abs, "(n_eff =", n / f_abs, ") | r^2 =", f_r2, "(n_eff =", n / f_r2, ")\n")

# ---- 3. Gráfico de la autocorrelación -----------------------------------------------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/09_acf.png", width = 1400, height = 900, res = 130)
par(mfrow = c(2, 2))
for (nm in names(series)) {
  a <- tab_acf[[nm]]
  plot(lags, a, type = "h", lwd = 3, ylim = range(c(a, -banda, banda)),
       xlab = "rezago (días)", ylab = "autocorrelación", main = nm)
  abline(h = c(-banda, banda), lty = 2, col = "red"); abline(h = 0)
}
dev.off()

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(tab_var = tab_var, dif = dif, tab_acf = tab_acf, banda = banda,
             efecto = c(abs = f_abs, r2 = f_r2)), "resultados/FaseD.rds")
cat("\nGuardado: resultados/FaseD.rds y figuras/09_acf.png\n")
