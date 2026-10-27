# =====================================================================
# TAREA 3 - Generadores de variables aleatorias y comprobaciones básicas
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_A.R y resultados/params_A.rds (de 02_mle.R)
# =====================================================================
set.seed(2026)
source("funciones_A.R")
p <- readRDS("resultados/params_A.rds")
k <- p$k_w; lam <- p$lam_w; mu <- p$mu_ln; sg <- p$sigma_ln

# Tamaño de la muestra sintética: unas 20 veces el de la muestra real (m = 1791)
m  <- 1791
N  <- 20 * m
cat("N =", N, "\n")

# ---- 1. Normal estándar por Box-Muller --------------------------------
z <- rnorm_bm(N)
cat("\n== Box-Muller ==\n")
cat("media (teórica 0):", mean(z), "\n")
cat("varianza (teórica 1):", mean((z - mean(z))^2), "\n")
cat("curtosis en exceso (teórica 0):",
    mean((z - mean(z))^4) / mean((z - mean(z))^2)^2 - 3, "\n")

# ---- 2. Weibull por transformada inversa ------------------------------
xw <- rweibull_ti(N, k, lam)
G  <- function(a) exp(lgamma(a))                       # Gamma(a) vía lgamma
mean_w <- lam * G(1 + 1 / k)
sd_w   <- lam * sqrt(G(1 + 2 / k) - G(1 + 1 / k)^2)
q_w    <- function(pp) lam * (-log(1 - pp))^(1 / k)   # cuantiles teóricos
cat("\n== Weibull (transformada inversa) ==\n")
cat("media simulada / teórica:", mean(xw), "/", mean_w, "\n")
cat("desv.est. simulada / teórica:", sd(xw), "/", sd_w, "\n")
cat("cuantiles 50,90,99 simulados:", quantile(xw, c(.5, .9, .99)), "\n")
cat("cuantiles 50,90,99 teóricos :", q_w(c(.5, .9, .99)), "\n")
# Verificación clave: F(X) debe ser Uniforme(0,1)
uu <- pweibull_propia(xw, k, lam)
cat("F(X): media (0.5) y varianza (1/12 = 0.08333):", mean(uu), mean((uu - mean(uu))^2), "\n")

# ---- 3. Lognormal vía Box-Muller + exp(mu + sigma Z) -------------------
xl <- rlnorm_bm(N, mu, sg)
cat("\n== Lognormal (Box-Muller) ==\n")
cat("media simulada / teórica:", mean(xl), "/", exp(mu + sg^2 / 2), "\n")
cat("mediana simulada / teórica:", median(xl), "/", exp(mu), "\n")
cat("media de ln X (mu) y desv.est. de ln X (sigma):", mean(log(xl)), sd(log(xl)), "\n")

# ---- 4. Signos ---------------------------------------------------------
s <- signos(N)
cat("\n== Signos ==\n")
cat("proporción de +1 (teórica 0.5):", mean(s == 1), "\n")
cat("media de los retornos simulados con signo (Weibull):", mean(s * xw), "\n")


