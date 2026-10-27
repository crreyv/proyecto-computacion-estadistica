# =====================================================================
# funciones_A.R - Generadores de la Fase A (SOLO definiciones, sin salida)
# Solo se usa runif() como generador base y funciones elementales.
# =====================================================================

# Weibull(k forma, lam escala) por TRANSFORMADA INVERSA
# F(x) = 1 - exp(-(x/lam)^k)  =>  F^{-1}(u) = lam * (-log(1 - u))^(1/k)
# Como 1 - U tiene la misma distribución que U, se usa -log(U).
rweibull_ti <- function(n, k, lam) {
  u <- runif(n)
  lam * (-log(u))^(1 / k)
}

# Normal estándar por BOX-MULLER (a partir de dos uniformes independientes)
rnorm_bm <- function(n) {
  u1 <- runif(n); u2 <- runif(n)
  sqrt(-2 * log(u1)) * cos(2 * pi * u2)
}

# Lognormal(mu, sigma): X = exp(mu + sigma * Z), Z ~ N(0,1)
rlnorm_bm <- function(n, mu, sigma) {
  exp(mu + sigma * rnorm_bm(n))
}

# Función de distribución acumulada de la Weibull (para verificaciones)
pweibull_propia <- function(x, k, lam) 1 - exp(-(x / lam)^k)

# Signos +/-1 con probabilidad 1/2 cada uno (supuesto del modelo, Fase A.4)
signos <- function(n) ifelse(runif(n) < 0.5, 1, -1)