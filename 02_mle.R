# =====================================================================
# TAREA 2 - Máxima verosimilitud: lognormal (forma cerrada) y Weibull (Newton)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Ejecutar desde la misma carpeta que 01_datos.R (ver instrucciones previas)
# =====================================================================
set.seed(2026)

# ---- 0. Datos (mismo procedimiento de la tarea 1) --------------------
archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]
gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct
x <- abs(r[r != 0])          # variable de trabajo
m <- length(x)
lx <- log(x)                 # logaritmos: se usan en ambas familias
cat("m =", m, "\n")

# ---- 1. LOGNORMAL: solución cerrada ----------------------------------
# ln X ~ N(mu, sigma^2)  =>  mu_hat = media(ln x),  sigma2_hat = media((ln x - mu_hat)^2)
mu_ln    <- mean(lx)
sigma_ln <- sqrt(mean((lx - mu_ln)^2))        # divisor m (MLE), no m-1

loglik_lognormal <- function(mu, sigma, x) {
  l <- log(x)
  sum(-l - log(sigma) - 0.5 * log(2 * pi) - (l - mu)^2 / (2 * sigma^2))
}
ll_ln  <- loglik_lognormal(mu_ln, sigma_ln, x)
aic_ln <- 2 * 2 - 2 * ll_ln                    # 2 parámetros

# ---- 2. WEIBULL: Newton para la forma k -------------------------------
# Ecuación de verosimilitud (perfil): g(k) = S1/S0 - 1/k - mean(ln x) = 0
#   S0 = sum(x^k),  S1 = sum(x^k ln x),  S2 = sum(x^k (ln x)^2)
# g'(k) = (S2*S0 - S1^2)/S0^2 + 1/k^2
g_weibull <- function(k, x) {
  l  <- log(x)
  xk <- x^k
  S0 <- sum(xk); S1 <- sum(xk * l); S2 <- sum(xk * l^2)
  c(g  = S1 / S0 - 1 / k - mean(l),
    gp = (S2 * S0 - S1^2) / S0^2 + 1 / k^2)
}

newton_weibull <- function(x, k0, tol = 1e-12, max_iter = 100) {
  k <- k0
  hist <- data.frame(iter = 0, k = k0, g = unname(g_weibull(k0, x)["g"]))
  for (it in 1:max_iter) {
    v  <- g_weibull(k, x)
    kn <- unname(k - v["g"] / v["gp"])      # paso de Newton (unname: evita nombres heredados)
    if (!is.finite(kn) || kn <= 0) kn <- k / 2   # protección: k debe ser > 0
    hist <- rbind(hist, data.frame(iter = it, k = unname(kn),
                                   g = unname(g_weibull(kn, x)["g"])))
    if (abs(kn - k) < tol) { k <- kn; break }
    k <- kn
  }
  list(k = unname(k), iter = it, hist = hist)
}

# Valor inicial: Var(ln X) = pi^2/(6 k^2) para la Weibull => k ~ 1.2 / sd(ln x)
# (aproximación mía para arrancar; no es parte de la guía)
k0 <- 1.2 / sqrt(mean((lx - mean(lx))^2))
nw <- newton_weibull(x, k0)
k_w   <- nw$k
lam_w <- (mean(x^k_w))^(1 / k_w)             # lambda_hat^k = media(x^k)

loglik_weibull <- function(k, lam, x) {
  m <- length(x)
  m * log(k) - m * k * log(lam) + (k - 1) * sum(log(x)) - sum((x / lam)^k)
}
ll_w  <- loglik_weibull(k_w, lam_w, x)
aic_w <- 2 * 2 - 2 * ll_w

cat("\n== Historial de Newton (Weibull) ==\n")
print(nw$hist, row.names = FALSE, digits = 10)

# ---- 3. Verificaciones propias (sin librerías) -------------------------
cat("\n== Verificaciones propias ==\n")
cat("g(k_hat) (debe ser ~0):", g_weibull(k_w, x)["g"], "\n")
# El máximo del perfil debe estar en k_hat: comparamos con puntos vecinos
perfil <- function(k) {
  lam <- (mean(x^k))^(1 / k)
  loglik_weibull(k, lam, x)
}
cat("perfil en k_hat-0.05, k_hat, k_hat+0.05:",
    perfil(k_w - 0.05), perfil(k_w), perfil(k_w + 0.05), "\n")
# Robustez al valor inicial
for (k_ini in c(0.3, 1, 3)) {
  cat("inicio", k_ini, "-> k =", newton_weibull(x, k_ini)$k, "\n")
}

# ---- 4. Tabla de resultados ---------------------------------------------
tabla <- data.frame(
  Familia   = c("Lognormal", "Weibull"),
  Parametro1 = c(mu_ln, k_w),        # lognormal: mu ; Weibull: forma k
  Parametro2 = c(sigma_ln, lam_w),   # lognormal: sigma ; Weibull: escala lambda
  logLik    = c(ll_ln, ll_w),
  AIC       = c(aic_ln, aic_w)
)
cat("\n== Resultados (oficiales: implementación propia) ==\n")
print(tabla, row.names = FALSE, digits = 8)
cat("Lognormal: (mu, sigma) | Weibull: (k forma, lambda escala)\n")
cat("Delta AIC (lognormal - Weibull):", aic_ln - aic_w, "\n")

# ---- 5. Comparación con librería (SOLO VERIFICACIÓN, no produce resultados) ----
cat("\n== Verificación con MASS::fitdistr (NO es el resultado oficial) ==\n")
res <- tryCatch({
  f_w  <- MASS::fitdistr(x, "weibull")
  f_ln <- MASS::fitdistr(x, "lognormal")
  print(rbind(
    Weibull_propio    = c(forma = k_w, escala = lam_w),
    Weibull_fitdistr  = unname(f_w$estimate)
  ))
  print(rbind(
    Lognormal_propio   = c(meanlog = mu_ln, sdlog = sigma_ln),
    Lognormal_fitdistr = unname(f_ln$estimate)
  ))
}, error = function(e) cat("MASS no disponible:", conditionMessage(e), "\n"))

# ---- 6. Guardar resultados para las tareas siguientes ----------------------
dir.create("resultados", showWarnings = FALSE)
saveRDS(list(mu_ln = mu_ln, sigma_ln = sigma_ln, k_w = k_w, lam_w = lam_w,
             ll_ln = ll_ln, ll_w = ll_w, aic_ln = aic_ln, aic_w = aic_w),
        "resultados/params_A.rds")
cat("\nGuardado: resultados/params_A.rds\n")


