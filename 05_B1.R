# =====================================================================
# B1 - Monte Carlo y muestreo por importancia: P(pérdida acumulada a 5 días > u)
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_A.R y resultados/params_A.rds
# Modelo (Fase A): pérdida diaria L_j = s_j * X_j, con X_j ~ Weibull(k, lam) iid
#   y s_j = +1 (pérdida) o -1 (ganancia) con prob 1/2, independiente de X_j.
#   Evento: S5 = L_1 + ... + L_5 > u
# =====================================================================
set.seed(2026)
source("funciones_A.R")
p <- readRDS("resultados/params_A.rds")
k <- p$k_w; lam <- p$lam_w

# ---- 1. Simulación del modelo (MC simple) en bloques para no saturar la memoria ----
sim_S5 <- function(n) {
  X  <- matrix(rweibull_ti(5 * n, k, lam), n, 5)
  Sg <- matrix(signos(5 * n), n, 5)
  rowSums(X * Sg)
}

# ---- 2. Simulación de referencia y elección del umbral u ----------------------------
bloque <- 2e5
Sref <- unlist(lapply(1:20, function(i) sim_S5(bloque)))     # 4.000.000 simulaciones
u0 <- quantile(Sref, 0.999, type = 7, names = FALSE)
u  <- round(u0, 2)                                           # umbral redondeado
cat("Cuantil 99.9% de S5 en la simulación de referencia:", u0, "\n")
cat("Umbral elegido u =", u, "\n")
rm(Sref)

# Referencia independiente (otra simulación) para VALIDAR la probabilidad en u
N_ref <- 4e6
conteo <- 0
for (i in 1:(N_ref / bloque)) conteo <- conteo + sum(sim_S5(bloque) > u)
p_ref  <- conteo / N_ref
se_ref <- sqrt(p_ref * (1 - p_ref) / N_ref)
cat(sprintf("Probabilidad de referencia p_ref = %.6f (EE = %.6f), N_ref = %d\n",
            p_ref, se_ref, N_ref))

# ---- 3. Muestreo por importancia -----------------------------------------------------
# Propuesta g: X_j ~ Weibull(k, c*lam) (misma forma, escala mayor) y
#              P(día de pérdida) = q en lugar de 1/2.
# Peso por observación: W = prod_j [ f(x_j)/g(x_j) ] * (0.5/q)^{#pérdidas} * (0.5/(1-q))^{#ganancias}
# con f(x)/g(x) = c^k * exp( -(x/lam)^k + (x/(c*lam))^k )  (se trabaja en logaritmos)
sim_IS <- function(n, c, q, u) {
  lp <- c * lam
  X    <- matrix(rweibull_ti(5 * n, k, lp), n, 5)
  perd <- matrix(runif(5 * n) < q, n, 5)                  # TRUE: día de pérdida
  S    <- rowSums(X * ifelse(perd, 1, -1))
  logwx <- rowSums(k * log(c) - (X / lam)^k + (X / lp)^k)
  nl    <- rowSums(perd)
  logws <- nl * log(0.5 / q) + (5 - nl) * log(0.5 / (1 - q))
  logw  <- logwx + logws
  y <- numeric(n)
  idx <- which(S > u)
  y[idx] <- exp(logw[idx])                                # indicador * peso
  list(y = y, logw = logw)
}

# ---- 4. Estudio piloto: cómo elegir la propuesta (y qué pasa si se elige mal) ----------
cat("\n== Piloto (n = 2e5 por combinación) ==\n")
cs <- c(0.7, 1, 1.5, 2, 3, 4, 8); qs <- c(0.5, 0.7, 0.85, 0.95)
piloto <- expand.grid(c = cs, q = qs)
piloto$estimacion <- NA; piloto$razon_var <- NA
for (i in 1:nrow(piloto)) {
  y <- sim_IS(2e5, piloto$c[i], piloto$q[i], u)$y
  piloto$estimacion[i] <- mean(y)
  v <- var(y)
  piloto$razon_var[i] <- if (v > 0) p_ref * (1 - p_ref) / v else NA   # Var_MC / Var_IS (por simulación)
}
print(piloto, digits = 4, row.names = FALSE)
cat("razon_var > 1: la propuesta reduce la varianza; < 1: la AUMENTA; NA: nunca alcanza el evento.\n")

# ---- 5. Estimaciones oficiales (misma N para ambos métodos) ----------------------------
N <- 1e5
c_is <- 2; q_is <- 0.85                                   # propuesta elegida (a partir del piloto)

y_mc <- as.numeric(sim_S5(N) > u)
p_mc <- mean(y_mc); se_mc <- sd(y_mc) / sqrt(N)

sis  <- sim_IS(N, c_is, q_is, u)
p_is <- mean(sis$y); se_is <- sd(sis$y) / sqrt(N)

cat("\n== Resultados B1 (N =", N, ") ==\n")
res <- data.frame(
  metodo = c("Monte Carlo simple", "Muestreo por importancia", "Referencia (4e6)"),
  estimacion = c(p_mc, p_is, p_ref),
  error_estandar = c(se_mc, se_is, se_ref)
)
print(res, digits = 5, row.names = FALSE)
cat("\nVarianza MC (empírica):", se_mc^2, " | Varianza IS (empírica):", se_is^2, "\n")
cat("Reducción de varianza (empírica, Var_MC/Var_IS):", se_mc^2 / se_is^2, "\n")
cat("Reducción de varianza (usando p_ref, más estable):",
    p_ref * (1 - p_ref) / var(sis$y), "\n")
cat("Tamaño de MC simple equivalente a la N del IS:", N * p_ref * (1 - p_ref) / var(sis$y), "\n")

# Validación frente a la referencia: diferencia en errores estándar
z_mc <- (p_mc - p_ref) / sqrt(se_mc^2 + se_ref^2)
z_is <- (p_is - p_ref) / sqrt(se_is^2 + se_ref^2)
cat("Diferencia con la referencia (en EE): MC =", z_mc, " IS =", z_is, "\n")
# Comprobación de los pesos: E_g[W] = 1
cat("Media de los pesos (debe ser ~1):", mean(exp(sis$logw)), "\n")

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(u = u, p_ref = p_ref, se_ref = se_ref, res = res, piloto = piloto,
             c_is = c_is, q_is = q_is, N = N), "resultados/B1.rds")





