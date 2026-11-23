# =====================================================================
# FASE C (parte 1) - Algoritmo EM para una mezcla de dos normales
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Parte del esqueleto de la docente (esqueleto_EM.R): se completan los bloques TODO,
# sin cambiar los nombres de las funciones ni la estructura de theta = list(pi, mu, sigma).
# =====================================================================
set.seed(2026)

# ---- 1. Utilidades (dadas por la docente) --------------------------------------------
dnorm_propia <- function(x, mu, sigma) exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))
rnorm_bm <- function(n) {
  u1 <- runif(n); u2 <- runif(n)
  sqrt(-2 * log(u1)) * cos(2 * pi * u2)
}
simular_mezcla <- function(n, pi1, mu, sigma) {
  comp <- ifelse(runif(n) < pi1, 1, 2)
  mu[comp] + sigma[comp] * rnorm_bm(n)
}

# ---- 2. Paso E -----------------------------------------------------------------------------
# gamma_i = P(componente 1 | x_i, theta) = pi f1(x_i) / ( pi f1(x_i) + (1 - pi) f2(x_i) )
paso_E <- function(x, theta) {
  a <- theta$pi * dnorm_propia(x, theta$mu[1], theta$sigma[1])         # aporte del componente 1
  b <- (1 - theta$pi) * dnorm_propia(x, theta$mu[2], theta$sigma[2])   # aporte del componente 2
  a / (a + b)
}

# ---- 3. Paso M -----------------------------------------------------------------------------
# Estimadores de máxima verosimilitud ponderados: pesos gamma (comp. 1) y 1 - gamma (comp. 2)
paso_M <- function(x, gamma) {
  w1 <- sum(gamma)            # "número efectivo" de datos del componente 1
  w2 <- sum(1 - gamma)        # idem, componente 2
  mu1 <- sum(gamma * x) / w1
  mu2 <- sum((1 - gamma) * x) / w2
  s1  <- sqrt(sum(gamma * (x - mu1)^2) / w1)
  s2  <- sqrt(sum((1 - gamma) * (x - mu2)^2) / w2)
  list(pi = w1 / length(x), mu = c(mu1, mu2), sigma = c(s1, s2))
}

# ---- 4. Log-verosimilitud de la mezcla --------------------------------------------------------
loglik <- function(x, theta) {
  sum(log(theta$pi * dnorm_propia(x, theta$mu[1], theta$sigma[1]) +
          (1 - theta$pi) * dnorm_propia(x, theta$mu[2], theta$sigma[2])))
}

# ---- 5. Bucle EM ------------------------------------------------------------------------------------
# Criterio de parada: cambio ABSOLUTO de la log-verosimilitud entre dos iteraciones < tol.
em_mezcla <- function(x, theta0, tol = 1e-8, max_iter = 2000) {
  theta <- theta0
  ll <- loglik(x, theta)                        # log-verosimilitud inicial
  for (it in 1:max_iter) {
    gamma <- paso_E(x, theta)                   # paso E con los parámetros actuales
    theta <- paso_M(x, gamma)                   # paso M con esas responsabilidades
    ll_nuevo <- loglik(x, theta)
    if (!is.finite(ll_nuevo))
      stop("Log-verosimilitud no finita en la iteración ", it, " (posible colapso de una desviación).")
    ll <- c(ll, ll_nuevo)                       # se registra la log-verosimilitud de cada iteración
    if (abs(ll_nuevo - ll[length(ll) - 1]) < tol) break
  }
  list(theta = theta, ll = ll, iter = it)
}

# ---- 6. Prueba con datos simulados de parámetros conocidos (de la docente) ------------------------
verificar_em <- function() {
  verdad <- list(pi = 0.7, mu = c(0, 0.1), sigma = c(0.6, 1.5))
  x <- simular_mezcla(5000, verdad$pi, verdad$mu, verdad$sigma)
  ini <- list(pi = 0.5, mu = c(-0.5, 0.5), sigma = c(0.5 * sd(x), 2 * sd(x)))
  r <- em_mezcla(x, ini)
  monotona <- all(diff(r$ll) >= -1e-9)
  cerca <- abs(r$theta$pi - 0.7) < 0.08 && all(abs(r$theta$sigma - c(0.6, 1.5)) < 0.2)
  cat("log-verosimilitud no decreciente:", if (monotona) "OK" else "FALLA", "\n")
  cat("parámetros cercanos a los verdaderos:", if (cerca) "OK" else "REVISAR", "\n")
  cat(sprintf("pi=%.3f  sigma=(%.3f, %.3f)  iteraciones=%d\n",
              r$theta$pi, r$theta$sigma[1], r$theta$sigma[2], r$iter))
  invisible(r)
}
cat("== Verificación con datos simulados (docente) ==\n")
verificar_em()

# ---- 7. EM sobre la serie real (retornos CON signo, incluidos los ceros) -------------------------
archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct; n <- length(r)

# Inicio de esta corrida base (los 10 inicios aleatorios vienen en el siguiente script):
# componentes con desviaciones distintas (si fueran idénticas, el EM no podría separarlos)
ini <- list(pi = 0.7, mu = c(mean(r), mean(r)), sigma = c(0.5 * sd(r), 2 * sd(r)))
res <- em_mezcla(r, ini)

cat("\n== EM sobre GBP/USD (tol = 1e-8) ==\n")
cat("iteraciones:", res$iter, "\n")
cat("log-verosimilitud inicial:", res$ll[1], " final:", res$ll[length(res$ll)], "\n")
cat("log-verosimilitud no decreciente (tolerancia 1e-9):", all(diff(res$ll) >= -1e-9), "\n")
cat("incremento mínimo entre iteraciones:", min(diff(res$ll)), "\n")
cat(sprintf("pi = %.4f | mu = (%.4f, %.4f) | sigma = (%.4f, %.4f)\n",
            res$theta$pi, res$theta$mu[1], res$theta$mu[2], res$theta$sigma[1], res$theta$sigma[2]))
cat("primeras 6 log-verosimilitudes:", head(res$ll, 6), "\n")

# ---- 8. Gráfico de convergencia ---------------------------------------------------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/06_convergencia_EM.png", width = 1400, height = 560, res = 130)
par(mfrow = c(1, 2))
plot(0:(length(res$ll) - 1), res$ll, type = "l", lwd = 2,
     xlab = "iteración", ylab = "log-verosimilitud", main = "Convergencia del EM")
inc <- pmax(diff(res$ll), 1e-16)
plot(1:length(inc), inc, type = "l", log = "y", lwd = 2,
     xlab = "iteración", ylab = "incremento de la log-verosimilitud (escala log)",
     main = "Incrementos por iteración")
dev.off()

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(res = res, ini = ini), "resultados/EM_base.rds")
cat("\nGuardado: resultados/EM_base.rds y figuras/06_convergencia_EM.png\n")
