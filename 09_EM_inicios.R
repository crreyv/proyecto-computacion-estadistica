# =====================================================================
# FASE C (parte 2) - Sensibilidad del EM: inicios aleatorios, tolerancias y casos degenerados
# Proyecto final - Computación Estadística 2026-II - GRUPO 2 - GBP/USD
# Requiere: funciones_EM.R
# =====================================================================
set.seed(2026)
source("funciones_EM.R")

archivo <- "retornos_diarios_12_series.csv"
candidatos <- c(file.path("data", archivo), archivo)
ruta <- candidatos[file.exists(candidatos)][1]
if (is.na(ruta)) stop("No encuentro el CSV. Carpeta actual: ", getwd())
datos <- read.csv(ruta, stringsAsFactors = FALSE, encoding = "UTF-8")
gbp <- datos[datos$serie == "GBP/USD", ]; gbp <- gbp[order(gbp$t), ]
r <- gbp$retorno_pct; n <- length(r); s <- sd(r)

# ---- Versión "diagnóstica" del bucle EM: no se detiene con error, informa el estado ----------
# (la lógica de los pasos E y M es la misma de funciones_EM.R)
em_diag <- function(x, theta0, tol = 1e-8, max_iter = 2000) {
  theta <- theta0; ll <- loglik(x, theta); estado <- "no convergió (max_iter)"
  for (it in 1:max_iter) {
    gamma <- paso_E(x, theta)
    theta <- paso_M(x, gamma)
    ll_nuevo <- loglik(x, theta)
    if (!is.finite(ll_nuevo) || any(!is.finite(unlist(theta)))) { estado <- "colapso (valor no finito)"; break }
    ll <- c(ll, ll_nuevo)
    if (abs(ll_nuevo - ll[length(ll) - 1]) < tol) { estado <- "convergió"; break }
  }
  list(theta = theta, ll = ll, iter = it, estado = estado)
}
# Etiquetado: el componente 1 es el de menor desviación (calma); si se intercambian, pi -> 1 - pi
ordenar <- function(theta) {
  if (all(is.finite(unlist(theta))) && theta$sigma[1] > theta$sigma[2])
    list(pi = 1 - theta$pi, mu = rev(theta$mu), sigma = rev(theta$sigma)) else theta
}
fila <- function(nombre, ini, out) {
  th <- ordenar(out$theta)
  data.frame(inicio = nombre, pi0 = ini$pi, sigma0_1 = ini$sigma[1], sigma0_2 = ini$sigma[2],
             estado = out$estado, iter = out$iter, ll_final = out$ll[length(out$ll)],
             pi_calma = th$pi, sigma_calma = th$sigma[1], sigma_turb = th$sigma[2],
             stringsAsFactors = FALSE)
}

# ---- 1. Doce inicios aleatorios ------------------------------------------------------------------
n_ini <- 12
inicios <- lapply(1:n_ini, function(k)
  list(pi = runif(1, 0.1, 0.9), mu = runif(2, -s, s), sigma = runif(2, 0.1 * s, 2 * s)))
salidas <- lapply(inicios, function(ini) em_diag(r, ini))
tab_ini <- do.call(rbind, lapply(1:n_ini, function(k) fila(paste0("I", k), inicios[[k]], salidas[[k]])))
cat("== 12 inicios aleatorios (tol = 1e-8); componente 1 = el de menor desviación ==\n")
print(tab_ini, digits = 5, row.names = FALSE)
ok <- tab_ini$estado == "convergió"
mejor <- max(tab_ini$ll_final[ok])
cat("\nInicios que convergieron:", sum(ok), "de", n_ini, "\n")
cat("Inicios con log-verosimilitud final a menos de 1e-3 de la mejor:",
    sum(ok & abs(tab_ini$ll_final - mejor) < 1e-3), "\n")
cat("Rango de pi_calma entre los convergidos:", range(tab_ini$pi_calma[ok]), "\n")
cat("Mejor log-verosimilitud:", mejor, "\n")

# ---- 2. Tolerancias de parada ---------------------------------------------------------------------
ini_base <- list(pi = 0.7, mu = c(mean(r), mean(r)), sigma = c(0.5 * s, 2 * s))
tols <- c(1e-4, 1e-6, 1e-8, 1e-10)
tab_tol <- do.call(rbind, lapply(tols, function(tl) {
  out <- em_diag(r, ini_base, tol = tl); th <- ordenar(out$theta)
  data.frame(tol = tl, iter = out$iter, ll_final = out$ll[length(out$ll)],
             pi_calma = th$pi, sigma_calma = th$sigma[1], sigma_turb = th$sigma[2])
}))
cat("\n== Efecto de la tolerancia de parada (mismo inicio) ==\n")
print(tab_tol, digits = 8, row.names = FALSE)

# ---- 3. Casos degenerados --------------------------------------------------------------------------
i_max <- which.max(abs(r))
deg <- list(
  "componentes identicos"          = list(pi = 0.5,  mu = c(0, 0),     sigma = c(s, s)),
  "estrecho en 0, sigma=0.03"      = list(pi = 0.05, mu = c(0, 0),     sigma = c(0.03, s)),
  "estrecho en 0, sigma=0.01"      = list(pi = 0.05, mu = c(0, 0),     sigma = c(0.01, s)),
  "estrecho en max|r|, sigma=0.01" = list(pi = 0.01, mu = c(r[i_max], 0), sigma = c(0.01, s)))
tab_deg <- do.call(rbind, lapply(names(deg), function(nm) {
  out <- em_diag(r, deg[[nm]], max_iter = 500)
  data.frame(inicio = nm, estado = out$estado, iter = out$iter,
             ll_ultima_finita = out$ll[length(out$ll)],
             pi_final = out$theta$pi, sigma_min_final = min(out$theta$sigma),
             stringsAsFactors = FALSE)
}))
cat("\n== Casos degenerados ==\n")
print(tab_deg, digits = 6, row.names = FALSE)
cat("Proporción de ceros exactos en los datos:", mean(r == 0), "\n")
cat("Log-verosimilitud de una sola normal (referencia):",
    sum(log(dnorm_propia(r, mean(r), sqrt(mean((r - mean(r))^2))))), "\n")

# ---- 4. Gráfico: trayectorias de log-verosimilitud de los 12 inicios -----------------------------------
dir.create("figuras", showWarnings = FALSE)
png("figuras/07_inicios_EM.png", width = 1400, height = 560, res = 130)
par(mfrow = c(1, 2))
mx <- max(sapply(salidas, function(o) length(o$ll)))
plot(NA, xlim = c(0, mx), ylim = c(-2140, mejor), xlab = "iteración",
     ylab = "log-verosimilitud", main = "12 inicios aleatorios")
for (o in salidas) lines(0:(length(o$ll) - 1), o$ll, col = "grey40")
plot(NA, xlim = c(0, mx), ylim = c(mejor - 0.05, mejor + 0.005), xlab = "iteración",
     ylab = "log-verosimilitud", main = "Zoom cerca del máximo")
for (o in salidas) lines(0:(length(o$ll) - 1), o$ll, col = "grey40")
abline(h = mejor, lty = 2, col = "red")
dev.off()

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(tab_ini = tab_ini, tab_tol = tab_tol, tab_deg = tab_deg), "resultados/EM_inicios.rds")
cat("\nGuardado: resultados/EM_inicios.rds y figuras/07_inicios_EM.png\n")
