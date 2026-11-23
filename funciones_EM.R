# =====================================================================
# funciones_EM.R - Funciones del EM (solo definiciones, sin salida)
# Son las mismas funciones de 08_EM.R (esqueleto de la docente completado),
# extraídas para reutilizarlas con source() en los demás scripts de la Fase C.
# theta = list(pi = prob. del componente 1, mu = c(mu1, mu2), sigma = c(s1, s2))
# =====================================================================
dnorm_propia <- function(x, mu, sigma) exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))

paso_E <- function(x, theta) {
  a <- theta$pi * dnorm_propia(x, theta$mu[1], theta$sigma[1])
  b <- (1 - theta$pi) * dnorm_propia(x, theta$mu[2], theta$sigma[2])
  a / (a + b)
}

paso_M <- function(x, gamma) {
  w1 <- sum(gamma); w2 <- sum(1 - gamma)
  mu1 <- sum(gamma * x) / w1
  mu2 <- sum((1 - gamma) * x) / w2
  s1  <- sqrt(sum(gamma * (x - mu1)^2) / w1)
  s2  <- sqrt(sum((1 - gamma) * (x - mu2)^2) / w2)
  list(pi = w1 / length(x), mu = c(mu1, mu2), sigma = c(s1, s2))
}

loglik <- function(x, theta) {
  sum(log(theta$pi * dnorm_propia(x, theta$mu[1], theta$sigma[1]) +
          (1 - theta$pi) * dnorm_propia(x, theta$mu[2], theta$sigma[2])))
}

em_mezcla <- function(x, theta0, tol = 1e-8, max_iter = 2000) {
  theta <- theta0
  ll <- loglik(x, theta)
  for (it in 1:max_iter) {
    gamma <- paso_E(x, theta)
    theta <- paso_M(x, gamma)
    ll_nuevo <- loglik(x, theta)
    if (!is.finite(ll_nuevo))
      stop("Log-verosimilitud no finita en la iteración ", it, " (posible colapso de una desviación).")
    ll <- c(ll, ll_nuevo)
    if (abs(ll_nuevo - ll[length(ll) - 1]) < tol) break
  }
  list(theta = theta, ll = ll, iter = it)
}
