# Proyecto final de Computación Estadística 2026-II
### Riesgo de pérdida en la serie de retornos GBP/USD

**Universidad Nacional de Colombia, sede de La Paz**
**Curso:** Computación Estadística, Grupo I (semestre 2026-II). **Docente:** Yoselin Nieto Gil
**Grupo de trabajo:** Grupo 2
**Integrantes:** 
- Adriana Jaimes
- Crisdailyn Rey Vargas
- John Rincón Toro
- Maria Vega Gutierrez

---

## 1. ¿De qué trata este proyecto?

Los mercados financieros tienen días tranquilos y días turbulentos, y las pérdidas grandes se concentran en pocos días. El proyecto responde una sola pregunta: ¿cuánto puede perderse en un día o en una semana al mantener un activo, y qué tan confiable es ese número?

Cada grupo trabaja con la serie de retornos diarios de un activo distinto y programa los métodos del curso desde cero y sin cajas negras: generación de variables aleatorias, integración por Monte Carlo, muestreo por importancia, remuestreo (jackknife y bootstrap) y modelos de mezcla con el algoritmo EM.

Nuestra serie es GBP/USD (tipo de cambio de la libra esterlina frente al dólar, según la denominación de la guía del proyecto):

| Característica | Valor |
|---|---|
| Observaciones | 1.866 |
| Período (según la guía) | 2 de enero de 1980 a 21 de mayo de 1987 |
| Variable | `retorno_pct`: retorno diario en porcentaje, $100\,(\ln P_t - \ln P_{t-1})$ |

El archivo de datos no trae fechas, solo el número de día `t`. El período se toma de la guía del proyecto.

## 2. Estado actual del trabajo

Este repositorio contiene el avance (primera entrega, 28 de octubre de 2026), que según la guía incluye:

| Parte | Contenido | Estado |
|---|---|---|
| Fase A | Simulación y generación de variables (completa) | Entregada en este avance |
| Fase B1 | Monte Carlo y muestreo por importancia | Entregada en este avance |
| Fase B2 | Bootstrap percentil | Entregada en este avance (preliminar) |
| Fase B2 | Jackknife, intervalo BCa e intervalos asintóticos | Pendiente (entrega final) |
| Fase C | Mezclas y algoritmo EM | Pendiente (entrega final) |
| Fase D | Interpretación y conclusiones | Pendiente (entrega final) |

## 3. Estructura del repositorio

```
proyecto-computacion-estadistica/
├── README.md                          Este archivo
├── run_all.R                          Ejecuta todo el análisis de principio a fin
├── funciones_A.R                      Generadores de variables aleatorias (solo definiciones)
├── 01_datos.R                         Carga de datos, ceros y variable de trabajo
├── 02_mle.R                           Máxima verosimilitud: lognormal y Weibull (Newton) y AIC
├── 03_generadores.R                   Comprobación de los generadores
├── 04_validacion.R                    QQ-plot, Kolmogorov-Smirnov a mano y cuantiles
├── 05_B1.R                            Monte Carlo y muestreo por importancia
├── 06_B2.R                            Bootstrap percentil
├── informe_avance.Rmd                 Informe en R Markdown (lee los resultados de los scripts)
├── data/
│   └── retornos_diarios_12_series.csv Datos entregados por la docente (12 series)
├── resultados/                        Archivos .rds generados por los scripts
└── figuras/                           Gráficos generados por los scripts
```

Los scripts buscan el archivo de datos en `data/` y, si no está allí, en la carpeta principal. Las carpetas `resultados/` y `figuras/` se crean solas al ejecutar el código.

## 4. Requisitos

- **R** (versión 4.x o superior recomendada) y, opcionalmente, **RStudio**.
- Para compilar el informe: los paquetes `rmarkdown` y `knitr`.
- Para el informe en PDF: una instalación de LaTeX (por ejemplo `tinytex::install_tinytex()`). El informe también puede generarse como HTML autónomo, sin LaTeX.
- El paquete `MASS` (incluido con las instalaciones habituales de R) se usa **únicamente para verificar** nuestros resultados con `fitdistr`. No produce ningún resultado oficial. Si no está disponible, el script lo avisa y sigue.

Todo lo demás se implementó con funciones base de R.

## 5. Cómo reproducir los resultados

La semilla aleatoria es fija: cada script ejecuta `set.seed(2026)` al inicio, y por eso los números, tablas y gráficos son los mismos en cada ejecución.

1. Abre R (o RStudio) y establece como carpeta de trabajo la carpeta del repositorio. En RStudio: *Session → Set Working Directory → Choose Directory*.
2. Verifica que el archivo `retornos_diarios_12_series.csv` esté en `data/` o en esa misma carpeta.
3. Ejecuta, en la consola de R:

   ```r
   source("run_all.R")
   ```

   Esto corre los seis scripts en orden. El script `05_B1.R` realiza simulaciones grandes y puede tardar un rato.
4. Compila el informe con *Knit* en RStudio, o con:

   ```r
   rmarkdown::render("informe_avance.Rmd")
   ```

   El informe lee los resultados guardados en `resultados/` y las figuras de `figuras/`, así que sus números coinciden con los del código. **El archivo `.Rmd` debe estar en la misma carpeta que `resultados/` y `figuras/`.**

También puedes ejecutar cada script por separado, en el orden numérico. Los scripts `03` a `06` dependen de los archivos que dejan los anteriores en `resultados/`.

Todas las rutas son **relativas**: no hay rutas locales del tipo `C:/Users/...`.

## 6. Resultados principales del avance (semilla 2026)

**Datos.** 75 de los 1.866 retornos valen exactamente cero (4,019 %). La variable de trabajo $|r|$ tiene 1.791 datos.

**Fase A: máxima verosimilitud y AIC.**

| Familia | Parámetros | Log-verosimilitud | AIC |
|---|---|---|---|
| Lognormal | $\hat\mu = -1{,}0035$, $\hat\sigma = 1{,}0411$ | −816,10 | 1636,21 |
| Weibull | $\hat k = 1{,}1367$, $\hat\lambda = 0{,}5997$ | −766,13 | 1536,26 |

Se elige la **Weibull** (diferencia de AIC de aproximadamente 100). La distancia de Kolmogorov-Smirnov (dos muestras) es 0,0221 para la Weibull y 0,0573 para la lognormal. La Weibull reproduce los cuantiles hasta el 99 %, mientras que la lognormal sobreestima mucho la cola (cuantil 99 % simulado cercano a 4,2 frente a 2,33 observado). La Weibull subestima algo el extremo de la cola.

**B1: pérdida acumulada de 5 días.** Umbral $u = 5{,}67$ (probabilidad de referencia $\approx 0{,}00098$).

| Método | Estimación | Error estándar |
|---|---|---|
| Monte Carlo simple ($N = 10^5$) | 0,000980 | $9{,}9\times10^{-5}$ |
| Muestreo por importancia ($N = 10^5$) | 0,000992 | $1{,}07\times10^{-5}$ |
| Referencia ($4\times10^6$) | 0,000981 | $1{,}6\times10^{-5}$ |

La **reducción de varianza** del muestreo por importancia es de aproximadamente **85 veces**.

**B2: bootstrap percentil (B = 5000).**

| Estadístico | Estimación | Error estándar | IC 95 % percentil |
|---|---|---|---|
| Desviación estándar | 0,7592 | 0,0200 | [0,7214; 0,7984] |
| VaR al 99 % | 1,9001 | 0,0973 | [1,7325; 2,1420] |

Estas cifras salen de la ejecución con la semilla fija 2026. Los valores exactos que obtengas deben coincidir con los de `resultados/` y con el informe.

## 7. Supuestos y limitaciones

- **Ceros.** Los ceros se excluyen solo en la Fase A, porque Weibull y lognormal solo admiten valores positivos. El modelo ajustado nunca produce un retorno exactamente cero, mientras que en los datos ocurre en el 4 % de los días.
- **Signos.** Asignar $+$ o $-$ con probabilidad 1/2, independiente de la magnitud, es un **supuesto del modelo**, no una propiedad demostrada de los datos.
- **Independencia.** Las simulaciones de B1 y el bootstrap ordinario de B2 suponen observaciones independientes. Los retornos reales muestran **volatilidad agrupada** (rachas de días turbulentos), así que la probabilidad real de rachas de pérdidas puede ser mayor y los intervalos pueden ser demasiado optimistas. Una mejora prevista es el bootstrap por bloques.
- **Cola.** La Weibull subestima el último 0,5 % de la cola, que es justo la zona del evento raro de B1.
- **Kolmogorov-Smirnov.** Como los parámetros se estimaron con los mismos datos, la distancia KS se interpreta como medida descriptiva, no como un contraste formal.
- **Valores cercanos a cero.** No hay datos entre 0 y 0,0198, aunque la Weibull esperaría unos 37. Es compatible con un redondeo de los precios originales, pero es solo una hipótesis que no podemos verificar con los archivos disponibles.

## 8. Fuente de los datos

Los datos provienen de los paquetes de R `Ecdat` (divisas e IBM) y `evir` (BMW y Siemens), tomados del repositorio Rdatasets:

- Arel-Bundock, V. *Rdatasets*. https://github.com/vincentarelbundock/Rdatasets

[Pendiente: completar la referencia del paquete `Ecdat` con `citation("Ecdat")` y la referencia original del conjunto `Garch` que aparece en `?Garch`.]

## 9. Próximos pasos

Para la entrega final: jackknife (sesgo y error estándar), intervalos bootstrap BCa, intervalos asintóticos, el algoritmo EM con mezclas de dos normales (Fase C), y la discusión financiera y las conclusiones (Fase D).

---

## Nota de transparencia sobre el uso de inteligencia artificial

Para la creación de los códigos en R de este proyecto (los scripts `01` a `06`, `funciones_A.R`, `run_all.R`) usamos Claude, de Anthropic (modelo Claude Sonnet 5.5, octubre de 2026). Claude también se usó para explicar los conceptos, las fórmulas y el razonamiento de cada método.

Todos los scripts fueron ejecutados por el grupo en R, y los resultados se contrastaron con valores de referencia (por ejemplo, `MASS::fitdistr` como verificación y las simulaciones de referencia de B1). Las cifras del informe se leen directamente de los archivos que genera el código.
