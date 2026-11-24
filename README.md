<<<<<<< HEAD
# Proyecto final de Computación Estadística 2026-II
### Riesgo de pérdida en la serie de retornos GBP/USD: simulación, remuestreo y mezclas
**Universidad Nacional de Colombia, sede de La Paz**
**Curso:** Computación Estadística, Grupo I (semestre 2026-II). **Docente:** Yoselin Nieto Gil
**Grupo de trabajo:** Grupo 2
**Integrantes:**
-Adriana Jaimes
-Crisdailyn Rey Vargas
-John Rincón Toro
-Maria Vega Gutierrez


---

## 1. ¿De qué trata este proyecto?

Los mercados financieros tienen días tranquilos y días turbulentos, y las pérdidas grandes se concentran en pocos días. El proyecto responde una sola pregunta: **¿cuánto puede perderse en un día o en una semana al mantener un activo, y qué tan confiable es ese número?**

Cada grupo trabaja con la serie de retornos diarios de un activo distinto y programa los métodos del curso **desde cero y sin cajas negras**: generación de variables aleatorias (Unidad 1), integración por Monte Carlo, muestreo por importancia, jackknife y bootstrap (Unidad 2), y modelos de mezcla con el algoritmo EM (Unidad 3).

Nuestra serie es **GBP/USD** (tipo de cambio de la libra esterlina frente al dólar, según la denominación de la guía del proyecto):

| Característica | Valor |
|---|---|
| Observaciones | 1.866 |
| Período (según la guía) | 2 de enero de 1980 a 21 de mayo de 1987 |
| Variable | `retorno_pct`: retorno diario en porcentaje, $100\,(\ln P_t - \ln P_{t-1})$ |

El archivo de datos no trae fechas, solo el número de día `t`. El período se toma de la guía del proyecto.

## 2. Estado del trabajo

| Parte | Contenido | Scripts |
|---|---|---|
| Datos | Carga, filtro de la serie, ceros y variable de trabajo | `01` |
| Fase A | Máxima verosimilitud (lognormal y Weibull), AIC, generadores, validación | `02`, `03`, `04` |
| Fase B1 | Monte Carlo y muestreo por importancia | `05` |
| Fase B2 | Bootstrap percentil, jackknife, BCa e intervalos asintóticos | `06`, `07` |
| Fase C | Algoritmo EM, inicios aleatorios, tolerancias, BIC, regímenes y cola implicada | `08`, `09`, `10` |
| Fase D | VaR por modelo y autocorrelación | `11` |
| Informes | Avance (28 de octubre) e informe final | `informe_avance.Rmd`, `informe_final.Rmd` |

## 3. Estructura del repositorio

```
proyecto-computacion-estadistica/
├── README.md                          Este archivo
├── run_all.R                          Ejecuta todo el análisis de principio a fin (scripts 01 a 11)
├── funciones_A.R                      Generadores de variables aleatorias (solo definiciones)
├── funciones_EM.R                     Pasos E y M del EM y log-verosimilitud (solo definiciones)
├── 01_datos.R                         Datos, ceros y variable de trabajo
├── 02_mle.R                           Máxima verosimilitud: lognormal y Weibull (Newton) y AIC
├── 03_generadores.R                   Comprobación de los generadores
├── 04_validacion.R                    QQ-plot, Kolmogorov-Smirnov a mano y cuantiles
├── 05_B1.R                            Monte Carlo y muestreo por importancia
├── 06_B2.R                            Bootstrap percentil
├── 07_B2_completo.R                   Jackknife, BCa e intervalos asintóticos
├── 08_EM.R                            Algoritmo EM (esqueleto de la docente completado) y convergencia
├── 09_EM_inicios.R                    12 inicios aleatorios, tolerancias y casos degenerados
├── 10_EM_interpretacion.R             BIC frente a una normal, regímenes y cola implicada
├── 11_FaseD.R                         VaR por modelo y autocorrelación
├── informe_avance.Rmd                 Informe del avance (Fase A, B1 y bootstrap percentil)
├── informe_final.Rmd                  Informe final (Fases A a D)
├── data/
│   └── retornos_diarios_12_series.csv Datos entregados por la docente (12 series)
├── resultados/                        Archivos .rds generados por los scripts
└── figuras/                           Gráficos .png generados por los scripts
```

Los scripts buscan el archivo de datos en `data/` y, si no está allí, en la carpeta principal. Las carpetas `resultados/` y `figuras/` se crean solas al ejecutar el código.

## 4. Requisitos

- **R** (versión 4.x o superior recomendada) y, opcionalmente, **RStudio**.
- Para compilar los informes: los paquetes `rmarkdown` y `knitr`.
- Para generar los informes en PDF: una instalación de LaTeX (por ejemplo `tinytex::install_tinytex()`). Los informes también pueden generarse como HTML autónomo, sin LaTeX, que es la opción por defecto.
- El paquete `MASS` (incluido en las instalaciones habituales de R) se usa **únicamente para verificar** nuestras estimaciones con `fitdistr`. Las funciones `pnorm` y `qnorm` de R base también se usan solo para comprobar nuestras implementaciones. Ninguna de ellas produce un resultado oficial.

Todo lo demás se implementó con funciones base de R (generador uniforme, `log`, `exp`, `sqrt`, ordenar, cuantiles, álgebra de matrices y gráficos).

## 5. Cómo reproducir los resultados

La semilla aleatoria es fija: cada script ejecuta `set.seed(2026)` al inicio, y por eso los números, tablas y gráficos son los mismos en cada ejecución.

1. Abre R (o RStudio) y establece como carpeta de trabajo la carpeta del repositorio. En RStudio: *Session → Set Working Directory → Choose Directory*. Comprueba con `getwd()` y `list.files()`.
2. Verifica que el archivo `retornos_diarios_12_series.csv` esté en `data/` o en esa misma carpeta.
3. Ejecuta, en la consola de R:

   ```r
   source("run_all.R")
   ```

   Esto corre los once scripts en orden. El script `05_B1.R` realiza simulaciones grandes y puede tardar varios minutos.
4. Compila el informe final con *Knit* en RStudio, o con:

   ```r
   rmarkdown::render("informe_final.Rmd")
   ```

   El informe lee los resultados guardados en `resultados/` y las figuras de `figuras/`, así que sus números coinciden con los del código. **El archivo `.Rmd` debe estar en la misma carpeta que `resultados/` y `figuras/`.**

También puedes ejecutar cada script por separado, en el orden numérico. Cada uno depende de los archivos que dejan los anteriores en `resultados/` (por ejemplo, `03` a `06` usan `params_A.rds` de `02`, y `11` usa los resultados de `06` y `10`).

Todas las rutas son **relativas**: no hay rutas locales del tipo `C:/Users/...`.

## 6. Qué hace cada script

### `01_datos.R`: datos, ceros y variable de trabajo
Lee el CSV, filtra **solo** la serie `GBP/USD`, la ordena por `t` y verifica que coincide con la guía (1.866 datos, `t` consecutivo, sin valores faltantes). Cuenta los retornos exactamente iguales a cero y calcula su porcentaje. Construye la variable de trabajo de la Fase A, $x = |r|$ sin los ceros exactos, y calcula estadísticos descriptivos. Guarda `figuras/01_exploracion.png`.

### `02_mle.R`: máxima verosimilitud y AIC
Estima por máxima verosimilitud, con implementación propia, dos familias para $x$:
- **Lognormal:** solución cerrada, $\hat\mu = \overline{\ln x}$ y $\hat\sigma^2 = \overline{(\ln x - \hat\mu)^2}$.
- **Weibull:** $\hat\lambda^k = \overline{x^k}$ y la forma $\hat k$ como raíz de una ecuación que se resuelve con el **método de Newton**. Esa ecuación tiene raíz única porque su derivada es siempre positiva.

Compara las familias con el **AIC**, verifica la convergencia desde varios valores iniciales y contrasta con `MASS::fitdistr` solo como verificación. Guarda `resultados/params_A.rds`.

### `funciones_A.R` y `03_generadores.R`: generación de variables aleatorias
Usan únicamente el generador uniforme base `runif`:
- **Weibull por transformada inversa:** $X = \lambda\,[-\ln U]^{1/k}$.
- **Normal estándar por Box-Muller:** $Z = \sqrt{-2\ln U_1}\cos(2\pi U_2)$.
- **Lognormal:** $X = \exp(\mu + \sigma Z)$.
- **Signos:** $+1$ o $-1$ con probabilidad 1/2.

`03_generadores.R` compara medias, desviaciones y cuantiles simulados con sus valores teóricos y comprueba que $F(X)$ es uniforme.

### `04_validacion.R`: validación de las muestras simuladas
Genera muestras sintéticas 20 veces más grandes que la real y las compara con los datos mediante: el gráfico de ajuste de densidades con **eje vertical logarítmico**, los gráficos **cuantil-cuantil**, la **distancia de Kolmogorov-Smirnov calculada a mano** (la función de error se programó por serie de Taylor) y los **cuantiles 50, 90 y 99**. Guarda `figuras/02_ajuste_log.png`, `figuras/03_qq.png` y `resultados/validacion_A.rds`.

### `05_B1.R`: Monte Carlo y muestreo por importancia
Estima la probabilidad de que la **pérdida acumulada de 5 días consecutivos supere un umbral $u$** bajo el modelo de la Fase A. El umbral se elige con una simulación de referencia de $4\times10^6$ réplicas para que la probabilidad sea cercana al 0,1 %. Se estima con Monte Carlo simple y con **muestreo por importancia**, con una propuesta diseñada por el grupo: la escala de la Weibull se multiplica por $c$ y la probabilidad de que un día sea de pérdida pasa de 1/2 a $q$. Incluye un **piloto** que compara muchas propuestas y muestra que una mala elección puede aumentar la varianza. Guarda `resultados/B1.rds`.

### `06_B2.R` y `07_B2_completo.R`: remuestreo
Sobre los **retornos originales, con los ceros**, estudian la **desviación estándar** y el **VaR al 99 %** diario (cuantil 99 % de la pérdida $-r$):
- `06`: bootstrap percentil con $B = 5000$ réplicas (remuestreo programado a mano, sin paquetes).
- `07`: **jackknife** (sesgo y error estándar), intervalo **BCa** (la aceleración sale de las réplicas jackknife de la muestra original) e **intervalos asintóticos** (método delta para la desviación estándar; estimación de la densidad con núcleo gaussiano para el VaR, con análisis de sensibilidad al ancho de banda). La normal estándar (cola y cuantil) se programó a mano.

Guardan `figuras/04_bootstrap.png`, `figuras/05_intervalos.png`, `resultados/B2.rds` y `resultados/B2_completo.rds`.

### `funciones_EM.R`, `08_EM.R`, `09_EM_inicios.R`, `10_EM_interpretacion.R`: Fase C
Parten del esqueleto de la docente (`esqueleto_EM.R`), con los bloques `TODO` completados sin cambiar los nombres de las funciones ni la estructura de `theta`:
- `08`: pasos E y M, bucle EM, prueba con datos simulados de parámetros conocidos (`verificar_em()`) y EM sobre la serie real (retornos **con signo**, incluidos los ceros). Registra la log-verosimilitud por iteración y comprueba que nunca disminuye. Guarda `figuras/06_convergencia_EM.png`.
- `09`: **12 inicios aleatorios**, **cuatro tolerancias de parada** y **casos degenerados** (componentes idénticos, desviación que tiende a cero). Guarda `figuras/07_inicios_EM.png`.
- `10`: comparación con una sola normal mediante **BIC y AIC**, regímenes de calma y turbulencia, probabilidad diaria de turbulencia y **cola implicada** $P(|r|>\text{percentil 99})$. Guarda `figuras/08_prob_turbulencia.png`.

### `11_FaseD.R`: VaR por modelo y dependencia temporal
Calcula el VaR a los niveles 95, 99, 99,5 y 99,9 % con el cuantil empírico, una normal única, la mezcla del EM y la Weibull de la Fase A con signos. Calcula, con función propia, la **autocorrelación** de $r$, $|r|$, $r^2$ y de la probabilidad de turbulencia, y un efecto de diseño aproximado. Guarda `figuras/09_acf.png`.

### `informe_avance.Rmd` e `informe_final.Rmd`
Informes en R Markdown que siguen la estructura de la guía del informe. Las cifras, tablas y figuras se leen de `resultados/` y `figuras/`, y el anexo incluye el código.

## 7. Resultados principales (semilla 2026)

**Datos.** 75 de los 1.866 retornos valen exactamente cero (4,019 %). La variable de trabajo $|r|$ tiene 1.791 datos. La curtosis en exceso de los retornos es 3,06, es decir, colas más pesadas que las de una normal.

**Fase A: máxima verosimilitud y AIC.**

| Familia | Parámetros | Log-verosimilitud | AIC |
|---|---|---|---|
| Lognormal | $\hat\mu = -1{,}0035$, $\hat\sigma = 1{,}0411$ | −816,10 | 1636,21 |
| Weibull | $\hat k = 1{,}1367$, $\hat\lambda = 0{,}5997$ | −766,13 | 1536,26 |

Se elige la **Weibull** (diferencia de AIC de aproximadamente 100). La distancia de Kolmogorov-Smirnov (dos muestras) es 0,0221 para la Weibull y 0,0573 para la lognormal. La Weibull reproduce los cuantiles hasta el 99 %, mientras que la lognormal sobreestima mucho la cola (cuantil 99 % simulado cercano a 4,2 frente a 2,33 observado). La Weibull subestima algo el último 0,5 % de la cola.

**B1: pérdida acumulada de 5 días.** Umbral $u = 5{,}67$ (probabilidad de referencia $\approx 0{,}00098$).

| Método | Estimación | Error estándar |
|---|---|---|
| Monte Carlo simple ($N = 10^5$) | 0,000980 | $9{,}9\times10^{-5}$ |
| Muestreo por importancia ($N = 10^5$) | 0,000992 | $1{,}07\times10^{-5}$ |
| Referencia ($4\times10^6$) | 0,000981 | $1{,}6\times10^{-5}$ |

La **reducción de varianza** del muestreo por importancia es de aproximadamente **85 veces**.

**B2: intervalos de confianza al 95 %** (retornos originales, $B = 5000$).

| Estadístico | Estimación | Percentil | BCa | Asintótico |
|---|---|---|---|---|
| Desviación estándar | 0,7592 | [0,7214; 0,7984] | [0,7250; 0,8029] | [0,7205; 0,7979] |
| VaR al 99 % | 1,9001 | [1,7325; 2,1420] | [1,7494; 2,1712] | [1,7324; 2,0677] |

Para la desviación estándar, los métodos coinciden. Para el VaR se separan: el jackknife subestima el error estándar (0,045 frente a 0,097 del bootstrap; el cuantil solo toma 3 valores al quitar una observación) y el intervalo asintótico es simétrico y depende del ancho de banda.

**Fase C: EM.** Mezcla de dos normales sobre los retornos con signo:

| Régimen | Proporción | Media | Desviación |
|---|---|---|---|
| Calma | 0,4542 | −0,0478 | 0,3675 |
| Turbulencia | 0,5458 | 0,0111 | 0,9703 |

El EM converge en 366 iteraciones (tolerancia $10^{-8}$; 458 con $10^{-10}$), con log-verosimilitud −2051,555 que nunca disminuye. Los 12 inicios aleatorios convergen a la misma solución. El BIC es 4140,8 para la mezcla y 4281,3 para una normal única (diferencia de 140,6). La probabilidad de que $|r|$ supere el percentil 99 de los datos es 0,0102 (observada), 0,0100 (mezcla), 0,0102 (Weibull) y 0,0026 (normal única). Un inicio con un componente muy estrecho en 0 colapsa con $\hat\pi = 0{,}0402 = 75/1866$, la proporción de ceros exactos.

**Fase D: VaR por modelo.**

| Nivel | Empírico | Normal | Mezcla | Weibull con signos |
|---|---|---|---|---|
| 95 % | 1,260 | 1,264 | 1,282 | 1,249 |
| 99 % | 1,900 | 1,781 | 2,017 | 1,991 |
| 99,5 % | 2,201 | 1,971 | 2,278 | 2,298 |
| 99,9 % | 3,002 | 2,361 | 2,808 | 2,992 |

La normal subestima el VaR en los niveles altos. Los retornos con signo no están autocorrelacionados, pero $|r|$ y $r^2$ sí (29 de 30 rezagos fuera de la banda de ruido blanco), lo que indica volatilidad agrupada.

Estas cifras salen de la ejecución con la semilla fija 2026. Los valores exactos que obtengas deben coincidir con los de `resultados/` y con el informe.

## 8. Cumplimiento de las reglas de programación

| Regla de la guía | Cómo se cumple |
|---|---|
| Sin funciones de ajuste de distribuciones (`fitdistr`) | La máxima verosimilitud es propia. `MASS::fitdistr` solo se usa para verificar, y se muestran ambos resultados lado a lado. |
| Sin generadores no uniformes (`rweibull`, `rlnorm`) | Todo se genera desde `runif` con transformada inversa y Box-Muller. |
| Sin paquetes de remuestreo ni de mezclas | El bootstrap y el jackknife se programaron a mano, y el EM se programó desde el esqueleto de la docente. |
| Función de error para la normal | Se construye con `erf` por serie de Taylor y fracción continua; `pnorm` y `qnorm` solo verifican. |
| Semilla fija | `set.seed(2026)` al inicio de cada script. |
| Ejecución de principio a fin sin intervención manual | `source("run_all.R")`. |
| Rutas relativas | No se usan rutas locales. |

## 9. Supuestos y limitaciones

- **Ceros.** Se excluyen solo en la Fase A, porque Weibull y lognormal solo admiten valores positivos. El modelo ajustado nunca produce un retorno exactamente cero, mientras que en los datos ocurre en el 4 % de los días.
- **Signos.** Asignar $+$ o $-$ con probabilidad 1/2, independiente de la magnitud, es un **supuesto del modelo**, no una propiedad demostrada de los datos.
- **Independencia.** Las simulaciones de B1, el bootstrap ordinario y la mezcla suponen observaciones independientes. La autocorrelación de $|r|$ lo contradice, así que los intervalos están probablemente subestimados. Un efecto de diseño aproximado de 3,3 sugiere que los errores estándar de estadísticos basados en $r^2$ podrían estar subestimados en un factor cercano a 1,8. Es una aproximación gruesa, no una corrección.
- **Cola.** La Weibull subestima el último 0,5 % de la cola, que es la zona del evento raro de B1.
- **Kolmogorov-Smirnov.** Como los parámetros se estimaron con los mismos datos, la distancia KS se interpreta como medida descriptiva, no como un contraste formal.
- **Jackknife para el VaR.** Un cuantil es un estadístico no suave: el jackknife no es consistente para él, y la aceleración del BCa en ese caso es poco confiable.
- **Mezcla.** La log-verosimilitud de una mezcla normal no está acotada: hay soluciones degeneradas (desviación que tiende a cero) que deben descartarse. La mezcla clasifica cada día solo por su retorno, sin mirar los días vecinos, y por eso no captura la persistencia de los regímenes.
- **Fechas.** El archivo no trae fechas. Para ligar los días de turbulencia con hechos del mercado solo se pudo acotar la fecha (el período tiene 1.927 días hábiles para 1.866 observaciones), así que la relación con los eventos es aproximada.
- **Valores cercanos a cero.** No hay datos entre 0 y 0,0198, aunque la Weibull esperaría unos 37. Es compatible con un redondeo de los precios originales, pero es solo una hipótesis que no podemos verificar con los archivos disponibles.

## 10. Fuente de los datos

Los datos provienen de los paquetes de R `Ecdat` (divisas e IBM) y `evir` (BMW y Siemens), tomados del repositorio Rdatasets:

- Arel-Bundock, V. *Rdatasets*. https://github.com/vincentarelbundock/Rdatasets

[Pendiente: completar la referencia del paquete `Ecdat` con `citation("Ecdat")` y la referencia original del conjunto `Garch` que aparece en `?Garch`.]

Las referencias de prensa usadas para contextualizar los días de turbulencia están en `informe_final.Rmd`.

---

## Nota de transparencia sobre el uso de inteligencia artificial

Para la creación de los **códigos en R** de este proyecto (los scripts `01` a `11`, `funciones_A.R`, `funciones_EM.R`, `run_all.R` y los borradores de `informe_avance.Rmd`, `informe_final.Rmd` y de este `README.md`) usamos **Claude**, de Anthropic (modelo Claude Sonnet 5.5, entre octubre y noviembre de 2026). Claude también se usó para explicar los conceptos, las fórmulas y el razonamiento de cada método, y para buscar y verificar las fuentes sobre los hechos del mercado citados en el informe.

Todos los scripts fueron ejecutados por el grupo en R, y los resultados se contrastaron con valores de referencia calculados de forma independiente (por ejemplo, `MASS::fitdistr` como verificación y las simulaciones de referencia de B1). Las cifras de los informes se leen directamente de los archivos que genera el código.

[Completar: qué parte del código y del análisis es producto propio del grupo.]
