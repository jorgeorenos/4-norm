# Proyecto 4-norm

## Objetivo y alcance actual

Construir un proyecto MATLAB para explorar los componentes del número de condición bajo la norma inducida 4. Se incluyen la generación de vectores unitarios, matrices Haar en SO(n), geometría en dimensión 2 y una estimación por potencia generalizada con multinicio de la norma inducida 4 de Q en SO(n).

Este archivo debe ubicarse en la raíz del proyecto `4-norm` y orientar el trabajo de Codex en todo el repositorio. Si la carpeta actual ya es la raíz del proyecto, no crear otra carpeta `4-norm` anidada.

El candidato numérico no se considera necesariamente x* global. No implementar todavía optimización de Q, minimización del número de condición ni evaluación de PQ o su inversa.

## Organización del repositorio

| Ruta relativa a `4-norm/` | Contenido |
| --- | --- |
| `AGENTS.md` | Instrucciones de trabajo para Codex. |
| `README.md` | Objetivo, estructura y forma de ejecutar el ejemplo. |
| `.gitignore` | Archivos temporales y resultados regenerables. |
| `docs/geometry_norm4.md` | Definiciones, métodos de generación e interpretación geométrica. |
| `docs/worklog.md` | Registro breve del trabajo y su verificación. |
| `src/helpers/haar_so.m` | Generación Haar de matrices de SO(n). |
| `scripts/example_2x2.m` | Ejemplo reproducible con las dos gráficas solicitadas. |
| `src/4-norm/norm_4.m` | Norma vectorial 4 escalada. |
| `src/4-norm/power_norm4_multiple_starts.m` | Iteración generalizada para múltiples inicios y matrices. |
| `src/4-norm/compute_induced_norm.m` | Estimaciones inducidas 4 por matriz con multinicio. |
| `scripts/example_norm4_2x2.m` | Estimación y comparación angular en SO(2). |
| `src/helpers/spherical_amplitude.m` | Amplitud de norma 4 en una dirección esférica de R³. |
| `src/helpers/angular_amplitude.m` | Amplitud de norma 4 en una dirección angular de R². |
| `scripts/example_norm4_3x3.m` | Estimación y comparación esférica en SO(3). |
| `scripts/check_norm4.m` | Comprobaciones pequeñas sin framework. |
| `docs/induced_norm4_algorithm.md` | Método, tolerancias y limitaciones. |

Las funciones sustantivas deben residir en `src/`; los scripts deben coordinar llamadas y visualizaciones. Mantener la documentación en español y los identificadores de código en inglés. Usar MATLAB base, sin dependencias de toolboxes adicionales ni lenguajes externos. Para estimar la norma inducida usar `compute_induced_norm(Q,@norm_4)`; la única salida es un vector columna de estimaciones para una entrada `n×n×M`, o un escalar para una matriz `n×n`. No devolver candidatos, diagnósticos ni historiales. La matriz Q permanece fija durante el multinicio; los handles equivalentes deben implementar matemáticamente la norma 4.

## Convenciones matemáticas

Trabajar con vectores y matrices reales en doble precisión. Cada vector es una columna. Para N vectores de dimensión n, usar una matriz X de tamaño n por N.

La norma vectorial es:

$$
\|x\|_4=\left(\sum_{j=1}^{n}|x_j|^4\right)^{1/4}.
$$

La esfera unitaria es la superficie S₄ = {x : ‖x‖₄ = 1}; la bola B₄ = {x : ‖x‖₄ ≤ 1} incluye el interior. El generador inicial debe producir puntos de la superficie. En dimensión 2, esta superficie es una curva cerrada.

No confundir la norma inducida de una matriz con la raíz cuarta de la suma de las cuartas potencias de sus entradas. Tampoco utilizar `norm(X,4)` para obtener las normas de las columnas de una matriz: calcularlas explícitamente por columna.

## Vectores para ejemplos y pruebas

Los scripts generan direcciones gaussianas con `randn(n,N)` o direcciones angulares con `[cos(theta); sin(theta)]` y normalizan por columnas con `norm_4(Z,1)`. El estimador genera internamente sus inicios gaussianos y los normaliza mediante el handle recibido. Esta construcción produce puntos sobre S₄, sin garantizar uniformidad respecto al área ni a la medida de cono. El espaciado angular no equivale a espaciado uniforme por longitud de arco; cerrar la curva solamente al graficar. No se necesita una función pública de generación de vectores.

## Función `haar_so`

Interfaz:

```matlab
Q = haar_so(n)
```

Aceptar cualquier entero positivo n y devolver una matriz real n por n distribuida según Haar en:

$$
SO(n)=\{Q:Q^TQ=I_n,\ \det(Q)=1\}.
$$

Implementar el siguiente método:

1. Generar `G = randn(n,n)` con entradas normales estándar independientes.
2. Obtener `[Q,R] = qr(G)` sin pivoteo.
3. Corregir los signos de las columnas con `d = sign(diag(R))` y `Q = Q * diag(d)`. Si alguna entrada de `diag(R)` es exactamente cero, regenerar G. Esta corrección corresponde a elegir la factorización con diagonal positiva en R y produce la distribución Haar en O(n).
4. Si `det(Q) < 0`, cambiar el signo de una columna fija, por ejemplo `Q(:,end) = -Q(:,end)`.
5. Devolver Q. Para n = 1, devolver exactamente 1, el único elemento de SO(1).

Explicar en la documentación que el paso final lleva las dos componentes de O(n) a SO(n) y conserva la invariancia Haar bajo multiplicación por elementos de SO(n). No basta con aplicar QR y corregir únicamente el determinante: también se necesita la corrección de los signos de R.

No usar `rand` en lugar de `randn` para G. No normalizar toda Q mediante un escalar para corregir el determinante. No exigir igualdad exacta en las comprobaciones de punto flotante.

Las funciones deben consumir el estado aleatorio existente: no llamar a `rng` ni reiniciar la semilla dentro de ellas. La reproducibilidad se controla desde los scripts.

## Script `scripts/example_2x2.m`

El script debe ser ejecutable desde la raíz del proyecto y localizar `src` a partir de su propia ubicación mediante `mfilename('fullpath')`, `fileparts` y `fullfile`, sin rutas absolutas dependientes de un equipo. Evitar `cd`, `savepath` y cambios persistentes a la configuración de MATLAB.

Secuencia requerida:

1. Definir una semilla explícita al comienzo, por ejemplo `rng(42,'twister')`.
2. Definir `n = 2` y `N = 1000` como parámetros visibles.
3. Generar `theta = (0:N-1)*(2*pi/N)`, `Z = [cos(theta); sin(theta)]` y `X = Z ./ norm_4(Z,1)`.
4. Generar una sola matriz `A = haar_so(n)` y mostrarla en la consola.
5. Calcular `Y = A*X`. No volver a normalizar Y: eso borraría la información de la transformación que se quiere observar.
6. Comprobar las normas de X, la ortogonalidad de A y su determinante; mostrar los errores numéricos.
7. Cerrar las curvas para el dibujo con `Xclosed = [X,X(:,1)]` y `Yclosed = [Y,Y(:,1)]`.
8. Crear dos gráficas, preferentemente en una sola figura con dos paneles mediante `tiledlayout(1,2)`.

### Primera gráfica

Mostrar únicamente el contorno de S₄ original, con ejes `x_1` y `x_2` y título «Esfera unitaria de la norma 4».

### Segunda gráfica

Superponer el contorno original en gris discontinuo y la imagen `A(S₄)` en naranja continuo. Usar ejes de coordenadas `u_1` y `u_2`, título «Esfera unitaria y su imagen bajo A» y leyenda que distinga «Original: norma 4 = 1» de «Transformada: y = Ax».

Usar `axis equal`, cuadrícula, líneas legibles y los mismos límites en ambos paneles. Calcular los límites a partir de X e Y, con margen, para evitar recortes y comparaciones engañosas. Usar los mismos puntos y la misma A en toda la visualización.

Explicar en comentarios y documentación:

- En SO(2), A es una rotación: conserva longitudes euclidianas, ángulos y área. La curva naranja es la curva original rotada, sin cizallamiento ni estiramiento euclidiano.
- Una rotación general no conserva la norma 4, porque su esfera unitaria no tiene simetría bajo todas las rotaciones. Por ello, las columnas de Y no tienen necesariamente norma 4 igual a uno.
- Algunas rotaciones especiales, como los múltiplos de π/2, sí hacen coincidir los contornos. No descartar muestras Haar ni repetir sorteos para conseguir una figura visualmente más llamativa.
- Estas gráficas ilustran geometría; el script geométrico no identifica x* ni calcula la norma inducida de A. La estimación se realiza por separado en `example_norm4_2x2.m`.

No usar `clear all`, `close all` ni modificar propiedades globales de figuras. La exportación de imágenes es opcional; si se incorpora, debe controlarse con una variable y guardar los resultados en una subcarpeta ignorada, por ejemplo `docs/figures/generated/`.

## README y documentación

El README debe comenzar con el nombre del proyecto y esta declaración:

> Este proyecto se creó para explorar los componentes del número de condición bajo la norma 4.

Aclarar que el interés matricial corresponde a la norma inducida 4. Puede incluirse, como motivación para matrices invertibles, la definición κ₄(A) = ‖A‖₄ ‖A⁻¹‖₄, sin implementar todavía el cálculo del número de condición.

Describir el alcance actual, las funciones, el formato n por N y las dos gráficas. Incluir la ejecución desde la raíz:

```matlab
run(fullfile('scripts','example_2x2.m'))
```

Documentar la versión de MATLAB efectivamente utilizada para verificar el proyecto; no afirmar compatibilidad verificada con versiones no probadas. Indicar que las dependencias previstas son únicamente MATLAB base.

En `docs/geometry_norm4.md`, explicar la normalización, la distinción entre esfera y bola, la distribución de los puntos generados y el algoritmo Haar. En `docs/worklog.md`, registrar fechas, cambios relevantes, decisiones y verificaciones efectivamente realizadas.

Referencias técnicas para documentar los métodos:

- Francesco Mezzadri, *How to generate random matrices from the classical compact groups*: https://arxiv.org/abs/math-ph/0609050
- MATLAB, `qr`: https://www.mathworks.com/help/matlab/ref/qr.html
- MATLAB, `randn`: https://www.mathworks.com/help/matlab/ref/double.randn.html
- MATLAB, `rng`: https://www.mathworks.com/help/matlab/ref/rng.html

## Git y archivos generados

Crear un `.gitignore` acotado para temporales de MATLAB, respaldos del editor, archivos del sistema operativo y `docs/figures/generated/`. Mantener versionados los archivos `.m`, `.md` y `.gitignore`. No ignorar indiscriminadamente todas las imágenes o todos los datos, pues podrían ser documentación o entradas del proyecto.

No crear repositorios remotos, publicar ni hacer push sin una solicitud del usuario.

## Criterios de aceptación

- Las funciones admiten dimensiones generales; solamente el modo angular y el ejemplo son específicos de dimensión 2.
- Verificar tamaños, finitud y normas de columnas para dimensiones representativas, por ejemplo n = 1, 2, 5 y 10.
- Usar como referencia tolerancias de `1e-12` para las normas vectoriales y `1e-12*max(1,n)` para `norm(Q'*Q-eye(n),'fro')` y `abs(det(Q)-1)` en esas dimensiones.
- Comprobar que restablecer la misma semilla y repetir las mismas llamadas reproduce los resultados.
- Revisar el manejo de dimensiones y cantidades inválidas.
- Ejecutar el ejemplo y revisar que existan ambos paneles, las curvas estén cerradas, la leyenda sea correcta y no haya recortes.
- Las comprobaciones numéricas de ortogonalidad y determinante no demuestran por sí solas la distribución Haar; esa propiedad depende del algoritmo documentado.
- No introducir un framework de pruebas ni dependencias adicionales. Bastan comprobaciones pequeñas de las propiedades matemáticas, el estimador y los ejemplos.
- Si MATLAB no está disponible, declarar explícitamente que el código no se ejecutó en MATLAB; no presentar una inspección estática como una prueba de ejecución.
- Al terminar una implementación, resumir los archivos creados, cómo ejecutar los ejemplos y qué se verificó. No declarar que la estimación multinicio certifique un máximo global.
