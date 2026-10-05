# 4-norm

Este proyecto se creó para explorar los componentes del número de condición bajo la norma 4.

El interés matricial corresponde a la norma inducida 4. Como motivación, para una matriz invertible se tiene $\(\kappa_4(A)=\lVert A\rVert_4\lVert A^{-1}\rVert_4\)$. El proyecto estima $\(\lVert Q\rVert_4\)$ para $\(Q\in SO(n)\)$, sin certificar el máximo global ni calcular todavía el número de condición o una Q óptima.

## Alcance actual

- `src/4-norm/norm_4.m`: evalúa la norma vectorial 4 con escalamiento; `norm_4(X,1)` evalúa las columnas de cada página.
- `src/4-norm/power_norm4_multiple_starts.m`: itera desde múltiples vectores iniciales para cada matriz.
- `src/4-norm/compute_induced_norm.m`: estima la norma inducida 4 únicamente con inicios aleatorios; acepta una matriz o un lote y devuelve únicamente un vector columna de estimaciones.
- `src/helpers/haar_so.m`: genera una matriz Haar en \(SO(n)\).
- `src/helpers/angular_amplitude.m` y `src/helpers/spherical_amplitude.m`: evalúan amplitudes para las referencias independientes bidimensional y tridimensional.
- `scripts/example_2x2.m`: dibuja el contorno de \(S_4\) en dimensión 2 y su imagen mediante una única rotación aleatoria.
- `scripts/example_norm4_2x2.m`: compara la estimación por potencia con una referencia angular refinada para la misma Q.
- `scripts/example_norm4_3x3.m`: compara la estimación por potencia con una referencia esférica refinada para la misma Q en dimensión 3 y visualiza la referencia independiente sobre \(S_4\).
- `scripts/example_norm4_9x9.m`: estima la norma inducida para las mismas 1000 matrices de \(SO(9)\) con distintos números de inicios aleatorios, y grafica nueve histogramas comparables. Para 1, 10 y 25 inicios realiza la iteración completa; a partir de 50 aplica la criba aproximada de cuatro iteraciones y tres finalistas.
- `scripts/example_norm4_100_starts.m`: mide únicamente el cálculo para 1000 matrices de \(SO(9)\) y 100 inicios por matriz; imprime las cantidades y el tiempo, sin crear figuras. Usa una criba aproximada: hasta cuatro iteraciones por inicio y refinamiento de los tres mejores de cada matriz.
- `scripts/check_norm4.m`: comprobaciones reproducibles sin framework.

Los ejemplos geométricos normalizan sus direcciones angulares directamente con `norm_4(X,1)`; los experimentos de inicios fijos usan `randn` y la misma normalización. No se necesitan generadores de vectores ni iteradores escalares separados.

Después de añadir `src` a la ruta de MATLAB, la interfaz mínima es:

```matlab
qNorm4 = compute_induced_norm(Q, @norm_4);
options = struct('NumRandomStarts', 50);
Qbatch = cat(3, haar_so(2), haar_so(2), haar_so(2));
qNorm4 = compute_induced_norm(Qbatch, @norm_4, options); % 3 por 1
```

`Q` debe ser double real de tamaño `n×n` o `n×n×M`; cada matriz debe pertenecer a SO(n). La salida tiene tamaño `M×1` y conserva el orden de las matrices; para una sola matriz es escalar. Cada matriz usa sus propios inicios aleatorios. `options` puede especificar parcialmente `NumRandomStarts` (50), `MaxIterations` (1000), `NormTolerance` (1e-10), `StationarityTolerance` (1e-8), `FeasibilityTolerance` (1e-12) y `MaxWorkingMemoryMB` (128). Opcionalmente, `ScreenIterations` (0 por defecto: desactivado) y `NumFinalists` (3) activan una criba aproximada: los inicios descartados podrían haber superado a los finalistas si hubieran seguido iterando. El handle debe implementar matemáticamente la norma 4. Ninguna modalidad certifica un máximo global; véase `docs/induced_norm4_algorithm.qmd`.

Desde la raíz del proyecto:

```matlab
run(fullfile('scripts','example_2x2.m'))
run(fullfile('scripts','example_norm4_2x2.m'))
run(fullfile('scripts','example_norm4_3x3.m'))
run(fullfile('scripts','example_norm4_100_starts.m'))
run(fullfile('scripts','check_norm4.m'))
```

Las dependencias previstas son únicamente MATLAB base. La interfaz ya no devuelve candidatos, diagnósticos ni historiales. Al alcanzar el límite de iteraciones conserva la mejor amplitud factible evaluada, sin reportar un estado de convergencia.

Verificado con MATLAB R2024b Update 9: comprobaciones numéricas, ejemplos 2D y 3D, comparación desde inicios fijos, convergencia multinicio y generación de figuras documentales. Se ejecutó también el ejemplo SO(9) completo con 1000 matrices y nueve cantidades de inicios, y los generadores de figuras documentales. Se comprobó la estructura de las figuras; no se realizó una inspección visual.

El multinicio procesa matrices e inicios por páginas mediante `pagemtimes`. Los inicios gaussianos se normalizan una sola vez mediante el handle, y cada trayectoria tiene su propia parada. Las páginas terminadas se retiran del bloque; con `@norm_4` también se compactan los espacios de inicios terminados cuando eso reduce suficientemente el trabajo, conservando el orden original de los resultados. `MaxWorkingMemoryMB` determina el tamaño inicial de bloque a partir de una estimación conservadora de los temporales; no limita la memoria total del proceso y se procesa como mínimo una matriz.

Todas las normalizaciones, amplitudes y normas de residuos se evalúan mediante `normHandle`. Con `@norm_4` se utiliza la interfaz por columnas `normHandle(X,1)`, conservando el escalamiento. Los handles equivalentes que solo aceptan vectores, como `@(x) norm_4(x)`, se evalúan por columna y pueden tardar más. La interfaz de una sola entrada de `norm_4` conserva el comportamiento para vectores fila y columna.

Para medir el caso de 1000 matrices de SO(9) con 25 inicios:

```matlab
run(fullfile('scripts','benchmark_norm4.m'))
```

El benchmark imprime la mediana de tres ejecuciones después de un calentamiento y excluye la generación de las matrices Q. Incluye validación, generación de inicios y estimación de las normas.

En MATLAB R2024b Update 9, la comparación con las mismas matrices y semillas pasó de 3.156 s a 0.671 s (medianas de tres ejecuciones, 1000 matrices de SO(9), 25 inicios). El benchmark independiente registró 0.682 s en una sesión nueva. La diferencia máxima de estimaciones frente a la implementación anterior fue 2.220e-16. Los tiempos dependen del equipo y del estado de la sesión.
