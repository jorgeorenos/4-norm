# 4-norm

Este proyecto se creó para explorar los componentes del número de condición bajo la norma 4.

El interés matricial corresponde a la norma inducida 4. Para una matriz invertible se tiene $\(\kappa_4(A)=\lVert A\rVert_4\lVert A^{-1}\rVert_4\)$. El proyecto estima $\(\lVert Q\rVert_4\)$ y $\(\kappa_4(Q)\)$ para $\(Q\in SO(n)\)$, y refina rotaciones canonicalizadas mediante Pattern Search para reducir $\(\kappa_4(PQ)\)$ con un factor de Cholesky fijo $P$. Las estimaciones multinicio y la búsqueda local no certifican máximos de norma ni un mínimo global.

## Alcance actual

- `src/4-norm/norm4_columns.m`: evalúa la norma vectorial 4 por columnas y páginas, con una ruta escalada para magnitudes extremas.
- `src/4-norm/power_norm4_multiple_starts.m`: itera desde múltiples vectores iniciales para cada matriz.
- `src/4-norm/compute_4_norm.m`: estima la norma inducida 4 únicamente con inicios aleatorios; acepta una matriz o un lote y devuelve únicamente un vector columna de estimaciones. No verifica la pertenencia a SO(n): quien llama la garantiza por construcción.
- `src/4-norm/compute_4_cond.m`: estima el número de condición κ₄(Q) = ‖Q‖₄‖Qᵀ‖₄ para Q en SO(n), usando la traspuesta como inversa exacta; acepta una matriz o un lote y devuelve únicamente un vector columna de estimaciones.
- `src/helpers/haar_so.m`: genera una matriz Haar en \(SO(n)\).
- `src/helpers/angular_amplitude.m` y `src/helpers/spherical_amplitude.m`: evalúan amplitudes para las referencias independientes bidimensional y tridimensional.
- `src/optimizer/create_norm4_start_banks.m`: crea bancos gaussianos forward/inverse reproducibles mediante un flujo local, sin modificar el RNG global.
- `src/optimizer/compute_cond4_pq_fixed_starts.m`: estima κ₄(PQ) para una matriz o lote de rotaciones usando los mismos bancos en cada comparación y la inversa $Q'P^{-1}$.
- `src/optimizer/pattern_search_cond4.m`: ejecuta Pattern Search oportunista sobre direcciones de Givens en SO(n), admite una reducción progresiva de trayectorias sin reiniciar su estado, canonicaliza las salidas y conserva diagnósticos por punto inicial.
- `scripts/example_2x2.m`: dibuja el contorno de \(S_4\) en dimensión 2 y su imagen mediante una única rotación aleatoria.
- `scripts/example_norm4_2x2.m`: compara la estimación por potencia con una referencia angular refinada para la misma Q.
- `scripts/example_norm4_3x3.m`: compara la estimación por potencia con una referencia esférica refinada para la misma Q en dimensión 3 y visualiza la referencia independiente sobre \(S_4\).
- `scripts/example_norm4_9x9.m`: estima la norma inducida para las mismas 1000 matrices de \(SO(9)\) con distintos números de inicios aleatorios, y grafica nueve histogramas comparables. Para 1, 10 y 25 inicios realiza la iteración completa; a partir de 50 aplica la criba aproximada de dos iteraciones y cuatro finalistas.
- `scripts/example_norm4_100_starts.m`: mide únicamente el cálculo para 1000 matrices de \(SO(9)\) y 100 inicios por matriz; imprime las cantidades y el tiempo, sin crear figuras. Usa una criba aproximada: hasta dos iteraciones por inicio y refinamiento de los cuatro mejores de cada matriz.
- `scripts/condition_number/example_cond4_100_starts.m`: análogo al anterior para el número de condición κ₄(Q), calculado con `compute_4_cond`.
- `scripts/condition_number/example_cond4_9x9.m`: estima κ₄(Q) para las mismas 1000 matrices de \(SO(9)\) con distintos números de inicios aleatorios, y grafica nueve histogramas comparables; aplica la misma criba a partir de 50 inicios.
- `scripts/check/check_norm4.m`, `scripts/check/check_norm4_packing.m` y `scripts/check/check_norm4_screening.m`: comprobaciones reproducibles del estimador de la norma, sin framework.
- `scripts/check_condition_number/check_cond4.m`: comprobaciones reproducibles del estimador del número de condición, sin framework.
- `scripts/optimizer/example_cond4_9x9_pattern_search.m`: genera una P sintética, selecciona 100 de 10000 rotaciones Haar y aplica un refinamiento progresivo 100→25→5 para reducir κ₄(PQ) con bancos independientes de selección y optimización.
- `scripts/optimizer/benchmark_cond4_9x9_pattern_search.m`: compara de forma reproducible la búsqueda completa con 100 trayectorias, la configuración de 50 trayectorias y el esquema progresivo 100→25→5 usando la misma P, rotaciones y bancos.
- `scripts/check_optimizer/check_pattern_search_cond4.m`: comprueba bancos reproducibles, valores conocidos, equivalencia entre lote y llamadas escalares, descenso, canonicalización y pertenencia a SO(2).

Los ejemplos geométricos normalizan sus direcciones angulares directamente con `norm4_columns(X)`; los experimentos de inicios fijos usan `randn` y la misma normalización. No se necesitan generadores de vectores ni iteradores escalares separados.

Después de añadir `src` a la ruta de MATLAB, la interfaz mínima es:

```matlab
qNorm4 = compute_4_norm(Q);
options = struct('NumRandomStarts', 50);
Qbatch = cat(3, haar_so(2), haar_so(2), haar_so(2));
qNorm4 = compute_4_norm(Qbatch, options); % 3 por 1
kappa4 = compute_4_cond(Qbatch, options); % 3 por 1
```

`Q` debe ser double real de tamaño `n×n` o `n×n×M`; no se verifica la pertenencia a SO(n): los ejemplos la garantizan construyendo las matrices con `haar_so`, y las cotas teóricas $1\leq\lVert Q\rVert_4\leq n^{1/4}$ solo aplican a entradas de SO(n). La salida tiene tamaño `M×1` y conserva el orden de las matrices; para una sola matriz es escalar. Cada matriz usa sus propios inicios aleatorios y el algoritmo emplea internamente `norm4_columns`. `options` puede especificar parcialmente `NumRandomStarts` (50), `MaxIterations` (1000), `NormTolerance` (1e-10), `StationarityTolerance` (1e-8), `FeasibilityTolerance` (1e-12) y `MaxWorkingMemoryMB` (128). Opcionalmente, `ScreenIterations` (0 por defecto: desactivado) y `NumFinalists` (3) activan una criba aproximada: los inicios descartados podrían haber superado a los finalistas si hubieran seguido iterando. Ninguna modalidad certifica un máximo global; véase `docs/induced_norm4_algorithm.qmd`.

`compute_4_cond(Q, options)` acepta las mismas entradas y opciones, y estima κ₄(Q) = ‖Q‖₄‖Qᵀ‖₄ para Q en SO(n), donde la traspuesta es la inversa exacta. Internamente evalúa un único lote intercalado Q, Qᵀ, Q, Qᵀ, ... con `compute_4_norm`, de modo que una llamada por lotes consume la misma secuencia aleatoria que las llamadas secuenciales por matriz. Cada factor es una cota inferior del máximo correspondiente, así que el producto subestima κ₄; en SO(n) el valor verdadero satisface $1\leq\kappa_4(Q)\leq n^{1/2}$.

Desde la raíz del proyecto:

```matlab
run(fullfile('scripts','example_2x2.m'))
run(fullfile('scripts','example_norm4_2x2.m'))
run(fullfile('scripts','example_norm4_3x3.m'))
run(fullfile('scripts','example_norm4_100_starts.m'))
run(fullfile('scripts','condition_number','example_cond4_100_starts.m'))
run(fullfile('scripts','condition_number','example_cond4_9x9.m'))
run(fullfile('scripts','check','check_norm4.m'))
run(fullfile('scripts','check_condition_number','check_cond4.m'))
run(fullfile('scripts','check_optimizer','check_pattern_search_cond4.m'))
run(fullfile('scripts','optimizer','example_cond4_9x9_pattern_search.m'))
run(fullfile('scripts','optimizer','benchmark_cond4_9x9_pattern_search.m'))
```

Las dependencias previstas son únicamente MATLAB base. Pattern Search está implementado en el repositorio y no requiere Global Optimization Toolbox. Las interfaces de estimación de norma conservan únicamente sus valores; el optimizador devuelve por separado un `struct` con candidatos y diagnósticos de la búsqueda. Al alcanzar un límite conserva el mejor candidato evaluado, pero no lo presenta como mínimo global.

Verificado con MATLAB R2024b Update 9: comprobaciones numéricas de la norma y del número de condición, ejemplos 2D y 3D, comparación desde inicios fijos, convergencia multinicio y generación de figuras documentales. Se ejecutaron también el ejemplo SO(9) completo con 1000 matrices y nueve cantidades de inicios, los dos ejemplos de número de condición en `scripts/condition_number/` (1000 matrices de SO(9)) y los generadores de figuras documentales. Las ejecuciones por lotes imprimen los resultados esperados pero el proceso termina con un fallo de salida conocido (`free(): chunks in smallbin corrupted`); no se certifica una salida limpia del proceso. Se comprobó la estructura de las figuras; no se realizó una inspección visual.

El multinicio procesa matrices e inicios por páginas mediante `pagemtimes`. Los inicios gaussianos se normalizan una sola vez con un kernel interno fijo de norma 4, y cada trayectoria tiene su propia parada. Las páginas terminadas se retiran del bloque y se compactan los espacios de inicios terminados cuando eso reduce suficientemente el trabajo, conservando el orden original de los resultados. `MaxWorkingMemoryMB` determina el tamaño inicial de bloque a partir de una estimación conservadora de los temporales; no limita la memoria total del proceso y se procesa como mínimo una matriz.

Para medir el caso de 1000 matrices de SO(9) con 25 inicios:

```matlab
run(fullfile('scripts','benchmark_norm4.m'))
```

El benchmark imprime la mediana de tres ejecuciones después de un calentamiento y excluye la generación de las matrices Q. Incluye validación, generación de inicios y estimación de las normas.

En MATLAB R2024b Update 9, la comparación con las mismas matrices y semillas pasó de 3.156 s a 0.671 s (medianas de tres ejecuciones, 1000 matrices de SO(9), 25 inicios). El benchmark independiente registró 0.682 s en una sesión nueva. La diferencia máxima de estimaciones frente a la implementación anterior fue 2.220e-16. Los tiempos dependen del equipo y del estado de la sesión.
