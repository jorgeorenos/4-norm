# 4-norm

Este proyecto se creó para explorar los componentes del número de condición bajo la norma 4.

El interés matricial corresponde a la norma inducida 4. Como motivación, para una matriz invertible se tiene $\(\kappa_4(A)=\lVert A\rVert_4\lVert A^{-1}\rVert_4\)$. El proyecto estima $\(\lVert Q\rVert_4\)$ para $\(Q\in SO(n)\)$, sin certificar el máximo global ni calcular todavía el número de condición o una Q óptima.

## Alcance actual

- `src/4-norm/generate_unit_l4_vectors.m`: genera una matriz real `n` por `N`; cada columna pertenece a la esfera \(S_4=\{x:\lVert x\rVert_4=1\}\).
- `src/4-norm/norm_4.m`: evalúa la norma vectorial 4 con escalamiento.
- `src/4-norm/power_norm4_single_start.m`: itera desde un vector inicial.
- `src/4-norm/compute_induced_norm.m`: estima la norma inducida 4 con inicios estructurados y aleatorios; devuelve un escalar, un vector candidato y diagnósticos opcionales.
- `src/helpers/haar_so.m`: genera una matriz Haar en \(SO(n)\).
- `src/helpers/angular_amplitude.m` y `src/helpers/spherical_amplitude.m`: evalúan amplitudes para las referencias independientes bidimensional y tridimensional.
- `scripts/example_2x2.m`: dibuja el contorno de \(S_4\) en dimensión 2 y su imagen mediante una única rotación aleatoria.
- `scripts/example_norm4_2x2.m`: compara la estimación por potencia con una referencia angular refinada para la misma Q.
- `scripts/example_norm4_3x3.m`: compara la estimación por potencia con una referencia esférica refinada para la misma Q en dimensión 3 y visualiza ambos candidatos sobre \(S_4\).
- `scripts/check_norm4.m`: comprobaciones reproducibles sin framework.

Después de añadir `src` a la ruta de MATLAB, la interfaz mínima es:

```matlab
qNorm4 = compute_induced_norm(Q, @norm_4);
[qNorm4, xBest, info] = compute_induced_norm(Q, @norm_4);
[qNorm4, xBest, info] = compute_induced_norm(Q, @norm_4, options);
```

`Q` debe ser double real de SO(n). `options` puede especificar parcialmente `NumRandomStarts` (50), `MaxIterations` (1000), `NormTolerance` (1e-10), `StationarityTolerance` (1e-8), `FeasibilityTolerance` (1e-12) y `StoreHistory` (false). El handle debe implementar matemáticamente la norma 4. La estimación no certifica un máximo global; véase `docs/induced_norm4_algorithm.md`.
