# Refinamiento de κ₄(PQ) mediante Pattern Search

## Problema

Para una matriz de covarianzas reducida simétrica definida positiva se fija

\[
P=\operatorname{chol}(\Sigma_e,\text{lower})
\]

y se busca una rotación canonicalizada \(Q\in SO(K)\) que reduzca

\[
\kappa_4(PQ)=\|PQ\|_4\,\|(PQ)^{-1}\|_4.
\]

Como \(QQ'=I\), la matriz de impacto conserva la covarianza:

\[
(PQ)(PQ)'=PP'=\Sigma_e.
\]

La inversa no es la traspuesta de \(PQ\). La implementación usa el sistema triangular

\[
(PQ)^{-1}=Q'P^{-1}
\]

mediante `P\eye(K)`, sin llamar a `inv`.

## Bancos fijos

`create_norm4_start_banks(K,S,seed)` crea dos arreglos gaussianos `K×S`:

- `banks.forward`, para estimar `norm(P*Q,4)`;
- `banks.inverse`, para estimar `norm((P*Q)^(-1),4)`.

La función usa un `RandStream` local y no modifica el RNG global. Cada banco se reutiliza para todas las matrices comparadas en una etapa. El ejemplo mantiene bancos independientes para selección y optimización; la validación independiente quedó fuera del alcance actual.

`compute_cond4_pq_fixed_starts` evalúa una matriz o un lote de rotaciones con los mismos bancos. Las opciones controlan la iteración de potencia y la criba, pero el número de inicios lo determina el número de columnas de cada banco.

## Pattern Search sobre SO(K)

`pattern_search_cond4` implementa una búsqueda por patrones determinista sin Global Optimization Toolbox. Desde cada matriz inicial sondea todas las rotaciones de Givens por la derecha,

\[
Q_{\text{trial}}=QG_{ij}(\pm\Delta),\qquad 1\le i<j\le K,
\]

de forma oportunista: acepta la primera mejora suficiente de cada punto inicial. Un sondeo exitoso expande la malla angular y uno fallido la contrae. El orden de las direcciones rota entre iteraciones. Como los pasos son productos de matrices de `SO(K)`, los candidatos permanecen ortogonales con determinante positivo salvo redondeo.

Los puntos iniciales y finales se canonicalizan hacia la identidad. La canonicalización se realiza fuera de los sondeos para no introducir discontinuidades discretas dentro de cada iteración. Debido a que un estimador con un número finito de inicios no es exactamente invariante bajo permutaciones con signos, las matrices finales canonicalizadas se vuelven a evaluar y nunca se devuelve para una trayectoria un valor peor que su punto inicial canonicalizado bajo los bancos de optimización.

El resultado incluye la mejor rotación canonicalizada, los valores iniciales y finales por trayectoria, los puntos finales antes y después de canonicalizar, tamaños finales de malla, iteraciones, evaluaciones y razones de terminación. Es una búsqueda local sobre una estimación multinicio; no certifica el mínimo global ni los máximos exactos de las normas inducidas.

## Ejecución

Desde la raíz del proyecto:

```matlab
run(fullfile('scripts','check_optimizer','check_pattern_search_cond4.m'))
run(fullfile('scripts','optimizer','example_cond4_9x9_pattern_search.m'))
```

El ejemplo genera una `P` sintética de dimensión 9 con una semilla exclusiva, selecciona 100 de 1000 rotaciones Haar mediante bancos de criba y refina las 100 con un segundo par de bancos. Al terminar deja disponibles `QStar`, `kappa4Star`, `impactMatrix` y `searchResult`.
