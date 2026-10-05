> Documento histórico de la implementación inicial. La interfaz vigente desde 2026-10-02 devuelve únicamente un vector columna para matrices `n×n×M`; se eliminaron las salidas adicionales y `StoreHistory`. Para el contrato actual y la ejecución, consultar `README.md` y `docs/induced_norm4_algorithm.qmd`. Los ejemplos de interfaz y diagnósticos que siguen están obsoletos.

# Implementar la norma inducida 4 en el proyecto MATLAB 4-norm

## Encargo para Codex

Implementa la siguiente etapa del proyecto existente: estimar la norma matricial inducida 4 de una matriz Q de SO(n), mediante un método de potencia generalizada con múltiples inicializaciones. Reutiliza las funciones de generación de vectores y matrices Haar que ya existan.

La interfaz obligatoria recibe Q como primer argumento y una función handle de norma vectorial como segundo argumento:

```matlab
qNorm4 = compute_induced_norm(Q, @norm_4);
```

El resultado principal debe ser un escalar con la estimación de la norma inducida 4. Permite obtener, opcionalmente, el mejor vector y diagnósticos:

```matlab
[qNorm4, xBest, info] = compute_induced_norm(Q, @norm_4);
[qNorm4, xBest, info] = compute_induced_norm(Q, @norm_4, options);
```

No cambiar el orden de los dos argumentos obligatorios. No devolver una estructura como primera salida. No pedir al usuario que genere los vectores iniciales para que la llamada mínima funcione.

Este encargo amplía explícitamente el alcance inicial del AGENTS.md: ahora se autoriza implementar la búsqueda de vectores candidatos y la estimación de la norma inducida 4. Actualiza las cláusulas del AGENTS.md que excluían esta etapa para que no contradigan el trabajo solicitado. Mantén las convenciones de directorios, MATLAB base y documentación en español. No implementar todavía la búsqueda de una Q óptima, la minimización del número de condición ni la evaluación de PQ o de su inversa.

## Archivos y responsabilidades

| Ruta | Responsabilidad |
| --- | --- |
| `src/4-norm/norm_4.m` | Norma 4 de un vector. |
| `src/4-norm/power_norm4_single_start.m` | Iteración desde un único vector inicial; recibe el handle de norma. |
| `src/4-norm/compute_induced_norm.m` | Validación de Q y del contrato del handle, inicializaciones, multinicio y selección del mejor resultado. |
| `scripts/example_norm4_2x2.m` | Ejemplo reproducible con matriz Haar 2 por 2, semilla y comparación angular. |
| `docs/induced_norm4_algorithm.md` | Formulación, algoritmo, criterios de parada y limitaciones. |
| `README.md`, `docs/worklog.md`, `AGENTS.md` | Actualizar alcance, instrucciones de ejecución y registro del trabajo. |

Preserva `scripts/example_2x2.m` y su propósito geométrico. No duplicar ni reemplazar innecesariamente las funciones existentes. Si un archivo ya contiene una implementación compatible, extenderlo en lugar de crear una versión paralela.

## Qué calcula cada función

Para un vector real x:

$$
\operatorname{norm\_4}(x)=\|x\|_4
=\left(\sum_{j=1}^n |x_j|^4\right)^{1/4}.
$$

Para una matriz Q:

$$
\|Q\|_4
=\max_{\|x\|_4=1}\|Qx\|_4.
$$

El segundo problema necesita una búsqueda numérica. No confundirlo con la norma 4 de las entradas de Q, ni con su norma de Frobenius, norma 2 o radio espectral. Una matriz ortogonal tiene norma 2 igual a uno, pero su norma inducida 4 puede ser mayor que uno.

Q permanece fija durante todo el multinicio. Haar se utiliza para generar la matriz antes de llamar al estimador; no se vuelve a generar Q dentro de la búsqueda.

## Contrato de `@norm_4` y uso de handles

Implementar:

```matlab
function value = norm_4(x)
```

- Aceptar vectores fila o columna reales, finitos, no vacíos, en doble precisión; también aceptar un escalar como vector de dimensión 1.
- Devolver un escalar real no negativo; para el vector cero, devolver cero.
- Rechazar matrices que no sean vectores, entradas complejas y entradas no finitas.
- Para evitar overflow o underflow innecesario, usar escalamiento: si `s = max(abs(x))` es positivo, calcular `s * sum((abs(x)/s).^4)^(1/4)`; si s es cero, devolver cero.

En `compute_induced_norm` y `power_norm4_single_start`, nombrar el segundo argumento `normHandle` y comprobar que es un `function_handle`. Todas las normalizaciones, evaluaciones de amplitud y medidas vectoriales indicadas en este documento deben usar `normHandle(vector)`; no llamar internamente a `norm_4` por su nombre ni repetir fórmulas de la norma en cada paso.

El handle debe implementar matemáticamente la norma vectorial 4. Esta etapa no ofrece un optimizador para normas arbitrarias: las potencias cúbicas y raíces cúbicas del algoritmo son específicas de p = 4. Explicar esta restricción en el help de las funciones. Aceptar handles equivalentes, por ejemplo `@(x) norm_4(x)`, sin inspeccionar su nombre mediante `func2str` para decidir si son válidos.

Validar que las evaluaciones efectivamente utilizadas del handle devuelvan escalares reales finitos y que sean positivas para vectores no nulos. La validación de tipo y de salida no demuestra por sí sola que un handle arbitrario represente la norma 4; documentar el contrato matemático. No sustituir silenciosamente un handle inválido por otra función.

## Opciones e inicializaciones

Usar una estructura opcional de opciones con valores predeterminados explícitos:

| Campo | Valor inicial propuesto | Significado |
| --- | ---: | --- |
| `NumRandomStarts` | 50 | Número de inicios aleatorios adicionales a los estructurados. |
| `MaxIterations` | 1000 | Máximo de iteraciones por inicio. |
| `NormTolerance` | 1e-10 | Tolerancia para el cambio relativo de la estimación. |
| `StationarityTolerance` | 1e-8 | Tolerancia del residuo relativo de estacionariedad. |
| `FeasibilityTolerance` | 1e-12 | Tolerancia para la norma unitaria de x. |
| `StoreHistory` | false | Guardar o no la trayectoria de cada inicio. |

Permitir sobrescribir campos sin exigir una estructura completa. Rechazar campos desconocidos, cantidades inválidas y tolerancias no positivas o no finitas. `NumRandomStarts` puede ser cero; `MaxIterations` debe ser entero positivo. Estos valores son ajustes iniciales, no una garantía de precisión global.

Validar que Q sea real, finita, cuadrada, no vacía y de doble precisión. En esta primera implementación, exigir SO(n): verificar ortogonalidad con `norm(Q'*Q-eye(n),'fro')` y determinante cercano a +1. Usar como tolerancia de referencia `1e-12*max(1,n)` para estas dos comprobaciones. Informar un error descriptivo si no se cumple el contrato. No proyectar ni modificar Q para hacerla ortogonal.

El uso de la norma matricial de Frobenius en esta comprobación de Q no reemplaza el handle vectorial de la optimización.

Generar dos grupos de inicios:

1. **Estructurados:** para cada i de 1 a n, tomar `z = Q(i,:).'`, que corresponde a Q transpuesta por el vector canónico e_i, y normalizar con `normHandle`.
2. **Aleatorios:** reutilizar `generate_unit_l4_vectors(n,NumRandomStarts,'random')` si ya existe. Después, normalizar cada columna mediante el handle recibido. Si no se solicita ningún inicio aleatorio, no llamar al generador con N = 0 si este no lo admite.

La inicialización estructurada satisface:

$$
x_i^{(0)}=\frac{Q^T e_i}{\|Q^T e_i\|_4},
\qquad
\|Qx_i^{(0)}\|_4=\frac{1}{\|Q^T e_i\|_4}.
$$

No es necesario que los inicios aleatorios sean uniformes respecto al área de superficie. Nunca cambiar signos para forzar x a ser no negativo: las Q Haar generalmente contienen entradas positivas y negativas.

No llamar a `rng`, `rng('shuffle')` ni reiniciar semillas dentro de ninguna función de `src/`. La función consumirá el estado aleatorio establecido por el script o por el llamador.

## Iteración de potencia generalizada

Interfaz interna sugerida:

```matlab
[value, xHat, runInfo] = power_norm4_single_start(Q, normHandle, x0, options);
```

1. Convertir x0 en columna, comprobar dimensión n y normalizarlo con el handle.
2. Evaluar y registrar su amplitud inicial `normHandle(Q*x)`.
3. Repetir hasta convergencia o hasta `MaxIterations`:

```matlab
y = Q*x;
v = Q'*(y.^3);
u = sign(v).*abs(v).^(1/3);
xNew = u / normHandle(u);
valueNew = normHandle(Q*xNew);
```

El producto correcto es Q por x; los vectores son columnas. `Q'*(y.^3)` eleva las componentes de y al cubo, no la matriz Q. La actualización conserva signos y realiza la raíz cúbica real. No usar `v.^(1/3)` sobre componentes negativas, pues puede producir valores complejos.

Para una matriz Q ortogonal y un x unitario, v y u no pueden ser nulos en aritmética exacta. Si aparece una salida nula o no finita, identificar la ejecución como fallo numérico; no dividir entre cero ni generar un nuevo inicio ocultamente.

En aritmética exacta, las amplitudes de esta iteración no disminuyen. Registrar descensos superiores a una tolerancia de redondeo razonable como anomalías numéricas; no ocultarlos declarando convergencia. Conservar el mejor candidato factible ya evaluado, incluida la inicialización.

No normalizar y antes de medir `normHandle(y)`, pues se perdería la amplificación que se busca maximizar. No introducir un paso de gradiente, un parámetro de aprendizaje o un solver externo dentro de esta primera versión de la potencia generalizada.

## Relación con Lagrange y criterios de parada

Definir:

$$
f(x)=\|Qx\|_4^4.
$$

En un punto estacionario unitario:

$$
Q^T(Qx)^{\odot3}=f(x)x^{\odot3}.
$$

Esta relación explica la actualización, pero no certifica un máximo global. La notación de potencias con el símbolo de producto por componentes indica operaciones componente a componente.

Después de cada actualización, comprobar en el **mismo xNew**:

```matlab
yNew = Q*xNew;
vNew = Q'*(yNew.^3);
lambdaNew = valueNew^4;
xCube = xNew.^3;
residual = vNew - lambdaNew*xCube;

denominator = max(normHandle(vNew) + ...
    abs(lambdaNew)*normHandle(xCube), realmin);
relativeResidual = normHandle(residual) / denominator;
feasibilityError = abs(normHandle(xNew) - 1);
relativeChange = abs(valueNew - valueOld) / ...
    max([1, abs(valueNew), abs(valueOld)]);
```

Declarar convergencia solamente cuando se cumplan conjuntamente las tolerancias de cambio de norma, residuo y factibilidad. Un cambio pequeño de la función no basta si el residuo sigue siendo grande.

Una ejecución que alcanza `MaxIterations` debe identificarse como tal, sin etiquetarla como convergente. Evaluar los diagnósticos del vector que realmente se devuelve; no mezclar el residuo del último iterado con un mejor vector guardado en otra iteración.

## Selección y resultados del multinicio

La primera salida es:

$$
\texttt{qNorm4}=\max_r\|Q\widehat x_r\|_4.
$$

Seleccionar el mejor candidato factible y finito evaluado entre todos los inicios, incluso si su ejecución terminó por límite de iteraciones. Conservar su estado de convergencia y emitir una advertencia descriptiva si el candidato seleccionado no convergió. No reemplazarlo silenciosamente por uno convergente de menor amplitud.

Si no existe ningún candidato válido, producir un error. No devolver NaN, cero ni un valor truncado a las cotas teóricas para ocultar un fallo.

`xBest` debe tener tamaño n por 1 y satisfacer la norma unitaria dentro de tolerancia. Debe verificarse que `qNorm4` coincida con `normHandle(Q*xBest)`.

`info` debe contener al menos:

- `method`: potencia generalizada con multinicio.
- `isEstimate = true`.
- `numStarts`, `numRandomStarts` y `bestStartIndex`.
- `bestConverged`, `bestTerminationReason` y `bestIterations`.
- `bestStationarityResidual` y `bestFeasibilityError`.
- `startValues`, `startConverged` y las razones de terminación por inicio.
- `theoreticalLowerBound = 1` y `theoreticalUpperBound = n^(1/4)`.
- Historiales solo cuando `StoreHistory` sea verdadero.

Para Q ortogonal se cumple `1 <= norm(Q,4) <= n^(1/4)` en sentido matemático de norma inducida; esta expresión no implica que MATLAB acepte `norm(Q,4)`. Verificar las cotas con tolerancia y reportar anomalías, sin forzar los resultados a entrar en el intervalo.

Llamar al resultado **estimación**. Ni el multinicio, ni un residuo pequeño, ni la coincidencia entre inicios prueban por sí solos que se encontró el máximo global. No declarar que xBest es necesariamente x*.

## Script `scripts/example_norm4_2x2.m`

Crear un script reproducible con este flujo, ajustando solamente la gestión de rutas y la presentación de resultados:

```matlab
% Localizar la raíz a partir de la ubicación del script y añadir src al path.
rng(42, 'twister');

n = 2;
Q = haar_so(n);
normHandle = @norm_4;

options = struct();
options.NumRandomStarts = 50;
options.MaxIterations = 1000;
options.StoreHistory = true;

[qNorm4, xBest, info] = compute_induced_norm(Q, normHandle, options);

disp(Q);
fprintf('Estimación de la norma inducida 4: %.12f\n', qNorm4);
fprintf('Norma 4 del vector candidato: %.12f\n', normHandle(xBest));
fprintf('Amplitud de Q*xBest: %.12f\n', normHandle(Q*xBest));
disp(xBest);
disp(info);
```

La semilla se define antes de generar tanto Q como los inicios aleatorios. No restaurarla o reiniciarla entre esas llamadas. El script debe funcionar desde la raíz usando:

```matlab
run(fullfile('scripts', 'example_norm4_2x2.m'))
```

Además, realizar una comparación angular independiente para esta misma Q:

1. Generar, por ejemplo, 10000 puntos ordenados mediante `generate_unit_l4_vectors(2,N,'angular')`.
2. Normalizarlos y evaluar cada amplitud mediante el handle, usando un bucle o `arrayfun` sobre las columnas; el contrato del handle es vectorial.
3. Encontrar los máximos locales de la muestra usando vecinos periódicos, sin `findpeaks` ni toolboxes. Manejar los casos planos o de empate sin errores.
4. Refinar los intervalos de esos candidatos con `fminbnd` aplicado al negativo de la amplitud angular. Parametrizar `x(theta) = [cos(theta);sin(theta)] / normHandle([cos(theta);sin(theta)])`. La función es periódica y puede evaluarse fuera de [0,2*pi) para refinar candidatos en los extremos.
5. Conservar también el mejor punto de la malla. Mostrar la diferencia entre la estimación de potencia y la referencia angular aproximada. No reemplazar con ella el resultado de la función general.
6. Dibujar una figura con la amplitud `normHandle(Q*x(theta))` frente al ángulo y marcar la amplitud del mejor candidato de potencia. Marcar su dirección con `atan2(xBest(2),xBest(1))`, ajustada a [0,2*pi).

La búsqueda angular solo pertenece a la validación de dimensión 2; no debe invocarse dentro de la función general ni presentarse como certificación exacta. Los puntos finitos de una malla pueden omitir picos, aunque el refinamiento reduce el error.

Usar figuras locales, etiquetas en español y rutas relativas. No usar `clear all`, `close all`, `savepath`, cambios de directorio ni propiedades gráficas globales. La exportación de gráficas puede ser opcional y debe ir a una carpeta de resultados ignorada por Git.

## Verificación requerida

Realizar comprobaciones pequeñas y significativas, sin añadir un framework ni dependencias externas:

| Caso | Resultado esperado |
| --- | --- |
| `norm_4([1;2])` | `17^(1/4)`. |
| Vector cero y versiones fila/columna | Cero y concordancia entre orientaciones. |
| `compute_induced_norm(eye(n),@norm_4)` | Aproximadamente 1; cubrir n = 1, 2 y 5. |
| Permutación con signos y determinante +1 | Aproximadamente 1. |
| `Q = [1,-1;1,1]/sqrt(2)` | Aproximadamente `2^(1/4)`, con tolerancia inicial de 1e-8. |
| Q Haar de dimensiones 2, 5 y 10 | Factibilidad, finitud, consistencia escalar/vector, diagnósticos y cotas ortogonales. |
| `@(x) norm_4(x)` | Misma interfaz y resultados que con `@norm_4`, al repetir el mismo estado aleatorio. |
| Entradas y opciones inválidas | Errores descriptivos. |
| `MaxIterations` insuficiente para satisfacer criterios | Estado de no convergencia correctamente comunicado. |

Comprobar reproducibilidad reiniciando la semilla **desde el llamador** y repitiendo la misma secuencia completa de generación de Q y evaluación. Comprobar también que dos llamadas consecutivas sin reiniciar la semilla consuman el estado aleatorio, sin afirmar que necesariamente producirán valores finales distintos.

La identidad es un caso especial en el que todos los vectores unitarios son maximizadores. No exigir un xBest único ni igualdad de vectores entre métodos cuando existan múltiples maximizadores.

Si MATLAB no está disponible, comunicar que la implementación no pudo ejecutarse en MATLAB. No presentar inspección estática o ejecución en otro lenguaje como prueba de ejecución MATLAB.

## Documentación y cierre de la implementación

Actualizar el README con la llamada mínima, las salidas opcionales y la ejecución del nuevo script. Mantener su declaración original sobre la exploración de los componentes del número de condición bajo la norma 4.

En `docs/induced_norm4_algorithm.md`, explicar la función del handle, la iteración, las inicializaciones estructuradas y aleatorias, las tolerancias y el carácter aproximado del resultado. Citar como referencia técnica:

Nicholas J. Higham (1992), *Estimating the matrix p-norm*, Numerische Mathematik 62, 539–555: https://nhigham.com/wp-content/uploads/2023/09/high92e.pdf

Las garantías de convergencia global que requieren matrices no negativas no deben trasladarse a las Q Haar, que generalmente contienen signos mezclados. Esta implementación es una estimación con multinicio, no un algoritmo certificado de optimización global.

Registrar en `docs/worklog.md` los archivos modificados y las comprobaciones efectivamente realizadas. No hacer commits, crear repositorios remotos ni publicar sin una solicitud del usuario.

Al terminar, informar cómo ejecutar el script, qué devuelve la función, qué se verificó y cualquier ejecución que no haya convergido.
