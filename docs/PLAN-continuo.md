# Plan: mapas rectangulares, formato `.obs` y modo off-lattice (continuo)

Estado: completado el 2026-10-22 (fases F0–F8). Versión objetivo: 2.2.0 (cambios aditivos; el modo red cuadrado
debe seguir produciendo el mismo PGM bit a bit con la misma semilla).

Referencias del proyecto consumidor (hyperGalton2D):
- formato: `~/granular/hyperGalton2D/docs/formats/obstacles-format.pdf` (versión 1, enmendada el
  2026-10-22, congelada);
- plan: `~/granular/hyperGalton2D/docs/plan/plan.pdf`. Hito M0a: extensión al continuo y exportador
  `.obs`, noviembre de 2026. Decisión D1: hiperuniformidad **isótropa** en banda periódica en x.

## Objetivo

1. **Cajas rectangulares** en los dos modos: red de *N*ₓ × *N*ᵧ sitios y caja continua *L*ₓ × *L*ᵧ.
   La campaña varía *L*ₓ/*a* de 20 a 320 con *L*ᵧ casi fijo, así que la caja no cuadrada es el caso
   normal.
2. **Salida en el formato `.obs` versión 1** en los dos modos. La visualización en PGM pasa a un script
   de Python que convierte `.obs` → PGM.
3. **Modo continuo** (`mode: continuum`): *N* puntos con posiciones reales. En el continuo E = 0 se
   alcanza en el régimen desordenado, por eso el algoritmo principal es minimización por gradiente
   (L-BFGS), como en la bibliografía (Uche, Stillinger y Torquato 2004; Zhang, Stillinger y
   Torquato 2015).

## Física y convenciones (según la especificación `.obs`)

### Vectores de onda y máscara

Para una caja [0, *L*ₓ) × [0, *L*ᵧ) periódica, **k** = 2π (*n*ₓ/*L*ₓ, *n*ᵧ/*L*ᵧ), (*n*ₓ, *n*ᵧ) ∈ ℤ².
La máscara es 0 < |**k**| ≤ *K*, con |**k**| la norma euclídea del vector físico.

**La máscara es un círculo en *k* físico, nunca en el espacio de índices.** Con *L*ₓ ≠ *L*ᵧ, la condición
*n*ₓ² + *n*ᵧ² ≤ *N*_K² es una elipse en *k* y da hiperuniformidad direccional, que D1 deja fuera de la
campaña. En índices, el círculo es una elipse de semiejes *K L*ₓ/2π y *K L*ᵧ/2π:
(*n*ₓ/*L*ₓ)² + (*n*ᵧ/*L*ᵧ)² ≤ (*K*/2π)².

- **Red:** *L*ₓ = *N*ₓ*a*, *L*ᵧ = *N*ᵧ*a*. La selección se hace con la clave entera exacta
  *q* = *n*ₓ²*N*ᵧ² + *n*ᵧ²*N*ₓ² (cabe en 64 bits para *N* ≤ 46340), que agrupa las capas sin errores de
  redondeo. Límite de Nyquist por eje: |*n*ₓ| < *N*ₓ/2, |*n*ᵧ| < *N*ᵧ/2.
- **Continuo:** sin límite de Nyquist; capas ordenadas por |**k**|² en double.
- **Test obligatorio:** para *L*ₓ ≠ *L*ᵧ, todos los modos de la máscara cumplen |**k**| ≤ *K* en *k* físico,
  ningún modo fuera la cumple, y el conjunto es invariante bajo *k*ₓ → −*k*ₓ y *k*ᵧ → −*k*ᵧ. Un test que
  use radios en índices debe fallar.

### Unidades de *K* y κ

- `K` en el `.obs`: **unidades reducidas absolutas, *d*⁻¹, con el 2π incluido.** No se normaliza por
  ninguna longitud de caja.
- κ = *K*/(2π√ρ_O) = *K a*/2π, con *a* = ρ_O^{−1/2}: el grupo adimensional para comparar entre formas de
  caja, densidades y generadores. Se escribe como `kappa`.
- Cuando el *K* sale de un χ objetivo, se escribe el punto medio entre la última capa incluida y la
  primera excluida. Así un lector que recalcule la máscara a partir de `K`, `Lx` y `Ly` obtiene el mismo
  conjunto aunque haya redondeo.

### χ y el conteo de modos

- `n_constrained_modes` = *M*(*K*): modos de la máscara en un semiespacio (uno de cada par ±**k**), sin
  **k** = 0. Es la cantidad sin ambigüedad y es **obligatoria** cuando hay `K`.
- `chi` es un dato derivado y va siempre con `chi_convention`:

| `chi_convention` | Fórmula | Generador |
|---|---|---|
| `point_pattern` | χ = *M*(*K*)/[2(*N* − 1)] | modo continuo |
| `binary_lattice` | χ = *M*(*K*)/(*n*_min − 1), *n*_min = min(ocupados, vacíos) | modo red (= 2.x) |

- **La validación cruzada red ↔ continuo se aparea por *K* (o κ) e igual ρ_O, nunca por χ.** Los χ que
  etiquetan la campaña (0,1, 0,3 y 0,45) son `point_pattern`. El control en red se genera con el
  **mismo *K*** y reporta el χ de su propia convención. Para densidad de red < 1/2, χ_binary ≈ 2 χ_point
  a igual *K* y *N*.
- Relación aproximada en el continuo (conteo de modos por área): χ_point ≈ πκ²/4, κ ≈ 2√(χ/π).
  Entonces χ = 0,1, 0,3 y 0,45 dan κ ≈ 0,357, 0,618 y 0,757.
- **Error en la documentación de la 2.x:** la ecuación del README (`docs/eq/chi.svg`) dice
  χ = *M*_K/[2(*n* − 1)], pero el código y la viñeta del README calculan *M*_K/(*n* − 1), que es
  `binary_lattice`. Se corrige la ecuación, no el código, para que la 2.x siga siendo reproducible.
  También se revisa el aviso "χ > 1/2", que se escribió suponiendo la convención de Torquato.

### Residuo

`residual_S` = ⟨*S*⟩_mask normalizado a 1 para un patrón sin correlación del mismo tipo: ⟨*S*⟩_mask en el
continuo (*S* → 1 para Poisson) y ⟨*S*⟩_mask/(1 − ρ) en la red (= `mask_residual` de 2.x). Se calcula
sobre las coordenadas escritas.

### Gradiente y fases (continuo)

Con ρ̂(**k**) = Σⱼ e^{−i**k**·**r**ⱼ}: ∂|ρ̂|²/∂**r**ⱼ = 2**k** Im[ρ̂\*(**k**) e^{−i**k**·**r**ⱼ}]. Una evaluación
de E y ∇E cuesta O(*N*·*M*). Las fases se obtienen por recurrencia,
e^{−i2π*n*x/*L*ₓ} = (e^{−i2π*x*/*L*ₓ})ⁿ, con un *n* máximo distinto por eje (*K L*ₓ/2π y *K L*ᵧ/2π).

## Requisitos del formato `.obs`

- **Unidades reducidas** (*d* = 1).
- **Sin solapamiento, como error duro:** distancia mínima ≥ 2*R*.
  - **Red:** el sitio (*i*, *j*) va en (*i a*, *j a*). Hace falta *a* ≥ 2*R*, porque dos sitios vecinos
    ocupados quedan a distancia *a*.
  - **Continuo:** los patrones *stealthy* de referencia tienen *d*_min/*a* ≈ 0,485 según el plan de
    hyperGalton2D, pero nada lo garantiza en cada realización. Se agrega un **término de núcleo duro**:
    E_core = Σ_{i<j} (σ − *r*ᵢⱼ)² θ(σ − *r*ᵢⱼ), con σ = `d_min_target` (por defecto 2*R*) y lista de celdas.
    Vale 0 en cualquier patrón que ya cumple la regla, así que solo modifica las realizaciones que
    fallarían. Se reporta cuántos pares tuvo que separar, porque el control RSA se ajusta a la
    distribución de distancias al vecino más próximo del patrón *stealthy*.
- **Contorno:** en producción es `periodic_x`. El patrón se genera en el toro y se declara `periodic_x`;
  `d_min` se mide con imagen mínima solo en x.
- **Encabezado coherente con lo escrito:** `rho_O`, `d_min`, el orden, el dominio, `residual_S` y
  `sha256_body` se calculan a partir de los números tal como quedan escritos. Todos los reales se
  escriben con la representación más corta que reproduce el double exacto (`std::format("{}")`).
- **Dominio:** 0 ≤ *x* < *L*ₓ estricto. Al plegar, `x + Lx` puede redondear a `Lx`; en ese caso se usa 0.
- **Orden:** *y* ascendente y después *x*.
- **SHA-256 del cuerpo**, desde la primera línea sin `#` inclusive (verificado con el ejemplo).
  Implementación propia en C++, con vectores NIST y comparación con `hashlib`.
- **Claves obligatorias para `class: stealthy`:** `K`, `n_constrained_modes`, `chi_convention` y, si se
  escribe, `chi`.
- **Claves recomendadas:** `chi`, `kappa`, `residual_S`, `lattice_constant` (en modo red; la regla 8 la
  necesita); `alpha_numvar` lo puede agregar el script de análisis, porque el checksum cubre solo el
  cuerpo.
- **Nombre de archivo:** `obs_file: auto` → `<class>_rho<value>_Lx<value>_Ly<value>_arr<index>.obs`.
  Opcional: `realizations: n` con semillas consecutivas.
- **Archivo `.sk` acompañante** (opcional): dos columnas *k* (en *d*⁻¹) y *S*(*k*).
- **El ejemplo de la especificación incumple la regla 4:** declara `rho_O: 0.288675135`, pero
  64/(16 × 13.8564065) = 0.2886751338, una diferencia de 1,2·10⁻⁹ > 10⁻⁹. Sigue igual después de la
  enmienda. Hay que corregirlo del lado de hyperGalton2D.

## Parámetros

### Máscara (una de las tres formas)

| Clave | Modo | Significado |
|---|---|---|
| `kappa` | los dos | κ = *K a*/2π; **la forma recomendada para la campaña**, porque aparea red y continuo |
| `chi` | los dos | χ objetivo en la convención del modo (`point_pattern` o `binary_lattice`); el programa avisa qué convención usó |
| `K` | solo red cuadrada | significado de 2.x (unidades de 2π/*N*), por compatibilidad; se rechaza con *N*ₓ ≠ *N*ᵧ y en el continuo |

Si se dan varias, gana `kappa`, después `chi`, y el programa avisa.

### Geometría y obstáculos

| Clave | Por defecto | Significado |
|---|---|---|
| `mode` | `lattice` | `lattice` o `continuum` (F0 ✅) |
| `Nx`, `Ny` | `N` | lados de la red en sitios |
| `lattice_constant` | 1 | *a*_red, en *d* |
| `Lx`, `Ly`, `rho_O` | — | caja y densidad en unidades reducidas (continuo; en red, alternativa a `Nx`, `Ny`, `lattice_constant`) |
| `R` | — | radio de los obstáculos |
| `d_min_target` | 2*R* | distancia que impone el núcleo duro (continuo) |
| `boundary` | `periodic_x` | lo que se declara en el `.obs` |

**Control en red a igual ρ_O y caja:** dados `Lx`, `Ly`, `rho_O` y la fracción de ocupación `density`,
el programa deriva *a*_red = √(`density`/ρ_O), *N*ₓ = round(*L*ₓ/*a*_red) y *N*ᵧ = round(*L*ᵧ/*a*_red), y
avisa del *L*ₓ, *L*ᵧ y ρ_O efectivos si hay que redondear.

**Límite de Nyquist en red:** *K* < π/*a*_red por eje, es decir κ < 1/(2√`density`) (≈ 0,91 para
`density` 0,3). Los tres χ de la campaña entran.

### Optimización y salida

| Clave | Por defecto | Significado |
|---|---|---|
| `minimizer` | `lbfgs` | `lbfgs`, `anneal` o `anneal+lbfgs` (continuo) |
| `lbfgs_memory`, `max_iterations`, `grad_tol` | 10, 10⁵, — | control de L-BFGS |
| `residual_tol` | ~10⁻²⁰ | parar cuando `residual_S` ≤ tol y no hay solapamientos |
| `mc_step` | 0,1 · *a* | paso inicial del MC continuo (adaptativo) |
| `init_file` | — | configuración inicial desde un `.obs` |
| `obs_file` | — | salida `.obs`; `auto` = nombre según la especificación |
| `arr_index`, `realizations` | 0, 1 | índice de realización; número de realizaciones |
| `sk_kmax` | 3*K* | alcance de *S*(*k*) en el continuo |

El PGM (`output_file`) se mantiene en modo red durante la serie 2.x por compatibilidad.

## Fases

**F0 — Refactor sin cambios de comportamiento.** ✅ Hecha (2026-10-22): PGM, `sk.dat` y
`annealing.dat` idénticos a 2.1.0 con 1 y 4 hilos (96 × 96 y 200 × 200); tests con sanitizers sin errores.

**F1 — Red rectangular *N*ₓ × *N*ᵧ con máscara circular en *k* físico.** ✅ Hecha (2026-10-22). Residuo de la red a igual *K* (128 × 128, densidad 0,3, `weight_power` 2), con 0,55 M / 12,8 M / 113 M propuestas: χ_point 0,1 → 0,0084 / 0,0065 / —; 0,3 → 0,078 / 0,064 / 0,054; 0,45 → 0,16 / 0,133 / 0,119.
- Tablas de fases por eje, Nyquist por eje, selección por la clave entera *q*, `choose_mask` por `kappa`
  o `chi`, DFT de *S*(*k*) rectangular con *k* en *d*⁻¹, PGM de ancho *N*ₓ y alto *N*ᵧ.
- `chi_convention` explícita en el código (`binary_lattice`), corrección de la ecuación del README y
  revisión del aviso χ > 1/2.
- Criterios: con *N*ₓ = *N*ᵧ y `K` de 2.x, salida idéntica a 2.1.0; con *N*ₓ ≠ *N*ᵧ, máscara circular en
  *k* físico (test de arriba), energía contra la versión directa, Parseval.

**F2 — Salida `.obs` y herramientas de Python.** ✅ Hecha (2026-10-22). El validador confirma que el ejemplo de la especificación incumple la regla 4. Con esto hyperGalton2D ya recibe mapas de red.
- C++: escritor `.obs` (claves obligatorias y de diagnóstico, formato exacto, orden, `d_min` con
  `periodic_x` o `walls_x`, SHA-256), nombre automático, `.sk`, `realizations`, derivación de la red a
  partir de `Lx`, `Ly`, `rho_O` y `density`.
- `scripts/obs_tools.py`:
  - `validate`: las ocho reglas de la especificación, más las claves de `class: stealthy`, más el
    recálculo de la máscara a partir de `K`, `Lx` y `Ly` (debe dar `n_constrained_modes`), más el aviso
    de hueco γ < 0;
  - `topgm`: en red, un píxel por sitio (debe reproducir el PGM de hyper2D); en el continuo, discos de
    radio *R* a `--px-per-d` píxeles por *d*, plegados en x si el contorno es periódico; P2 o P5;
    opción de orientación.
- Tests: vectores NIST; el `.obs` de hyper2D pasa `validate`; `.obs` → PGM igual al PGM de hyper2D.

**F3 — Núcleo continuo.** ✅ Hecha (2026-10-22): `include/continuum.hpp`, `src/continuum.cpp`, solo como
biblioteca (todavía no se usa desde el archivo de parámetros). Tests: fases por recurrencia con error
1,0·10⁻¹² hasta *n* = 600 (menor que `std::polar` en double, 1,5·10⁻¹²), E contra la versión directa,
gradientes espectral y de núcleo duro contra diferencias finitas (errores 2·10⁻⁷), suma nula de los
gradientes, máscara circular en caja 60 × 17, χ_point por capas, RSA. Costo de E + ∇E con un hilo
(caja de aspecto 4, `weight_power` 2): *N* = 1000 → 2,3 ms; 4900 → 43 ms (κ = 0,618) y 61 ms (κ = 0,757);
10 000 → 166 ms. Unos 1,5 ns por par punto·modo, ∝ *N*².
- `ContinuumConfig`, máscara circular en *k* físico con caja rectangular, fases por recurrencia por eje,
  ρ̂, E, E_w, ∇E y versiones directas. Término de núcleo duro con lista de celdas y su gradiente.
- Tests: recurrencia contra una referencia en `long double` (error < 4·*n*·2π·ε, el límite de escribir el ángulo en double; no peor que `std::polar`), E contra la directa, ∇E contra diferencias finitas
  (las dos energías), Σⱼ ∇ⱼE = 0, máscara circular.

**F4 — Minimizador y driver.** ✅ Hecha (2026-10-22). L-BFGS propio (Wolfe fuerte, Nocedal y Wright) en
lugar de liblbfgs, para no agregar una dependencia ni descargar código; mismas iteraciones que L-BFGS-B
de SciPy en los problemas de prueba. La parada es sobre el residuo sin ponderar. 4900 puntos en
140 × 35, 4 hilos, *R* = 0,1: χ_point 0,3 → 1853 iteraciones / 57 s (`weight_power` 0) y 253 / 15 s
(`weight_power` 2), residuo ≤ 10⁻²⁰; χ_point 0,45 → 1371 / 29 s y 564 / 13 s. *S* < 6·10⁻²⁰ dentro de la
máscara en todas las direcciones, α = 0,99 y 1,04. Con χ = 0,3 el núcleo duro actúa (pares en σ); con
χ = 0,45 el patrón mantiene ~0,65 *a* por sí solo. Pendiente para F5: `init_file`.
- L-BFGS detrás de `minimize(f, x, opciones, callback)`; propuesta: incluir liblbfgs (MIT).
- Parada por `residual_tol` sin solapamientos, `grad_tol`, `max_iterations`. Reporte: E, `residual_S`,
  residuo |**k**| ≤ *K*/4, máx *S*(*k*) en la máscara, `d_min`, pares separados por el núcleo duro.
- Tests: cuadrática, Rosenbrock; *N* = 200, χ_point = 0,2 debe llegar a residuo < 10⁻²⁰ en segundos.

**F5 — `init_file` y MC continuo.** ✅ Hecha (2026-10-22). `init_file` arranca de un `.obs` (red o continuo),
con verificación de formato, conteo y `sha256_body`, y toma su caja; un `Lx` distinto se rechaza. MC de
desplazamientos de un punto (`minimizer: anneal` y `anneal+lbfgs`): ρ̂ incremental (deriva 10⁻¹⁵), paso
adaptativo, núcleo duro como restricción, mismo resultado con 1 y 4 hilos. Solo MC: residuo 7·10⁻⁵ con
400 puntos en 0,3 s. Como primera etapa redujo L-BFGS de 274 a 46 iteraciones con 400 puntos, pero con
4900 (χ = 0,3, `weight_power` 2) tardó 45 s contra 15 s de L-BFGS solo: `lbfgs` sigue por defecto.

**F6 — Paralelismo y rendimiento.** ✅ Hecha (2026-10-22). Solo C++ estándar (`std::thread`, sin OpenMP ni
directivas). Máscara en tramos de igual *n*ₓ; ρ̂ en 32 grupos fijos de puntos con sumas parciales sumadas en
orden; gradiente factorizado por tramo y en bloques de 32 puntos con el punto como bucle interno. Resultado
idéntico bit a bit con cualquier cantidad de hilos. E + ∇E con un hilo: 0,7–0,95 ns por par punto·modo (antes
2,9); escala ×4,5 con 6 núcleos. *N* = 4900: 13,7 / 4,1 / 3,2 ms (1 / 4 / 6 hilos); 20 000: 171 / 52 / 39 ms;
40 000: 640 / 182 / 142 ms. Corridas completas (4900 puntos, `weight_power` 2, 6 hilos, 4 semillas):
χ = 0,3 → 1,0–3,8 s; χ = 0,45 → 3,1–3,4 s (antes 15 y 13 s con 4 hilos). Corregido de paso: la continuación
de L-BFGS se cortaba si había un par justo en σ. El MC continuo conserva su núcleo anterior.
- ρ̂ en paralelo sobre modos; ∇E sobre partículas. Tiempos hasta el *N* más grande de la campaña
  (*N* = (*L*ₓ/*a*)(*L*ᵧ/*a*), con *L*ₓ/*a* = 320).

**F7 — Script de análisis.** ✅ Hecha (2026-10-22): `scripts/analyze_obs.py` (y `test_hyper2D.py` le pasa los
`.obs`). ρ̂ por producto de matrices de las tablas por eje; *S*(*k*) radial y por sectores; σ²(*R*) espectral
exacta en el toro con cola de Parseval, y conteo directo (único método si *y* no es periódico); α; *g*(*r*);
`--write-alpha` sin tocar el cuerpo. Continuo χ = 0,3: ⟨*S*⟩_mask = 4,9·10⁻²¹ igual en los tres sectores,
α = 0,998, conteo directo dentro del 1,7 %. Control en red 160 × 40: ⟨*S*⟩ = 3,8·10⁻³, α = 0,985. Tests en
`ctest` (`tests/test_analysis.py`).
- `test_hyper2D.py` lee `.obs`; cajas rectangulares; `periodic_x`; *k* en *d*⁻¹.
- σ²(*R*) espectral exacta con cola cerrada por Parseval; conteo directo como control; *g*(*r*);
  opción para escribir `alpha_numvar`.

**F8 — Documentación y versión.** ✅ Hecha (2026-10-22): versión 2.2.0 (CMake, CHANGELOG con la sección
[2.2.0] reorganizada), `examples/continuum.in` y `examples/lattice_control.in` apareados por κ, README
revisado (versión, requisitos con SciPy, tests, inicio rápido, tiempos del MC con los núcleos de F6,
referencias de Uche et al. 2004 y Nocedal y Wright), `hyper2D --help` con los modos. Remedido el MC como
primera etapa con los núcleos nuevos: no conviene (4900 puntos: 5,1–9,7 s contra 1,8–3,9 s de L-BFGS solo).
- README, ejemplos de parámetros (uno por modo, apareados por κ), figuras, CHANGELOG, 2.2.0.

## Validación física

- Continuo, χ_point = 0,1, 0,3, 0,45: *S*(*k*) = 0 en la máscara, α ≈ 1, hueco de exclusión en *g*(*r*).
- **Isotropía en cajas no cuadradas** (*L*ₓ/*L*ᵧ = 1, 2, 4 y el rango de la campaña): *S*(*k*) y σ²(*R*)
  sin dependencia de la dirección. Comparar con el atajo de índices, que debe mostrar anisotropía.
- **Red contra continuo a igual *K* (κ) e igual ρ_O**, nunca a igual χ. ✅ Hecha el 2026-10-30, 16 realizaciones
  por familia a κ = 0,618 y 0,757: [`validacion-cruzada.md`](validacion-cruzada.md). Las dos familias son
  hiperuniformes (α ≈ 0,996, H = *S*(0)/*S*(*k*_pico) ≤ 1,3·10⁻⁵ ≪ 10⁻⁴). A χ = 0,45 coinciden salvo por el
  residuo de la red, que queda en 0,0114, **por encima del umbral 0,011**. A χ = 0,3 **no coinciden** fuera de la
  máscara (la red tiene un pico *S* = 1,37 en 1,08 *K* y Λ un 4 % mayor), y no por el residuo: los patrones de red
  relajados a *S* = 0 en el continuo conservan esa estructura. Es la dependencia con el protocolo de los
  estados fundamentales degenerados. Controles de Λ: triangular 0,50837 (0,508347), cuadrada 0,45765 (0,4576).
- Con núcleo duro: residuo y fracción de realizaciones afectadas en función de 2*R*/*a*.

## Control en red *stealthy* (análisis del 2026-10-22)

Con densidad de ocupación ρ = 0,3, la red a igual *K* no es *stealthy* para χ_point = 0,3 y 0,45, y
enfriar 200 veces más lento solo baja el residuo un 30 %. **El residuo lo fija ρ, no el enfriamiento:**
con κ = 0,618 (χ_point = 0,3), ~4900 obstáculos y el mismo enfriamiento, `residual_S` ≈ 0,2 ρ:

| ρ | 0,30 | 0,15 | 0,08 | 0,04 | 0,02 | 0,01 |
|---|---|---|---|---|---|---|
| `residual_S` | 0,066 | 0,030 | 0,016 | 0,0076 | 0,0040 | 0,0019 |

- Mismo costo por propuesta (~2 min cada corrida): el número de modos depende de κ y *N*, no de ρ.
- Explicación: redondear un patrón de puntos a una red de espaciado *a*_red da *S* ≈ (*K a*_red)²/12 ∝ κ²ρ
  en el borde de la máscara; el recocido elimina la mayor parte.
- Sin exclusión, la distancia mínima es un espaciado de red (1,6 % de los puntos con un vecino a un sitio,
  ρ = 0,02). La opción **`min_distance`** (implementada) prohíbe pares a menos de σ: con σ hasta 0,5 *a* el
  residuo no cambia. En el `.obs`, σ = 2*R*/*a*_red garantiza la regla 7.
- α = 0,96–1,00 con ρ ≤ 0,04.

**Decidido (2026-10-22):** el control en red usa red fina, ρ ≤ 0,02, con exclusión 2*R*/*a*_red (implícita
al dar `R`). *a*_red = *a*√ρ.

χ_point = 0,45 (κ = 0,757), con exclusión 0,4 *a*: `residual_S` = 0,0113 con ρ = 0,02 y 0,0064 con ρ = 0,01
(≈ 0,6 ρ, tres veces más que con χ = 0,3); α = 0,95–0,99. El patrón mantiene por sí solo
*d*_min = 0,58–0,67 *a*, así que la exclusión no actúa. Para ~0,004 con χ = 0,45 hace falta ρ ≈ 0,006.

**Umbral aceptado (2026-10-22): `residual_S` ≤ 0,011** para el control en red. Con ρ = 0,02 cubre los tres χ
de la campaña (0,0040 con χ = 0,3; 0,011 con χ = 0,45).

**Revisión (2026-10-30):** con 16 realizaciones, χ = 0,45 da `residual_S` = 0,01142 ± 0,00003 (todas entre 0,0112
y 0,0116): el umbral 0,011 no se cumple con ρ = 0,02. Hay que subirlo a 0,012 o usar ρ ≈ 0,01 para χ = 0,45
([`validacion-cruzada.md`](validacion-cruzada.md)).

Pendiente: el costo de *S*(*k*) y de las tablas en redes de 10⁶–10⁷ sitios (la DFT de *S*(*k*) es
O(*N*ₓ*N*ᵧ(*N*ₓ + *N*ᵧ)); pasar a FFT si hace falta).

## Riesgos

- **El control en red a igual *K* puede no ser *stealthy*.** Prueba rápida con 2.x (red 128 × 128,
  densidad 0,3, enfriamiento corto): χ_binary 0,2, 0,6 y 0,9 (≈ χ_point 0,1, 0,3 y 0,45) dan
  `residual_S` ≈ 0,008, 0,078 y 0,16. Con χ_point = 0,3 y 0,45 la red queda lejos de E = 0. Hay que
  medirlo con enfriamiento largo y `weight_power` al principio de F1, y decidir con hyperGalton2D si el
  control en red cubre los tres χ o solo el más bajo.
- Costo O(*N*·*M*) ∝ χ*N*²: con *L*ₓ/*a* = 320 el *N* puede llegar a ~10⁴–10⁵ según *L*ᵧ/*a*. Por encima
  de ~2·10⁴ hará falta NUFFT (FINUFFT) o más hilos; se mide en F6.
- χ cerca de 1/2: L-BFGS puede estancarse; reinicios o `anneal+lbfgs`.
- Cajas muy alargadas: el espaciado en *k* del eje corto es grueso (2π/*L*ᵧ), así que el círculo se
  resuelve peor en esa dirección y χ se discretiza más.
- Un patrón generado en el toro y declarado `periodic_x`: correcto para la estadística de volumen, pero
  las paredes en y cortan un patrón correlacionado a través del borde.

## Decisiones

- Máscara: círculo en *k* físico; `K` en *d*⁻¹; entrada recomendada `kappa`.
- χ: dato derivado, siempre con `chi_convention`; apareamiento por *K*/κ.
- Optimizador: L-BFGS por defecto, MC continuo opcional; liblbfgs propuesta (se cierra en F4).
- Salida: `.obs` v1 en los dos modos; PGM de red conservado en 2.x; PGM del continuo con `obs_tools.py`.
- Versión: 2.2.0 si el modo red cuadrado queda idéntico.
