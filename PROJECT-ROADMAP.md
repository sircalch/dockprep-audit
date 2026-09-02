# DockPrep Audit — estado, decisiones y hoja de ruta

**Última actualización:** 20 de agosto de 2026 (America/Hermosillo)  
**Responsable científico:** Andrés Monreal Hernández  
**Estado:** prototipo funcional + cohorte piloto provisional verificada; aún no existe un benchmark de docking terminado.

## 1. Visión del proyecto

DockPrep Audit será una herramienta abierta, local y reproducible para detectar y documentar decisiones estructurales que pueden afectar la preparación de receptores para docking molecular. La inspiración general es el tipo de utilidad científica que ofrece AMDock, pero DockPrep Audit no copiará su interfaz, código ni marca y tendrá un propósito diferente: auditoría, trazabilidad y evaluación de riesgo metodológico antes de modificar una estructura.

El programa no declarará automáticamente que una estructura es biológicamente correcta ni aplicará reparaciones silenciosas. Su función será hacer visibles decisiones como selección de conformaciones alternativas, tratamiento de aguas, metales, residuos no estándar, anotaciones incompletas y duplicados.

## 2. Producto científico previsto

El resultado principal se plantea como un **artículo de software y métodos con validación original**. Podrá defenderse como contribución original si el benchmark demuestra que las características auditadas predicen diferencias medibles entre políticas de preparación o resultados de redocking.

Pregunta de trabajo:

> ¿Las características estructurales detectables antes de preparar un receptor predicen cuándo diferentes políticas de preparación cambian la reproducibilidad del docking molecular?

La herramienta por sí sola constituye software científico. Los resultados originales deberán proceder de un protocolo predefinido, controles, comparaciones pareadas, métricas cuantitativas e interpretación prudente.

## 3. Separación respecto a trabajos anteriores

El nuevo estudio debe permanecer diferenciado de:

1. El manuscrito de conformaciones alternativas y recuperación de poses sobre la cohorte previa.
2. El manuscrito del flujo de preparación y auditoría de 102 receptores DUD-E enviado a *Journal of Molecular Modeling*.

Se sincronizaron los registros públicos completos de ambos proyectos. El archivo `benchmark/previous-study-exclusions.csv` contiene **132 registros de exclusión** obtenidos de:

- `sircalch/dude-receptor-prep-audit`, registro de 102 objetivos DUD-E.
- `sircalch/docking-reference-audit`, registro de candidatos y expansiones del estudio anterior.

Durante la selección se detectó que `2P16` ya pertenecía al registro anterior y fue retirado. La exclusión actual se basa en identificadores PDB. Antes del artículo habrá que decidir si también se excluirán blancos biológicos repetidos aunque usen otra estructura PDB.

## 4. Lo que ya funciona

### Motor de auditoría v0.1.0

El paquete Python ya puede:

- Leer archivos PDB sin subirlos a servicios externos.
- Detectar conformaciones alternativas (`altLoc`).
- Detectar aguas, metales y residuos no estándar.
- Detectar ocupancias o elementos faltantes.
- Detectar registros atómicos duplicados.
- Asignar severidad y estado de revisión.
- Generar informes JSON y HTML portátiles.
- Ejecutarse desde línea de comandos con `dockprep-audit audit`.

La prueba automática incluida pasa correctamente.

### Inventario reproducible

El script `scripts/build_pilot_inventory.py`:

- Descarga PDB desde RCSB.
- Conserva URL de origen.
- Calcula SHA-256.
- Recupera título, método experimental y resolución.
- Registra heterogrupos.
- Ejecuta la auditoría.
- Excluye automáticamente identificadores de estudios anteriores.

### Verificación espacial

El script `scripts/verify_pilot_eligibility.py`:

- Comprueba una instancia explícita del ligando: componente, cadena y número de residuo.
- Identifica las cadenas del receptor situadas a 6 Å del ligando.
- Cuenta aguas a 4 Å del ligando.
- Calcula la distancia mínima entre ligando y conformaciones alternativas.
- Calcula la distancia mínima entre ligando y metales.
- Impide considerar relevante un hallazgo que aparece lejos del sitio de unión.

Este control descartó varios candidatos que parecían adecuados a escala de estructura completa pero no tenían el hallazgo cerca del ligando.

## 5. Cohorte piloto provisional actual

Los 12 casos pasan la verificación automática de ligando, cadena y relevancia espacial. Siguen siendo **provisionales** hasta revisar ensamblaje biológico y congelar el manifiesto antes del docking.

| PDB | Estrato | Ligando declarado | Instancia | Cadena receptora próxima | Evidencia local |
|---|---|---|---|---|---|
| 1M17 | Conformación alternativa | AQ4 | A:999 | A | altLoc a 3.105 Å |
| 1T46 | Conformación alternativa | STI | A:3 | A | altLoc a 2.841 Å |
| 4RJ3 | Conformación alternativa | 3QS | A:302 | A | altLoc a 4.328 Å |
| 5A2S | Metal/cofactor | OTF | A:1 | A | metal a 2.006 Å |
| 1CBX | Metal/cofactor | BZS | A:500 | A | metal a 2.305 Å |
| 4EXS | Metal/cofactor | X8Z | A:301 | A | metal a 2.274 Å |
| 1OHR | Política de aguas | 1UN | A:201 | A; B | 5 aguas a 4 Å |
| 3FNU | Política de aguas | 006 | A:329 | A | 4 aguas a 4 Å |
| 4GID | Política de aguas | 0GH | A:501 | A | 14 aguas a 4 Å |
| 1A28 | Control de bajo riesgo | STR | A:1 | A | sin altLoc/metal local; 1 agua cercana |
| 1QKT | Control de bajo riesgo | EST | A:600 | A | sin altLoc/metal local; 3 aguas cercanas |
| 1RBP | Control de bajo riesgo | RTL | A:183 | A | sin altLoc/metal local; 2 aguas cercanas |

Nota: “agua cercana” no significa “agua funcional” ni obliga a conservarla. Solo establece proximidad geométrica para una revisión posterior.

Nota de sustitución (2026-08-21): `1CPS` fue reemplazado por `1CBX` (carboxipeptidasa A + L-benzilsuccinato) tras fallar de forma irreparable en el caso de humo técnico por un problema geométrico real del depósito PDB en el residuo `TYR A:204` (ver sección 15). `1CBX` ya había pasado la verificación espacial en una ronda de candidatos previa y fue reverificado contra los 12 casos actuales (12/12 pass) y contra el pipeline técnico completo (RMSD 0.682 Å). Detalle en `benchmark/PILOT_SELECTION.md`.

## 6. Entorno técnico actual

Se creó un entorno aislado en `.venv` con Python 3.12.

Instalado correctamente:

- Meeko 0.7.1
- Gemmi 0.7.5
- RDKit 2026.3.5
- NumPy 2.5.2
- Pillow 12.3.0

Resuelto el 21 de agosto de 2026:

- **AutoDock Vina:** se descargaron `vina_1.2.7_win.exe` y `vina_split_1.2.7_win.exe` desde el release oficial `ccsb-scripps/AutoDock-Vina` v1.2.7 en GitHub, se guardaron en `tools/vina/` (excluido de git, checksums registrados en `tools/vina/MANIFEST.md`) y se verificó que `vina_1.2.7_win.exe --version` responde `AutoDock Vina v1.2.7`.
- **ProDy:** no fue necesario instalarlo ni instalar Visual C++ Build Tools. Meeko 0.7.1 expone `Polymer.from_pdb_string()` como ruta de preparación de receptor que no depende de ProDy (solo la CLI `mk_prepare_receptor` lo exige). Se instaló `scipy` (wheel precompilado, sin compilador) porque `meeko.PDBQTReceptor` lo requiere. Las clases clave (`Polymer`, `ResidueChemTemplates`, `MoleculePreparation`, `PDBQTWriterLegacy`, `PDBQTMolecule`, `RDKitMolCreate`) importan correctamente. La extracción reproducible de receptor/ligando (pendiente 5 de la sección 15) deberá escribirse contra la API de Meeko directamente, no contra su CLI.
- **dimorphite-dl 2.0.2:** instalado para la protonación de ligandos a pH fisiológico (sección 15, pendiente #9). Fija `rdkit<2026`, así que bajó nuestra versión de RDKit de 2026.3.5 a **2025.9.6** — verificado que no cambia ningún resultado del pipeline (mismo RMSD exacto en `1CBX` antes y después). `scrubber` (la herramienta que recomienda la guía oficial de Vina) no instala en Python 3.12 (`ModuleNotFoundError: 'imp'`); `dimorphite-dl` es su motor de protonación subyacente y sí instala limpio.

Pendiente:

- No están instalados Smina ni Open Babel (no urgente: la ruta principal usa Vina + Meeko).

No se ha ejecutado docking en este proyecto.

## 7. Decisiones ya tomadas

- El programa será abierto y local-first.
- El prototipo se desarrolla en Python.
- El primer formato auditado es PDB.
- Los reportes primarios serán JSON y HTML.
- El artículo combinará software, método y validación original.
- La primera validación será un piloto de 12 casos antes de ampliar a 60.
- Ningún hallazgo se tratará como biológicamente relevante solo por aparecer en el archivo completo.
- Las políticas de preparación se compararán sin declarar de antemano una opción universalmente correcta.
- Los resultados del nuevo artículo no reutilizarán como nuevos los resultados de los manuscritos anteriores.

## 8. Decisiones pendientes antes del docking

1. ~~Revisar el ensamblaje biológico y confirmar la cadena receptora de cada caso.~~ **Hecho 2026-08-21** (sección 15; `3FNU` resuelto por decisión explícita del usuario).
2. ~~Confirmar visualmente la identidad y pose del ligando declarado.~~ **Hecho 2026-08-21** (sección 15, vía métricas de validación wwPDB).
3. ~~Decidir si la independencia se exigirá solo por PDB o también por blanco biológico.~~ **Hecho 2026-08-21 (12 casos), extendido a los 60 el 2026-08-24**: independencia por PDB ID (ya cumplida); el traslape de blanco biológico se documenta explícitamente en vez de reabrir la cohorte — 6/12 casos originales, y 14/48 casos nuevos genuinos (3 más eran falso positivo por el péptido coactivador NCOA1 compartido, verificado directamente contra el UniProt real de la estructura excluida). Total: **20/60 casos (33%)** comparten blanco real con el estudio DUD-E previo. Ver `benchmark/PILOT_SELECTION.md` § Biological-target independence.
4. ~~Definir cómo tratar receptores multiméricos, especialmente `1OHR`.~~ **Hecho 2026-08-21 en la práctica, marcado explícitamente 2026-08-24**: `1OHR` ya corre como dímero A;B desde la verificación inicial (concuerda geométricamente, sin necesitar override); el mecanismo general para receptores multiméricos es `benchmark/receptor_chain_overrides.csv`, usado explícitamente en `3FNU` (2026-08-21) y `1SN5` (2026-08-23) cuando el hallazgo geométrico de 6 Å no coincide con el ensamblaje biológico real declarado en `REMARK 350`.
5. ~~Definir una política explícita para aminoácidos con altLoc A/B.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
6. ~~Definir qué aguas se conservan: ninguna, todas las cercanas o una selección predeclarada.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
7. ~~Definir el tratamiento de metales sin distorsionar coordinación o cargas.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
8. ~~Elegir la ruta de protonación y preparación del receptor.~~ **Parcialmente hecho 2026-08-21** (tautómeros de histidina coordinante de metal corregidos geométricamente; resto de residuos usa el default de Meeko — ver sección 15).
9. ~~Elegir la ruta de preparación de ligandos y el estado de protonación.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
10. ~~Congelar tamaño y centro de caja antes de ejecutar Vina.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
11. ~~Definir semillas, exhaustividad, número de modos y repeticiones.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).
12. ~~Registrar criterios de éxito, exclusión y fallo antes de observar resultados.~~ **Hecho 2026-08-21** (ver detalle abajo en la sección 15).

## 9. Plan experimental propuesto

### Fase A — congelar el piloto — **COMPLETADA 2026-08-21**

- ~~Completar la revisión estructural de los 12 casos.~~ Hecho (ensamblaje/cadena y ligando, sección 15).
- ~~Sustituir cualquier caso ambiguo antes de ejecutar docking.~~ Hecho (`1CPS` → `1CBX`).
- ~~Cambiar el estado de `provisional` a `frozen` en un manifiesto versionado.~~ Hecho: `benchmark/pilot_cases.csv` en `frozen`; manifiesto consolidado en `benchmark/pilot_manifest_frozen.csv`.
- ~~Guardar URL, checksum, cadena, ligando, caja y motivo de inclusión.~~ Hecho para URL/checksum/cadena/ligando/motivo (columnas `receptor_chains_basis` y `receptor_chains_override_reason` documentan la base de cada decisión de cadena). **La caja de docking queda deliberadamente fuera de este congelamiento** — es una decisión de Fase B/C (pendientes #10 de la sección 8), no de la composición de la cohorte.
- ~~Crear un SHA-256 del manifiesto congelado.~~ Hecho: `benchmark/pilot_manifest_frozen.sha256.txt` (`4d11807e0120db50c739b58143f6397159f3f674313f5b330a093e19d3829723`).

Generado por `scripts/freeze_pilot_manifest.py` (consolida `pilot_cases.csv` + `cohort_candidates.csv` + `eligibility_results.csv`; no introduce juicio científico nuevo, solo empaqueta lo ya auditado). **A partir de aquí, ningún caso puede añadirse, quitarse o reclasificarse de estrato sin una entrada nueva y fechada en `benchmark/PILOT_SELECTION.md` § Substitutions.** Las decisiones de política de preparación (Fase B: altLoc, aguas, metales, protonación, caja, semillas) siguen abiertas y no afectan la composición ya congelada de la cohorte.

### Fase B — preparar condiciones comparables

Se proponen tres políticas:

1. **Referencia conservadora:** conserva la información fuente y aplica una selección explícita cuando exista altLoc.
2. **Simplificada convencional:** aplica una política documentada de eliminación de solvente/heterogrupos.
3. **Flujo abierto predeterminado:** utiliza una herramienta abierta nombrada, con versión y parámetros registrados.

El motor, caja, semilla, exhaustividad y ligando deben mantenerse comparables dentro de cada caso.

### Fase C — piloto de redocking

- Extraer receptor y ligando de referencia.
- Preparar PDBQT con manifiestos por etapa.
- Ejecutar AutoDock Vina 1.2.7.
- Retener todas las poses solicitadas, afinidades, logs y parámetros.
- No descartar fallos de preparación o ejecución después de observarlos.

**Decisión 2026-08-22: solo 2 políticas reales, no 3.** Al ir a implementar la Fase C se confirmó lo ya anotado en la política de aguas: las Políticas 2 ("simplificada convencional") y 3 ("flujo abierto predeterminado") son **computacionalmente idénticas** en este pipeline (mismo Meeko, mismo Vina, `water_policy=none` en ambas, sin ninguna otra diferencia implementada) — correr "3 políticas" habría dado el mismo número dos veces bajo etiquetas distintas. El usuario decidió fusionarlas y correr **2 políticas reales**: `conservative_water` (Política 1, aguas puente) y `simplified_no_water` (Políticas 2/3 fusionadas, sin agua).

**Corrida lanzada 2026-08-22:** `scripts/run_phase_c.py` — 12 casos × 2 políticas × 3 semillas predeclaradas (42, 123, 2024), exhaustividad 32, umbral de éxito 2.0 Å. Resultados por corrida en `benchmark/phase-c/raw_runs/`, resumen (mediana por caso×política) en `benchmark/phase-c/phase_c_summary.csv`.

**Resultado de la Fase C, 2026-08-22 — 72/72 corridas completadas, 0 fallas de preparación.** Esto valida que las decisiones de la Fase B (altLoc, metales, His, protonación, caja) sostienen la ejecución real con repeticiones, no solo el diagnóstico técnico de un caso aislado.

| PDB | Estrato | `conservative_water` (mediana) | `simplified_no_water` (mediana) | Diferencia entre políticas |
|---|---|---|---|---|
| 1M17 | altLoc | sampling_fail (8.315 Å) | sampling_fail (7.969 Å) | Ninguna — falla en ambas |
| 1T46 | altLoc | success (0.831 Å) | success (0.830 Å) | Ninguna |
| 4RJ3 | altLoc | success (0.460 Å) | success (0.337 Å) | Ninguna (éxito en ambas) |
| 5A2S | metal | success (1.571 Å) | success (1.571 Å) | Idéntica (sin agua puente) |
| 1CBX | metal | success (1.416 Å) | success (1.416 Å) | Idéntica (sin agua puente) |
| 4EXS | metal | sampling_fail (2.670 Å) | sampling_fail (2.670 Å) | Idéntica (sin agua puente) |
| 1OHR | agua | success (0.831 Å) | success (0.794 Å) | Ninguna (éxito en ambas) |
| 3FNU | agua | sampling_fail (4.504 Å) | sampling_fail (4.513 Å)* | Ninguna — falla en ambas |
| 4GID | agua | inestable (2/3 `scoring_fail`, 1/3 `sampling_fail`) | inestable (2/3 `sampling_fail`, 1/3 `scoring_fail`) | Ninguna consistente |
| 1A28 | control | success (0.356 Å) | success (0.356 Å) | Idéntica |
| 1QKT | control | inestable (1/3 éxito) | inestable (2/3 éxito) | Leve, dentro del ruido de semilla |
| 1RBP | control | success (0.704 Å) | success (0.704 Å) | Idéntica |

*`3FNU` semilla 123 dio 7.919 Å en `simplified_no_water` frente a ~4.50 Å en las otras dos combinaciones — varianza real entre semillas, no un patrón de política.

**Lectura honesta:** con esta implementación concreta de "agua puente" (átomo de oxígeno rígido con carga fija de TIP3P, sin término de desolvatación), **la política de agua no cambió el resultado de forma sistemática en ningún caso**. Donde había agua puente detectada, incluirla o quitarla no mejoró ni empeoró consistentemente el redocking; las diferencias que aparecen (`1OHR`, `4RJ3`) son pequeñas y van en ambas direcciones, compatibles con ruido de semilla más que con un efecto real de la política. Esto es un hallazgo válido y reportable, pero debe declararse con su limitación: Vina, con un punto de carga rígido, puede no capturar el efecto de esa agua de la forma en que lo haría un método con energía de solvatación explícita (ej. MM/GBSA) o un scoring especializado.

**3 de 12 casos fallan consistentemente en ambas políticas con las 3 semillas coincidiendo** (`1M17`, `4EXS`, `3FNU`) — no es ruido, es una limitación real de ese receptor/ligando/caja específico, independiente de la política de agua.

**`4GID` y `1QKT` son intrínsecamente inestables entre semillas** (cambian de categoría de éxito con la misma política, solo cambiando la semilla) — confirma que la decisión de correr 3 semillas (pendiente #11) era necesaria: una sola corrida habría dado una impresión falsa de certeza en estos 2 casos.

### Fase D — métricas

Métricas primarias:

- RMSD de átomos pesados de la pose respecto al ligando cristalográfico.
- Éxito de redocking con umbral predefinido de RMSD ≤ 2 Å.
- Diferencia pareada de RMSD entre políticas.
- Concordancia de poses y rankings entre políticas.

Métricas secundarias:

- Frecuencia y coocurrencia de hallazgos.
- Acuerdo entre auditoría automática y revisión humana ciega.
- Tiempo de revisión/preparación.
- Fallos de preparación por condición.
- Sensibilidad a semilla o repetición cuando corresponda.

**Agregado por estrato, 2026-08-22 (sobre los resultados de la Fase C) — fracción de éxito promedio (3 semillas) y RMSD medio, por política:**

| Estrato | `conservative_water` (éxito / RMSD medio) | `simplified_no_water` (éxito / RMSD medio) |
|---|---|---|
| Conformación alternativa | 0.67 / 3.202 Å | 0.67 / 3.045 Å |
| Metal/cofactor | 0.67 / 1.886 Å | 0.67 / 1.886 Å |
| Política de aguas | 0.33 / 2.672 Å | 0.33 / 2.650 Å |
| Control de bajo riesgo | 0.78 / 2.435 Å | 0.89 / 0.993 Å |

**Lectura honesta — este agregado NO es un hallazgo confiable todavía, y hay que decir por qué:** en 3 de los 4 estratos las políticas son prácticamente indistinguibles (misma fracción de éxito, RMSD medio a menos de 0.2 Å de diferencia). La única brecha aparente (`low_risk_control`: 0.78 vs 0.89) está enteramente causada por **un solo caso**, `1QKT`, cuya diferencia pareada es +4.325 Å — no porque la política de agua cambiara algo en `1QKT` (que no tiene agua puente funcionalmente relevante), sino porque `1QKT` es intrínsecamente inestable entre semillas (ya documentado en la Fase C) y la semilla 123 dio un resultado distinto en cada política por azar. Con solo 3 casos por estrato, un solo caso inestable puede mover el promedio del estrato entero — esto es precisamente la limitación de tamaño de muestra que la sección 9 (Fase E) ya anticipaba. **No se debe reportar "el control de bajo riesgo favorece la política sin agua" como conclusión** — es ruido de un caso, no una señal del estrato.

**Conclusión de la Fase D con los 12 casos (revisada 2026-08-22, ver investigación abajo):** la política de agua no muestra diferencia sistemática en la categoría binaria de éxito (umbral 2.0 Å) en ningún estrato, y la varianza entre semillas domina sobre la varianza entre políticas en los casos límite (`4GID`, `1QKT`). **Pero esta conclusión, mirada solo con la categoría binaria, esconde un efecto real** — ver la investigación de los 3 casos con fallo consistente inmediatamente abajo, que corrige esta lectura para `1M17`.

**Investigación de los 3 casos con fallo consistente (`1M17`, `4EXS`, `3FNU`), 2026-08-22:** se revisaron las 9 poses completas (no solo la mejor reportada) en las 6 corridas (3 semillas × 2 políticas) de cada caso.

- **`4EXS`**: mínimo RMSD consistente entre **2.187 y 2.252 Å** en las 6 corridas — un "casi acierto" robusto y reproducible, justo arriba del umbral. Sin efecto de la política de agua (no tiene agua puente detectada; las corridas con y sin agua dan arrays de RMSD idénticos byte a byte en varias semillas).
- **`3FNU`**: mínimo RMSD consistente entre **2.896 y 2.968 Å** en las 6 corridas — también un casi acierto reproducible, sin efecto de la política de agua.
- **`1M17`**: aquí la política de agua **sí tiene un efecto real y grande**, invisible en la categoría binaria: con agua puente, la mejor pose está entre **2.17 y 2.38 Å** (las 3 semillas); sin agua, entre **4.29 y 4.78 Å** (las 3 semillas). Es una diferencia consistente de ~2 Å completos, reproducible en las 3 semillas — el agua ayuda genuinamente a que Vina muestree cerca de la pose correcta en este caso, aunque ninguna de las dos políticas cruce el umbral de 2.0 Å para clasificarse como `success`.

**Conclusión corregida:** afirmar que "la política de agua no tuvo efecto" habría sido una lectura incompleta, causada por reducir el resultado a una categoría binaria. Con RMSD continuo, **al menos 1 de 12 casos (`1M17`) muestra un efecto de política real y grande**; los otros 2 casos de fallo consistente (`4EXS`, `3FNU`) sí son casos genuinamente insensibles a la política de agua, no solo "escondidos" por el umbral. Esto refuerza reportar siempre el RMSD continuo pareado junto con la categoría de éxito (ya listado como métrica primaria en la sección 9), no solo la tasa de éxito binaria — un umbral único puede ocultar exactamente el tipo de señal que este proyecto busca detectar.

**Análisis final con RMSD continuo, mejor-de-9 poses, mediana de 3 semillas por caso (2026-08-22):** se recalculó la diferencia pareada (`conservative_water` − `simplified_no_water`) usando la **mejor de las 9 poses** (capacidad de muestreo) en vez de solo la pose #1 (que depende del ranking de Vina, más ruidoso), para separar el efecto de la política del ruido de la función de puntuación.

| PDB | Estrato | Mejor pose, con agua (mediana) | Mejor pose, sin agua (mediana) | Diferencia |
|---|---|---|---|---|
| `1M17` | altLoc | 2.358 Å | 4.737 Å | **−2.379 Å** |
| `4GID` | agua | 1.907 Å | 2.283 Å | −0.376 Å |
| `4RJ3` | altLoc | 0.460 Å | 0.337 Å | +0.123 Å |
| `1OHR` | agua | 0.831 Å | 0.794 Å | +0.037 Å |
| `3FNU` | agua | 2.898 Å | 2.948 Å | −0.050 Å |
| `1QKT` | control | 1.448 Å | 1.459 Å | −0.011 Å |
| resto (6 casos) | — | — | — | 0.000 Å |

**Agregado por estrato (diferencia media, mejor-de-9):** conformación alternativa −0.752 Å (enteramente por `1M17`; los otros 2 casos ≈0), metal/cofactor 0.000 Å, política de aguas −0.130 Å (sobre todo `4GID`), control de bajo riesgo −0.004 Å.

**Hallazgo adicional al usar mejor-de-9 en vez de solo la pose #1:** `1QKT` y `4GID`, que parecían "inestables entre semillas" en el análisis de la Fase C (basado en la categoría de la pose #1), resultan **mucho más estables** al mirar la mejor pose muestreada (`1QKT`: 1.448 vs 1.459 Å entre políticas, prácticamente idéntico). Esto separa dos fenómenos que el análisis anterior mezclaba: la **inestabilidad real de muestreo** (¿existe una pose cercana a la correcta?) y el **ruido de la función de puntuación** (¿cuál pose pone Vina en el puesto #1?) — en `1QKT`/`4GID` el problema es lo segundo, no lo primero.

**Conclusión final, defendible y cuantificada:** de los 12 casos, **solo `1M17` muestra un efecto de política de agua grande y reproducible** (−2.38 Å, consistente en las 3 semillas). El resto de la variación observada en corridas anteriores era ruido de puntuación, no de preparación. Con n=12 y 1 solo caso con efecto claro, esto no permite generalizar por estrato todavía.

**Explicación mecanística de por qué `1M17` es distinto, verificada 2026-08-22:** la agua puente de `1M17` (`A:10`) está a **2.78 Å del átomo N3 del ligando `AQ4`** (el análogo de erlotinib) — dentro de la distancia clásica de puente de hidrógeno (2.7–3.2 Å). N3 es el nitrógeno del anillo de quinazolina, el farmacóforo central que ancla los inhibidores tipo erlotinib/gefitinib a la región bisagra (*hinge*) del dominio quinasa de EGFR. Esta agua también contacta el receptor (`CYS 751 SG`, ver sección de política de aguas más arriba). Es decir: **no es una agua puente genérica cualquiera — es mecanísticamente análoga a la famosa agua de la flap de la proteasa del VIH** (mismo patrón: puente literal entre el farmacóforo del inhibidor y el receptor), solo que en el sitio de unión de una quinasa en vez de una proteasa aspártica. Esto explica de forma coherente por qué eliminarla degrada tanto la pose muestreada: se está quitando un punto de anclaje real del ligando, no un adorno estructural.

**Hipótesis probada y refutada, 2026-08-22:** se planteó que el efecto debería correlacionar con si la agua puente forma un contacto a distancia de puente de hidrógeno (≤3.2 Å) con un heteroátomo (N/O) del ligando, no solo con "hay agua cerca". Al aplicar este criterio más estricto a los 12 casos, **9 de 12 califican** (`1M17`, `1T46`, `4RJ3`, `5A2S`, `1OHR`, `3FNU`, `4GID`, `1A28`, `1QKT`) — incluyendo dos controles de bajo riesgo que no mostraron ningún efecto real en el docking. **La hipótesis no se sostiene**: el contacto a distancia de H-bond con un heteroátomo del ligando es común y no distingue por sí solo cuáles casos son sensibles a la política de agua. Solo `1M17` mostró el efecto grande pese a que 8 casos más cumplen el mismo criterio geométrico. Se reporta esto honestamente en vez de forzar una narrativa más limpia de la que los datos sostienen: la explicación mecanística de `1M17` (agua puente a N3 de la quinazolina, contacto con Cys751) sigue siendo válida como descripción de *ese caso*, pero no generaliza a una regla predictiva simple para los otros 11. Identificar qué hace a `1M17` distinto de los otros 8 casos con contacto similar requeriría, como mínimo, examinar la geometría de la cavidad (¿la ausencia del agua deja un hueco que Vina explora mal?) — no se investigó más a fondo por ahora.

### Fase D sobre los 60 casos — 2026-08-24

**Fase C de los 60 completada, 2026-08-24:** 360/360 corridas (60 casos × 2 políticas × 3 semillas, exhaustividad 32), **0 fallas de preparación**. Los 72 resultados de los 12 casos originales se reutilizaron sin recomputar (parámetros idénticos); solo se corrieron 288 nuevas para los 48 casos de Fase E. Script nuevo y reproducible: `scripts/run_phase_c_60.py`. Se generalizó también el análisis de agregado por estrato en `scripts/run_phase_d.py` (antes hecho a mano) y se validó contra los 12 casos originales antes de usarlo — reproduce exactamente los números ya documentados arriba.

**Agregado por estrato, n=15 por estrato — fracción de éxito (mejor-de-9 poses) y RMSD medio, por política:**

| Estrato | `conservative_water` (éxito / RMSD medio) | `simplified_no_water` (éxito / RMSD medio) | Diferencia (éxito) |
|---|---|---|---|
| Conformación alternativa | 0.778 / 1.286 Å | 0.822 / 1.545 Å | −0.044 |
| Metal/cofactor | 0.689 / 2.171 Å | 0.667 / 2.224 Å | +0.022 |
| **Política de aguas** | **0.644 / 2.963 Å** | **0.556 / 3.123 Å** | **+0.088** |
| Control de bajo riesgo | 0.933 / 1.459 Å | 0.933 / 1.470 Å | 0.000 |

**Con n=15 (vs. n=3 del piloto), el patrón cambia de "un solo caso atípico" a una señal distribuida y en la dirección esperada:** el estrato `water_policy` es el único con una diferencia de éxito notable entre políticas (+8.8 puntos porcentuales a favor de conservar agua puente), y es también el estrato con el RMSD medio más alto en términos absolutos (más difícil en general). Se recalculó caso por caso (diferencia `simplified_no_water` − `conservative_water` en RMSD de mejor pose; positivo = el agua ayuda):

| Estrato | Casos donde el agua ayuda claramente (Δ>0.3 Å) | Casos donde perjudica (Δ<−0.3 Å) | Proporción que ayuda |
|---|---|---|---|
| **Política de aguas** | `1CVZ` (+0.47), `1EPP` (+0.47), `1WBK` (+0.50), `4GID` (+0.38), `6ASH` (+1.42) | `1EED` (−0.31), `2F25` (−0.37) | **5/15 (33%)** |
| Conformación alternativa | `1M17` (+2.38 en pose top1, el mismo efecto documentado en el piloto), `2I4H` | `5E0J` | 2/15 (13%) |
| Metal/cofactor | `1KJO` | `3LXE` | 1/15 (7%) |
| Control de bajo riesgo | `3BEJ` | (ninguno) | 1/15 (7%) |

**Lectura, cuantificada y con las mismas salvedades que en el piloto de 12 (n pequeño, no es un test estadístico formal):** la proporción de casos con un beneficio claro de retener agua puente es ~3× a ~5× mayor en el estrato que el motor de auditoría marca por presencia de agua puente relevante que en los otros tres estratos combinados. Esto es **más evidencia a favor de la pregunta central del proyecto** que lo que permitía concluir el piloto de 12 (donde solo `1M17` mostraba el efecto y no se podía generalizar por estrato) — pero sigue sin ser una relación determinista: 8 de 15 casos del propio estrato `water_policy` no muestran ningún efecto, y 2 muestran el efecto contrario (agua perjudica). La hipótesis geométrica simple ya refutada arriba (contacto H-bond ≤3.2 Å) sigue sin explicar por qué unos casos sí y otros no dentro del mismo estrato — sigue siendo trabajo pendiente identificar qué distingue a los 5 casos que sí responden.

Archivos: `benchmark/phase-c-60/phase_c_summary.csv`, `benchmark/phase-c-60/phase_d_by_case.csv` (120 filas, caso × política), `benchmark/phase-c-60/phase_d_by_stratum.csv` (8 filas, estrato × política).

### Fase E — ampliación

Solo si el piloto es ejecutable y las métricas son interpretables:

- Ampliar a 60 complejos, 15 por estrato.
- Congelar una cohorte independiente antes de nuevos resultados.
- Añadir análisis pareados e intervalos de confianza.
- Separar análisis confirmatorio de exploraciones posteriores.

**Búsqueda de candidatos completa, 2026-08-22 — 48/48 casos nuevos encontrados y verificados. Meta de 60 casos (15 por estrato) cumplida exactamente.** Búsqueda vía la API de RCSB (texto completo + palabras clave de `struct_keywords`), auditoría con `dockprep_audit` y filtro geométrico (altLoc/metal/agua-puente ≤6 Å del ligando, igual criterio que los 12 originales), con tope de 2 casos por proteína (UniProt) para evitar que un solo blanco domine un estrato — decisión tomada explícitamente tras encontrar que una búsqueda inicial de "kinase inhibitor" daba 6 de 13 candidatos como la misma proteína (PKA).

| Estrato | Nuevos encontrados | Total (con los 3 originales) | Meta |
|---|---|---|---|
| Conformación alternativa | 12 | **15 ✓** | 15 |
| Metal/cofactor | 12 | **15 ✓** | 15 |
| Política de aguas | 12 | **15 ✓** | 15 |
| Control de bajo riesgo | 12 | **15 ✓** | 15 |
| **Total** | **48** | **60 ✓** | **60** |

Los 48 candidatos se verificaron con el **mismo motor formal** usado para los 12 originales (`scripts/verify_pilot_eligibility.py`): **48/48 pasan sin necesitar revisión manual** (ligando encontrado, instancia única, cadena de receptor geométricamente relevante — algunos casos correctamente detectan que necesitan más de una cadena, ej. dímeros de proteasas tipo pepsina/renina, sin necesitar anulación manual como `3FNU`). Ninguno está en el registro de exclusión de estudios previos (132 registros).

**Nota de control de calidad al cerrar la búsqueda:** dos candidatos encontrados vía la búsqueda "HIV-2 protease inhibitor" (`1WBK`, `1WBM`) resultaron ser, según su propio encabezado PDB, proteasa del **VIH-1** (no VIH-2) — la búsqueda de texto completo de RCSB había hecho coincidir un término mencionado de pasada en el artículo, no el blanco real. Se detectó al leer el `TITLE`/`COMPND` del archivo antes de aceptarlo, evitando duplicar sin darse cuenta el blanco de `1OHR`; se usó solo `1WBK`, respetando el límite de 2 casos por proteína (ya son 2 entradas de proteasa del VIH-1 en el estrato de agua: `1OHR` y `1WBK`).

Archivos: `benchmark/expansion_cases.csv` (los 48 casos, con estrato/ligando/cadena/resseq), PDBs copiados a `benchmark/pilot-inventory/raw-pdb/`, resultados de elegibilidad en `benchmark/expansion-eligibility/`.

**Cohorte de 60 congelada, 2026-08-22.** `expansion_cases.csv` se fusionó con `pilot_cases.csv` (60 filas, estado `frozen`), se regeneró la elegibilidad completa (60/60 pass) y se regeneró el manifiesto congelado: `benchmark/pilot_manifest_frozen.csv`, checksum `72346f3e5fee4284599d0f8b72cc45a0974c5864e22382b5c405f2c18322c66a`.

**Nota de robustez encontrada al congelar:** la API GraphQL de RCSB (`data.rcsb.org`) estuvo caída de forma sostenida durante esta sesión (confirmado con `curl` directo, no solo desde el script) — el chequeo de independencia por blanco biológico (columnas `target_uniprot`/`shares_target_*`) no pudo completarse para los 60 casos. Se corrigió `scripts/freeze_pilot_manifest.py` para que este fallo externo **no bloquee el congelamiento** (reintentos con espera creciente, y si persiste, congela igual con esas columnas vacías, avisando explícitamente) — la proveniencia central (URL, checksum, cadena, ligando) no depende de esa API y no debía quedar rehén de una caída ajena. **Pendiente real:** volver a correr `scripts/freeze_pilot_manifest.py` cuando la API se recupere, para rellenar la independencia biológica de los 48 casos nuevos (los 12 originales ya se verificaron manualmente en la sección de arriba: 6 de 12 comparten blanco con el estudio DUD-E previo).

**Pendiente antes de correr la Fase C sobre los 60 casos:** (1) correr el pipeline técnico completo (extracción, protonación, Vina) sobre los 48 casos nuevos — un trabajo de cómputo considerable, dado que algunos casos ya tardaron varios minutos por corrida en el piloto de 12; (2) revisión visual de ensamblaje/cadena para los casos multi-cadena detectados automáticamente, antes de dar por buena la extracción de receptor; (3) backfill de independencia biológica mencionado arriba.

**Smoke test técnico de los 48 casos nuevos, 2026-08-22 (`scripts/run_expansion_smoke.py`, política por defecto: sin agua, exhaustividad 8, semilla 42):** 36/48 corrieron sin error en la primera pasada; 12 fallaron. Diagnóstico caso por caso:

- **Bug real corregido — altLoc del ligando nunca se resolvía (`3P0M`, `1SN5`):** `extract_ligand_atoms()` en `scripts/smoke_redock_case.py` solo filtraba por componente/cadena/resseq, sin resolver `altLoc` — a diferencia de la extracción del receptor, que ya usaba `choose_altloc_conformers()`. `3P0M` (ligando `4SB`, altLoc A/B a 0.54/0.46) y `1SN5` (ligando `T3`, altLoc A/B/C a 0.40/0.40/0.20) exponían ambos/los tres conformeros a la vez, dando coordenadas duplicadas que rompían el emparejamiento exacto de átomos (`build_reference_to_pdbqt_index_map`, tolerancia 0.01 Å) más adelante ("No ligand.pdbqt atom within 0.01 A of reference atom"). Corregido aplicando la misma política ya establecida para el receptor (mayor ocupancia, empate en `A`). Reverificado: `3P0M` corre limpio (`scoring_fail`, RMSD top-1 8.197 Å / mejor pose 1.925 Å) y `1SN5` corre limpio (`sampling_fail`, RMSD 6.6–8.6 Å) — ya no son fallos de preparación, son resultados de docking legítimos.
- **Fallo de red, no de código (`3MWU`):** timeout contra `data.rcsb.org` durante la caída sostenida de esa sesión. Reintentado tras la recuperación de la API: corre limpio (`success`, RMSD top-1 1.128 Å).
- **Límite genuino de Meeko — residuo no estándar sin plantilla de "padding" (`1PSO`, `1WKR`, `1XDH`):** los tres son complejos de una aspártico-proteasa tipo pepsina con **pepstatina**, cuyo residuo de estatina (`STA`) rompe el enlace peptídico normal C-N. Meeko sí reconoce `STA` como parte del polímero (llega como `ATOM`, no `HETATM`, y nuestro filtro de residuos no estándar ya lo deja pasar) pero su lógica de *padding* de enlaces truncados no tiene plantilla para el patrón de conectividad de `STA` (`RuntimeError: Expected 2 paddings ... but got 1`). No es un bug de este proyecto — es una limitación conocida de la cobertura de plantillas de Meeko para péptidomiméticos, misma categoría que el problema geométrico irreparable de `1CPS`. Nota: esto concentra 3 de los 12 fallos en el estrato `water_policy` (pepstatina es clásica para estudiar agua puente), donde además `1EPQ` y `1ENT` fallan por una razón distinta (ver abajo) — 5 de 12 casos de ese estrato con problema técnico, pendiente de decisión de sustitución.
- **Fallos genuinos de química del ligando, sin variante de protonación válida (`6FTF`, `7ORS`, `1EPQ`, `1THL`, `1ENT`, `3FLI`):** en cada caso, todas las variantes de protonación propuestas por dimorphite-dl fallan al construir el mol de RDKit (`AssignBondOrdersFromTemplate`/sanitización, típicamente valencia de nitrógeno excedida) o dimorphite-dl no devuelve ninguna variante. Misma categoría que el fallo de `1CPS`/`7CI`: complejidad química real del ligando (nucleósido-análogo con ciano+purina en `6FTF`; guanidinas/heterociclos aromáticos con N cargado en los demás), no un bug de extracción.

**Resuelto 2026-08-23 — los 9 casos sustituidos, misma cohorte re-congelada.** Decisión del usuario: sustituir los 9 (no documentar como exclusión). Un agente en segundo plano buscó y verificó candidatos de reemplazo (elegibilidad formal + smoke test técnico completo, mismo estándar que `1CPS`→`1CBX`); el primer candidato para `1EPQ` (`1MEM`, catepsina K + inhibidor vinil-sulfona) también falló el smoke test por la misma razón de fondo (ligando covalentemente unido al esqueleto, mismo límite de plantillas de Meeko que los casos `STA`), así que se descartó por `1EED`.

| Caso sustituido | Motivo | Reemplazo | Ligando (cadena/resseq) |
|---|---|---|---|
| `6FTF` | química de ligando | `3PNA` | CMP / A / 250 |
| `7ORS` | química de ligando | `5K8S` | CMP / A / 501 |
| `1THL` | química de ligando | `1KJO` | THR / A / 1317 |
| `1PSO` | STA/pepstatina | `1EPP` | 1Z1 / E / 333 |
| `1WKR` | STA/pepstatina | `1PPM` | 0P1 / E / 324 |
| `1XDH` | STA/pepstatina | `1HRN` | 03D / A / 391 |
| `1EPQ` | química de ligando (`1MEM` también falló, ligando covalente) | `1EED` | 0EO / P / 327 |
| `1ENT` | química de ligando | `1WBM` | BLL / B / 1100 |
| `3FLI` | química de ligando | `3BEJ` | MUF / A / 473 |

Los 9 reemplazos pasaron elegibilidad formal (`verify_pilot_eligibility.py`, 60/60 pass) y el smoke test técnico completo (extracción + protonación + redocking Vina sin excepción — resultados `success`/`scoring_fail`/`sampling_fail`, todos resultados de docking legítimos, no fallos de preparación). Detalle completo, incluyendo la nota sobre el falso positivo del tope de 2-por-blanco (`Q15788`/NCOA1, un péptido coactivador compartido, no un blanco real) en `benchmark/PILOT_SELECTION.md`, sección "Substitutions".

**Manifiesto re-congelado, 2026-08-23:** `pilot_cases.csv`, `expansion_cases.csv`, `pilot-inventory/cohort_candidates.csv` y `pilot-eligibility/eligibility_results.csv` regenerados con los 9 reemplazos; `pilot_manifest_frozen.csv` re-generado con `freeze_pilot_manifest.py`. Nuevo checksum SHA-256: `9313e6a7498c07a436fda25c28ae56f0a10065f89edff8c57bb099578b06dbba`. **Bono:** la API de RCSB ya estaba disponible en este re-congelamiento, así que el backfill de independencia por blanco biológico (`target_uniprot`) que quedó pendiente el 2026-08-22 se completó automáticamente para los 60 casos — ya no es un pendiente.

**Bloqueo de entorno encontrado y resuelto durante este trabajo:** Windows Smart App Control bloqueó `RDKitMolStandardize....dll` a la 12:28:57 a.m. del 2026-08-23 (justo al cambiar la fecha del sistema), rompiendo cualquier ejecución de `smoke_redock_case.py` que dependiera de `dimorphite_dl`/`rdkit.Chem.MolStandardize`. No se intentó sortear — se reportó al usuario, quien lo desactivó manualmente. Documentado por si reaparece en una sesión futura tras otro cambio de fecha o actualización de Windows.

**Verificación wwPDB + ensamblaje biológico de los 48 casos nuevos, 2026-08-23 (pendientes #1 y #2 de la sección 8, ahora aplicados también a los 48):** un agente en segundo plano hizo la misma verificación cuantitativa usada en los 12 originales (RSCC/RSR, outliers de Mogul, choques, ocupancia vía `data.rcsb.org/graphql`) más el chequeo de ensamblaje biológico (`REMARK 350` de cada PDB crudo vs. `receptor_chains` ya asignado geométricamente). Resultados en `benchmark/expansion-validation/wwpdb_validation.csv`:

- **Identidad/pose: 37/48 limpios** (RSCC ≥ 0.888), **11/48 sin RSCC/RSR disponible** — no solo por ser pre-1998 (el precedente de los 12 originales), sino porque el depósito rutinario de factores de estructura no fue universal ni siquiera después de 1998: se confirmó directamente que 6 casos posteriores a 1998 (`1E6U`/`1E7S` 2000, `1O86` 2002, `1UZE` 2004, `1CVZ` 1999, `1G5Y` 2000) tampoco tienen archivo de factores de estructura en RCSB (404 en `files.rcsb.org/download/{id}-sf.cif`), junto con 5 pre-1998 esperados (`1DTH`, `1EED`, `1EPP`, `1PPM`, `1HRN`). **Ningún caso con RSCC < 0.8** (el más débil, `1WBM`, 0.817 — mismo rango que `1M17` en el piloto original). Sin discrepancias de identidad química (ligando declarado vs. nombre real vs. título de la estructura).
- **Ensamblaje/cadena: 47/48 de acuerdo**, incluyendo dos casos resueltos con cálculo directo en vez de solo comparar letras de cadena (`1E6U`/`1E7S`: la "segunda copia" del `REMARK 350` resultó ser una réplica por simetría cristalográfica del mismo monómero, medida a 12+ Å del ligando — irrelevante; `4G9L`/`4JA1`: había que emparejar con el biomolecule correcto del `REMARK 350`, no por defecto el primero).
- **1 caso con override real, mismo patrón que `3FNU` — `1SN5` (transtiretina, `low_risk_control`, ligando T3):** `REMARK 350` declara explícitamente `AUTHOR DETERMINED BIOLOGICAL UNIT: TETRAMERIC`, pero la regla geométrica de 6 Å solo capturó A;C porque B y D quedan a 6.61/6.86 Å — pasando el corte por poco, no como una copia cristalográfica lejana (contraste: los casos de `1E6U`/`1E7S` estaban a >12 Å). El sitio de unión de T3/T4 en transtiretina está documentado en la literatura como el canal central entre los dos dímeros del tetrámero. **Decisión del usuario 2026-08-23: aplicar el override A;B;C;D** (mismo criterio que `3FNU`). Aplicado en `benchmark/receptor_chain_overrides.csv`, re-verificado (60/60 pass, `receptor_chains_geometric` conserva el hallazgo crudo A;C para transparencia) y re-probado con el smoke test técnico — corre limpio (`sampling_fail`, sin error de preparación).

**Manifiesto re-congelado, 2026-08-23 (segunda vez el mismo día):** nuevo checksum SHA-256 `ab37330651f5787bc25b4a7a301f7b7c4a98bfdd580139aa027898d6e5d4abfe`.

**Nota de consistencia (2026-08-27):** esta línea decía "Fase C multi-semilla sobre los 60 casos" como pendiente -- quedó desactualizada: esa Fase C se corrió y completó más adelante la misma sesión (ver "Fase D sobre los 60 casos" arriba en esta sección: 360/360 corridas, ambas políticas × 3 semillas × exhaustividad 32, checksum de resultados en `benchmark/phase-c-60/`). Lo único que sigue pendiente de lo mencionado originalmente aquí es la Fase B como *políticas alternativas adicionales* más allá de las 2 ya comparadas (conservador/simplificado) — no hay trabajo de Fase C bloqueado.

## 10. Figuras previstas para el artículo

1. Resumen gráfico: estructura PDB → auditoría → decisión explícita → docking.
2. Diagrama de elegibilidad, exclusiones y selección de cohorte.
3. Arquitectura de DockPrep Audit y flujo de procedencia.
4. Frecuencia y coocurrencia de hallazgos estructurales.
5. Heatmap caso × hallazgo × política de preparación.
6. Distribución de RMSD y éxito de redocking por política.
7. Casos 3D: altLoc, agua y metal en el sitio de unión.
8. Concordancia entre auditoría automática y revisión humana.
9. Árbol de decisión para documentar preparación reproducible.
10. Capturas de la interfaz gráfica y exportación del reporte.

Las figuras deben responder preguntas; no se añadirán imágenes únicamente para aumentar su número.

**Figuras 1, 3, 4, 5, 6, 7, 9 completadas, 2026-08-24** (`scripts/build_figures.py`, `scripts/build_diagrams.py`, `scripts/build_site_figures.py` → `benchmark/figures/`, PNG 300dpi + PDF vectorial). Paleta categórica validada colorblind-safe (skill de dataviz del proyecto, azul/naranja para las 2 políticas), tipografía Arial (registrada desde `C:/Windows/Fonts/`, convención de la mayoría de revistas de modelado molecular). Instalado `matplotlib`+`seaborn` (no estaban en el entorno).

**Figura adicional de control de calidad, 2026-08-27** (`scripts/build_qc_figure.py` → `fig_qc_rscc_by_year.png`): año de depósito vs. disponibilidad de RSCC para los 48 casos de Fase E, visualizando directamente el hallazgo de la sección 9 ("la brecha de validación no es un corte limpio pre-1998") en vez de dejarlo solo en prosa. Se agregó al redactar el manuscrito ([manuscript/draft.md](manuscript/draft.md)) al notar que el dato ya existía (`benchmark/expansion-validation/wwpdb_validation.csv`) pero no tenía figura propia.

- **Figura 2 (embudo de selección) omitida a propósito:** los conteos de candidatos crudos provienen de ~20 lotes de búsqueda ad hoc nunca deduplicados de forma confiable entre sí; publicar un número que parezca exacto sin serlo se juzgó peor que no tener la figura.
- **Figura 7 (sitio de unión 3D), tres iteraciones:**
  - **v1 (2026-08-24):** se intentó primero `pymol-open-source` (PyPI); el wheel resultó estar pensado para distribución conda -- depende de 5 DLLs nativas (GLEW, FreeType, libpng16, libxml2, netCDF) no incluidas en el paquete ni presentes en el sistema, y no hay conda instalado. Instalar una distribución conda completa solo para esto se juzgó un cambio de entorno más grande de lo pedido; se presentaron 3 opciones al usuario y se eligió construir la figura con matplotlib 3D directamente desde las coordenadas reales del PDB, con enlaces inferidos por distancia de covalencia (≤1.75 Å). Archivo conservado como `fig_binding_site_examples_matplotlib_v1.png/.pdf`.
  - **v2 (2026-08-26):** al comparar contra artículos Q1 de ejemplo que el usuario compartió, matplotlib 3D se veía plano frente a un render molecular real. Se probó `pyvista` (VTK vía pip, sin el problema de DLLs nativas de PyMOL -- instala limpio, renderiza headless sin problema en este sistema) para esferas/cilindros con sombreado real y oclusión ambiental (`plotter.enable_ssao`). Mismas coordenadas y distancias que v1, mismos 3 casos. Nota: intentar primero rehacer las 3 figuras *conceptuales* (resumen gráfico, arquitectura, árbol de decisión) en Canva no funcionó -- el tipo de diseño "infographic" de Canva está optimizado para redes sociales (fabricó estadísticas falsas, desordenó el orden de los pasos en 2 de 4 candidatos); se descartó y esas 3 se rehicieron en matplotlib con un estilo más austero (cajas blancas con acento de color delgado en el borde, en vez de rellenos pastel).
  - **v3 (2026-08-27, actual):** las etiquetas de texto embebidas en la escena 3D de PyVista (`add_point_labels`) seguían tapando parte de la estructura o solapándose entre sí en el panel de metal (4 coordinadores convergiendo en un punto) pese a varias rondas de ajuste de fuente/offset/ángulo de cámara. Se rediseñó el pipeline: cada panel se renderiza limpio (sin texto embebido) con PyVista, las coordenadas 3D de interés se proyectan a píxeles 2D exactos (`vtkCoordinate.GetComputedDisplayValue`), y las etiquetas se dibujan como overlay 2D con matplotlib (`ax.annotate`) con línea guía delgada hacia el átomo real -- estilo "callout" estándar de figuras científicas. Se probó primero un buscador automático de espacio libre (`find_clear_spot`, prueba posiciones en anillos crecientes y mide blancura de píxeles) pero resultó frágil cuando 4-5 anclas caían muy juntas (panel de metal): encontraba huecos técnicamente blancos pero mutuamente solapados. Se resolvió alejando la cámara de los paneles de agua/metal (más margen real) y fijando las posiciones de etiqueta a mano usando las coordenadas 2D proyectadas reales (impresas con un script auxiliar), en vez de seguir ajustando el buscador automático a ciegas. Script final: `scripts/build_site_figures_pv.py`.
- **Figura 8 (concordancia auditoría automática vs. revisión humana) y figura 10 (capturas de GUI) siguen sin aplicar:** la "revisión humana" de este proyecto fue la verificación wwPDB/ensamblaje ya hecha (sección 9, Fase E), no una anotación manual independiente en formato de matriz de confusión; la GUI (v0.3.0) no se ha construido.

## 11. Interfaz gráfica prevista

La interfaz tipo aplicación se construirá después de estabilizar el motor. Funciones previstas:

- Arrastrar o seleccionar un PDB/mmCIF.
- Vista general de severidades y hallazgos.
- Visor molecular con resaltado de altLoc, aguas, metales y ligando.
- Selección explícita de cadena, ligando y política.
- Explicación de por qué cada decisión puede importar.
- Exportación de reporte HTML, JSON, imágenes y manifiesto reproducible.
- Ejecución opcional de preparación y Vina, separada de la auditoría.
- Historial local de proyectos sin subir estructuras.

Se debe evitar presentar una corrección automática como verdad biológica. La interfaz siempre mostrará las decisiones y sus consecuencias.

## 12. Versiones planeadas

### v0.1.0 — prototipo actual

- Auditoría PDB.
- CLI.
- Reportes JSON/HTML.
- Inventario y verificación espacial del piloto.

### v0.2.0 — piloto reproducible

- Manifiesto congelado.
- Extracción de receptor/ligando.
- Preparación registrada.
- Vina integrado mediante binario versionado.
- Resultados del piloto y primeras figuras.

### v0.3.0 — interfaz gráfica

- Carga de estructuras.
- Panel visual de hallazgos.
- Visor molecular y exportación.
- Pruebas de usabilidad.

### v1.0.0 — versión asociada al artículo

- Cohorte completa o validación final acordada.
- Documentación, tutoriales y ejemplos.
- Entorno reproducible.
- Archivo Zenodo con DOI.
- Lanzamiento GitHub.
- Manuscrito y material suplementario.

## 13. Publicación prevista

**Decisión, 2026-08-27: dos artículos, no uno.** Patrón estándar para este tipo de proyecto (herramienta + validación original) — un artículo de resultados/métodos en una revista de dominio, y un artículo de software separado en JOSS del mismo repositorio, citándose mutuamente. No es redundante: JOSS explícitamente espera evidencia de "necesidad de investigación" (statement of need), y un artículo previo que use la herramienta es justo esa evidencia.

**Artículo 1 — resultados y métodos, envío más próximo:**
- Destino principal: *Journal of Molecular Modeling* (mismo destino que el manuscrito previo de auditoría DUD-E del autor).
- Alternativa: revista de quimioinformática/modelado si el benchmark de 60 casos termina pesando más que la herramienta en sí.
- Base de datos ya lista: los 60 casos con Fase C completa (360 corridas, 2 políticas × 3 semillas × exhaustividad 32) y el hallazgo de la señal distribuida de política de agua (sección 9). Pendiente: solo escribir el manuscrito — los datos y las 7 figuras ya están congelados.
- Envío vía Editorial Manager (sistema de Springer); requiere carta de presentación, declaración de disponibilidad de datos/código (Zenodo + GitHub), declaración de conflicto de interés.

**Artículo 2 — software, unos meses después (cuando v0.3.0 esté lista):**
- Destino: JOSS (*Journal of Open Source Software*).
- JOSS revisa el **repositorio en GitHub en abierto**, no solo el texto (~250–1000 palabras): exige README con instalación + ejemplo de uso, suite de pruebas (`tests/`, todavía no existe), documentación mínima de API, licencia OSI (ya existe `LICENSE`), y evidencia de uso en investigación real -- citar el Artículo 1 cumple ese último punto.
- No es obligatorio tener GUI para JOSS, pero si de todas formas se construye v0.3.0, mejora la revisión.
- Enviar primero Artículo 1 (o al menos subir un preprint citable) antes de preparar JOSS, para tener la referencia de "statement of need" lista.

Material público previsto (compartido por ambos artículos):

- Repositorio GitHub separado.
- Lanzamientos versionados.
- Zenodo con DOI de versión y DOI conceptual.
- CITATION.cff, licencia, instrucciones y archivo de entorno.
- Datos derivados, tablas, scripts de figuras y manifiestos.
- No redistribuir coordenadas PDB si la práctica adecuada es recuperarlas desde sus fuentes; conservar identificadores, URLs y checksums.

**Repositorio GitHub creado (2026-09-01):** https://github.com/sircalch/dockprep-audit (público, cuenta sircalch). Commit inicial con código, manuscrito, roadmap, y las tablas/manifiestos derivados del benchmark (~16 MB, 8620 archivos) -- coordenadas PDB originales, `*.pdbqt`, `vina.log`/`receptor_prep.log` y las carpetas `raw-pdb/`/`raw_graphql/` explícitamente excluidos vía `.gitignore`, siguiendo la política de no-redistribución ya declarada en la sección Data availability del manuscrito (identificadores PDB + checksums del manifiesto congelado bastan para regenerarlos). `manuscript/draft.md` actualizado con la URL real, reemplazando el placeholder `[DATA:]`. Pendiente aún: DOI de Zenodo (sección Data availability del manuscrito, segundo placeholder `[DATA:]`).

## 14. Herramientas futuras relacionadas

Después de DockPrep Audit podrían desarrollarse, sin fragmentar artificialmente publicaciones:

- **DockRepro:** ejecución y reporte reproducible de docking.
- **PoseCompare:** comparación de poses, RMSD e interacciones entre motores.
- **PocketGuide:** selección y documentación de sitio/caja.
- **DockBench:** ejecución de benchmarks comparativos.
- **Bioassay Dataset Auditor:** detección de duplicados, fugas y errores en conjuntos QSAR/ML.
- **claimtestR Studio:** interfaz visual para afirmaciones estadísticas ejecutables.

Primero debe consolidarse DockPrep Audit; las demás herramientas solo deberían separarse si resuelven preguntas claramente diferentes.

## 15. Pendientes inmediatos al retomar

Orden recomendado:

1. ~~Descargar el binario oficial `vina_1.2.7_win.exe` dentro de una carpeta controlada y registrar URL, tamaño y SHA-256.~~ **Hecho 2026-08-21** (ver sección 6).
2. ~~Confirmar si Meeko puede cubrir la ruta de receptor requerida sin ProDy o definir una alternativa reproducible con Gemmi.~~ **Hecho 2026-08-21**: Meeko cubre la ruta vía `Polymer.from_pdb_string()`, sin ProDy.
3. Revisar visualmente los 12 complejos y documentar ensamblaje/cadena. *(pendiente, es revisión humana)*
4. Congelar `benchmark/pilot_cases.csv` y su checksum. *(pendiente, depende del paso 3)*
5. ~~Implementar extracción reproducible de receptor y ligando (contra la API de Meeko, no su CLI).~~ **Hecho 2026-08-21**: `scripts/smoke_redock_case.py` extrae receptor (ATOM, cadena declarada, altLoc en blanco/"A" con la bandera limpiada) y ligando (HETATM del componente/cadena/resseq declarados), asigna órdenes de enlace del ligando desde el SMILES del componente químico de RCSB (`AllChem.AssignBondOrdersFromTemplate`), y prepara ambos PDBQT con Meeko sin ProDy.
6. ~~Preparar y ejecutar **un solo caso de humo**, sin incorporarlo aún como resultado científico.~~ **Hecho 2026-08-21**: caso 1A28 (control de bajo riesgo) preparado y redocado con Vina 1.2.7.
7. ~~Verificar RMSD y conservación de correspondencia atómica.~~ **Hecho 2026-08-21**: RMSD de átomos pesados pose 1 vs. cristal = **0.383 Å** (23/23 átomos, correspondencia de índice 1:1 verificada), muy por debajo del umbral de éxito planeado (≤ 2 Å).
8. Solo después ejecutar los 12 casos. *(pendiente — ver nota de robustez abajo)*

**Diagnóstico técnico ampliado, 2026-08-21 (NO son resultados científicos):** se corrió `scripts/smoke_redock_case.py` sobre los 12 casos con una única política técnica arbitraria (altLoc en blanco/"A", solo registros ATOM de la cadena declarada, sin metales ni aguas, caja centrada en el ligando con padding de 20 Å, exhaustividad 8, semilla 42, una sola repetición). Esto **no** corresponde a ninguna de las tres políticas de la Fase B ni tiene controles, y no debe citarse como hallazgo del artículo. Sirve solo para medir qué tan robusto es el pipeline técnico.

Se encontró y corrigió un bug de correspondencia de átomos: el escritor de árbol de torsión de Meeko reordena los átomos del ligando en el PDBQT respecto al PDB de entrada (confirmado: es un no-op para un ligando rígido como progesterona, pero reordena uno flexible como retinol), y el remark `REMARK SMILES IDX` que parecía ofrecer la corrección resultó referirse al orden atómico de la cadena SMILES, no al de entrada — dio una correspondencia incorrecta. La solución verificada (distancia 0.000 Å) es emparejar átomos por coincidencia exacta de coordenadas contra `ligand.pdbqt`, que Meeko nunca desplaza respecto al PDB original.

RMSD de átomos pesados de la pose 1 vs. cristal, con esta política técnica única:

`scripts/smoke_redock_case.py` ahora calcula el RMSD de las 9 poses (no solo la top-1) y clasifica cada caso con el esquema de 3 vías estándar en benchmarks de redocking (usado p. ej. en evaluaciones del conjunto Astex Diverse): **success** (pose #1 ≤ 2 Å), **scoring_fail** (alguna pose ≤ 2 Å pero no la #1 — falla de la función de puntuación de Vina, no del sitio de unión) y **sampling_fail** (ninguna pose ≤ 2 Å — puede indicar que la política de preparación rompió el sitio de unión real). También registra qué residuos de receptor fueron ignorados por Meeko (`receptor_residues_ignored`), en línea con el principio de no ocultar decisiones de preparación.

| PDB | Estrato | RMSD #1 (Å) | Mejor RMSD (rank) | Resultado | Residuos ignorados |
|---|---|---|---|---|---|
| 1T46 | Conformación alternativa | 0.331 | 0.331 (#1) | **success** | ninguno |
| 5A2S | Metal/cofactor | 0.569 | 0.569 (#1) | **success** | 10 |
| 1CBX | Metal/cofactor | 0.682 | 0.682 (#1) | **success** | ninguno |
| 1OHR | Política de aguas | 0.816 | 0.816 (#1) | **success** | 15 |
| 1QKT | Control de bajo riesgo | 0.866 | 0.866 (#1) | **success** | ninguno |
| 1A28 | Control de bajo riesgo | 0.383 | 0.383 (#1) | **success** | 5 |
| 1RBP | Control de bajo riesgo | 0.694 | 0.694 (#1) | **success** | 1 |
| 1M17 | Conformación alternativa | 5.921 | 1.335 (#4) | scoring_fail | ninguno |
| 4RJ3 | Conformación alternativa | 7.236 | 0.282 (#5) | scoring_fail | ninguno |
| 4EXS | Metal/cofactor | 4.865 | 2.427 (#5) | sampling_fail (límite) | ninguno |
| 3FNU | Política de aguas | 7.812 | 2.974 (#3) | sampling_fail (límite) | ninguno |
| 4GID | Política de aguas | 9.831 | 6.926 (#2) | **sampling_fail** | ninguno |
| ~~1CPS~~ | Metal/cofactor | — | — | sustituido por 1CBX el 2026-08-21 | — |

Lectura con esta clasificación (sigue sin ser un resultado científico — política técnica única, sin controles):
- Los 3 controles de bajo riesgo y 4 de los 5 casos altLoc/metal/agua fueron `success` directo.
- `1M17` y `4RJ3` (conformación alternativa) son `scoring_fail`: Vina sí muestreó una pose casi nativa (1.34 Å y 0.28 Å) pero no la puso en el top-1 — comportamiento documentado de la función de puntuación de Vina, no un problema de nuestra preparación.
- `4EXS` y `3FNU` quedan justo en el límite (2.43 Å y 2.97 Å) — no concluyente con exhaustividad 8 y una sola semilla.
- `4GID` es el único `sampling_fail` claro (mejor pose a 6.93 Å): ninguna de las 9 poses se acercó al sitio real. Es el caso que más amerita revisión visual antes de congelar.

Fuentes usadas para diseñar este diagnóstico: guía oficial de Vina para docking con zinc ([autodock-vina.readthedocs.io/en/latest/docking_zinc.html](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html)); AutoDock4Zn como forcefield especializado para metaloenzimas ([PubMed 24931227](https://pubmed.ncbi.nlm.nih.gov/24931227/)); brecha documentada entre tasa de muestreo (~93%) y tasa de acierto en el top-1 (~35–40%) de Vina en el conjunto Astex Diverse ([PMC6102569](https://pmc.ncbi.nlm.nih.gov/articles/PMC6102569/)).

**Recomendación para la Fase B (metales):** no usar Vina "de fábrica" para los casos metal/cofactor; adoptar la guía oficial de zinc de AutoDock Vina o el forcefield AutoDock4Zn en vez de improvisar una solución propia, ya que la limitación de Vina estándar con geometría tetraédrica del zinc está bien documentada.

**Revisión estructural de `4GID`, 2026-08-21 (el único `sampling_fail` claro):** no se pudo usar el visor 3D de RCSB en esta sesión (el panel del navegador no renderizó), así que la revisión se hizo directamente sobre las coordenadas del PDB depositado. Hallazgo: `4GID` tiene un **segundo ligando cocristalizado**, `HET LPD A 502` (L-prolinamida, un fragmento de 8 átomos), a ~10 Å del ligando principal declarado (`0GH`), en la misma cadena. Nuestra extracción solo conserva el componente declarado y los registros `ATOM` del receptor, así que `LPD` desaparece por completo — y se verificó que **4 de sus 8 átomos caen dentro de la caja de docking** calculada (bbox de 0GH + padding de 20 Å). Es decir, la caja de búsqueda incluye una cavidad que en el cristal real está ocupada por otra molécula pero que en nuestra preparación queda vacía, dándole a Vina espacio artificial para explorar poses que no existen en el complejo real. Se descartaron otras dos hipótesis: la unidad biológica es monomérica (cadena A sola es correcta, confirmado por `REMARK 350`) y la ocupancia de cadena A es normal (10/3043 átomos <1.0, sin anomalías).

Esto **no se corrigió automáticamente** (sería una reparación silenciosa); queda como decisión explícita pendiente de Fase B: ¿la caja debe recortarse para no invadir el espacio de co-ligandos secundarios cercanos, o deben tratarse como parte del contexto del receptor? Aplica en principio a cualquier caso con más de un HETATM no acuoso/no metálico cerca del sitio — vale la pena revisar los otros 11 casos por la misma razón antes de congelar el manifiesto.

**Generalización a los 12 casos, 2026-08-21:** se buscaron heterogrupos no acuosos/no metálicos distintos del ligando declarado, en la misma cadena, dentro de 15 Å del ligando. Aparte de `4GID` (`LPD`, ver arriba), aparecieron dos casos más:

- **`1T46`**: dos iones fosfato (`PO4`) a 13.08 Å y 7.83 Å — casi con certeza artefactos del buffer de cristalización, no relevantes biológicamente. Bajo riesgo.
- **`4RJ3`**: **`ALY A:33` a solo 3.49 Å del ligando** — y esto NO es un ligando externo. `ALY` es N6-acetil-lisina, un **residuo de aminoácido modificado que forma parte de la cadena principal de la proteína** (confirmado: aparece en `SEQRES` entre `LEU 32` y `LYS 34`, con átomos de backbone N/CA/C/O completos). El PDB lo registra como `HETATM` por ser no estándar. **Nuestro script de extracción del receptor solo conserva registros `ATOM`, así que borra por completo este residuo del backbone** — no es una decisión de política, es un hueco real en la cadena que probablemente explica por qué Meeko reportó fallas de plantilla en residuos vecinos por enlaces inter-residuo rotos en otros casos (ver `1A28`, sección de arriba: el mismo mecanismo — quitar un residuo intermedio deja huérfano el enlace peptídico de sus vecinos).

**Esto sí es un bug de extracción a corregir, no una decisión científica pendiente**: la extracción de receptor debe incluir residuos modificados que sean parte continua de la cadena polimérica (identificables por continuidad de `SEQRES`/numeración y presencia de átomos de backbone), no descartarlos junto con aguas/ligandos reales. La decisión científica aparte (¿mapear `ALY` a `LYS` estándar vía `--set_template`, dejarlo como residuo modificado, o eliminarlo explícitamente vía `--delete_residues`?) sigue siendo tuya — Meeko soporta las tres rutas de forma explícita y documentada, ninguna es automática.

**Corregido, 2026-08-21:** `extract_receptor_atoms()` en `scripts/smoke_redock_case.py` ahora detecta residuos HETATM con backbone completo (N, CA, C) continuos con la cadena y los incluye en la extracción del receptor (función nueva `find_polymer_hetatm_residues`), en vez de descartarlos junto con aguas/ligandos. Verificado en `4RJ3`: `ALY A:33` ahora aparece completo y correctamente parametrizado en `receptor.pdbqt` (Meeko sí tiene una plantilla nativa para este residuo modificado — cero advertencias). Se re-corrieron los 12 casos: sin regresiones; los residuos ignorados que quedan en `1A28`, `5A2S`, `1OHR`, `1RBP` son distintos a este mecanismo (probablemente huecos reales de la cadena o incompletitud del depósito, no residuos HETATM continuos) y siguen documentados como antes.

**Observación adicional:** el RMSD de `4GID` cambió de 9.831 Å a 6.998 Å entre corridas con la misma semilla (42) y sin cambios relacionados en su extracción — Vina no es completamente determinista con `--cpu 0` (auto-detección de hilos) pese a fijar la semilla. Relevante para la Fase C / pendiente #11 (sección 8): considerar `--cpu 1` para reproducibilidad estricta, o aceptar la varianza y usar repeticiones con múltiples semillas como ya está planeado.

**Decisión sobre `ALY A:33` en `4RJ3`, 2026-08-21: se conserva como residuo modificado (no se mapea a LYS estándar).** `4RJ3` es CDK2 (`COMPND: CYCLIN-DEPENDENT KINASE 2`, `ENGINEERED: YES`), y el residuo 33 en la numeración de CDK2 es la **lisina catalítica invariante del motivo VAIK** (forma el par iónico con Glu51 que coordina los fosfatos del ATP). Que aparezca acetilada (`ALY`), en una estructura marcada explícitamente como modificada genéticamente, y a 3.49 Å del inhibidor, indica que la acetilación de esta lisina es casi con certeza el objeto deliberado del estudio (incorporación genética de acetil-lisina para estudiar su efecto sobre la unión de inhibidores), no un artefacto de modelado. Mapearlo a LYS estándar sería una reparación silenciosa que redocaría contra un bolsillo electrostáticamente distinto al que existió en el cristal real, invalidando la propia verificación de RMSD. Se mantiene como `ALY` (Meeko ya lo parametriza nativamente sin errores). **Nota para la revisión visual pendiente de `4RJ3`:** esta modificación puede ser una segunda característica estructural relevante del caso, además del hallazgo de conformación alternativa que motivó su inclusión en ese estrato — vale la pena tenerlo presente al revisar el caso.

**Revisión de ensamblaje biológico y cadena de los 12 casos, 2026-08-21 (pendiente #1 de la sección 8):** no se pudo usar el visor 3D (el panel del navegador no renderizó en esta sesión); la confirmación se hizo cruzando `REMARK 350` (unidad biológica según autor y según software) contra las cadenas usadas en `benchmark/pilot-eligibility/eligibility_results.csv`.

- **11/12 confirmados correctos:** `1M17`, `1T46`, `4RJ3`, `1CBX`, `1RBP` (monoméricos, cadena A); `5A2S`, `4EXS` (dos copias monoméricas independientes en la unidad asimétrica, no un dímero real — B es redundante); `4GID` (4 copias monoméricas independientes); `1A28`, `1QKT` (dominio de unión a hormona de receptor nuclear — la unidad biológica completa es dimérica, pero cada monómero tiene su propio bolsillo de ligando independiente, por lo que modelar solo A es correcto); `1OHR` (dímero real tipo proteasa del VIH, ya usa correctamente A;B).
- **`3FNU` requiere decisión explícita, sin resolver:** `REMARK 350` confirma (autor Y software de acuerdo) que A+B es un dímero real, igual que C+D — HAP (Histo-Aspartic Protease de *P. falciparum*) es de la misma familia que la proteasa del VIH (mismo caso que `1OHR`, que sí usa A;B). Pero en este depósito específico: (1) el ligando `006` aparece en 4 copias independientes, una por cadena (A, B, C, D), no una compartida entre A y B; (2) el residuo catalítico His32 de la cadena B está a 26.7 Å del ligando de la cadena A (nada cerca, a diferencia de 1OHR donde ambas cadenas quedan dentro de 6 Å). La geometría de este caso concreto sugiere un bolsillo autocontenido en la cadena A, aunque la arquitectura biológica de la familia sea dimérica. **Decisión pendiente del usuario:** ¿extraer el receptor de `3FNU` solo con cadena A (como usa actualmente `pilot-eligibility`) o con A+B (como exige la arquitectura biológica real)? Ambas posturas son defendibles; no se resolvió unilateralmente.

**Resuelto 2026-08-21: se usa A+B, como el dímero real.** Se implementó un mecanismo explícito de anulación (`benchmark/receptor_chain_overrides.csv`, aplicado por `scripts/verify_pilot_eligibility.py`) que documenta la razón y conserva visible el hallazgo geométrico original (`receptor_chains_geometric=A`) junto al valor aplicado (`receptor_chains=A;B`) — nada se sobrescribe en silencio. Se re-corrieron los 12 casos (12/12 pass) y el caso de humo técnico de `3FNU` con ambas cadenas: pipeline limpio, sin residuos ignorados, resultado `sampling_fail` (mejor RMSD 3.071 Å, rank 5) — similar al valor con solo A (2.974 Å); el objetivo de este cambio era la corrección biológica del receptor, no mejorar el RMSD de este caso concreto.

**Confirmación visual/cuantitativa de identidad y pose del ligando, 12 casos, 2026-08-21 (pendiente #2 de la sección 8):** no se pudo usar el visor 3D de forma sostenida; se usaron las métricas de validación wwPDB vía la API de datos de RCSB (`RSCC`, `RSR`, outliers geométricos Mogul, choques intermoleculares, ocupancia) para cada instancia de ligando declarada — el estándar cuantitativo cristalográfico, más riguroso que una inspección visual.

- **Identidad química confirmada en 12/12**: el nombre/fórmula real en RCSB coincide con lo esperado para cada caso (AQ4≈análogo de erlotinib en 1M17, STI=imatinib en 1T46, 3QS en 4RJ3, OTF en 5A2S, BZS=L-benzilsuccinato en 1CBX, X8Z=L-captopril en 4EXS, 1UN en 1OHR, 006 en 3FNU, 0GH en 4GID, STR=progesterona en 1A28, EST=estradiol en 1QKT, RTL=retinol en 1RBP).
- **Ocupancia 1.0 en los 9 casos con datos** — instancia única, completamente modelada, consistente con el propio motor de auditoría.
- **RSCC/RSR disponibles y aceptables en 9/12** (rango 0.864–0.947 de RSCC, todos por encima del umbral típico ~0.8 de fiabilidad).
- **Limitación real en 3/12** (`1CBX` 1988, `1OHR` 1997, `1RBP` 1994): depósitos anteriores a la exigencia rutinaria de subir factores de estructura; RCSB no puede calcular RSCC/RSR. Se confía en las coordenadas de los autores originales sin cruce automático moderno — mencionar como limitación en el manuscrito, no es corregible.
- **`1M17` es el caso más débil** (RSCC 0.866, 12 outliers de enlace + 13 de ángulo, 11 choques intermoleculares) — sigue por encima del umbral aceptable, pero con menos margen que el resto.
- **`3FNU`**: las 4 copias cristalográficas (A, B, C, D) dan RSCC consistente entre sí (0.864–0.869), buena señal de reproducibilidad de la pose entre copias independientes.
- **`4GID`**: se confirma numéricamente lo ya encontrado antes — el ligando principal `0GH` tiene RSCC 0.94 (bien soportado), mientras que `LPD` (el co-ligando ignorado por la extracción) tiene RSCC más bajo (0.65–0.76), consistente con ser una molécula secundaria de menor prioridad en el refinamiento, no el objeto principal del estudio.

El pipeline técnico corre de extremo a extremo en 11/12 casos. `1CPS` falla de forma reproducible durante la percepción de enlaces de RDKit en el residuo `TYR A:204`: su oxígeno de carbonilo queda simultáneamente a ~1.23 Å de C, ~1.28 Å de CA y ~1.69 Å de N del mismo residuo — geométricamente incompatible con un carbonilo real. Es un problema de geometría del propio depósito PDB (posible error de refinamiento en una estructura antigua), no un error del script. Es precisamente el tipo de hallazgo que DockPrep Audit debe exponer, no corregir en silencio; no se aplicó ninguna heurística de reparación.

Patrón observado (no concluyente, sin controles ni repeticiones): los 3 controles de bajo riesgo reprodujeron la pose casi exactamente (< 1 Å); los casos de conformación alternativa, metal y agua tuvieron resultados mixtos, con varios RMSD altos. Esto es consistente con la hipótesis de trabajo del proyecto, pero **no puede usarse como evidencia** hasta la Fase B/C con manifiesto congelado, comparación de políticas, semillas múltiples y controles apropiados.

**Antes de correr los 12 casos como resultado científico**, hay que decidir: (a) ~~qué hacer con `1CPS`~~ **resuelto 2026-08-21: sustituido por `1CBX`**, (b) si el manejo de errores por caso restante (residuos con átomos faltantes, gaps de cadena) se resuelve caso por caso o mediante una política uniforme documentada en el manifiesto.

**Política de altLoc, decidida 2026-08-21 (pendiente #5 de la sección 8):** por cada residuo con conformaciones alternativas, se conserva la de mayor ocupancia reportada; en empate exacto se usa `A` como desempate declarado (no arbitrario). Esta regla **es una sola y se aplica igual en las tres políticas de preparación de la Fase B** — es una decisión de fidelidad estructural (qué conformación existió realmente en el cristal), no de filosofía de preparación, así que no forma parte del eje de comparación entre políticas; ese eje sigue siendo aguas/metales/protonación, como ya estaba planeado en la sección 9.

Implementada en `choose_altloc_conformers()` dentro de `scripts/smoke_redock_case.py`, con registro transparente por residuo (ocupancias vistas, confórmero elegido, base de la decisión) en `altloc_decisions.json` por caso — nunca una elección silenciosa. Verificación empírica en los 3 casos del estrato "conformación alternativa": en `1T46` la ocupancia real decide (`A`=0.65 vs `B`=0.35, `basis=occupancy`); en `1M17` y `4RJ3` todos los residuos con altLoc están empatados exactos 0.5/0.5 (`basis=tie_break_A`). Es decir, para esta cohorte, "quedarse con A" y "quedarse con la mayor ocupancia" dan el mismo resultado, pero la regla implementada es la generalizable (ocupancia primero), no la coincidencia alfabética. Se re-corrieron los 12 casos: sin regresiones, resultados de RMSD/outcome idénticos a la corrida anterior.

**Política de aguas, decidida 2026-08-21 (pendiente #6 de la sección 8):** a diferencia de altLoc, esta política **sí varía entre las tres políticas de preparación de la Fase B** — es el eje de comparación que ya estaba planeado en la sección 9.

- **Política 1 (referencia conservadora):** se conserva únicamente el agua "puente" — oxígeno a ≤3.0 Å de un átomo del ligando **y** ≤3.0 Å de un átomo del receptor. Es el criterio "loose water" usado en estudios de docking (estilo GOLD/Verdonk), más estricto que el filtro de 4 Å a cualquier átomo del ligando que el motor de auditoría ya usa solo para *señalar* candidatos para revisión humana, no para decidir qué conservar en la preparación.
- **Política 2 (simplificada convencional):** se eliminan todas las aguas sin excepción.
- **Política 3 (flujo abierto predeterminado):** coincide con la Política 2 en este pipeline (Meeko + Vina) — la extracción de receptor ya descarta aguas por defecto sin intervención adicional, igual que la mayoría de tutoriales/pipelines abiertos estándar. No se fuerza una diferencia artificial solo por tener tres casillas distintas.

**Validación empírica:** se aplicó el criterio de agua puente a los 12 casos. En `1OHR` identificó exactamente la famosa "agua de la flap" de la proteasa del VIH (`A:303`, contacta `ILE 50` N/H/O — el puente Ile50↔inhibidor descrito ampliamente en la literatura de este sistema), la mejor confirmación posible de que la regla funciona. Encontró 3 aguas puente en `1OHR`, 1 en `4GID`, 1 en `1QKT` (control de bajo riesgo, consistente con lo ya anotado en la sección 5), 1 cada uno en `1M17`/`1T46`/`4RJ3` (aguas puente existen fuera del estrato "agua" también, sin contradecir la clasificación por estrato dominante), y **ninguna en `3FNU`** — el caso que motivó el estrato "política de aguas" no tiene ninguna agua puente bajo este criterio más estricto (el hallazgo original de 4 Å laxo sigue siendo válido como señal de auditoría, pero no sobrevive como base para retener agua en la preparación real). Se decidió no reabrir la composición de la cohorte por esto; queda anotado como limitación conocida de `3FNU` bajo la Política 1.

Implementado en `scripts/smoke_redock_case.py`: `find_bridging_waters()` (detección) y `append_bridging_waters_to_pdbqt()` (añade las aguas retenidas al receptor preparado como átomos de oxígeno rígidos, tipo AutoDock `OA`, carga parcial -0.834 e de TIP3P — Meeko no tiene plantilla nativa para `HOH` porque el agua no es un residuo polimérico, así que no puede pasar por su pipeline de plantillas químicas; añadirla directamente al PDBQT es la práctica estándar para tratar unas pocas aguas ordenadas como parte de un receptor rígido de AutoDock/Vina). Nuevo flag `--water-policy {none,conservative}` (por defecto `none`, sin cambios de comportamiento respecto a las corridas anteriores). Verificado extremo a extremo en `1OHR`: con `--water-policy conservative` se añaden las 3 aguas puente al `receptor.pdbqt` y Vina corre sin error (RMSD top-1 0.791 Å, `success`, comparable al resultado sin agua).

**Política de metales, decidida 2026-08-21 (pendiente #7 de la sección 8):** a diferencia del agua, para los metales **no hay una alternativa defendible de "eliminarlos"** — un ion de zinc catalítico en una metaloenzima no es una decisión de filosofía de preparación, es parte estructural obligatoria del sitio activo. Por eso, igual que con altLoc, esta regla **es una sola y se aplica igual en las tres políticas de la Fase B**, no es un eje de comparación.

**Hallazgo de bug al investigar esto:** Meeko ya tiene plantillas químicas nativas para iones metálicos monoatómicos con la carga formal correcta (`ZN`→`[Zn+2]`, `MG`→`[Mg+2]`, `CA`→`[Ca+2]`, `MN`→`[Mn+2]`, `FE`→`[Fe+3]`; confirmado en `meeko/data/residue_chem_templates.json`). El problema no era que Meeko no supiera manejar metales — era que **nuestra propia extracción los descartaba**, exactamente el mismo mecanismo que el bug de `ALY` (residuo HETATM sin átomos de backbone N/CA/C, filtrado por el criterio "solo ATOM"). Corregido con `find_meeko_supported_metal_residues()` en `scripts/smoke_redock_case.py`, que identifica estos residuos y los incluye junto con los HETATM poliméricos.

**No se implementó el pseudo-átomo tetraédrico de AutoDock4Zn** (`zinc_pseudo.py` + `AD4Zn.dat` + `autogrid4`, documentado en la guía oficial de Vina para zinc, https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html) **como método principal**, porque exigiría cambiar la función de puntuación a `ad4` solo para los 3 casos de metal, introduciendo una variable de confusión: cualquier diferencia de RMSD podría deberse al motor de puntuación distinto, no a la política de preparación que es lo que realmente se quiere comparar. Se optó por mantener la función de puntuación de Vina uniforme en los 12 casos, representando el metal como un átomo rígido no enlazado con su carga formal correcta pero sin geometría de coordinación explícita — una limitación documentada, no oculta.

**Resultado empírico honesto (política técnica única, sin controles):** incluir el zinc con carga correcta pero sin geometría de coordinación **no mejoró uniformemente los resultados** frente a omitirlo — `5A2S` se mantuvo prácticamente igual (0.569→0.615 Å, success), `1CBX` empeoró pero sigue siendo éxito (0.682→1.391 Å), y `4EXS` empeoró (mejor pose 2.427→4.323 Å, sigue `sampling_fail`). Esto es consistente con lo documentado en la literatura: la función de puntuación estándar de Vina/AutoDock4 sin la extensión AD4Zn tiene dificultades reales con la coordinación tetraédrica del zinc. **Recomendación registrada para Fase C:** si el manuscrito necesita resultados más rigurosos para el estrato metal/cofactor, correr el flujo AD4Zn completo como análisis de sensibilidad declarado por separado, no como parte de la comparación principal de políticas.

Se re-corrieron los 12 casos con ambas correcciones (altLoc + metal): 12/12 sin errores.

**Protonación del receptor, decidida parcialmente 2026-08-21 (pendiente #8 de la sección 8):** Meeko no tiene información experimental de hidrógenos, así que para residuos ambiguos (HIS entre HID/HIE/HIP; ASP/ASH; GLU/GLH; LYS/LYN) elige silenciosamente la primera opción de una lista de prioridad fija (confirmado en `residue_chem_templates.json`: `HIS -> [HIE, HID, HIP, ...]`, gana la primera coincidencia). Para la mayoría de residuos esto es el default razonable y estándar del campo (Asp/Glu desprotonados, Lys protonada, His neutra) — se deja igual, compartido por las tres políticas de la Fase B (igual que altLoc/metales, es fidelidad estructural, no filosofía de preparación).

**Excepción verificada y corregida — histidinas coordinantes de metal:** el tautómero de una histidina que coordina un metal (Zn) no puede decidirse a ciegas — depende de cuál nitrógeno (ND1 o NE2) está físicamente coordinando el ion en el cristal (esa posición debe quedar sin protonar). Se verificó geométricamente (distancia ≤2.6 Å) en los 3 casos de metal: **8 de 12 histidinas coordinantes (contando ambas cadenas) requerían `HID`, no el `HIE` que Meeko asigna por defecto** — un error sistemático real que el flujo estándar habría cometido en silencio. Implementado en `assign_metal_coordinating_histidine_tautomers()` (`scripts/smoke_redock_case.py`): determina el tautómero correcto por geometría y reescribe el nombre del residuo (`HIS`→`HIE`/`HID`) antes de pasarlo a Meeko, con registro transparente en `his_tautomer_decisions.json` por caso. Verificado extremo a extremo: los 12 casos corren sin error; `5A2S` y `1CBX` sin cambio de RMSD (coincidía con el default en 3 de 4 residuos), `4EXS` cambió levemente (4.846→4.89 Å, sigue `sampling_fail`).

**Protonación del ligando, resuelta 2026-08-21 (pendiente #9 de la sección 8):** `scrubber` (la herramienta que recomienda la guía oficial de Vina para zinc) no instala en Python 3.12 (`ModuleNotFoundError: No module named 'imp'`, módulo eliminado en 3.12; su `setup.py` es anterior a ese cambio). En vez de forzar esa instalación o improvisar reglas propias, se instaló **`dimorphite-dl`** — el motor de protonación subyacente, activamente mantenido, más liviano, y el mismo enfoque metodológico (asignación de estado de ionización dominante a un pH dado) que usa `scrubber` internamente. Instala limpio en Python 3.12 (aunque fija `rdkit<2026`, bajando nuestra versión de 2026.3.5 a 2025.9.6 — verificado que no cambia ningún resultado: `1CBX` dio exactamente el mismo RMSD, 1.391 Å, antes y después del downgrade).

**Regla implementada:** cada ligando se protona a pH fisiológico (7.4) con `dimorphite_dl.protonate_smiles()` sobre el SMILES canónico de RCSB, antes de usarlo como plantilla para asignar órdenes de enlace sobre las coordenadas cristalográficas. Verificado el impacto real: `BZS` (ácido L-bencilsuccínico, `1CBX`) está depositado neutro (ambos -COOH, carga formal 0) pero a pH 7.4 se asigna correctamente como dicarboxilato (-COO⁻ ambos, carga -2) — el RMSD cambió de 1.391 a 1.551 Å (sigue `success`), reflejando el cambio químico real, no ruido.

**Caso límite encontrado y resuelto — `1OHR` (nelfinavir):** tiene un fenol y una amina terciaria con pKa cercanos a 7.4, así que `dimorphite-dl` devolvió **4 variantes ambiguas** en vez de una sola. La primera falló al construir la molécula (`Explicit valence for atom # 31 O, 2`, un error interno de la propia herramienta). Se implementó un mecanismo de reintento documentado: probar cada variante en orden hasta que una construya y sanee correctamente, registrando qué falló y por qué (`ligand_protonation.json` por caso) — **no** es un re-ranking por resultado de docking, solo un filtro de validez química. La variante que funcionó (fenol neutro, amina protonada `[NH+]`) coincide con la química esperada de nelfinavir en literatura farmacológica.

Se re-corrieron los 12 casos: **12/12 sin errores**, cada uno con su registro de protonación transparente.

**Independencia por blanco biológico, resuelta 2026-08-21 (pendiente #3 de la sección 8):** el registro de exclusión de estudios previos (132 registros) está indexado solo por PDB ID, así que no detectaba si un PDB *distinto* correspondía al *mismo* blanco biológico ya usado en un manuscrito previo. Se verificó comparando accesos UniProt (inequívoco, independiente de qué entrada PDB se depositó): **6 de los 12 casos comparten UniProt con un registro excluido**, los 6 provenientes del mismo estudio previo (`dude_receptor_prep_audit`, la auditoría de 102 blancos DUD-E) — `1A28`↔`3KBA` (receptor de progesterona), `1M17`↔`2RGP` (EGFR), `1QKT`↔`1SJ0` (receptor de estrógeno), `1T46`↔`3G0E` (c-KIT), `4GID`↔`3L5D` (BACE1), `4RJ3`↔`1H00` (CDK2).

**Decisión: se conservan los 12 casos, se documenta el traslape en vez de reabrir la cohorte.** Razones: (1) son entradas PDB distintas, no datos reciclados; (2) la pregunta de aquel estudio (prevalencia de reglas de auditoría en 102 blancos DUD-E) es categóricamente distinta a la de este (¿la política de preparación cambia la reproducibilidad del redocking?); (3) exigir blancos 100% nuevos para un piloto estratificado por riesgo estructural es casi inviable — la riqueza estructural visible (altLoc/agua/metal) se concentra justo en los blancos más estudiados del campo. Reemplazar 6 de 12 casos ya congelados y verificados se juzgó peor que la transparencia. **Debe declararse como limitación explícita en el manuscrito.** Detalle completo en `benchmark/PILOT_SELECTION.md` § Biological-target independence; `benchmark/pilot_manifest_frozen.csv` ahora incluye las columnas `target_uniprot`, `shares_target_with_prior_pdb`, `shares_target_prior_study` (checksum actualizado: `716cba92f1863153758dc454dde235fd3813351210e48545f1e39091f4c0fb35` — la composición de la cohorte no cambió, solo se añadió esta trazabilidad).

**Caja de docking, decidida 2026-08-21 (pendiente #10 de la sección 8):** se reemplazó la heurística arbitraria usada durante la validación del pipeline (bbox del ligando + padding de 20 Å) por una regla con respaldo bibliográfico: **caja cúbica de tamaño 2.9 × radio de giro (Rg) del ligando**, centrada en su centro geométrico — el tamaño que Feinstein & Brylinski (2015, [PMC4468813](https://pmc.ncbi.nlm.nih.gov/articles/PMC4468813/)) encontraron que maximiza la recuperación de sitio/pose (RMSD 4.9→4.0 Å, recuperación de residuos de sitio 0.78→0.92 frente a una caja ad hoc). Implementado en `compute_ligand_box()` (`scripts/smoke_redock_case.py`); el padding antiguo queda disponible solo como parámetro opcional para reproducir corridas previas, no como política.

**Resultado honesto al re-correr los 12 casos — mixto, no una mejora uniforme:**

| PDB | Antes (bbox+20Å) | Ahora (2.9×Rg) | Cambio |
|---|---|---|---|
| `4RJ3` | scoring_fail (5.316 Å) | **success (0.306 Å)** | mejora |
| `1QKT` | success (0.866 Å) | scoring_fail (6.243 Å) | empeora |
| `1M17` | scoring_fail (5.921 Å) | sampling_fail (7.942 Å) | empeora |
| `4EXS` | sampling_fail (mejor 2.744 Å) | sampling_fail (mejor 2.482 Å, rank 1) | mejora, sigue fallando |
| `5A2S` | success (0.615 Å) | success (1.635 Å) | peor pero sigue éxito |
| resto (6 casos) | sin cambio de categoría | | |

El total de `success` se mantuvo en 7/12, pero la composición cambió: se gana `4RJ3`, se pierde `1QKT`. Esto es evidencia real de que **la caja interactúa con el ranking de Vina de forma no trivial por caso** — no se intentó ajustar más allá de aplicar la regla decidida; no se buscó la combinación que diera "mejores números", eso sería ajustar la política a los resultados, exactamente lo que el proyecto prohíbe.

**Semillas, exhaustividad, número de modos y repeticiones, decidido 2026-08-21 (pendiente #11 de la sección 8):** partiendo del hallazgo ya documentado de que Vina no es completamente determinista con `--cpu 0` pese a fijar la semilla (`4GID` varió 9.831→6.998 Å entre corridas idénticas), la Fase C real **no puede reportar una sola corrida por caso** — eso ocultaría varianza real del motor de búsqueda, no solo de la política de preparación.

- **Exhaustividad:** subir de 8 (usado solo para las pruebas de humo de esta sesión, por velocidad) a **32** para la Fase C real. Costo verificado: ~22 s por caso en `1CBX` (15 átomos pesados) — perfectamente viable para 12 casos × 3 políticas × repeticiones.
- **Número de modos:** mantener el default de Vina (9) — ya es la base de la clasificación de 3 vías (`success`/`scoring_fail`/`sampling_fail`) implementada en la sección 15.
- **Repeticiones y semillas:** cada combinación caso×política se correrá con **3 semillas predeclaradas fijas: 42, 123, 2024** (no aleatorias en el momento, para que la corrida sea reproducible) — no una sola corrida. El resumen por caso×política será la **mediana** de las 3 repeticiones (más robusta que la media frente a una corrida atípica); las 3 corridas completas quedan siempre registradas, nunca colapsadas en silencio a un solo número.

**Criterios de éxito, exclusión y fallo, decididos 2026-08-21 (pendiente #12 de la sección 8) — predeclarados antes de observar resultados de la Fase C:**

- **Éxito (primario):** RMSD de átomos pesados de la pose #1 vs. la pose cristalográfica ≤ 2.0 Å — umbral estándar en la literatura de redocking.
- **Clasificación secundaria (diagnóstica, ya implementada):** `success` (pose #1 ≤ 2 Å), `scoring_fail` (alguna de las 9 poses ≤ 2 Å pero no la #1), `sampling_fail` (ninguna ≤ 2 Å).
- **Distinción crítica — fallo de preparación vs. fallo de docking, nunca deben confundirse:**
  - *Fallo de preparación* (técnico, descalifica el caso): la estructura no puede convertirse en un receptor o ligando válido para **alguna** de las tres políticas por un problema de calidad de datos del propio depósito (ejemplo real: `1CPS`, geometría de carbonilo imposible en `TYR A:204`). Un caso con fallo de preparación se sustituye **antes** de correr docking, documentando la razón técnica exacta — igual que ya se hizo con `1CPS`→`1CBX`. Nunca se decide por el resultado de RMSD.
  - *Fallo de docking* (`scoring_fail`/`sampling_fail`): **esto NO es una razón de exclusión — es el dato mismo que el estudio busca medir.** Excluir casos con RMSD alto porque "salieron mal" sería exactamente el sesgo de selección que el proyecto prohíbe explícitamente (sección 7: "ningún hallazgo se tratará como... resultados del nuevo artículo no reutilizarán..."). Un `sampling_fail` en la política "simplificada" y un `success` en la "conservadora" para el mismo caso es, potencialmente, el resultado central del artículo.
- **Comparación predeclarada:** para cada caso, comparar las 3 políticas de forma pareada (mismo caso, distinta política) en RMSD de la pose #1 y en categoría de 3 vías; agregar por estrato de riesgo estructural para responder la pregunta de investigación central (sección 1). Separar siempre el análisis confirmatorio (los 12/60 casos predeclarados) de cualquier exploración posterior, como ya establece la sección 9, Fase E.

## 16. Otros pendientes del programa de investigación

### claimtestR

- `claimtestR` 0.1.2 está publicado en GitHub y Zenodo.
- La comprobación tipo CRAN terminó con 0 errores, 0 advertencias y 1 nota esperada de envío nuevo.
- El envío definitivo a CRAN quedó pendiente porque el portal estaba cerrado del 5 al 19 de agosto de 2026. Debe comprobarse que reabrió, cargar el tarball y confirmar el correo del mantenedor.

### Perfiles académicos

- Añadir los DOI y artículos aceptados a ResearchGate cuando corresponda.
- Comprobar indexación en Google Scholar después de publicación; no se pueden fabricar ni garantizar citas.
- Mantener ORCID, GitHub, Zenodo y metadatos de autor consistentes.

## 17. Limitaciones a declarar en el manuscrito (piloto de 12 casos, 2026-08-22)

Redactado tras completar las Fases A–D del piloto. Ninguna de estas limitaciones invalida el trabajo; ocultarlas sí lo haría.

1. **Modelo de agua sin desolvatación.** La política "agua puente" representa el agua retenida como un átomo de oxígeno rígido con carga fija de TIP3P (−0.834 e), sin término de energía de desolvatación ni reoptimización de su posición durante el docking. La función de puntuación de Vina no está diseñada para modelar correctamente el costo/beneficio termodinámico de desplazar o conservar una molécula de agua explícita. El resultado nulo de la Fase D (ninguna diferencia sistemática entre políticas de agua) debe interpretarse **con esta limitación explícita**: no se puede distinguir si el efecto de la política de agua es genuinamente nulo o si el operacionalismo elegido es demasiado débil para detectarlo. Antes de ampliar a 60 casos con este mismo modelo, vale la pena decidir si se invierte en una representación más rigurosa (p. ej. energía de desolvatación explícita, o un método de puntuación distinto para el subconjunto de casos con agua).
2. **Metales sin geometría de coordinación explícita.** Los iones metálicos se representan con su carga formal correcta pero sin pseudo-átomos de coordinación tetraédrica (AD4Zn). Documentado con literatura de respaldo en la sección 15; el resultado en los 3 casos de metal es consistente con esta limitación conocida del campo, no con un error de nuestro pipeline.
3. **Traslape de blanco biológico con un estudio previo.** 20 de los 60 casos (33%) comparten blanco biológico real (UniProt) con un registro ya usado en la auditoría DUD-E de 102 blancos — misma familia de estudios del autor, distinta pregunta de investigación, distinta entrada PDB en todos los casos. (6 de los 12 originales; 14 de los 48 nuevos — 3 candidatos adicionales resultaron ser un falso positivo por un péptido coactivador compartido, no un traslape real de blanco, verificado directamente contra el UniProt de la estructura excluida antes de descartarlos de la cuenta.) Detalle completo en `benchmark/PILOT_SELECTION.md` § Biological-target independence.
4. ~~**Tamaño de muestra por estrato (n=3).**~~ **Resuelto 2026-08-24 mediante la ampliación a 60 casos (15 por estrato, Fase E).** Con n=15, el agregado por estrato ya no depende de un solo caso inestable (ver sección 9, "Fase D sobre los 60 casos") — la señal de política de agua es ahora distribuida (5/15 casos con efecto claro) en vez de un único caso atípico. Las conclusiones agregadas siguen sin ser un test estadístico formal (n sigue siendo pequeño para inferencia frecuentista estricta), pero ya no son solo "exploratorias" en el sentido de que un caso ruidoso pueda dominar el promedio.
5. **Cobertura de la validación estructural en depósitos antiguos.** 3 de 12 casos (`1CBX` 1988, `1OHR` 1997, `1RBP` 1994) son anteriores a la exigencia rutinaria de depositar factores de estructura, así que RCSB no puede calcular RSCC/RSR para confirmar cuantitativamente el ajuste de sus ligandos a la densidad electrónica; se confía en las coordenadas de los autores originales.
6. **No determinismo de Vina.** Incluso con semilla fija, Vina no es completamente reproducible bajo `--cpu 0` (multihilo); mitigado corriendo 3 semillas predeclaradas por caso×política, pero no eliminado — vale la pena declarar el rango observado, no solo la mediana.

## 18. Criterio de finalización

DockPrep Audit estará listo para un artículo cuando:

- El software tenga pruebas y documentación suficiente.
- La cohorte y las decisiones se hayan congelado antes de observar resultados.
- El benchmark sea reproducible desde fuentes públicas.
- Los resultados incluyan controles, fallos y sensibilidad, no solo casos exitosos.
- Las figuras se generen por código a partir de datos congelados.
- Las conclusiones no excedan lo medido.
- El repositorio y el archivo Zenodo correspondan exactamente a la versión descrita en el manuscrito.

## 19. Fase F — expansión hacia Q1 (decidida 2026-08-27)

**Contexto:** el borrador del Artículo 1 (`manuscript/draft.md`) se terminó de punta a punta (título → declaraciones) el 2026-08-27. El usuario compartió una revisión crítica externa (ChatGPT) del manuscrito completo. Evaluación propia de esa crítica, punto por punto, antes de actuar:

- **Válido y ya corregido en el manuscrito (2026-08-27):** el lenguaje de Discusión/Conclusión sobre-interpretaba lo que realmente se manipuló experimentalmente. La política de agua (conservador/simplificado) sí se aplicó a los 60 casos de los 4 estratos (Table 3/S3 tiene RMSD bajo ambas políticas para *todos* los estratos) — eso es un diseño válido: "una sola manipulación, aplicada uniformemente, ¿su efecto se concentra donde la auditoría lo predice?". Pero el texto sonaba a que también se habían probado políticas alternativas de altLoc/metal, cuando esas nunca variaron (altLoc siempre por mayor ocupancia, metal siempre conservado). Corregido en Discusión y Conclusión para ser precisos: la política de agua es la única manipulación probada; los otros 2 estratos son grupo de comparación bajo esa misma manipulación, no estratos con su propia política probada.
- **Válido, no corregido aún (requiere trabajo de software):** la regla `WATERS_PRESENT` del motor (cualquier HOH en toda la estructura, 59/60 casos) no tiene capacidad de triage real — el predictor que sí funciona es el criterio site-local de agua puente (≤3.0 Å ligando y receptor), que nunca se expuso como *finding* propio del motor de auditoría, solo como criterio de estratificación de la cohorte. Esto es una brecha real entre lo que el software dice detectar y lo que el paper valida.
- **Correcciones metodológicas técnicas válidas, evaluar caso por caso:** RMSD sin corrección de simetría (puede inflar RMSD en ligandos con grupos químicamente equivalentes — ácidos carboxílicos, anillos fenilo simétricos); diagnóstico más profundo de 1EPP/1PPM (~9 Å bajo ambas políticas); verificar que la fórmula $L=2.9R_g$ da espacio suficiente para ligandos pequeños/muy elongados.
- **Sugerencias de ampliación de alcance, no correcciones:** probar políticas alternativas de altLoc/metal; cohorte de validación externa independiente (discovery/validation split); evaluar el sistema como predictor formal (odds ratio, sensibilidad/especificidad, ROC/AUC, regresión logística); segunda representación de agua explícita (con H's, no solo carga TIP3P rígida) validando los casos positivos más fuertes (6ASH, 1WBK, 1CVZ, 4GID).

**Decisión del usuario 2026-08-27: Camino B — invertir el trabajo adicional necesario para apuntar a una revista Q1**, en vez de enviar la versión actual (Camino A, "corregir y enviar a JMM pronto") o solo el arreglo mínimo de redacción. Esto es un proyecto nuevo de varias semanas, no una revisión del manuscrito — se ejecuta como una fase más, con las mismas reglas de todo el proyecto (nada se congela sin decisión explícita, nada se descarta sin justificar, todo resultado —incluidos los negativos— se reporta).

**Plan de trabajo (orden tentativo, sujeto a revisión conforme avance):**

1. **F1 — Findings site-local reales en el motor de auditoría.** Agregar `SITE_BRIDGING_WATER_PRESENT` (agua ≤3.0 Å ligando *y* receptor — el mismo criterio ya usado para estratificar, ahora expuesto como salida real del motor v0.1.0), y evaluar si `SITE_ALTLOC_PRESENT`/`SITE_METAL_PRESENT` (versión site-local de los findings existentes) tienen sentido agregar también. Esto resuelve la brecha software-vs-paper directamente. Requiere: cambios a `dockprep_audit`, re-correr la auditoría sobre los 60 casos, regenerar `fig_finding_frequency.png`.
2. **F2 — Comparación real de política de altLoc.** Definir una segunda política defendible (ej. conformador alternativo específico vs. mayor ocupancia) y correr Fase C completa (3 semillas) con esa segunda política sobre los 60 casos, igual que ya se hizo para agua.
3. **F3 — Comparación real de política de metal.** Definir una segunda política defendible (ej. remover metal vs. conservarlo, o AD4Zn con coordinación explícita vs. carga formal simple) y correr la misma comparación.
4. **F4 — RMSD corregido por simetría.** Recalcular RMSD de las 360+ corridas ya existentes (y las nuevas de F2/F3) usando correspondencia por automorfismo molecular (ej. `rdkit.Chem.rdMolAlign.GetBestRMS`) en vez de matching por coordenada exacta únicamente. Mantener los valores actuales en material suplementario para trazabilidad, usar los corregidos como análisis principal.
5. **F5 — Diagnóstico secundario de 1EPP/1PPM.** Investigar por qué ambos fallan consistentemente (~9 Å) bajo cualquier política: ¿caja insuficiente, protonación, flexibilidad del ligando, mala referencia multiconformacional?
6. **F6 — Verificación de la fórmula de caja.** Confirmar que $L=2.9R_g$ reproduce fielmente el protocolo de Feinstein & Brylinski 2015 y no genera cajas insuficientes para ligandos pequeños o muy elongados en la cohorte actual.
7. **F7 — Cohorte de validación externa.** Construir una segunda cohorte independiente (60-100 PDB nuevos, no vistos durante el diseño de la hipótesis), enriquecida deliberadamente en agua puente site-local, y repetir el análisis sin ajustar ninguna definición después de ver resultados — el patrón discovery/validation que le da más peso a una conclusión Q1.
8. **F8 — Evaluación como predictor formal.** Con F1 (finding real) y F7 (cohorte de validación) listos, evaluar `SITE_BRIDGING_WATER_PRESENT` como predictor binario de sensibilidad a política (ej. $|\Delta \text{RMSD}| > 0.3$ Å, definido *antes* de correr F7): sensibilidad, especificidad, VPP/VPN, y si el tamaño de muestra lo permite, curva ROC/AUC con intervalo de confianza.
9. **F9 — Validación con segunda representación de agua (opcional, subconjunto).** Para los casos positivos más fuertes (6ASH, 1WBK, 1CVZ, 4GID), probar si el efecto se mantiene con una representación de agua explícita más rica (con hidrógenos, orientación de dipolo) en vez de solo la carga TIP3P rígida actual — sin necesitar FEP/MD completo.

**No se toca todavía:** el manuscrito actual (`manuscript/draft.md`) queda como está (con las correcciones de lenguaje del 2026-08-27 ya aplicadas) hasta que los resultados de F1-F9 estén listos para integrarse; no se reescribe la Discusión/Resultados de forma especulativa antes de tener los datos nuevos.

**F1 — completado (2026-08-28).** `SITE_BRIDGING_WATER_PRESENT`, `SITE_ALTLOC_PRESENT`, `SITE_METAL_PRESENT` agregados a `dockprep_audit.audit_pdb()` (parámetro opcional `ligand`, y `receptor_chains` para restringir qué cadenas cuentan como receptor). Versión del motor subida a 0.2.0 (`pyproject.toml`, `__init__.py`). 3 tests nuevos en `tests/test_audit.py`, todos pasan. Validado contra los 60 casos reales: para 1M17 y 1OHR, el nuevo finding reproduce exactamente el `bridging_waters_kept` ya usado en las corridas de Fase C.

**Hallazgo real al validar F1 (no un bug):** al correr `SITE_BRIDGING_WATER_PRESENT` (Eq. 1: agua ≤3.0 Å a ligando *y* receptor) contra las 60 fuentes, solo recuperó 12/15 casos del estrato `water_policy`, no 15/15. Causa: `scripts/verify_pilot_eligibility.py` (selección de cohorte) usó un criterio más laxo -- agua a ≤4.0 Å del ligando *solamente*, sin exigir contacto con el receptor -- que nunca fue el mismo que Eq. 1. Tres casos (1HRN, 1PPM, 3FNU) pasaron el filtro laxo pero no tienen ninguna agua que cumpla Eq. 1: confirmado contra `benchmark/phase-c-60/raw_runs/`, los tres tienen `bridging_waters_kept: 0` bajo la política "conservadora" -- para ellos, conservador y simplificado son el mismo input físico, así que cualquier diferencia de RMSD es solo el no-determinismo de Vina (Limitación 6), no un efecto de agua.

**Decisión del usuario 2026-08-28:** recalcular Tabla 3/Tabla 4 excluyendo esos 3 casos (n=12 para el estrato water_policy), en vez de solo documentarlo como limitación o posponerlo. Script nuevo: `scripts/recompute_water_stratum_excl_zero_bridge.py`. Resultado: el efecto se **fortalece** al excluir los casos triviales (brecha de éxito best-of-9 pasa de 0.088 a 0.111; brecha de RMSD de 0.160 Å a 0.196 Å), aunque el p-valor del Wilcoxon empeora levemente (0.28→0.38, n=15→12) por pérdida de potencia estadística, no por debilitamiento del efecto. `manuscript/draft.md` actualizado en: Abstract, Tabla 1 (nuevos finding codes), Sección 2.1 (versión motor), Sección 2.3 (párrafo nuevo explicando la discrepancia de criterios), Tabla 3 (fila water_policy a n=12 + fila de referencia n=15), Tabla 4 (título/nota), párrafo 5/2/8 → 5/2/5, test de Wilcoxon, párrafo 1EPP/1PPM en 3.4 (ahora con explicación causal), y Limitación 4.

**Verificado 2026-08-28:** se corrió la misma comprobación con `SITE_ALTLOC_PRESENT`/`SITE_METAL_PRESENT` contra los estratos `alternate_location` y `metal_or_cofactor` -- 15/15 confirmados en ambos, sin discrepancia. La causa raíz era específica del criterio de agua (selección con 4.0 Å solo-ligando vs. Eq. 1 con 3.0 Å simétrico ligando+receptor); altLoc y metal usaron el mismo umbral de 6 Å tanto en la selección de cohorte como en el finding site-local nuevo, así que no hay un problema equivalente ahí. F1 se da por cerrado.

**Segundo bug encontrado durante F2 (2026-08-28):** al construir la lista de casos para F2 (comparación real de política altLoc), `SITE_ALTLOC_PRESENT` recuperaba 17/60 en vez de 15/60. Causa: `altlocs = [a for a in atoms if a["altloc"]]` en `audit.py` contaba altLoc de **cualquier** átomo (ATOM+HETATM), incluyendo aguas desordenadas y el propio ligando -- no solo conformadores reales del receptor. Para 1GS4 el "altLoc site-local" era `HOH A2026` (agua); para 1SN5 era `T3 C601` (el propio ligando). **Fix:** excluir aguas de `altlocs` (afecta tanto `ALTLOC_PRESENT` como `SITE_ALTLOC_PRESENT`) y excluir además los átomos del ligando declarado específicamente de `SITE_ALTLOC_PRESENT`. Tras el fix, `SITE_ALTLOC_PRESENT` recupera exactamente 15/15 del estrato `alternate_location`, igual de limpio que `SITE_METAL_PRESENT`. 3 tests siguen pasando. `ALTLOC_PRESENT` whole-structure no cambió (sigue 22/60): ningún caso dependía únicamente de una agua/ligando para ese finding.

**Nota abierta, no resuelta:** al verificar que 1GS4/1SN5 (excluidos correctamente del grupo primario de F2) fueran comparaciones triviales como en F1, se confirmó para 1GS4 (receptor byte-idéntico entre políticas) pero **no** para 1SN5: su archivo tiene una segunda copia del ligando T3 (cadena D, resseq 602) que `extract_receptor_atoms()` retiene como parte del receptor (promovida de HETATM a ATOM) independientemente de la instancia de ligando declarada (C601) y del radio site-local de 6 Å -- y el altLoc de esa segunda copia sí se resuelve distinto entre políticas, produciendo una diferencia real de RMSD (~2.4 Å en la semilla 42). Por qué esta segunda copia se retiene como "receptor" en absoluto no se investigó todavía -- posible candidato para F5/diagnóstico adicional, no bloquea F2 ya que 1SN5 de cualquier forma queda fuera del grupo primario n=15.

**F2 — completado (2026-08-28).** 102 corridas nuevas (15 casos SITE_ALTLOC_PRESENT x 2 políticas x 3 semillas, agua fija en "conservador"), 0 fallos de preparación. Resultado (`benchmark/phase-f2-altloc/phase_f2_summary.csv`, `phase_f2_by_case.csv`):

| Política | Succ(top1) | RMSD(top1) | Succ(best-9) | RMSD(best-9) | Unstable |
|---|---|---|---|---|---|
| highest_occupancy | 0.556 | 2.662 | 0.800 | 1.249 | 3/15 |
| lowest_occupancy | 0.600 | 2.649 | 0.867 | 1.427 | 0/15 |

Wilcoxon pareado (best-of-9, n=15): **W=27.0, p=0.0637**, diff. media pareada (high-low) = -0.178 Å -- más cerca de significancia convencional que el hallazgo de agua (p=0.28-0.38), pese a un n menor. Concentración: 4/15 ayudan (mayor ocupancia mejor) / 1/15 perjudica / 10/15 sin efecto -- mismo patrón cualitativo que agua (minoría con efecto real). Los dos efectos más grandes van en direcciones opuestas: 2I4H (+2.275 Å, favorece mayor ocupancia) vs 5E0J (-1.706 Å, favorece menor ocupancia), un tira y afloja real, no un sesgo direccional limpio.

**Implicación para el manuscrito:** esto responde directamente a la crítica original ("solo se varió la política de agua") -- ahora hay DOS manipulaciones de política reales y probadas. El lenguaje actual de Discusión/Limitaciones que dice "esto requeriría una comparación dedicada que este estudio no corrió" (citando F2/F3) queda desactualizado para F2 y debe reescribirse. Pendiente: decidir cómo integrar (¿estrato agua+altLoc combinados en Tabla 3, o tabla nueva paralela?; ¿mencionar en Abstract?).

**F2 integrado al manuscrito (2026-08-28):** nueva §3.6 "Alternate-location-policy effect" (Tabla 5 + Tabla 6), paralela a §3.5. Abstract, Introducción, Discusión, Conclusiones y §2.3 actualizados para reflejar dos manipulaciones probadas. 2 referencias desactualizadas encontradas y corregidas en la revisión de consistencia posterior (Limitación 7 y caption Fig. S1, que seguían diciendo que altLoc "no recibió manipulación dedicada").

**F3 — completado (2026-08-28).** 90 corridas nuevas (15 casos SITE_METAL_PRESENT x 2 políticas x 3 semillas, agua fija en "conservador", altLoc fija en "highest_occupancy"), 0 fallos de preparación. Segunda política: `remove` -- elimina el ion metálico por completo; los tautómeros de histidina coordinante ya no se fuerzan geométricamente y vuelven al default de Meeko (forzar un tautómero por un metal ausente no tiene justificación física).

Resultado (`benchmark/phase-f3-metal/phase_f3_summary.csv`, `phase_f3_by_case.csv`):

| Política | Succ(top1) | RMSD(top1) | Succ(best-9) | RMSD(best-9) | Unstable |
|---|---|---|---|---|---|
| retain | 0.533 | 2.995 | 0.689 | 2.151 | 2/15 |
| remove | 0.533 | 3.038 | 0.800 | 1.953 | 2/15 |

Wilcoxon pareado (best-of-9, n=15): W=41.5, **p=0.293** (el más lejano de significancia de las tres comparaciones). Concentración: 1/15 ayuda / 2/15 perjudica / 12/15 sin efecto -- el más débil y menos concentrado de los tres, pero contiene el efecto individual más grande de todo el estudio: **1DTH** (metaloproteasa de zinc), retain=4.205 Å vs remove=1.475 Å (Δ=-2.730 Å).

**Investigación de 1DTH (antes de integrar al manuscrito, a pedido del usuario):** el top-1 es casi idéntico entre políticas (~7.8-7.9 Å, siempre falla) -- el efecto viene enteramente de que la mejor de las 9 poses cruza el umbral de 2.0 Å solo bajo `remove` (1.47-1.82 Å) vs `retain` (2.26-4.31 Å). Es un efecto de sampling, no de scoring: remover el metal le permite a Vina muestrear una pose casi nativa que bajo `retain` nunca aparece entre las 9. **Confound real encontrado:** remover el metal también revierte el tautómero de las 3 histidinas coordinantes (A:142/146/152, forzadas a HID bajo `retain` por coordinar Zn vía NE2) al default de Meeko -- así que la manipulación "metal" no aísla el átomo solo, cambia metal + protonación de histidinas juntos. Verificado que esto NO es un confound trivial/universal: 1CBX, 5A2S, 4EXS y 1KJO también tienen histidinas coordinantes que revierten tautómero bajo `remove` (2-4 cada uno) pero no muestran efecto comparable -- así que el confound existe ampliamente pero solo se traduce en un efecto grande en 1DTH específicamente.

**Integrado al manuscrito:** nueva §3.7 "Metal-policy effect" (Tabla 7 + Tabla 8) con la nota de confound explícita, nueva Limitación 8 dedicada al confound, y §2.3/Discusión/Conclusiones/Abstract/Introducción actualizados para las TRES manipulaciones. El three-way (altLoc, metal, agua) que la crítica original pedía queda cerrado.

**Pendiente real, no resuelto:** aislar el efecto del átomo de metal del cambio de tautómero de histidina requeriría una tercera política (metal removido, tautómero aún forzado geométricamente) -- no corrida. Candidato para trabajo futuro si se quiere un artículo de seguimiento, no bloquea la publicación actual (el confound está documentado con transparencia, no oculto).

**F4 — completado (2026-08-28).** RMSD corregido por simetría (`rdkit.Chem.rdmolops.GetSubstructMatches(mol, mol, useChirality=True)` sobre el grafo pesado del ligando, mínimo sobre automorfismos, SIN realineación de cuerpo rígido) recomputado sobre las 540 corridas ya existentes de agua/F2/F3 -- sin docking nuevo, solo post-procesamiento de archivos ya en disco (`scripts/compute_symmetry_rmsd.py`, `scripts/compute_symmetry_rmsd_aggregate.py`, `scripts/build_symmetry_corrected_tables.py`).

**Hallazgo real, no trivial:** la corrección **invierte cuál de las tres comparaciones es la más significativa**:
- **Agua:** p mejora levemente (0.28-0.38 → 0.21, n=11 tras excluir también 1OHR), sin cambio cualitativo.
- **AltLoc:** p **empeora** (0.064 → 0.21) -- el caso 1T46 (Δ+0.570 Å) resultó ser casi enteramente un artefacto de etiquetado simétrico en su ligando; al corregirlo, el efecto ya no se acerca a significancia.
- **Metal:** p **mejora dramáticamente** (0.293 → 0.0084, ¡ahora significativo!) -- el caso 4G9L (Δ+0.320 Å, "retener ayuda") caía por debajo del umbral tras la corrección (+0.072 Å), eliminando el único caso que iba en contra de la tendencia; sin él, el patrón queda 100% consistente (1DTH y 1KJO, ambos "retener perjudica").

**Verificado, no es un bug mío:** confirmado caso por caso (1DTH/1KJO prácticamente sin cambio = efectos reales; 4G9L y 1T46 colapsan = artefactos de simetría). 1 caso (1OHR) excluido de la corrección por completo: el SMILES de protonación elegido resuelve a 40 átomos pesados vs. 44 en la estructura depositada -- discrepancia de datos de RCSB para el componente 1UN, no un bug del pipeline.

**Decisión del usuario 2026-08-28:** promover los valores corregidos a análisis principal (Tablas 3-8 del manuscrito), mantener los valores originales (exact-index) en Table S3 para trazabilidad, per el plan original de F4. Manuscrito actualizado en: §2.4 (nueva metodología + exclusión de 1OHR documentada), Tablas 3-8 completas reescritas, Abstract, Introducción, Discusión (3 párrafos), Limitaciones (ítems 2/4/7 actualizados + nuevo ítem 9 sobre alcance del método), Conclusiones, captions de Fig. 3/4/S1 (marcadas como no corregidas, por decisión de no reconstruir esas figuras en esta sesión). Roadmap: 3 tests de `dockprep_audit` sin afectar (F4 no toca ese motor).

**Pendiente, no bloqueante:** Fig. 3, Fig. 4 y Fig. S1 siguen construidas sobre RMSD exact-index (sin corregir) -- decisión explícita de no reconstruirlas esta sesión dado que la corrección de agua específicamente es modesta (ningún caso >0.5 Å) y no se espera que cambie el patrón visual; queda como candidato de trabajo futuro si se quiere consistencia visual completa con las Tablas 3-8.

**F5 — completado (2026-08-28).** Diagnóstico de por qué 1EPP/1PPM fallan consistentemente (~9 Å) bajo cualquier política (`scripts/diagnose_1epp_1ppm.py`), revisando las 4 hipótesis del roadmap:
- **Caja insuficiente:** descartada como explicación general -- 1PPM cabe cómodo (12.9 Å en caja de 17.5 Å); 1EPP cabe justo (16.5 en 17.2); pero 1QRP tiene un ligando que literalmente EXCEDE su caja (20.9 en 18.9 Å) y aun así dockea bien (1.8 Å) -- una caja ajustada no predice fallo en este cohorte.
- **Protonación:** descartada -- ambos casos tuvieron éxito en el primer intento de variante SMILES, sin fallback.
- **Flexibilidad del ligando:** factor real pero no suficiente. 1EPP es el ligando MÁS flexible de los 60 casos (23 enlaces rotables, el máximo del cohorte); 1PPM está en el tercio superior (15). Hallazgo más amplio, no buscado deliberadamente: **el estrato `water_policy` completo es ~2x más flexible en promedio que los otros tres** (media 14.0 enlaces rotables, mediana 15, vs. 5.1/7.9/6.2 para altLoc/metal/low-risk) -- probablemente porque el agua puente es una característica bien documentada de sitios activos de aspartato-proteasas, y esa familia biológica trae ligandos peptidomiméticos grandes y flexibles. Pero 3 casos igual de flexibles en el mismo estrato (1QRP, 1WBK, 1WBM, 15-21 enlaces rotables) dockean bien -- la flexibilidad eleva el riesgo sin ser determinística.
- **Referencia multiconformacional de mala calidad:** no se puede descartar -- ambas estructuras (1994, 1992) están entre los 11 casos sin archivo de factores de estructura disponible (sin RSCC/RSR verificable), consistente con la Sección 3.2 ya documentada.

**Integrado al manuscrito:** nuevo párrafo extenso en §3.4 con las 4 hipótesis evaluadas; nueva **Limitación 10** documentando el desbalance de flexibilidad por estrato (aclarando que NO invalida la comparación pareada de agua, ya que compara el mismo ligando consigo mismo bajo ambas políticas, pero sí afecta comparar el RMSD absoluto del estrato `water_policy` contra los otros tres estratos cara a cara); Discusión actualizada con referencia cruzada.

**F6 — completado (2026-08-28).** Verificación de la fórmula de caja (`scripts/verify_box_formula.py`).

- **Fidelidad a la fuente:** confirmado vía fetch en vivo al paper (PMC4468813) -- la fórmula $L=2.9\times R_g$ reproduce textualmente "the dimensions of the search space are 2.9 times larger than the radius of gyration", validado en ligandos de 6-100 átomos pesados con caja cúbica. El paper NO menciona ninguna advertencia explícita para ligandos elongados o muy pequeños.
- **Chequeo empírico en los 60 casos:** **13/60 (22%)** tienen la extensión del ligando (bounding box eje-alineado) excediendo la caja generada, hasta 2.0 Å en el peor caso (1QRP). Concentrado en los ligandos más elongados (9 de los 10 con mayor ratio de elongación tienen margen cero o negativo), NO en los más pequeños (todos los casos con <15 átomos pesados tienen margen positivo cómodo, 0.6-1.6 Å) -- confirma que la preocupación a priori era correcta: es un problema de forma (Rg es un escalar isotrópico), no de tamaño.
- **Costo real, no catastrófico:** los 13 casos con caja insuficiente dockean modestamente peor (éxito 0.692 vs 0.804, RMSD medio 2.57 vs 1.69 Å) pero NO son la causa dominante de los peores fallos del estudio -- solo 2/13 (1SN5, 1EED) llegan a los ~9 Å de las fallas más severas, y ambos ya tienen explicación más específica documentada (F2's second-ligand-copy issue para 1SN5; F5's flexibility finding para 1EED).

**Integrado al manuscrito:** nuevo párrafo en §2.3 (junto a Eq. 2) citando la fuente textualmente y notando el supuesto de isotropía; nuevo párrafo en §3.3 con el chequeo empírico completo; nueva **Limitación 11** proponiendo una caja consciente de la forma (ejes principales en vez de Rg escalar) como refinamiento futuro, no implementado porque cambiaría retroactivamente las cajas de los 60 casos ya congelados.

**F7 — completado (2026-08-28/29).** Cohorte de validación externa, protocolo discovery/validation con todas las definiciones fijadas ANTES de screening (`scripts/f7_screen_candidates.py`, `f7_freeze_cohort.py`, `f7_run_redocking.py`, `f7_substitute.py`, `f7_aggregate.py`).

- **Construcción:** búsqueda RCSB (710 candidatos: X-ray, ≤2.2 Å, cadena única, con datos de afinidad de unión caracterizados), excluyendo los 60 casos + registro de 132 PDBs de estudios previos (`previous-study-exclusions.csv`). Screening local con `SITE_BRIDGING_WATER_PRESENT` (mismo finding real del motor): 453/710 (64%) calificaron. Selección determinística por orden de PDB ID (sin ver ningún resultado de docking) con tope de 2 casos por component_id de ligando → 80 casos congelados con checksums.
- **Redocking:** 480 corridas primera pasada (agua conservador/simplificado × 3 semillas), 384 OK (80% limpio, consistente con la cohorte principal). 16 fallos técnicos (mismas categorías que Tabla S2: Meeko/dimorphite-dl) → **una sola ronda de sustitución** (mismo protocolo ya establecido), 16 reemplazos tomados del pool ya calificado en el mismo orden fijo. 4 de los 16 sustitutos también fallaron y NO se sustituyeron de nuevo (regla: una sola ronda) → **n=76 final** (dentro del rango 60-100 planeado).
- **Resultado (RMSD corregido por simetría, mismo método que el resto del estudio): el hallazgo de agua NO replica.** Wilcoxon pareado n=76: **W=1371.0, p=0.78** (vs. p=0.21 en n=11 de la cohorte principal). Sin dirección consistente (simplificado tiene incluso mejor éxito: 0.803 vs 0.789). Concentración: 10 ayuda / 11 perjudica / 55 sin efecto -- casi perfectamente balanceado (la cohorte principal: 5 ayuda / 1 perjudica / 5 sin efecto en n=11).
- **Verificado que no es un bug:** tasas de éxito 79-80%, más altas que la mayoría de estratos de la cohorte principal -- el pipeline funciona correctamente, simplemente no hay señal de agua en esta muestra.
- **Hipótesis no confirmada:** F5 ya había encontrado que el estrato `water_policy` original está dominado por peptidomiméticos de aspartato-proteasas excepcionalmente flexibles (14 enlaces rotables en promedio); la cohorte F7 no tiene ese sesgo químico. Plausible pero no probado.

**Decisión del usuario 2026-08-29:** reportar con honestidad total, suavizar el framing central del paper en vez de esconder o reencuadrar el hallazgo como "más matizado". Integrado en: nueva **§3.9 "External validation cohort"** (Tabla 9 + análisis completo), Abstract, Introducción (contribución ahora "cuádruple", no triple), Discusión (párrafo 1 reescrito por completo -- el "hallazgo central" ahora es la no-replicación en sí misma, no el efecto de agua), Conclusiones, y nueva **Limitación 12** (diferencias de curaduría entre las dos cohortes como explicación alternativa no descartada). AltLoc y metal NO fueron validados externamente -- el texto ahora dice explícitamente que su significancia dentro-de-cohorte no debe leerse como más confiable solo porque agua falló la validación.

**F8 — completado (2026-08-29).** Evaluación predictiva formal de `SITE_BRIDGING_WATER_PRESENT` (`scripts/f8_predictive_evaluation.py`), usando F1 (finding real) + F7 (cohorte de validación).

- **Diseño:** dado que F7 es 100% positivo por construcción (se filtró específicamente por este finding), no puede dar una matriz de confusión completa por sí sola. En cambio, usé los **60 casos de la cohorte principal** (todos pasaron por la misma manipulación agua conservador/simplificado, Tabla 3), recalculando `SITE_BRIDGING_WATER_PRESENT` para los 60 (no solo el estrato water_policy): 36/60 positivos. Umbral de sensibilidad (|ΔRMSD|>0.3 Å) ya estaba pre-especificado desde Tabla 4/6/8 -- no se ajustó para esta evaluación.
- **Matriz de confusión (n=59, excluye 1OHR):** TP=9, FP=26, FN=1, TN=23. **Sensibilidad=0.900, VPN=0.958** (excelentes) pero **Especificidad=0.469, VPP=0.257** (débiles) -- el flag casi nunca deja pasar un caso realmente sensible, pero dispara muchas falsas alarmas.
- **Hallazgo clave -- la magnitud SÍ predice, aunque la dirección no:** diferencia de |ΔRMSD| media entre flagged/unflagged: 0.329 Å vs 0.023 Å (14x), Mann-Whitney U=801.5, **p<0.0001** -- altamente significativo. Y esto **SÍ replica en F7**: media/mediana de F7 (100% flagged) = 0.372/0.097 Å, muy cercano a la media/mediana del grupo flagged de la cohorte principal (0.329/0.110 Å), NO al grupo unflagged (0.023/0.000 Å). El VPP también replica de cerca (F7: 0.276 vs principal: 0.257).
- **Conclusión matizada:** el flag predice QUE habrá sensibilidad de mayor magnitud si la hay, replicado en datos independientes -- pero NO predice EN QUÉ DIRECCIÓN (eso es justo lo que F7 ya había refutado). Esto reconcilia F7 y F8: no es que el flag no sirva para nada, es que sirve para una cosa más estrecha de lo que el diseño original del estudio afirmaba.

**Integrado al manuscrito:** nueva **§3.10 "Predictive evaluation"** (Tabla 10 + análisis completo); Abstract, Introducción, Discusión (párrafo 1 y 3) y Conclusiones actualizados para reflejar la distinción magnitud-vs-dirección en vez de un "no replica" plano.

Con F1-F8 completos, queda solo **F9** (segunda representación de agua explícita, opcional) de la lista original de Fase F.

**F9 — completado (2026-09-01).** Robustez del efecto de agua frente a una segunda representación explícita, más rica que el oxígeno rígido sin H (`scripts/smoke_redock_case.py` -- nueva política `conservative_oriented`; `scripts/f9_oriented_water.py`, `scripts/f9_aggregate.py`).

- **Decisión del usuario 2026-09-01:** proceder con orientación de H's hacia el contacto de puente-H más cercano (geometría TIP3P idealizada) y correr solo en los 4 casos con el efecto más fuerte de la cohorte principal (6ASH, 1WBK, 1CVZ, 4GID; Tabla 4), no en toda la cohorte.
- **Método:** cada agua puente retenida se modela ahora como un cuerpo rígido O+H+H (enlace O-H=0.9572 Å, ángulo H-O-H=104.52°, cargas TIP3P completas O=-0.834e/H=+0.417e cada uno), con los dos H's orientados hacia los mismos dos átomos (ligando + receptor) que ya calificaron a esa agua como "puente" en `find_bridging_waters` -- sin introducir un criterio de orientación nuevo y separado. Se agregó `nearest_ligand_atom`/`nearest_protein_atom` a los metadatos de cada agua puente para soportar esto.
- **Bugs encontrados y corregidos durante la implementación (antes de correr nada costoso):** (1) desalineación de columnas en el formato PDBQT de los átomos H (campo de nombre de átomo de 3 vs 4 caracteres) que Vina rechazaba con "Coordinate is not valid" -- detectado en un smoke test de un solo caso antes del run completo; (2) una condición de carrera transitoria de visibilidad de archivo (`receptor_raw.pdb` no visible de inmediato a un subproceso hijo en este checkout sincronizado con OneDrive) -- mitigada con una pequeña espera (`_wait_for_file`), sin efecto en los resultados científicos de ningún run anterior.
- **12/12 corridas de Vina completadas sin error** (4 casos × 3 seeds, exhaustiveness=32, igual que la cohorte principal).
- **Resultado (RMSD mejor-de-9 corregida por simetría, mediana de 3 seeds):**

  | Caso | Δ oxígeno rígido (bare-O) | Δ agua orientada (F9) |
  |---|---|---|
  | 6ASH | +1.417 Å | **+2.581 Å** |
  | 1WBK | +0.497 Å | **+0.548 Å** |
  | 1CVZ | +0.480 Å | **+0.484 Å** |
  | 4GID | +0.485 Å | +0.251 Å (cae bajo el umbral de 0.3 Å, pero mantiene dirección) |

  Los 4/4 casos mantienen la dirección "agua ayuda"; 3/4 muestran un efecto igual o mayor con la representación más rica; 1/4 (4GID) muestra un efecto menor que ya no cruza el umbral de 0.3 Å pero sigue siendo positivo.
- **Conclusión:** en este subconjunto específico de los 4 casos más fuertes, el modelo de oxígeno rígido sin H fue, si acaso, conservador respecto al efecto -- no una fuente de inflación artificial. No es una prueba estadística formal (n=4, chequeo dirigido, no una re-corrida de toda la cohorte).

**Integrado al manuscrito:** nuevo párrafo en **§3.5** (con tabla comparativa inline) justo antes de §3.6, cita de datos a `benchmark/symmetry-rmsd/f9_oriented_water_symmetry_rmsd.csv`; **Limitación 1** actualizada para referenciar este chequeo sin eliminar la limitación general (no se re-corrió toda la cohorte, y el chequeo no aborda el término de desolvatación faltante).

**Con esto, Fase F completa: F1-F9, los 9 items originales de PROJECT-ROADMAP.md §19.**

**Post-Fase F, mantenimiento (2026-09-01):**

- **Proyecto reubicado** de `...OneDrive\Documents\ChatGPT\Zenodo y cuentas\dockprep-audit` a `C:\Users\Andre\Proyectos doctorado\dockprep-audit`, para unificar con los demás proyectos de investigación. El `.venv` no se movió (paquetes de PyMOL sin usar por el pipeline quedaron marcados "solo en la nube" por OneDrive, imposibles de copiar localmente) -- se recreó limpio con las mismas versiones exactas (`rdkit==2025.9.6`, `meeko==0.7.1`, `dimorphite_dl==2.0.2`, `scipy==1.18.0`, más `gemmi==0.7.5`, dependencia real de meeko no declarada en su metadata de PyPI). Verificado con un smoke test de redocking completo (6ASH) reproduciendo exactamente el mismo resultado numérico que antes de la mudanza. No había historial de git dentro de `dockprep-audit` (pertenecía a la carpeta padre, 0 commits) -- nada que migrar.
- **`README.md` actualizado** de "v0.1" a v0.2.0: documenta los tres findings site-local (`SITE_BRIDGING_WATER_PRESENT`, `SITE_ALTLOC_PRESENT`, `SITE_METAL_PRESENT`) y referencia el manuscrito/benchmark companion.
- **Nuevos tests unitarios** (`tests/test_smoke_redock_case.py`, 4 tests, ninguno existía antes para `smoke_redock_case.py`): geometría de `_oriented_hydrogen_positions` (longitud de enlace, ángulo H-O-H, caso degenerado sin colapsar) y un test de regresión específico para el bug de alineación de columnas PDBQT encontrado durante F9 -- verificado deliberadamente reintroduciendo el bug (falla) y restaurándolo (pasa) antes de dejarlo en el repo. 7/7 tests pasan en total (`python -m unittest discover -s tests`).

