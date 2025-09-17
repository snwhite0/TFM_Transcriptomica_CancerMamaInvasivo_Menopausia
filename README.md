# TFM_Transcriptomica_CancerMamaInvasivo_Menopausia
# Análisis Transcriptómico de Cáncer de Mama Invasivo según Estado Menopáusico

## Descripción
Este proyecto implementa un **flujo de trabajo de expresión génica diferencial (DEA)** utilizando **edgeR** en muestras de cáncer de mama invasivo. El objetivo es identificar **genes diferencialmente expresados (DEGs)** considerando el efecto del **estado menopáusico** y los **subtipos moleculares** (Basal, HER2, LumA, LumB, Normal).

### Objetivos
1. **Estado menopáusico (cohorte completa):**  
   Identificar genes cuya expresión varía según el estado menopáusico (pre-, peri-, post-), ajustando por el subtipo tumoral: ~ subtipo_tumoral + estado_menopausico
2. **Estado menopáusico dentro de cada subtipo:**  
Analizar la expresión diferencial por estado menopáusico en cada subtipo molecular, usando modelos independientes: ~ estado_menopausico

## Flujo del Script
1. **Instalación y carga de librerías**  
2. **Carga y preprocesamiento de datos**  
- Metadatos clínicos de pacientes y sus respectivas muestras: limpieza y organización  
- Datos de expresión génica: conteos de RNA-seq  
3. **Análisis exploratorio de datos**  
- PCA para visualización de la variabilidad global  
4. **Análisis de expresión diferencial (DEA)**  
- Preparación y filtrado de datos  
- DEA para los objetivos 1 y 2  
5. **Visualización de resultados**  
- Volcano plots, Upset plots, Heatmaps  
- Boxplots para DEGs significativos (test no-paramétrico de Kruskal-Wallis) 
6. **Análisis funcional de DEGs**  
- Enriquecimiento por sobre-representación (ORA) con EnrichR  
- Bubble plots para rutas (KEEG) y ontologías (GO) significativas
7. **Análisis de supervivencia**  
- OS, DSS, DFS, PFS  
- Curvas de Kaplan-Meier y test de Log-rank

## Requisitos
- R >= 4.2  
- Paquetes: `edgeR`, `enrichR`, `ggplot2`, `dplyr`, `stringr`, `survival`, `survminer`, `readr` ...  
- Otros paquetes de tidyverse según scripts
  
## Uso
1. Configurar rutas de datos y metadatos en los scripts.  
2. Ejecutar los scripts en orden:  
   - DEA objetivo 1 → cohorte completa  
   - DEA objetivo 2 → por subtipo  
   - Análisis funcional (EnrichR)  
   - Análisis de supervivencia  
3. Los resultados se guardan en subcarpetas dentro de `resultados/`.

## Estructura del repositorio
├── TFM_script_definitivo.html # Archivo HTML con el código completo

├── scripts/ # Scripts de R por tipo de análisis

├── data/ # Datos de entrada (raw data)

├── README.md # Este archivo
