### Preparación de librerías y paquetes

## Configuración inicial
# Crear función para instalar paquetes CRAN si faltan
instalar_si_falta = function(libreria) {
  if (!requireNamespace(libreria, quietly = TRUE)) {
    install.packages(libreria)
  }
}

# Crear función para instalar paquetes Bioconductor si faltan
instalar_bioconductor_si_falta = function(libreria) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace(libreria, quietly = TRUE)) {
    BiocManager::install(libreria, update = FALSE, ask = FALSE)
  }
}


## Ejecutar instalación y carga de librerías y paquetes necesarios
paquetes_cran = c("data.table", "tidyverse", "viridis", "ggrepel","UpSetR", "pheatmap", "RColorBrewer", "enrichR", "survival")   # Paquetes CRAN
paquetes_bioc = c("edgeR", "limma", "topGO", "AnnotationDbi", "org.Hs.eg.db", "GO.db", "Rgraphviz", "dorothea", "viper")         # Paquetes Bioconductor

# Instalar paquetes si faltan
invisible(sapply(paquetes_cran, instalar_si_falta))
invisible(sapply(paquetes_bioc, instalar_bioconductor_si_falta))
cat("¡Listo! Las librerías necesarias han sido instaladas y cargadas correctamente.\n")

# Cargar paquetes silenciosamente
suppressPackageStartupMessages({
  suppressWarnings({
    sapply(c(paquetes_cran, paquetes_bioc), require, character.only = TRUE)
    library(grid)
    library(gridExtra)
  })
})


## Registro de versiones de paquetes y librerías instalados
# Función para obtener versiones y origen
info_paquetes = function(lista, origen) {
  data.frame(
    Paquete = lista,
    Version = sapply(lista, function(p) as.character(packageVersion(p))),
    Origen = origen,
    stringsAsFactors = FALSE
  )
}

# Generar tabla con versiones y guardar en archivo CSV
df_cran = info_paquetes(paquetes_cran, "CRAN")
df_bioc = info_paquetes(paquetes_bioc, "Bioconductor")
df_otros = info_paquetes(c("grid", "gridExtra"), "Base/Extra")
df_total = rbind(df_cran, df_bioc, df_otros) # Combinar todo

write.csv(df_total, file = "versiones_paquetes_cargados.csv", row.names = FALSE)
cat("Se ha generado el archivo 'versiones_paquetes_cargados.csv'.\n")
    
