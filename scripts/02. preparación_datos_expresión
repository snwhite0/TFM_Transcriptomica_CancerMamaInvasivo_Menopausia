# 2.1. Cargar matriz de recuentos
expr_matriz_tumor = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
cat("Dimensión de la matriz de datos de expresión INICIAL: ",dim(expr_matriz_tumor),"\n")
print(expr_matriz_tumor[1:5,1:5])

# 2.2. Extraer matriz de conteos y redondear valores
conteo_tumor = as.matrix(expr_matriz_tumor[, -c("Hugo_Symbol", "Entrez_Gene_Id")]) # Extraer matriz de conteos

# Normalización de nombres de genes: Unificando Entrez_ID con HUGO_Symbols
entrez_ids = expr_matriz_tumor$Entrez_Gene_Id
entrez_ids = as.character(expr_matriz_tumor$Entrez_Gene_Id)
hugo_oficial = mapIds(org.Hs.eg.db,
                      keys = entrez_ids,
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")

hugo_oficial[is.na(hugo_oficial)] = entrez_ids[is.na(hugo_oficial)] 
nombres_filas = hugo_oficial 
duplicados = duplicated(nombres_filas) | duplicated(nombres_filas, fromLast = TRUE) 
nombres_filas[duplicados] = paste0(hugo_oficial[duplicados], "_", entrez_ids[duplicados]) # Para duplicados, añadir Entrez ID para hacer nombres únicos
rownames(conteo_tumor) = nombres_filas # Asignar como rownames

cat("Dimensión de la matriz de conteos: ",dim(conteo_tumor),"\n")
cat("Matriz de conteos sin redondear:"); print(conteo_tumor[1:5,1:3])

# Redondear valor de conteos
conteo_tumor_redondeado = round(conteo_tumor)
cat("Matriz de conteos redondeada:"); print(conteo_tumor_redondeado[1:5,1:3])

cat("¡Listo! Preprocesamiento de los datos de expresión génica completado.\n")
