### Preparación de datos de expresión génica
# Carga de la matriz de recuentos
expr_matriz_tumor = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
cat("Dimensión de la matriz de datos de expresión INICIAL: ",dim(expr_matriz_tumor),"\n")
print(expr_matriz_tumor[1:5,1:5])

# Extraer matriz de conteos
conteo_tumor = as.matrix(expr_matriz_tumor[, -c("Hugo_Symbol", "Entrez_Gene_Id")]) # Extraer matriz de conteos

# Normalización de nombres de genes: unificar Entrez_ID con HUGO_Symbols
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

# Redondear valores de conteos
conteo_tumor_redondeado = round(conteo_tumor)
