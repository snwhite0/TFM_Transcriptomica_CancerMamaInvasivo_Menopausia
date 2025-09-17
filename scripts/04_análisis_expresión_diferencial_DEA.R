### Análisis de Expresión Diferencial (DEA)

## PREPARACIÓN PREVIA 
# Garantizar el filtrado y orden correcto de los datos de interés
muestras_comunes = intersect(colnames(dge_tumor_filtrado_norm), pca_df_tumor$sample_id) 
dge_tumor_filtrado_norm_2 = dge_tumor_filtrado_norm[, muestras_comunes]
pca_df_tumor_filtrado = pca_df_tumor[pca_df_tumor$sample_id %in% muestras_comunes, ]
stopifnot(all(colnames(dge_tumor_filtrado_norm_2) == pca_df_tumor_filtrado$sample_id))

cat("Dimensión de la matriz de conteos preprocesada: ", dim(dge_tumor_filtrado_norm), "\n")
cat("Dimensión de la matriz de conteos preprocesada filtrada: ", dim(dge_tumor_filtrado_norm_2), "\n")
cat("Dimensión de la tabla de metadatos: ", dim(pca_df_tumor), "\n")
cat("Dimensión de la tabla de metadatos filtrados: ", dim(pca_df_tumor_filtrado), "\n")

# Renombrar variables para mayor claridad y guardar en archivos CSV
dge = dge_tumor_filtrado_norm_2
metadatos = pca_df_tumor_filtrado

output_dir = "../datos/procesados/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(dge, file.path(output_dir, "dge"), row.names = FALSE)
cat(sprintf("Archivo 'dge' guardado en: %s\n", output_dir))
write.csv(metadatos, file.path(output_dir, "metadatos"), row.names = FALSE)
cat(sprintf("Archivo 'metadatos' guardado en: %s\n", output_dir))

# Limpieza de factores y definición de los niveles posibles
metadatos$MENOPAUSE_STATUS = as.factor(trimws(metadatos$MENOPAUSE_STATUS))
metadatos$SUBTYPE = as.factor(trimws(metadatos$SUBTYPE))
metadatos$MENOPAUSE_STATUS = factor(trimws(metadatos$MENOPAUSE_STATUS),
                                    levels = c("Premenopausia", "Perimenopausia", "Postmenopausia"))
metadatos$SUBTYPE = factor(trimws(metadatos$SUBTYPE),
                           levels = c("BRCA_Basal", "BRCA_Her2", "BRCA_LumA", "BRCA_LumB", "BRCA_Normal"))

# Verificar distribución de los datos
cat("VERIFICACION DE LA DISTRIBUCIÓN DE LAS MUESTRAS:"); print(table(metadatos$MENOPAUSE_STATUS)); print(table(metadatos$SUBTYPE))

# Eliminar subtipo 'normal' del estudio por contener un número de muestras insuficiente
metadatos = metadatos[metadatos$SUBTYPE != "BRCA_Normal", ]
metadatos$SUBTYPE = droplevels(metadatos$SUBTYPE) # Asegurar que el factor SUBTYPE no incluye niveles vacíos
dge = dge[, colnames(dge) %in% metadatos$sample_id]
write.csv(dge, file.path(output_dir, "dge_sin_normal"), row.names = FALSE)
cat(sprintf("Archivo 'dge' guardado en: %s\n", output_dir))
write.csv(metadatos, file.path(output_dir, "metadatos_sin_normal"), row.names = FALSE)
cat(sprintf("Archivo 'metadatos' guardado en: %s\n", output_dir))


## DEA para alcanzar el Objetivo 1
# PREPARAR VARIABLES Y DISEÑO
metadatos$MENOPAUSE_STATUS = relevel(metadatos$MENOPAUSE_STATUS, ref = "Premenopausia") # Definir el nivel de referencia para los contrastes de menopausia
metadatos_obj1 = metadatos # Crear copia de los metadatos para el análisis del Objetivo 1
design_obj1 = model.matrix(~ MENOPAUSE_STATUS + SUBTYPE, data = metadatos_obj1) # Construir el modelo de diseño para el GLM sin interacciones
dge_obj1 = estimateDisp(dge, design_obj1) # Estimar la dispersión para el modelo de GLM
fit_obj1 = glmQLFit(dge_obj1, design_obj1) # Ajustar el modelo GLM usando la dispersión estimada

# DEFINIR CONTRASTES
contrastes_obj1 = makeContrasts(
  Post_vs_Pre = MENOPAUSE_STATUSPostmenopausia,
  Post_vs_Peri = MENOPAUSE_STATUSPostmenopausia - MENOPAUSE_STATUSPerimenopausia,
  Peri_vs_Pre = MENOPAUSE_STATUSPerimenopausia,
  levels = design_obj1
)

# APLICAR CONTRASTES Y FILTRAR RESULTADOS
resultados_completos_obj1 = list()
ids_genes_significativos_obj1 = list()

p_valor_umbral = 0.05    # FDR (p-valor ajustado)
logfc_umbral = 0.5       # Log2 Fold Change mínimo

for (nombre_contraste in colnames(contrastes_obj1)) {
  cat(paste0("\nContraste: ", nombre_contraste))
  qlf_obj1 = glmQLFTest(fit_obj1, contrast = contrastes_obj1[, nombre_contraste])    # Test de razón de verosimilitud para el contraste actual
  resultados_contraste = topTags(qlf_obj1, n = Inf, adjust.method = "BH")$table      # Obtener tabla completa de resultados ajustados por BH
  resultados_completos_obj1[[nombre_contraste]] = resultados_contraste               # Guardar tabla completa
  genes_significativos = subset(resultados_contraste, FDR < p_valor_umbral & abs(logFC) > logfc_umbral)  # Filtrar genes significativos según FDR y logFC
  ids_genes_significativos_obj1[[nombre_contraste]] = rownames(genes_significativos)   # Guardar únicamente los nombres de los genes significativos

# GUARDAR RESULTADOS EN ARCHIVOS CSV
dir_resultados_obj1 = "../resultados/5.DEGs_obj1_Estado_Menopausico_indep_Subtipo"
if (!dir.exists(dir_resultados_obj1)) dir.create(dir_resultados_obj1, recursive = TRUE)

for (contraste in names(resultados_completos_obj1)) { # Guardar tablas completas
  archivo_csv = file.path(dir_resultados_obj1, paste0("DEGs_", contraste, ".csv"))
  write.csv(resultados_completos_obj1[[contraste]], file = archivo_csv, row.names = TRUE)
  cat("Resultados guardados en:", archivo_csv, "\n")
}

resumen_significativos_obj1 = data.frame( # Guardar resumen DEGs significativos por contraste
  contraste = names(ids_genes_significativos_obj1),
  num_genes_significativos = sapply(ids_genes_significativos_obj1, length)
)

archivo_resumen = file.path(dir_resultados_obj1, "Resumen_Num_Genes_Significativos.csv")
write.csv(resumen_significativos_obj1, file = archivo_resumen, row.names = FALSE)
cat("Resumen de DEGs guardado en:", archivo_resumen, "\n")

# GENERAR VOLCANO PLOTS
make_volcano = function(res_df, contraste, dir_salida, p_cutoff=0.05, logfc_cutoff=0.5) {
  res_df$significance = "NS"
  res_df$significance[res_df$FDR < p_cutoff & res_df$logFC > logfc_cutoff] = "Up"
  res_df$significance[res_df$FDR < p_cutoff & res_df$logFC < -logfc_cutoff] = "Down"
  
  top_genes = rownames(head(res_df[order(res_df$FDR), ], 15))
  
  p = ggplot(res_df, aes(x=logFC, y=-log10(FDR), color=significance)) +
    geom_point(alpha=0.7, size=1.5) +
    scale_color_manual(values=c("Down"="blue", "Up"="red", "NS"="grey")) +
    geom_text_repel(data=subset(res_df, rownames(res_df) %in% top_genes),
                    aes(label=rownames(subset(res_df, rownames(res_df) %in% top_genes))),
                    size=3, max.overlaps=10) +
    theme_minimal(base_size=14) +
    labs(title=paste("Volcano plot -", contraste),
         x="log2 Fold Change", y="-log10(FDR)") +
    theme(legend.position="right")
  
  ggsave(filename=file.path(dir_salida, paste0("Volcano_", contraste, ".pdf")),
         plot=p, width=7, height=6)
}

dir_volcanos_obj1 = file.path(dir_resultados_obj1, "volcano_plots")
if (!dir.exists(dir_volcanos_obj1)) dir.create(dir_volcanos_obj1)

for (contraste in names(resultados_completos_obj1)) { 
  res_df = resultados_completos_obj1[[contraste]]
  make_volcano(res_df, contraste, dir_volcanos_obj1, p_valor_umbral, logfc_umbral)
}


## DEA para alcanzar el Objetivo 2
# PREPARAR VARIABLES Y DISEÑO
metadatos_obj2 = metadatos
design_obj2 = model.matrix(~ 0 + MENOPAUSE_STATUS:SUBTYPE, data = metadatos_obj2)
colnames(design_obj2) = gsub("MENOPAUSE_STATUS", "", colnames(design_obj2))
colnames(design_obj2) = gsub(":", "_", colnames(design_obj2))
colnames(design_obj2) = gsub("_SUBTYPEBRCA_", "_", colnames(design_obj2))
dge_obj2 = estimateDisp(dge, design_obj2)
fit_obj2 = glmQLFit(dge_obj2, design_obj2)

# DEFINIR CONTRASTES
contrastes_obj2 = makeContrasts(
  LumA_Post_vs_Pre  = Postmenopausia_LumA - Premenopausia_LumA,
  LumA_Post_vs_Peri = Postmenopausia_LumA - Perimenopausia_LumA,
  LumA_Peri_vs_Pre  = Perimenopausia_LumA - Premenopausia_LumA,

  LumB_Post_vs_Pre  = Postmenopausia_LumB - Premenopausia_LumB,
  LumB_Post_vs_Peri = Postmenopausia_LumB - Perimenopausia_LumB,
  LumB_Peri_vs_Pre  = Perimenopausia_LumB - Premenopausia_LumB,

  Basal_Post_vs_Pre  = Postmenopausia_Basal - Premenopausia_Basal,
  Basal_Post_vs_Peri = Postmenopausia_Basal - Perimenopausia_Basal,
  Basal_Peri_vs_Pre  = Perimenopausia_Basal - Premenopausia_Basal,

  Her2_Post_vs_Pre  = Postmenopausia_Her2 - Premenopausia_Her2,
  Her2_Post_vs_Peri = Postmenopausia_Her2 - Perimenopausia_Her2,
  Her2_Peri_vs_Pre  = Perimenopausia_Her2 - Premenopausia_Her2,

  levels = design_obj2
)

# APLICAR CONTRASTES Y FILTRAR RESULTADOS
resultados_completos_obj2 = list()
ids_genes_significativos_obj2 = list()

for (nombre_contraste in colnames(contrastes_obj2)) {
  cat(paste0("\nContraste: ", nombre_contraste))
  qlf_obj2 = glmQLFTest(fit_obj2, contrast = contrastes_obj2[, nombre_contraste])
  resultados_contraste_obj2 = topTags(qlf_obj2, n = Inf, sort.by = "PValue", adjust.method = "BH", p.value = 1)
  df_res_obj2 = as.data.frame(resultados_contraste_obj2$table)
  genes_signif_obj2 = df_res_obj2 %>% filter(FDR < p_valor_umbral, abs(logFC) > logfc_umbral)
  resultados_completos_obj2[[nombre_contraste]] = df_res_obj2
  ids_genes_significativos_obj2[[nombre_contraste]] = rownames(genes_signif_obj2)

# GUARDAR RESULTADOS EN ARCHIVOS CSV
dir_resultados_obj2 = "../resultados/5.DEGs_obj2_Estado_Menopausico_x_Subtipo"
if (!dir.exists(dir_resultados_obj2)) dir.create(dir_resultados_obj2, recursive = TRUE)

for (contraste in names(resultados_completos_obj2)) { 
  archivo_csv = file.path(dir_resultados_obj2, paste0("DEGs_", contraste, ".csv"))
  write.csv(resultados_completos_obj2[[contraste]], file = archivo_csv, row.names = TRUE)
  cat("Resultados guardados en:", archivo_csv, "\n")
}

# Guardar resumen DEGs significativos por contraste
resumen_significativos_obj2 = data.frame(
  contraste = names(ids_genes_significativos_obj2),
  num_genes_significativos = sapply(ids_genes_significativos_obj2, length)
)

archivo_resumen2 = file.path(dir_resultados_obj2, "Resumen_Num_Genes_Significativos.csv")
write.csv(resumen_significativos_obj2, file = archivo_resumen2, row.names = FALSE)
cat("Resumen de DEGs guardado en:", archivo_resumen2, "\n")

# GENERAR VOLCANO PLOTS
dir_volcanos_obj2 = file.path(dir_resultados_obj2, "volcano_plots")
if (!dir.exists(dir_volcanos_obj2)) dir.create(dir_volcanos_obj2)

for (contraste in names(resultados_completos_obj2)) {
  res_df = resultados_completos_obj2[[contraste]]
  make_volcano(res_df, contraste, dir_volcanos_obj2, p_valor_umbral, logfc_umbral)
}
