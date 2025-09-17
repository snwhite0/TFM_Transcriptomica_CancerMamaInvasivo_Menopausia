### Análisis de Componentes Principales (PCA)
# Crear y modificar matriz de expresión específicamente para PCA
dge_tumor = DGEList(counts = conteo_tumor_redondeado)
filtrado_tumor = filterByExpr(dge_tumor)
dge_tumor_filtrado = dge_tumor[filtrado_tumor, ]
dge_tumor_filtrado_norm = calcNormFactors(dge_tumor_filtrado)
expr_matriz_tumor = cpm(dge_tumor_filtrado_norm, log = TRUE)
varianzas_tumor = apply(expr_matriz_tumor, 1, var)
n_top_genes_pca = 500 # Valor seleccionado aleatoriamente
top_genes_pca_tumor = names(sort(varianzas_tumor, decreasing = TRUE)[1:n_top_genes_pca])
pca_resultado_tumor = prcomp(t(expr_matriz_tumor[top_genes_pca_tumor, ]), scale. = TRUE)
cat(capture.output(summary(pca_resultado_tumor))[1:9], sep = "\n")

# Crear data frame para ggplot
pca_df_tumor = data.frame(
  PC1 = pca_resultado_tumor$x[, 1],
  PC2 = pca_resultado_tumor$x[, 2],
  PC3 = pca_resultado_tumor$x[, 3],
  sample_id = rownames(pca_resultado_tumor$x)
)
pca_df_tumor = left_join(pca_df_tumor, datos, by = c("sample_id" = "SAMPLE_ID"))
rownames(pca_df_tumor) = pca_df_tumor$sample_id
  
# Eliminar filas con valores NA/vacías
pca_df_tumor = pca_df_tumor %>%
  filter(
    !is.na(SEX),
    !is.na(AGE),
    !is.na(MENOPAUSE_STATUS),
    !is.na(CANCER_TYPE_DETAILED),
    !is.na(AJCC_PATHOLOGIC_TUMOR_STAGE),
    AJCC_PATHOLOGIC_TUMOR_STAGE != "",
    !is.na(SUBTYPE),
    SUBTYPE != ""
  )

# Función para generar y guardar PCA plot 
resultados_dir = "../resultados/3.Analisis_Exploratorio/"
pca_plots_dir = file.path(resultados_dir, "pca_plots")
if (!dir.exists(pca_plots_dir)) dir.create(pca_plots_dir, recursive = TRUE)

crear_y_guardar_pca_plot = function(df, x_pc, y_pc, color_var, pca_resultado, titulo, filename) {
  var_x = round(summary(pca_resultado)$importance[2, as.numeric(sub("PC", "", x_pc))] * 100, 2)
  var_y = round(summary(pca_resultado)$importance[2, as.numeric(sub("PC", "", y_pc))] * 100, 2)

  xlab_text = paste0(x_pc, " (", var_x, "%)")
  ylab_text = paste0(y_pc, " (", var_y, "%)")

  color_var_is_numeric = is.numeric(df[[color_var]])

  p = ggplot(df, aes(x = !!sym(x_pc), y = !!sym(y_pc), colour = !!sym(color_var))) +
    geom_point(size = 2, alpha = 0.8) +
    (if (color_var_is_numeric) scale_colour_viridis(discrete = FALSE) else scale_colour_viridis(discrete = TRUE)) +
    xlab(xlab_text) + ylab(ylab_text) +
    ggtitle(titulo) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "right"
    )

  print(p)
  ggsave(filename = file.path(pca_plots_dir, filename), plot = p, width = 6, height = 8, units = "in", dpi = 300)
}

variables = list( # Listado de variables a analizar mediante PCA
  AGE = "Edad",
  MENOPAUSE_STATUS = "Estado Menopáusico",
  CANCER_TYPE_DETAILED = "Histología del Tumor",
  AJCC_PATHOLOGIC_TUMOR_STAGE = "Estadío del Cáncer",
  SUBTYPE = "Subtipo Molecular"
)

pcs_combos = list(
  c("PC1", "PC2"),
  c("PC1", "PC3"),
  c("PC2", "PC3")
)

for (var in names(variables)) {
  for (pcs in pcs_combos) {
    titulo = paste0("PCA (", pcs[1], " vs ", pcs[2], ") por ", variables[[var]])
    filename = paste0("pca_", tolower(var), "_", tolower(pcs[1]), "_", tolower(pcs[2]), ".pdf")
    crear_y_guardar_pca_plot(pca_df_tumor, pcs[1], pcs[2], var, pca_resultado_tumor, titulo, filename)
  }
}

ggsave(filename = file.path(pca_plots_dir, "pca_sexo.pdf"), plot = pca_plot_sexo, width = 6, height = 8, units = "in", dpi = 300)

