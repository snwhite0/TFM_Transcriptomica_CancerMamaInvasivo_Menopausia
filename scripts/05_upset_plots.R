### Upset Plots
resultados_dir = "../resultados/6.Superposicion_DEGs"
upset_plots_dir = file.path(resultados_dir, "upset_plots")

if (!dir.exists(upset_plots_dir)) {
  dir.create(upset_plots_dir, recursive = TRUE)
}

## Upset Plots para los resultados del DEA de cada objetivo (total = 2)
# Crear función para generar y guardar los plots
crear_y_guardar_upset_plot = function(lista_genes_significativos, nombre_objetivo) { 
  if (length(lista_genes_significativos) > 1 && sum(sapply(lista_genes_significativos, length) > 0) >= 2) {
    genes_for_upset = lista_genes_significativos[sapply(lista_genes_significativos, length) > 0]
    
    upset_plot = upset(fromList(genes_for_upset),
                       nsets = length(genes_for_upset),
                       keep.order = FALSE,
                       order.by = "freq",
                       decreasing = TRUE,
                       mainbar.y.label = "Tamaño intersección",
                       sets.x.label = "Tamaño del grupo",
                       query.legend = "bottom",
                       main.bar.color = "#E65113",  # naranja fuerte
                       matrix.color = "#5C4033",    # marrón chocolate
                       sets.bar.color = "#FBCBA4",  # melocotón claro 
                       #shade.color = "#fcebd5",    # beige claro
                      # number.angles = 45,
                      # point.size = 4,
                      # matrix.dot.alpha = 0.6,
                      # line.size = 0.8,
                      # shade.alpha = 0.5,
                      # set_size.numbers_size = 6,
                      # text.scale = 1.2
    )
    
    print(upset_plot)

    archivo_pdf_upset = file.path(upset_plots_dir, paste0("upset_plot_", nombre_objetivo, ".pdf"))
    pdf(file = archivo_pdf_upset, width = 12, height = 10)
    print(upset_plot)
    grid.text(paste("Superposición DEGs -", nombre_objetivo),
              x = 0.5, y = 0.98, just = "top", gp = gpar(fontsize = 20, fontface = "bold"))
    dev.off()

    cat(paste0("Upset Plot para ", nombre_objetivo, " guardado en: ", archivo_pdf_upset, "\n"))

  } else {
    cat(paste0("Advertencia: No se pudo generar el Upset Plot para ", nombre_objetivo,
               ". Asegúrese de que la lista contenga al menos dos conjuntos no vacíos.\n"))
  }
}

# Ejecutar función para cada objetivo
cat("Upset Plot: DEGs Compartidos y Específicos en Obj1"); crear_y_guardar_upset_plot(ids_genes_significativos_obj1, "obj1")
cat("Upset Plot: DEGs Compartidos y Específicos en Obj2"); crear_y_guardar_upset_plot(ids_genes_significativos_obj2, "obj2")

## Upset Plots para los resultados del DEA por subtipo molecular del objetivo 2 (total = 4)
target_subtypes = c("Basal", "Her2","LumA", "LumB")

for (subtipo in target_subtypes) {
  cat("Upset Plot: DEGs compartidos en el subtipo", subtipo)

  listas_subtipo = ids_genes_significativos_obj2[grep(paste0("^", subtipo, "_"), names(ids_genes_significativos_obj2))] # Filtrar resultados por subtipo
  listas_subtipo = listas_subtipo[sapply(listas_subtipo, length) > 0]   # Eliminar listas vacías
  
  if (length(listas_subtipo) >= 2) {
    upset_obj2 = upset(fromList(listas_subtipo),
                      nsets = length(listas_subtipo),
                      keep.order = TRUE,
                      order.by = "freq",
                      decreasing = TRUE,
                      mainbar.y.label = "Tamaño intersección",
                      sets.x.label = "Tamaño del grupo",
                      query.legend = "bottom",
                      main.bar.color = "#E65113",  # naranja fuerte
                      matrix.color = "#5C4033",    # marrón chocolate
                      sets.bar.color = "#FBCBA4",  # melocotón claro 
                      #shade.color = "#fcebd5",    # beige claro
                      number.angles = 45,
                      point.size = 5,
                      matrix.dot.alpha = 0.6,
                      line.size = 0.8,
                      shade.alpha = 0.5,
                      set_size.numbers_size = 6,
                      text.scale = 1.2
    )
    
    print(upset_obj2)

    pdf(file = file.path(upset_plots_dir, paste0("upset_plot_obj2_", subtipo, ".pdf")),
    width = 12, height = 9)
    
    print(upset_obj2)
    dev.off()
  }
}
