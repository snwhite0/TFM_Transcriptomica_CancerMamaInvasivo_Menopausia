### Heat Maps

## Configuración inicial
hm_dir = "../resultados/7.Perfiles_Expresion_DEGs/heat_maps"
dir_obj1 = file.path(hm_dir, "obj1")
dir_obj2 = file.path(hm_dir, "obj2")
dir.create(hm_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_obj1, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_obj2, showWarnings = FALSE, recursive = TRUE)

# Normalizar datos de expresión (log-CPM) y guardar en archivos CSV
logCPM_obj1 = cpm(dge_obj1, log = TRUE)
logCPM_obj2 = cpm(dge_obj2, log = TRUE)

write.csv(logCPM_obj1, file.path(output_dir, "logCPM_obj1"), row.names = FALSE)
cat(sprintf("Archivo 'logCPM_obj1' guardado en: %s\n", output_dir))
write.csv(logCPM_obj2, file.path(output_dir, "logCPM_obj2"), row.names = FALSE)
cat(sprintf("Archivo 'logCPM_obj2' guardado en: %s\n", output_dir))

# Definir paleta de colores para las anotaciones
cols_anno = list(
  MENOPAUSE_STATUS = c(
    Premenopausia   = "#E65113",   # Naranja 
    Perimenopausia  = "#B87333",   # Cobre rojizo 
    Postmenopausia  = "#FFB347"    # Ocre 
  ),
  SUBTYPE = c(
    Basal   = "#F7C6CE",  # Rosa pastel
    Her2    = "#8FD8A0",  # Verde pastel 
    LumA    = "#CBA0FF",  # Morado pastel
    LumB    = "#A0C8FF"   # Azul claro pastel
  )
)

# Homogeneizar nombres de subtipos
metadatos_obj1$SUBTYPE = gsub("BRCA_", "", metadatos_obj1$SUBTYPE)
metadatos_obj2$SUBTYPE = gsub("BRCA_", "", metadatos_obj2$SUBTYPE)

# Generar función principal para crear los heat maps
crear_y_guardar_heat_maps = function(genes, matriz_expr, metadatos, titulo, ruta_salida,
                                   anotaciones_a_incluir = c("MENOPAUSE_STATUS", "SUBTYPE"), # Eliminar 'SUBTYPE' para obj1
                                   cols_anno, cluster_cols = FALSE) {
  
  genes_presentes = genes[genes %in% rownames(matriz_expr)]   # Filtrar los genes presentes en la matriz de expresión
  if (length(genes_presentes) == 0) {
    cat("AVISO (", titulo, "): No hay genes presentes en la matriz de expresión. Omitiendo heatmap.\n")
    return(NULL)
  }
  
  expresion = matriz_expr[genes_presentes, , drop = FALSE]   # Extraer la expresión de los genes filtrados y manejar NAs
  expresion_filtrada_na = expresion[complete.cases(expresion), , drop = FALSE]   # Manejar NAs y asegurar que 'expresion' sigue siendo una matriz
  if (NROW(expresion_filtrada_na) == 0) {
    cat("AVISO (", titulo, "): Después de manejar NAs, no quedan genes válidos. Omitiendo heatmap.\n")
    return(NULL)
  }
  expresion = expresion_filtrada_na 

  muestras_comunes = intersect(colnames(expresion), metadatos$sample_id)   # Alinear las muestras entre la expresión y los metadatos
  if (length(muestras_comunes) == 0) {
    cat("AVISO (", titulo, "): No hay muestras comunes entre la matriz de expresión y los metadatos. Omitiendo heatmap.\n")
    return(NULL)
  }
  
  expresion = expresion[, muestras_comunes, drop = FALSE]
  if (NROW(expresion) < 2) {
    cat("AVISO (", titulo, "): Menos de 2 genes para clusterizar (N = ", NROW(expresion), "). Omitiendo clustering de filas y heatmap.\n")
    cluster_rows_final = FALSE # Desactivar el clustering de filas
    show_rows = TRUE 
  } else {
    cluster_rows_final = TRUE
    show_rows = NROW(expresion) < 100 # Decide si mostrar los nombres de los genes según la cantidad total de DEGs para el heatmap
  }
  
  if (cluster_cols && NCOL(expresion) < 2) {
    cat("AVISO (", titulo, "): 'cluster_cols' es TRUE pero menos de 2 muestras para clusterizar columnas (N = ", NCOL(expresion), "). Desactivando clustering de columnas.\n")
    cluster_cols = FALSE # Se fuerza a FALSE si no hay suficientes columnas
  }

  ann = metadatos[match(muestras_comunes, metadatos$sample_id), ]   # Preparar el dataframe de anotación para las columnas del heatmap
  rownames(ann) = ann$sample_id
  anotaciones_presentes = anotaciones_a_incluir[anotaciones_a_incluir %in% colnames(ann)]
  if (length(anotaciones_presentes) == 0) {
    cat("AVISO (", titulo, "): Ninguna columna de anotación especificada fue encontrada. El heatmap se generará sin anotaciones de columna.\n")
    ann = NA 
  } else {
    ann = ann[, anotaciones_presentes, drop = FALSE]

    if(is.vector(ann) && ncol(expresion) > 1) {  # Asegurar que 'ann' es un data.frame si se reduce a una sola columna/fila para evitar problemas de ordenamiento
        ann_df = data.frame(ann)
        colnames(ann_df) = anotaciones_presentes
        ann = ann_df
    } else if (is.matrix(ann) && !is.data.frame(ann)) {
        ann = as.data.frame(ann)
    }

    if ("MENOPAUSE_STATUS" %in% colnames(ann)) { # Ordenar las muestras por estado menopáusico si la columna existe
      orden_menopausia = names(cols_anno$MENOPAUSE_STATUS) # Usar el orden de la leyenda
      ann$MENOPAUSE_STATUS = factor(as.character(ann$MENOPAUSE_STATUS), levels = orden_menopausia)
    }

    if ("SUBTYPE" %in% colnames(ann)) {  # Ordenar las muestras por subtipo si la columna existe
      orden_subtipo = names(cols_anno$SUBTYPE) # Usar el orden de la leyenda
      ann$SUBTYPE = factor(as.character(ann$SUBTYPE), levels = orden_subtipo)
    }

    if ("SUBTYPE" %in% colnames(ann) && "MENOPAUSE_STATUS" %in% colnames(ann)) { # Ordenar el dataframe de anotación completo. Por SUBTYPE primero, luego por MENOPAUSE_STATUS.
        ann = ann[order(ann$SUBTYPE, ann$MENOPAUSE_STATUS), , drop = FALSE]
    } else if ("SUBTYPE" %in% colnames(ann)) {
        ann = ann[order(ann$SUBTYPE), , drop = FALSE]
    } else if ("MENOPAUSE_STATUS" %in% colnames(ann)) {
        ann = ann[order(ann$MENOPAUSE_STATUS), , drop = FALSE]
    }
    
    expresion = expresion[, rownames(ann), drop = FALSE]  # Reordenar la matriz de expresión para que coincida con el orden de las anotaciones
}
  
    hm = pheatmap(expresion, scale = "row", # Diseñar el heat map
         cluster_rows = cluster_rows_final, 
         cluster_cols = cluster_cols,       
         show_rownames = show_rows, 
         show_colnames = FALSE,
         fontsize = 8, fontsize_row = 8, border_color = NA,
         annotation_col = ann, annotation_colors = cols_anno,
         color = colorRampPalette(rev(brewer.pal(5,"RdBu")))(50),
         main = titulo)
    
    invisible({
      pdf(file = ruta_salida, width = 12, height = 10)
      print(hm) 
      dev.off()
  })
}


## Ejecutar función para crear heat maps para los resultados del DEA del Objetivo 1
# Heat map individual para cada contraste
ruta_contrastes_obj1_dir = file.path(dir_obj1, "contrastes_indiv")
dir.create(ruta_contrastes_obj1_dir, recursive = TRUE, showWarnings = FALSE)

for (nombre_contraste in names(ids_genes_significativos_obj1)) {
  cat("Generando heatmap para el contraste", nombre_contraste, "...\n")
  genes_contraste = ids_genes_significativos_obj1[[nombre_contraste]] 
  
  if (!is.null(genes_contraste) && length(genes_contraste) > 0) {
    crear_y_guardar_heat_maps(genes_contraste, logCPM_obj1, metadatos_obj1,
                            titulo = paste("Obj1: Contraste", nombre_contraste),
                            ruta_salida = file.path(ruta_contrastes_obj1_dir, paste0(nombre_contraste, ".pdf")), 
                            cols_anno = cols_anno)
  } else {
    cat("AVISO: El contraste '", nombre_contraste, "' (Obj1) no tiene genes significativos o no fue encontrado. Omitiendo heatmap.\n")
  }
}

# Heat maps para DEGs comunes entre N contrastes (N >= 2)
genes_comunes_a_todos_los_3_contrastes = Reduce(intersect, ids_genes_significativos_obj1)
all_obj1_contrasts_names = names(ids_genes_significativos_obj1)

ruta_comunes_inter_N_contr_obj1_dir = file.path(dir_obj1, "comunes_inter_N_contrasts")
dir.create(ruta_comunes_inter_N_contr_obj1_dir, recursive = TRUE, showWarnings = FALSE)

min_contrasts_for_intersection = 2
num_obj1_contrasts = length(all_obj1_contrasts_names)

if (num_obj1_contrasts < min_contrasts_for_intersection) {
  cat(sprintf("AVISO: Menos de %d contrastes para Obj1, omitiendo heatmaps para 7.3.1.\n", min_contrasts_for_intersection))
} else {
  for (n_contrasts_in_combo in min_contrasts_for_intersection:num_obj1_contrasts) {
    combinations_of_obj1_contrasts = combn(all_obj1_contrasts_names, n_contrasts_in_combo, simplify = FALSE)
    
    cat(sprintf("Generando heatmap de DEGs comunes entre %d contrastes ...\n", n_contrasts_in_combo))

    for (combo_contrasts in combinations_of_obj1_contrasts) {
      list_of_gene_sets_for_combo = ids_genes_significativos_obj1[combo_contrasts]
      common_genes_in_combo = Reduce(intersect, list_of_gene_sets_for_combo)
      
      if (length(common_genes_in_combo) > 0) {
        combo_name_file = paste(sort(combo_contrasts), collapse = "_")
        combo_name_title = paste(sort(combo_contrasts), collapse = " y ")
        
        title = paste("DEGs Comunes:", combo_name_title)
        filename = file.path(ruta_comunes_inter_N_contr_obj1_dir, paste0("DEGs_comunes_", combo_name_file, ".pdf"))
        
        crear_y_guardar_heat_maps(common_genes_in_combo, logCPM_obj1, metadatos_obj1,
                                titulo = title,
                                cols_anno = cols_anno,
                                ruta_salida = filename)
      } else {
        cat(sprintf("AVISO: No se encontraron DEGs comunes para la combinaci\u00f3n de contrastes: '%s'. Omitiendo heatmap.\n", paste(combo_contrasts, collapse = ", ")))
      }
    }
  }
}


## Ejecutar función para crear heat maps para los resultados del DEA del Objetivo 2
for (subtipo in target_subtypes) {
  cat("Generando heatmaps para el subtipo", subtipo, "...\n")

  dir_subtipo = file.path(dir_obj2, subtipo)
  dir.create(dir_subtipo, showWarnings = FALSE, recursive = TRUE)

  # Heat maps para cada contraste individual de cada subtipo
  contrastes_subtipo = grep(paste0("^", subtipo, "_"), names(ids_genes_significativos_obj2), value = TRUE)

  for (nombre_contraste in contrastes_subtipo) {
    genes = ids_genes_significativos_obj2[[nombre_contraste]]
    if (is.null(genes) || length(genes) == 0) next

    muestras = metadatos_obj2[metadatos_obj2$SUBTYPE == subtipo, ]$sample_id
    muestras_validas = intersect(muestras, colnames(logCPM_obj2))
    if (length(muestras_validas) == 0) next

    expr = logCPM_obj2[, muestras_validas, drop = FALSE]
    meta = metadatos_obj2[metadatos_obj2$sample_id %in% muestras_validas, , drop = FALSE]

    archivo_salida = file.path(dir_subtipo, paste0(nombre_contraste, ".pdf"))
    crear_y_guardar_heat_maps(
      genes, expr, meta,
      titulo = paste("Contraste", nombre_contraste),
      cols_anno = cols_anno,
      ruta_salida = archivo_salida
    )
  }
  
  # Heat maps para DEGs comunes entre N contrastes (N >= 2) de un subtipo dado
  genes_list = lapply(contrastes_subtipo, function(x) ids_genes_significativos_obj2[[x]])
  names(genes_list) = contrastes_subtipo
  genes_list = Filter(function(x) !is.null(x) && length(x) > 0, genes_list)

  if (length(genes_list) >= 1) {
    for (n_contrastes in 2:length(genes_list)) {
      combos = combn(names(genes_list), n_contrastes, simplify = FALSE)
      
      for (combo in combos) {
        genes_comunes = Reduce(intersect, genes_list[combo])
        
        if (length(genes_comunes) > 0) {
          muestras_validas = intersect(
            metadatos_obj2[metadatos_obj2$SUBTYPE == subtipo, ]$sample_id,
            colnames(logCPM_obj2)
          )
  
          if (length(muestras_validas) > 0) {
            expr = logCPM_obj2[, muestras_validas, drop = FALSE]
            meta = metadatos_obj2[metadatos_obj2$sample_id %in% muestras_validas, , drop = FALSE]
  
            nombre_combo = paste(combo, collapse = "_")
            crear_y_guardar_heat_maps(
              genes_comunes, expr, meta,
              titulo = paste("DEGs comunes en", nombre_combo),
              cols_anno = cols_anno,
              ruta_salida = file.path(dir_subtipo, paste0("DEGs_comunes_", nombre_combo, ".pdf"))
            )
          }
        } else {
          cat("  ⚠ No hay DEGs comunes para la combinación:", paste(combo, collapse = ", "), "\n")
        }
      }
    }
  } else {
    cat("  ⚠ No hay DEGs para generar heatmaps comunes en", subtipo, "\n")
  }

}

# Determinar y guardar los DEGs comunes entre contrastes para cada subtipo en un archivo CSV
genes_dir = file.path(dir_obj2, "genes_comunes_listas")
dir.create(genes_dir, showWarnings = FALSE, recursive = TRUE)

for (subtipo in target_subtypes) {
  contrastes_subtipo = grep(paste0("^", subtipo, "_"), names(ids_genes_significativos_obj2), value = TRUE)

  genes_list = lapply(contrastes_subtipo, function(x) ids_genes_significativos_obj2[[x]])
  names(genes_list) = contrastes_subtipo
  genes_list = Filter(function(x) !is.null(x) && length(x) > 0, genes_list)

  if (length(genes_list) >= 2) {
    for (n_contrastes in 2:length(genes_list)) {
      combos = combn(names(genes_list), n_contrastes, simplify = FALSE)

      for (combo in combos) {
        genes_comunes = Reduce(intersect, genes_list[combo])
        if (length(genes_comunes) > 0) {
          nombre_combo = paste(combo, collapse = "_")
          archivo_genes = file.path(genes_dir, paste0("DEGs_comunes_", subtipo, "_", nombre_combo, ".txt"))
          write.table(genes_comunes, archivo_genes, quote = FALSE, row.names = FALSE, col.names = FALSE)
          cat("Guardados", length(genes_comunes), "genes comunes para:", nombre_combo, "en", archivo_genes, "\n")
        }
      }
    }
  } else {
    cat("⚠ No hay suficientes contrastes con DEGs para generar combinaciones en", subtipo, "\n")
  }
}                   
