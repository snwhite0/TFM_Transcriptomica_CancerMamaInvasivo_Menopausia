### Test estadístico no-paramétrico Kruskal-Wallis

## Configuración inicial
dir_boxplots = file.path("..", "resultados", "8.Significancia_DEGs", "boxplots")
dir.create(dir_boxplots, recursive = TRUE, showWarnings = FALSE)

# Inicializar un data.frame vacío para almacenar los resultados del test estadístico
all_kruskal_results = data.frame(
  Gene = character(),
  Subtype = character(),
  Contrast = character(), 
  Context = character(),
  P_Value = numeric(),
  Significance_Threshold = numeric(),
  Is_Significant = logical(),
  Reason_Not_Plotted = character(),
  stringsAsFactors = FALSE
)

# Generar función para crear y guardar box plots de DEGs significativos según el umbral de p-valor establecido
crear_y_guardar_box_plot = function(gen_id, expresion_matrix, metadatos_df,
                                         ruta_salida_dir,
                                         p_threshold = 0.05,
                                         contexto_analisis = "Global",
                                         contraste_actual = NA_character_,
                                         subtipo_actual = NA_character_,
                                         show_plot_in_rmarkdown = TRUE) {

  reason_not_plotted = NA_character_
  
  # 1. Validar que el gen exista en la matriz
  if (!gen_id %in% rownames(expresion_matrix)) {
    reason_not_plotted = "Gen no encontrado en matriz de expresión."
    all_kruskal_results <<- bind_rows(all_kruskal_results,
                                            data.frame(
                                              Gene = gen_id, Subtype = subtipo_actual, Contrast = contraste_actual,
                                              Context = contexto_analisis, 
                                              P_Value = NA_real_,
                                              Significance_Threshold = p_threshold, Is_Significant = FALSE,
                                              Reason_Not_Plotted = reason_not_plotted, stringsAsFactors = FALSE
                                            ))
    return(invisible(NULL))
  }
  
  # 2. Preparar datos de expresión para el gen
  expresion_gen_df = data.frame(
    sample_id = colnames(expresion_matrix),
    exp = as.numeric(expresion_matrix[gen_id, ])
  )

  datos_plot = expresion_gen_df %>% inner_join(metadatos_df, by = "sample_id")

  # 3. Filtrar por subtipo si corresponde
  if (!is.na(subtipo_actual) && "SUBTYPE" %in% colnames(datos_plot) && subtipo_actual != "Todos" && subtipo_actual != "NA_character_") {
    datos_plot = datos_plot %>% filter(SUBTYPE == subtipo_actual)
  }
  
  # Validar que haya datos tras el filtrado
  if (nrow(datos_plot) == 0) {
    reason_not_plotted = "No hay datos de muestras para el gen o subtipo."
    all_kruskal_results <<- bind_rows(all_kruskal_results,
                                            data.frame(
                                              Gene = gen_id, Subtype = subtipo_actual, Contrast = contraste_actual,
                                              Context = contexto_analisis, 
                                              P_Value = NA_real_,
                                              Significance_Threshold = p_threshold, Is_Significant = FALSE,
                                              Reason_Not_Plotted = reason_not_plotted, stringsAsFactors = FALSE
                                            ))
    return(invisible(NULL))
  }
  
  # 4. Reordenar factores y validar grupos
  orden_menopausia = c("Premenopausia", "Perimenopausia", "Postmenopausia")
  datos_plot = datos_plot %>%
    mutate(MENOPAUSE_STATUS = factor(MENOPAUSE_STATUS, levels = intersect(orden_menopausia, unique(MENOPAUSE_STATUS)))) %>%
    drop_na(MENOPAUSE_STATUS)

  grupos_presentes = unique(datos_plot$MENOPAUSE_STATUS)
  if (length(grupos_presentes) < 2 || any(table(datos_plot$MENOPAUSE_STATUS) < 2)) {
    reason_not_plotted = "Menos de 2 grupos o menos de 2 muestras por grupo para test estadístico."
    all_kruskal_results <<- bind_rows(all_kruskal_results,
                                            data.frame(
                                              Gene = gen_id, Subtype = subtipo_actual, Contrast = contraste_actual,
                                              Context = contexto_analisis,
                                              P_Value = NA_real_,
                                              Significance_Threshold = p_threshold, Is_Significant = FALSE,
                                              Reason_Not_Plotted = reason_not_plotted, stringsAsFactors = FALSE
                                            ))
    return(invisible(NULL))
  }
  
  # 5. Aplicar test Kruskal-Wallis
  test_result = tryCatch({
    kruskal.test(exp ~ MENOPAUSE_STATUS, data = datos_plot)
  }, error = function(e) {
    message(sprintf("ERROR Kruskal-Wallis: Gen %s, Subtipo %s, Contraste %s, Contexto %s. Detalle: %s",
                    gen_id, subtipo_actual, contraste_actual, contexto_analisis, e$message))
    return(NULL)
  })

  if (is.null(test_result)) {
    reason_not_plotted = "Error durante el test Kruskal-Wallis."
    all_kruskal_results <<- bind_rows(all_kruskal_results,
                                            data.frame(
                                              Gene = gen_id, Subtype = subtipo_actual, Contrast = contraste_actual,
                                              Context = contexto_analisis,
                                              P_Value = NA_real_,
                                              Significance_Threshold = p_threshold, Is_Significant = FALSE,
                                              Reason_Not_Plotted = reason_not_plotted, stringsAsFactors = FALSE
                                            ))
    return(invisible(NULL))
  }

  global_p_value = test_result$p.value
  test_label = "kruskal.test"
  
  # 6. Evaluar significancia
  is_significant_for_plot = !is.na(global_p_value) && global_p_value < p_threshold

  all_kruskal_results <<- bind_rows(all_kruskal_results,
                                          data.frame(
                                            Gene = gen_id, Subtype = subtipo_actual, Contrast = contraste_actual,
                                            Context = contexto_analisis,
                                            P_Value = global_p_value,
                                            Significance_Threshold = p_threshold, Is_Significant = is_significant_for_plot,
                                            Reason_Not_Plotted = if(is_significant_for_plot) NA_character_ else "P-value >= threshold",
                                            stringsAsFactors = FALSE
                                          ))
  if (!is_significant_for_plot) {
    message(sprintf("  P-valor para gen %s (Subtipo %s, Contraste %s): %.4f >= umbral %.4f. No se generará el boxplot.\n",
                    gen_id, subtipo_actual, contraste_actual, global_p_value, p_threshold))
    return(invisible(NULL)) # No mostrar plot si no es significativo
  }
  
  # 7. Generar boxplot con ggplot2
  colores_menopausia_map = cols_anno$MENOPAUSE_STATUS
  colores_actuales = colores_menopausia_map[names(colores_menopausia_map) %in% unique(datos_plot$MENOPAUSE_STATUS)]
  
  p = ggplot(datos_plot, aes(x = MENOPAUSE_STATUS, y = exp, fill = MENOPAUSE_STATUS)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6, size = 1) +
    scale_fill_manual(values = colores_actuales) +
    labs(
      title = paste0("Expresión de ", gen_id, "\n(Subtipo: ", subtipo_actual, ", Contraste: ", contraste_actual, ")"),
      x = "Estado Menopáusico",
      y = "Expresión (logCPM)",
      subtitle = paste0("P-valor (", test_label, "): ", sprintf("%.3e", global_p_value))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "darkred"),
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
  
  # 8. Agregar estrellas de significancia
  signif_stars = ""
  if (!is.na(global_p_value)) {
    if (global_p_value <= 0.0001) {
      signif_stars = "****"
    } else if (global_p_value <= 0.001) {
      signif_stars = "***"
    } else if (global_p_value <= 0.01) {
      signif_stars = "**"
    } else if (global_p_value <= 0.05) {
      signif_stars = "*"
    } else {
      signif_stars = "ns"
    }
  }

  max_y = max(datos_plot$exp, na.rm = TRUE)
  min_y = min(datos_plot$exp, na.rm = TRUE)
  y_range = max_y - min_y
  y_position_signif = max_y + (y_range * 0.15)

  if (is_significant_for_plot) {
    p = p + annotate("text",
                     x = length(unique(datos_plot$MENOPAUSE_STATUS)) / 2 + 0.5,
                     y = y_position_signif,
                     label = signif_stars,
                     size = 6, color = "black", fontface = "bold")
  }

  # 9. Guardar boxplot en archivo PDF
  nombre_archivo_limpio = gsub("[^A-Za-z0-9_.-]", "_", paste0("boxplot_", gen_id, ".pdf"))
  ruta_completa_salida = file.path(ruta_salida_dir, nombre_archivo_limpio)

  dir.create(ruta_salida_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(ruta_completa_salida, plot = p, width = 8, height = 6, units = "in", dpi = 300)
  message(sprintf("Boxplot guardado en: %s", ruta_completa_salida))
  
  # 10. Retornar plot en RMarkdown
  if (show_plot_in_rmarkdown) {
    return(p) # Retorna el objeto ggplot
  } else {
    return(invisible(NULL))
  }
}


## Aplicar análisis estadístico Kruskal-Wallis para los resultados del DEA del Obj1
ruta_base_contr_indiv_obj1 = file.path(dir_boxplots, "obj1", "contrastes_indiv")
dir.create(ruta_base_contr_indiv_obj1, recursive = TRUE, showWarnings = FALSE)

# A. APLICAR A LOS DEGs RESULTANTES DE CADA CONTRASTE INDIVIDUAL
for (nombre_contraste in names(ids_genes_significativos_obj1)) {
  
  # Crear directorio de salida para este contraste
  ruta_boxplots_contraste_obj1 = file.path(ruta_base_contr_indiv_obj1, nombre_contraste)
  dir.create(ruta_boxplots_contraste_obj1, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("Generando boxplots para el contraste '%s' ...\n", nombre_contraste))
  
  genes_contraste = ids_genes_significativos_obj1[[nombre_contraste]]

  if (!is.null(genes_contraste) && length(genes_contraste) > 0) {
    
    # Procesar cada gen significativo en este contraste
    for (gen in genes_contraste) {
      p_plot = crear_y_guardar_box_plot(gen_id = gen,
                                              expresion_matrix = logCPM_obj1,
                                              metadatos_df = metadatos_obj1,
                                              ruta_salida_dir = ruta_boxplots_contraste_obj1,
                                              p_threshold = 0.01,
                                              contexto_analisis = paste0("Obj1_ContrasteIndiv_", nombre_contraste),
                                              contraste_actual = nombre_contraste,
                                              subtipo_actual = "Todos",
                                              show_plot_in_rmarkdown = TRUE)
      if (!is.null(p_plot)) {
        print(p_plot) # Mostrar en RMarkdown
        cat("\n")
      }
    }
  
    # Guardar resultados estadísticos para el contraste dado
    results_for_this_context = all_kruskal_results %>%
      filter(Context == paste0("Obj1_ContrasteIndiv_", nombre_contraste))

    if (nrow(results_for_this_context) > 0) {
      write.csv(results_for_this_context,
                file.path(ruta_boxplots_contraste_obj1, paste0("kruskal_wallis_results_", nombre_contraste, ".csv")),
                row.names = FALSE)
      message(sprintf("  Resultados de Kruskal-Wallis guardados en: %s\n", file.path(ruta_boxplots_contraste_obj1, paste0("kruskal_wallis_results_", nombre_contraste, ".csv"))))
    } else {
      cat(sprintf("  No hay resultados de Kruskal-Wallis para guardar en Obj1_ContrasteIndiv_%s.\n", nombre_contraste))
    }
  } else {
    cat(sprintf("AVISO: El contraste '%s' (Obj1) no tiene genes significativos. Omitiendo boxplots para 8.3.2.\n", nombre_contraste))
  }
}

# B. APLICAR A LOS DEGs COMUNES ENTRE N CONTRASTES
ruta_base_comunes_inter_N_contr_obj1 = file.path(dir_boxplots, "obj1", "comunes_inter_N_contrasts")
dir.create(ruta_base_comunes_inter_N_contr_obj1, recursive = TRUE, showWarnings = FALSE)

min_contrasts_for_intersection = 2
num_obj1_contrasts = length(all_obj1_contrasts_names)

if (num_obj1_contrasts < min_contrasts_for_intersection) {
  cat(sprintf("AVISO: Menos de %d contrastes para Obj1, omitiendo boxplots para 8.3.1.\n", min_contrasts_for_intersection))
} else {
  
  # Iterar por número de contrastes en cada combinación
  for (n_contrasts_in_combo in min_contrasts_for_intersection:num_obj1_contrasts) {
    combinations_of_obj1_contrasts = combn(all_obj1_contrasts_names, n_contrasts_in_combo, simplify = FALSE)
    
    cat(sprintf("Generando boxplots de DEGs significativos entre %d contrastes...\n", n_contrasts_in_combo))

    for (combo_contrasts in combinations_of_obj1_contrasts) {
      
      # Genes comunes entre los contrastes de esta combinación
      list_of_gene_sets_for_combo = ids_genes_significativos_obj1[combo_contrasts]
      common_genes_in_combo = Reduce(intersect, list_of_gene_sets_for_combo)
      
      if (length(common_genes_in_combo) > 0) {
        combo_name_file = paste(sort(combo_contrasts), collapse = "_")
        combo_name_title = paste(sort(combo_contrasts), collapse = " y ")
        
        ruta_boxplots_comunes_obj1_combo = file.path(ruta_base_comunes_inter_N_contr_obj1, combo_name_file)
        dir.create(ruta_boxplots_comunes_obj1_combo, recursive = TRUE, showWarnings = FALSE)
        
        cat(sprintf("Procesando contrastes '%s'...\n", combo_name_title))
        
        # Procesar cada gen común en esta combinación
        for (gen in common_genes_in_combo) {
          p_plot = crear_y_guardar_box_plot(gen_id = gen,
                                                  expresion_matrix = logCPM_obj1,
                                                  metadatos_df = metadatos_obj1,
                                                  ruta_salida_dir = ruta_boxplots_comunes_obj1_combo, 
                                                  p_threshold = 0.05,
                                                  contexto_analisis = paste0("Obj1_Comunes_N_Contr_", combo_name_file),
                                                  contraste_actual = combo_name_title,
                                                  subtipo_actual = "Todos", 
                                                  show_plot_in_rmarkdown = TRUE)
          if (!is.null(p_plot)) {
            print(p_plot)
            cat("\n")
          }
        }
        
        # Guardar resultados estadísticos
        results_for_this_context = all_kruskal_results %>%
          filter(Context == paste0("Obj1_Comunes_N_Contr_", combo_name_file))

        if (nrow(results_for_this_context) > 0) {
          write.csv(results_for_this_context,
                    file.path(ruta_boxplots_comunes_obj1_combo, paste0("kruskal_wallis_results_", combo_name_file, ".csv")),
                    row.names = FALSE)
          message(sprintf("  Resultados de Kruskal-Wallis guardados en: %s\n", file.path(ruta_boxplots_comunes_obj1_combo, paste0("kruskal_wallis_results_", combo_name_file, ".csv"))))
        } else {
          cat(sprintf("  No hay resultados de Kruskal-Wallis para guardar en Obj1_Comunes_N_Contr_%s.\n", combo_name_file))
        }
      } else {
        cat(sprintf("AVISO: No se encontraron DEGs comunes para la combinaci\u00f3n de contrastes: '%s'. Omitiendo boxplots.\n", paste(combo_contrasts, collapse = ", ")))
      }
    }
  }
}


## Aplicar análisis estadístico Kruskal-Wallis para los resultados del DEA del Obj2
dir_individuales_bp = file.path(dir_boxplots, "obj2")
dir.create(dir_individuales_bp, showWarnings = FALSE)

# Iterar por cada subtipo definido en target_subtypes
for (subtipo in target_subtypes) {
  cat("Generando boxplots para el subtipo:", subtipo, "...\n")
  
  # Crear directorio para el subtipo dado
  dir_subtipo_bp = file.path(dir_individuales_bp, subtipo)
  dir.create(dir_subtipo_bp, showWarnings = FALSE, recursive = TRUE)

  # C. APLICAR A LOS DEGs RESULTANTES DE CADA CONTRASTE POR SUBTIPO
  contrastes_subtipo = grep(paste0("^", subtipo, "_"), names(ids_genes_significativos_obj2), value = TRUE) # Seleccionar contrastes que pertenecen a un subtipo dado
  
  for (nombre_contraste in contrastes_subtipo) {
    genes = ids_genes_significativos_obj2[[nombre_contraste]]
    if (is.null(genes) || length(genes) == 0) {
      cat("  No hay DEGs para el contraste", nombre_contraste, ". Omitiendo\n")
      next
    }

    # Directorio de salida
    ruta_salida_bp_contraste = file.path(dir_subtipo_bp, gsub("_", "-", nombre_contraste))
    dir.create(ruta_salida_bp_contraste, recursive = TRUE, showWarnings = FALSE)

    # Generar boxplots
    for (gen_id in genes) {
      p_plot = crear_y_guardar_box_plot(
        gen_id = gen_id,
        expresion_matrix = logCPM_obj2,
        metadatos_df = metadatos_obj2,
        ruta_salida_dir = ruta_salida_bp_contraste,
        p_threshold = 0.01,
        contexto_analisis = paste0("Obj2_ContrasteIndividual_", nombre_contraste),
        contraste_actual = nombre_contraste,
        subtipo_actual = subtipo,
        show_plot_in_rmarkdown = TRUE
      )

      if (!is.null(p_plot)) {
        print(p_plot)
        cat("\n")
      }
    }

    # Guardar resultados estadísticos
    context_filter_string = paste0("Obj2_ContrasteIndividual_", nombre_contraste)
    results_for_this_context = all_kruskal_results %>% filter(Context == context_filter_string)

    if (nrow(results_for_this_context) > 0) {
      write.csv(results_for_this_context,
                file.path(ruta_salida_bp_contraste, paste0("kruskal_wallis_results_", gsub("_", "-", nombre_contraste), ".csv")),
                row.names = FALSE)
      message(sprintf("  Resultados de Kruskal-Wallis guardados en: %s\n", file.path(ruta_salida_bp_contraste, paste0("kruskal_wallis_results_", gsub("_", "-", nombre_contraste), ".csv"))))
    } else {
      cat(sprintf("  No hay resultados de Kruskal-Wallis para guardar en %s.\n", context_filter_string))
    }
  }

  # D. APLICAR A LOS DEGs COMUNES ENTRE N CONTRASTES DEL MISMO SUBTIPO  
  # Guardar lista de genes por contraste filtrando vacíos
  genes_list = lapply(contrastes_subtipo, function(x) ids_genes_significativos_obj2[[x]])
  genes_list = Filter(function(x) !is.null(x) && length(x) > 0, genes_list)

  if (length(genes_list) == 0) {
    cat("  No hay DEGs en ningún contraste para", subtipo, "\n")
    next
  }

  # Analizar intersecciones para combinaciones de 2 hasta n contrastes
  for (n in 2:length(genes_list)) {
    cat("\n Analizando genes comunes en combinaciones de", n, "contrastes para", subtipo, "...\n")
    
    combos = combn(contrastes_subtipo, n, simplify = FALSE)
    for (combo in combos) {
      
      # Genes comunes en una combinación dada de contrastes
      genes_comunes_combo = Reduce(intersect, lapply(combo, function(x) ids_genes_significativos_obj2[[x]]))
      
      if (length(genes_comunes_combo) == 0) {
        cat("    No hay genes comunes para la combinación:", paste(combo, collapse = " & "), "\n")
      } else {
        cat("    Encontrados", length(genes_comunes_combo), "genes comunes para la combinación:", paste(combo, collapse = " & "), "\n")
        
        # Directorio de salida
        ruta_salida_bp_comunes_combo = file.path(dir_subtipo_bp, paste0("genes_comunes_", paste(combo, collapse = "_")))
        dir.create(ruta_salida_bp_comunes_combo, recursive = TRUE, showWarnings = FALSE)
        
        # Generar boxplots
        for (gen_id in genes_comunes_combo) {
          p_plot = crear_y_guardar_box_plot(
            gen_id = gen_id,
            expresion_matrix = logCPM_obj2,
            metadatos_df = metadatos_obj2,
            ruta_salida_dir = ruta_salida_bp_comunes_combo,
            p_threshold = 0.05,
            contexto_analisis = paste0("Obj2_IntraSubtipo_Comunes_", subtipo, "_combo"),
            contraste_actual = "DEGS_Comunes",
            subtipo_actual = subtipo,
            show_plot_in_rmarkdown = TRUE
          )

          if (!is.null(p_plot)) {
            print(p_plot)
            cat("\n")
          }
        }

        # Guardar resultados estadísticos
        context_filter_string = paste0("Obj2_IntraSubtipo_Comunes_", subtipo, "_combo")
        results_for_this_context = all_kruskal_results %>% filter(Context == context_filter_string)
        
        if (nrow(results_for_this_context) > 0) {
          write.csv(results_for_this_context,
                    file.path(ruta_salida_bp_comunes_combo, paste0("kruskal_wallis_results_comunes_", paste(combo, collapse = "_"), ".csv")),
                    row.names = FALSE)
          message(sprintf("    Resultados de Kruskal-Wallis guardados en: %s\n",
                          file.path(ruta_salida_bp_comunes_combo, paste0("kruskal_wallis_results_comunes_", paste(combo, collapse = "_"), ".csv"))))
        } else {
          cat(sprintf("    No hay resultados de Kruskal-Wallis para guardar en %s.\n", context_filter_string))
        }
      }
    }
  }
}
