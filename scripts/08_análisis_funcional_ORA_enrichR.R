### Análisis de Enriquecimiento Funcional (ORA) con EnrichR
# El análisis se aplica únicamente a los grupos de DEGs que han superado los filtros de significancia del test de Kruskal-Wallis, 
es decir, aquellos para los que se han generado boxplots.

## Configuración inicial del análisis funcional de DEGs
dir_analisis_funcional = file.path("..", "resultados", "9.Analisis_Funcional_DEGs")
dir.create(dir_analisis_funcional, recursive = TRUE, showWarnings = FALSE)

# Definir rutas a los boxplots generados previamente en Obj1 y Obj2
ruta_boxplots_obj1 = "~/Documentos/TFM_Sarasua_Naia/resultados/8.Significancia_DEGs/boxplots/obj1/"
ruta_boxplots_obj2 = "~/Documentos/TFM_Sarasua_Naia/resultados/8.Significancia_DEGs/boxplots/obj2/"

# Parámetros generales para el Over-Representation Analysis (ORA)
min_genes_ora = 1   # Número mínimo de DEGs requeridos para ejecutar el ORA
top_n = 3           # Número máximo de términos significativos a mostrar en los bubbleplots

# Generar función para extraer listas de genes a partir de subcarpetas con boxplots
#    - Busca archivos PDF de boxplots en cada subcarpeta
#    - Extrae los identificadores de genes a partir de los nombres de archivo
#    - Devuelve una lista nombrada con los genes por subcarpeta

obtener_genes_boxplots = function(ruta_boxplots) {
  subcarpetas = list.dirs(ruta_boxplots, recursive = TRUE, full.names = TRUE)
  lista_genes = list()
  
  for (subcarpeta in subcarpetas) {
    archivos_pdf = list.files(subcarpeta, pattern = "boxplot_.*\\.pdf$", full.names = TRUE)
    if (length(archivos_pdf) == 0) next
    
    nombre_lista = gsub(ruta_boxplots, "", subcarpeta)
    nombre_lista = gsub("/", "_", nombre_lista)
    
    genes_subcarpeta = sub("boxplot_(.*)\\.pdf$", "\\1", basename(archivos_pdf))
    lista_genes[[nombre_lista]] = unique(genes_subcarpeta)
  }
  
  return(lista_genes)
}

# Generar función para crear y guardar bubbleplots
#    - Lee resultados de ORA desde un archivo CSV
#    - Filtra términos con p < 0.05
#    - Selecciona los top N términos más significativos
#    - Genera bubbleplots (tamaño = nº genes, color = -log10(p))
#    - Guarda los plots en PDF

crear_y_guardar_bubble_plot = function(csv_file, top_n = 3, guardar = FALSE, output_dir_ORA = NULL) {
  ora_res = read_csv(csv_file, show_col_types = FALSE)
  ora_sig = ora_res %>% 
    filter(p_value < 0.05) %>% 
    arrange(p_value) %>% 
    slice_head(n = top_n) %>% 
    mutate(term_name = str_wrap(term_name, width = 40)) %>% 
    mutate(log_p = -log10(p_value))
  
  if (nrow(ora_sig) == 0) return(NULL)
  
  p = ggplot(ora_sig, aes(x = source, y = reorder(term_name, log_p))) +
    geom_point(aes(size = intersection_size, color = log_p), alpha = 0.8) +
    scale_color_gradientn(colors = c("yellow", "orange", "darkorange", "brown")) +
    labs(title = paste("Top", top_n, "GO -", tools::file_path_sans_ext(basename(csv_file))),
         x = "Fuente de la ontología", y = NULL,
         size = "Genes", color = "-log10(p)") +
    theme_minimal(base_size = 16)
  
  if (guardar && !is.null(output_dir_ORA)) {
    output_file_ORA = file.path(output_dir_ORA, paste0(tools::file_path_sans_ext(basename(csv_file)), "_bubbleplot.pdf"))
    ggsave(output_file_ORA, plot = p, width = 10, height = 6)
  }
  
  return(p)
}

# Generar función para ejecutar ORA con Enrichr
#    - Recibe una lista de genes y bases de datos de referencia (GO, KEGG, etc.)
#    - Filtra términos con p < 0.05
#    - Estandariza nombres de columnas para evitar errores (term_name, p_value, intersection, etc.)
#    - Calcula el tamaño de la intersección (nº genes asociados a cada término)
#    - Selecciona los top N términos significativos por base de datos
#    - Guarda resultados en CSV
#    - Genera y exporta un bubble plot resumen

ejec_enrichr_ora = function(gene_list,
                            analysis_name,
                            dir_analisis_funcional,
                            objeto = NULL,
                            databases = c("GO_Biological_Process_2025", 
                                          "GO_Molecular_Function_2025", 
                                          "GO_Cellular_Component_2025",
                                          "KEGG_2021_Human"
                                          #"Reactome Pathways 2024",
                                          #"GTEx Aging Signatures 2021",
                                          #"GTEx Tissue Expression Up",
                                          #"GTEx Tissue Expression Down",
                                          #"ClinVar 2025",
                                          #"OMIM Disease"
                                          )) {
  
  cat(sprintf("Ejecutando ORA con Enrichr para: %s\n", analysis_name))
  
  gene_list_symbols = gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list_symbols) < min_genes_ora) {
    cat("Genes insuficientes para análisis.\n")
    return(NULL)
  }
  
  enrich_results = enrichr(genes = gene_list_symbols, databases = databases)
  
  output_dir_ORA = file.path(dir_analisis_funcional, "enrichr_ora", objeto)
  dir.create(output_dir_ORA, recursive = TRUE, showWarnings = FALSE)
  
  todos_resultados = data.frame()
  
  db_nombres_cortos = c(   # Diccionario con nombres más claros y concisos para los bubbleplots
    "GO_Biological_Process_2025"   = "BP",
    "GO_Molecular_Function_2025"   = "MF",
    "GO_Cellular_Component_2025"   = "CC",
    "KEGG_2021_Human"              = "KEGG"
  )

  for (db in databases) {
    if (is.null(enrich_results[[db]]) || nrow(enrich_results[[db]]) == 0) next
    
    res_df = enrich_results[[db]]
    if ("Adjusted.P.value" %in% colnames(res_df)) {
      res_df = dplyr::rename(res_df, p_value = Adjusted.P.value)
    } else if ("P.value" %in% colnames(res_df)) {
      res_df = dplyr::rename(res_df, p_value = P.value)
    } else {
      res_df$p_value = NA
    }
    
    # Filtrado por p < 0.05
    res_df = dplyr::filter(res_df, p_value < 0.05)
    
    if (nrow(res_df) == 0) 
      next
    
    # Renombrar columnas estándar
    if ("Term" %in% colnames(res_df)) names(res_df)[names(res_df) == "Term"] = "term_name"
    if ("Genes" %in% colnames(res_df)) names(res_df)[names(res_df) == "Genes"] = "intersection"
    res_df$source = db_nombres_cortos[[db]]

    if ("intersection" %in% colnames(res_df)) {
      res_df$intersection_size = sapply(strsplit(res_df$intersection, ";"), length)
    }
    
    # Seleccionar top n
    res_df = res_df %>% arrange(p_value) %>% slice_head(n = top_n)
    
    todos_resultados = rbind(todos_resultados, res_df)
  }
  
  if (nrow(todos_resultados) == 0) {
    cat("Sin términos significativos para: ", analysis_name, "\n")
    return(NULL)
  }
  
  # Guardar CSV filtrado
  archivo_csv = file.path(output_dir_ORA, paste0(analysis_name, ".csv"))
  write.csv(todos_resultados, archivo_csv, row.names = FALSE)
  
  # Generar y guardar bubble plot
  todos_resultados$log_p = -log10(todos_resultados$p_value)
  todos_resultados$term_name_wrapped = stringr::str_wrap(todos_resultados$term_name, width = 40)
  
  p = ggplot(todos_resultados, aes(x = source, y = reorder(term_name_wrapped, log_p))) +
    geom_point(aes(size = intersection_size, color = log_p), alpha = 0.8) +
    scale_color_gradientn(colors = c("yellow", "orange", "darkorange", "brown")) +
    labs(title = paste("Top", top_n, "términos ORA -", analysis_name),
         x = "Fuente de ontología", y = NULL,
         size = "Genes", color = "-log10(p)") +
    theme_minimal(base_size = 16)
  
  archivo_plot = file.path(output_dir_ORA, paste0(analysis_name, "_bubbleplot.pdf"))
  ggsave(archivo_plot, plot = p, width = 10, height = 6)
  
  return(p)
}                            


## Ejecución del ORA para los resultados significativos del Objetivo 1
# - Se extraen los genes asociados a cada subcarpeta de boxplots
# - Para cada lista de genes con tamaño >= min_genes_ora:
#       · Se prepara un nombre corto de análisis
#       · Se ejecuta la función ejec_enrichr_ora()
#       · Resultados (CSV y bubbleplots) se guardan en la carpeta del análisis

lista_genes_obj1 = obtener_genes_boxplots(ruta_boxplots_obj1)

for (nombre_lista in names(lista_genes_obj1)) {
  genes = lista_genes_obj1[[nombre_lista]]
  
  if (length(genes) >= min_genes_ora) {
    nombre_corto = sub(".*(obj1.*)", "\\1", nombre_lista)
    
    cat(sprintf("Ejecutando ORA para subcarpeta OBJ1: %s (%d genes)\n", nombre_corto, length(genes)))
    ejec_enrichr_ora(
      gene_list = genes,
      analysis_name = nombre_corto,
      dir_analisis_funcional = dir_analisis_funcional,
      objeto = "obj1"
    )
  }
}


## Ejecución del ORA para los resultados significativos del Objetivo 2
# - Proceso equivalente al de Obj1, pero aplicado a las listas de genes obtenidas de los boxplots de Obj2

lista_genes_obj2 = obtener_genes_boxplots(ruta_boxplots_obj2)

for (nombre_lista in names(lista_genes_obj2)) {
  genes = lista_genes_obj2[[nombre_lista]]
  
  if (length(genes) >= min_genes_ora) {
    nombre_corto = sub(".*(obj2.*)", "\\1", nombre_lista)
    
    cat(sprintf("Ejecutando ORA para subcarpeta OBJ2: %s (%d genes)\n", nombre_corto, length(genes)))
    ejec_enrichr_ora(
      gene_list = genes,
      analysis_name = nombre_corto,
      dir_analisis_funcional = dir_analisis_funcional,
      objeto = "obj2"
    )
  }
}
