### Análisis de Supervivencia

## Configuración inicial
outdir_supervivencia = "../resultados/10.Analisis_Supervivencia/obj1"
dir.create(outdir_supervivencia, recursive = TRUE, showWarnings = FALSE)
outdir_supervivencia = "../resultados/10.Analisis_Supervivencia/obj2"
dir.create(outdir_supervivencia, recursive = TRUE, showWarnings = FALSE)

# Preparación de metadatos
# - Crear variables binarias de eventos de supervivencia:
#     · OS_EVENT: Overall Survival (1=deceased, 0=alive)
#     · DSS_EVENT: Disease-Specific Survival (1=dead with tumor, 0=alive/tumor-free)
#     · DFS_EVENT: Disease-Free Survival (1=recurred/progressed, 0=disease-free)
#     · PFS_EVENT: Progression-Free Survival (1=progression, 0=censored)

metadatos_obj1$OS_EVENT = ifelse(substr(metadatos_obj1$OS_STATUS, 1, 1) == "1", 1, 0)
metadatos_obj1$DSS_EVENT = ifelse(metadatos_obj1$DSS_STATUS == "", NA,
                            ifelse(substr(metadatos_obj1$DFS_STATUS,1,1)=="1",1,0))
metadatos_obj1$DFS_EVENT = ifelse(metadatos_obj1$DFS_STATUS == "", NA,
                            ifelse(substr(metadatos_obj1$DFS_STATUS,1,1)=="1",1,0))
metadatos_obj1$PFS_EVENT = ifelse(substr(metadatos_obj1$PFS_STATUS, 1, 1) == "1", 1, 0) 

metadatos_obj2$OS_EVENT = ifelse(substr(metadatos_obj2$OS_STATUS, 1, 1) == "1", 1, 0) 
metadatos_obj2$DSS_EVENT = ifelse(metadatos_obj2$DSS_STATUS == "", NA,
                            ifelse(substr(metadatos_obj2$DFS_STATUS,1,1)=="1",1,0)) 
metadatos_obj2$DFS_EVENT = ifelse(metadatos_obj2$DFS_STATUS == "", NA,
                            ifelse(substr(metadatos_obj2$DFS_STATUS,1,1)=="1",1,0)) 
metadatos_obj2$PFS_EVENT = ifelse(substr(metadatos_obj2$PFS_STATUS, 1, 1) == "1", 1, 0) 

# Definir lista de tipos de supervivencia con columnas de tiempo y evento
supervivencias = list(
  OS  = list(time = "OS_MONTHS",  event = "OS_EVENT"),
  DSS = list(time = "DSS_MONTHS", event = "DSS_EVENT"),
  DFS = list(time = "DFS_MONTHS", event = "DFS_EVENT"),
  PFS = list(time = "PFS_MONTHS", event = "PFS_EVENT")
)

# Obtener listado de DEGs a analizar a partir de los boxplots de Obj1/2
genes_listado = obtener_genes_boxplots(ruta_boxplots_obj1)
genes_a_analizar_superv_obj1 = unique(unlist(genes_listado))
genes_listado = obtener_genes_boxplots(ruta_boxplots_obj2)
genes_a_analizar_superv_obj2 = unique(unlist(genes_listado))

# Generar función de análisis de supervivencia Kaplan-Meier por mediana
analisis_supervivencia = function(df, time_col, event_col, gen_expr_col, tipo, carpeta, subtitulo="") {
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Filtrar valores NA
  expr_vec = as.numeric(df[[gen_expr_col]])
  tiempo = df[[time_col]]
  evento = df[[event_col]]
  
  idx_validos = !is.na(expr_vec) & !is.na(tiempo) & !is.na(evento)
  if(sum(idx_validos) < 2) return(NULL)
  
  expr_vec = expr_vec[idx_validos]
  tiempo = tiempo[idx_validos]
  evento = evento[idx_validos]
  
  # 2. Dividir muestras en grupos "ALTO" y "BAJO" según mediana
  grupo = ifelse(expr_vec >= median(expr_vec, na.rm = TRUE), "ALTO", "BAJO")
  grupo = factor(grupo, levels = c("ALTO", "BAJO"))
  
  N_ALTO = sum(grupo == "ALTO", na.rm = TRUE)
  N_BAJO = sum(grupo == "BAJO", na.rm = TRUE)
  
  # 3. Ajustar modelo Kaplan-Meier y test log-rank
  surv_obj = Surv(time = tiempo, event = evento)
  fit = survfit(surv_obj ~ grupo)
  lr = survdiff(surv_obj ~ grupo)
  pval = 1 - pchisq(lr$chisq, df=1)
  
  # 4. Generar y guardar plot PDF
  fname = sprintf("%s/%s%s.pdf", carpeta, tipo, ifelse(subtitulo=="","",paste0("_",subtitulo)))
  pdf(fname)
  plot(fit, col = c("#E65113","#FF8C00"), lwd=3,
       main = sprintf("%s %s (p=%.4f)", tipo, subtitulo, pval),
       xlab="Meses", ylab="Probabilidad de supervivencia")
  legend("bottomleft",
         legend = c(sprintf("ALTO (n=%d)", N_ALTO), sprintf("BAJO (n=%d)", N_BAJO)),
         col = c("#E65113","#FF8C00"), lwd=3)
  dev.off()
  
  # 5. Retornar resumen con conteos y p-valor
  return(data.frame(
    Tipo = tipo,
    Menopausia = ifelse(subtitulo=="","Global",subtitulo),
    N_total = length(expr_vec),
    N_ALTO = N_ALTO,
    N_BAJO = N_BAJO,
    pval_logrank = pval,
    stringsAsFactors = FALSE
  ))
}


## Ejecución deL análisis de supervivencia para los resultados significativos del Objetivo 1
resumen_final_obj1 = list()

# Iterar sobre cada gen significativo
for (gen in genes_a_analizar_superv_obj1) {
  if (!gen %in% rownames(logCPM_obj1)) next
  
  # Guardar datos de expresión en los metadatos
  metadatos_obj1[[paste0(gen,"_expr")]] = logCPM_obj1[gen, match(metadatos_obj1$sample_id, colnames(logCPM_obj1))]
  
  carpeta_gen = file.path(outdir_supervivencia, gen)
  dir.create(carpeta_gen, recursive = TRUE, showWarnings = FALSE)
  
  # ANÁLISIS GLOBAL
  res_global = do.call(rbind, lapply(names(supervivencias), function(tipo) {
    analisis_supervivencia(metadatos_obj1,
                           time_col=supervivencias[[tipo]]$time,
                           event_col=supervivencias[[tipo]]$event,
                           gen_expr_col=paste0(gen,"_expr"),
                           tipo=tipo,
                           carpeta=file.path(carpeta_gen,"supervivencia_global"))
  }))
  
  # ANÁLISIS ESTRATIFICADO POR MENOPAUSIA
  res_menop = do.call(rbind, lapply(unique(metadatos_obj1$MENOPAUSE_STATUS), function(estado) {
    datos_filtrados = subset(metadatos_obj1, MENOPAUSE_STATUS==estado)
    if (nrow(datos_filtrados) < 5) return(NULL)
    do.call(rbind, lapply(names(supervivencias), function(tipo) {
      analisis_supervivencia(datos_filtrados,
                             time_col=supervivencias[[tipo]]$time,
                             event_col=supervivencias[[tipo]]$event,
                             gen_expr_col=paste0(gen,"_expr"),
                             tipo=tipo,
                             carpeta=file.path(carpeta_gen,"supervivencia_menopausia"),
                             subtitulo=estado)
    }))
  }))
  
  # Guardar resultados combinados en CSV por gen
  resumen_final_obj1[[gen]] = rbind(res_global, res_menop)
  write.csv(resumen_final_obj1[[gen]],
            file.path(carpeta_gen, sprintf("resumen_supervivencia_%s.csv", gen)),
            row.names=FALSE)
}


## Ejecución deL análisis de supervivencia para los resultados significativos del Objetivo 2
resumen_final_obj2 = list()

# Iterar sobre cada gen significativo
for (gen in genes_a_analizar_superv_obj2) {
  if (!gen %in% rownames(logCPM_obj2)) next
  
  # Guardar datos de expresión en los metadatos
  metadatos_obj2[[paste0(gen,"_expr")]] = logCPM_obj2[gen, match(metadatos_obj2$sample_id, colnames(logCPM_obj2))]
  
  carpeta_gen = file.path(outdir_supervivencia, gen)
  dir.create(carpeta_gen, recursive = TRUE, showWarnings = FALSE)
  
  # ANÁLISIS GLOBAL
  res_global = do.call(rbind, lapply(names(supervivencias), function(tipo) {
    analisis_supervivencia(metadatos_obj2,
                           time_col=supervivencias[[tipo]]$time,
                           event_col=supervivencias[[tipo]]$event,
                           gen_expr_col=paste0(gen,"_expr"),
                           tipo=tipo,
                           carpeta=file.path(carpeta_gen,"supervivencia_global"))
  }))
  
  # ANÁLISIS ESTRATIFICADO POR SUBTIPO
  res_subtipo = do.call(rbind, lapply(unique(metadatos_obj2$SUBTYPE), function(subtipo) {
    datos_filtrados = subset(metadatos_obj2, SUBTYPE==subtipo)
    if (nrow(datos_filtrados) < 5) return(NULL)
    do.call(rbind, lapply(names(supervivencias), function(tipo) {
      analisis_supervivencia(datos_filtrados,
                             time_col=supervivencias[[tipo]]$time,
                             event_col=supervivencias[[tipo]]$event,
                             gen_expr_col=paste0(gen,"_expr"),
                             tipo=tipo,
                             carpeta=file.path(carpeta_gen,"supervivencia_subtipo"),
                             subtitulo=subtipo)
    }))
  }))
  
  # ANALISIS ESTRATIFICADO POR MENOPAUSIA Y SUBTIPO
  res_menop_subtipo = do.call(rbind, lapply(unique(metadatos_obj2$MENOPAUSE_STATUS), function(menop) {
    do.call(rbind, lapply(unique(metadatos_obj2$SUBTYPE), function(subtipo) {
      datos_filtrados = subset(metadatos_obj2, MENOPAUSE_STATUS==menop & SUBTYPE==subtipo)
      if (nrow(datos_filtrados) < 5) return(NULL)
      do.call(rbind, lapply(names(supervivencias), function(tipo) {
        analisis_supervivencia(datos_filtrados,
                               time_col=supervivencias[[tipo]]$time,
                               event_col=supervivencias[[tipo]]$event,
                               gen_expr_col=paste0(gen,"_expr"),
                               tipo=tipo,
                               carpeta=file.path(carpeta_gen,"supervivencia_menopausia_subtipo"),
                               subtitulo=paste(menop, subtipo, sep="_"))
      }))
    }))
  }))
  
  # Guardar resultados combinados en CSV por gen
  resumen_final_obj2[[gen]] = rbind(res_global, res_subtipo, res_menop_subtipo)
  write.csv(resumen_final_obj2[[gen]],
            file.path(carpeta_gen, sprintf("resumen_supervivencia_%s.csv", gen)),
            row.names=FALSE)
}
