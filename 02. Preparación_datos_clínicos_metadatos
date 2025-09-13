# 1.1. Cargar archivos
tcga2018_clinicos_pacientes = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip=4)
tcga2018_clinicos_muestras = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt", skip=4)
tcga2015_clinicos_pacientes = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pub2015/data_clinical_patient.txt", skip=4)
tcga2015_clinicos_muestras = data.table::fread("~/Documentos/TFM_Sarasua_Naia/datos/crudos/brca_tcga_pub2015/data_clinical_sample.txt", skip=4)

# 1.2. Combinar datos
tcga2018_todos_clinicos = merge(tcga2018_clinicos_muestras, tcga2018_clinicos_pacientes, by = "PATIENT_ID", all = TRUE)
cat("Dimensión de la tabla de datos clínicos/metadatos INICIAL:", dim(tcga2018_todos_clinicos), "\n\n")

datos = merge(tcga2015_clinicos_pacientes[, c("PATIENT_ID","MENOPAUSE_STATUS")],
               tcga2018_todos_clinicos,
               by = "PATIENT_ID", all = FALSE) # Las nuevas columnas se fusionan solo con las filas de las muestras presentes en ambos conjuntos de datos.

datos = merge(tcga2015_clinicos_muestras[, c("PATIENT_ID","ER_STATUS_BY_IHC","ER_STATUS_IHC_PERCENT_POSITIVE","PR_STATUS_BY_IHC",
                "PR_STATUS_IHC_PERCENT_POSITIVE","IHC_HER2","HER2_IHC_PERCENT_POSITIVE")],
              datos,
              by = "PATIENT_ID", all = FALSE)

# 1.3. Seleccionar variables (columnas) de interés para el estudio
#nombres_con_indices = paste0(colnames(datos), " [", 1:ncol(datos), "]")
#print(nombres_con_indices)
datos = datos[, c(1,2,3,4,5,6,7,8,9,10,12,27,30,31,32,37,51,55,56,57,58,59,60,61,62)] 

# 1.4. Asignar a las filas vacías para la variable 'MENOPAUSE_STATUS' valores basándose en la edad
valores_distintos_menopause_antes = table(datos$MENOPAUSE_STATUS)
datos$MENOPAUSE_STATUS[datos$MENOPAUSE_STATUS == "[Not Available]"] = NA
datos$MENOPAUSE_STATUS[datos$MENOPAUSE_STATUS == "Indeterminate (neither Pre or Postmenopausal)"] = NA

asignar_estado_menopausico = function(AGE, MENOPAUSE_STATUS) {
  if (is.na(MENOPAUSE_STATUS)) {  # Solo imputar si el estado menopáusico es NA
    if (is.na(AGE)) {
      return(NA)  # Si la edad también es NA, no se puede imputar
    } else if (AGE >= 37 & AGE <= 45) {
      return("Premenopausia")
    } else if (AGE >= 46 & AGE <= 52) {
      return("Perimenopausia")
    } else if (AGE >= 53 & AGE <= 70) {
      return("Postmenopausia")
    } else {
      return(NA)  # Si la edad está fuera de los rangos, mantener NA
    }
  } else {
    return(MENOPAUSE_STATUS)  # Si el valor no es NA, mantener el valor original
  }
}

datos$MENOPAUSE_STATUS = mapply(asignar_estado_menopausico, datos$AGE, datos$MENOPAUSE_STATUS)

# Estandarizar nombres de los valores de la variable "MENOPAUSE_STATUS'
datos$MENOPAUSE_STATUS[datos$MENOPAUSE_STATUS == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"] = "Premenopausia"
datos$MENOPAUSE_STATUS[datos$MENOPAUSE_STATUS == "Peri (6-12 months since last menstrual period)"] = "Perimenopausia"
datos$MENOPAUSE_STATUS[datos$MENOPAUSE_STATUS == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"] = "Postmenopausia"

# 1.5. Agrupar valores de la variable 'AJCC_PATHOLOGIC_TUMOR_STAGE'
datos$AJCC_PATHOLOGIC_TUMOR_STAGE[datos$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE I", "STAGE IA", "STAGE IB")] = "STAGE I"
datos$AJCC_PATHOLOGIC_TUMOR_STAGE[datos$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II", "STAGE IIA", "STAGE IIB")] = "STAGE II"
datos$AJCC_PATHOLOGIC_TUMOR_STAGE[datos$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC")] = "STAGE III"

# Observar distribución de las muestras según 'MENOPAUSE_STATUS' y 'SUBTYPE'
cat("Dimensión de la tabla de datos clínicos/metadatos ANTES de limpieza:", dim(datos), "\n") 
cat("VERIFICACION DE LA DISTRIBUCIÓN DE LAS MUESTRAS ANTES DE LIMPIEZA:"); print(table(datos$MENOPAUSE_STATUS)); print(table(datos$MENOPAUSE_STATUS))

# 1.6 Eliminar filas que correspondan a tipos de cáncer que no son de interés para el estudio por contener una cantidad de muestras insuficiente. 
datos = filter(datos, CANCER_TYPE_DETAILED != "Metaplastic Breast Cancer")
datos = filter(datos, CANCER_TYPE_DETAILED != "Breast Invasive Mixed Mucinous Carcinoma")
datos = filter(datos, CANCER_TYPE_DETAILED != "Invasive Breast Carcinoma")
datos = filter(datos, SEX != "Male") 

# 1.7. Eliminar filas con NA/vacíos en variables clave del estudio
datos = datos[!is.na(datos$SEX), ]
datos = datos[!is.na(datos$AGE), ]
datos = datos[!is.na(datos$MENOPAUSE_STATUS), ]
datos = datos[!is.na(datos$CANCER_TYPE_DETAILED), ]
datos = datos[!is.na(datos$AJCC_PATHOLOGIC_TUMOR_STAGE), ]
datos = datos[datos$AJCC_PATHOLOGIC_TUMOR_STAGE != "", ]
datos = datos[!is.na(datos$SUBTYPE), ]
datos = datos[datos$SUBTYPE != "", ]

# Observar distribución de las muestras según 'MENOPAUSE_STATUS' y 'SUBTYPE'
cat("Dimensión de la tabla de datos clínicos/metadatos DESPUÉS de limpieza:", dim(datos), "\n")
cat("VERIFICACION DE LA DISTRIBUCIÓN DE LAS MUESTRAS DESPUÉS DE LIMPIEZA:"); print(table(datos$MENOPAUSE_STATUS)); print(table(datos$MENOPAUSE_STATUS))

cat("¡Listo! Preprocesamiento de los datos clínicos/metadatos completado.\n")
