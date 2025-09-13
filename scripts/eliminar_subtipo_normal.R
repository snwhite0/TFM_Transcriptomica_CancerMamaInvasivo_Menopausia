metadatos = metadatos[metadatos$SUBTYPE != "BRCA_Normal", ]
metadatos$SUBTYPE = droplevels(metadatos$SUBTYPE) # Asegurar que el factor SUBTYPE no incluye niveles vacíos

dge = dge[, colnames(dge) %in% metadatos$sample_id]

write.csv(dge, file.path(output_dir, "dge_sin_normal"), row.names = FALSE)
cat(sprintf("Archivo 'dge' guardado en: %s\n", output_dir))

write.csv(metadatos, file.path(output_dir, "metadatos_sin_normal"), row.names = FALSE)
cat(sprintf("Archivo 'metadatos' guardado en: %s\n", output_dir))

cat("¡Listo! El subtipo 'normal' ya no forma parte de este estudio.")
