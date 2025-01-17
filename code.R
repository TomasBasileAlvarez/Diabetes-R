# PROYECTO DE MAESTRÍA - MARTÍN ROBERTO BASILE ALVAREZ








#### Cargar paquetes ####
require(pacman)
pacman::p_load(haven, tidyverse, ggpubr, lmtest, nortest, gtools, data.table, caret, glmnet, survival, flextable, blandr, BlandAltmanLeh, corrplot, ggplot2,
               rms, bestNormalize, flexsurv, pROC, timeROC, fmsb, factoextra, gridExtra,  nhanesA, wesanderson,forestmodel, ggedit,dummy,lmtest,kableExtra, 
               FactoMineR, fpc, NbClust, skimr, ggimage, glmnet, gtsummary,  ggsci, survminer, cluster, dplyr, ggplotify, UpSetR, nortest, viridis, officer, survival,
               magrittr, tidytidbits, Epi, gt, lubridate, glue, ggtext, fmsb, ggplotify, tidyr, openxlsx)


dir = "martin/"
####---- 0. Dataset magament ----####
library(readxl)
library(openxlsx)

CAIPADI_MASTER <- read_xlsx(paste0(dir,"1. Base CAIPaDi 2013-2019 master.xlsx"), guess_max = 3200)

####---- 1. Nombres de las columnas (Diccionario) ----####
# Crear un data frame más simple con información de las columnas
column_info <- data.frame(
  Variable = names(CAIPADI_MASTER),
  Tipo = unlist(lapply(CAIPADI_MASTER, function(x) {
    if (is.numeric(x)) {
      return("numérica")
    } else if (is.character(x)) {
      return("texto")
    } else if (inherits(x, "Date") || inherits(x, "POSIXct")) {
      return("fecha")
    } else if (is.logical(x)) {
      return("lógica")
    } else {
      return("otro")
    }
  }))
)

# Mostrar las primeras 10 filas del diccionario
print(head(column_info, 10))



####---- 1.1 Imputar ----####



  
CAIPADI_MASTER <- CAIPADI_MASTER %>%
      group_by(CLASE_V1, SEXO) %>%
      mutate(across(all_of("NUT_V1_Circunferencia de cintura"), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))






####---- 1.2 Tomar Datos para meter al Non insulin based model ----####
Vi = "V7" #Elegir cualquier numero de visita para genera la tabla que luego puede ser ingresada en la calculadora de omar.
# Sin embargo, las visitas 2 y 3 no tienen circunferencia de cintura. 


CAIPADI_Vi_MODEL <- CAIPADI_MASTER[, c(
paste("NUT_",Vi,"_IMC",sep = ""),
"EDAD", 
paste("ENF_",Vi,"_HbA1c",sep = ""), 
paste("ENF_",Vi,"_Glucosa",sep=""),
paste("ENF_",Vi,"_Colesterol HDL", sep = ""),
paste("ENF_",Vi,"_Triglicê©ridos",sep=""),
paste("NUT_",Vi,"_Talla m",sep = ""),
paste("NUT_",Vi,"_Circunferencia de cintura",sep = ""),
"SEXO","ANIOS_DX_A_PRIMERA_VISITA")]


CAIPADI_Vi_MODEL <- CAIPADI_Vi_MODEL %>%
  rename("bmi3" = paste("NUT_",Vi,"_IMC",sep = ""),
         "age4" = "ANIOS_DX_A_PRIMERA_VISITA",
         "hba1c3" = paste("ENF_",Vi,"_HbA1c",sep = ""),
         "glucose1" = paste("ENF_",Vi,"_Glucosa",sep=""),
         "tg" = paste("ENF_",Vi,"_Triglicê©ridos",sep=""),
         "hdl" =paste("ENF_",Vi,"_Colesterol HDL", sep = ""),
         "age" = "EDAD",
         "waist" = paste("NUT_",Vi,"_Circunferencia de cintura",sep = ""),
         "height" = paste("NUT_",Vi,"_Talla m",sep = ""),
         "sex_bin" = "SEXO")



CAIPADI_Vi_MODEL <- CAIPADI_Vi_MODEL %>%
  mutate(height = height * 100)

CAIPADI_Vi_MODEL <- CAIPADI_Vi_MODEL %>%
  mutate(sex = ifelse(sex_bin == 1, "male", "female"))

CAIPADI_Vi_MODEL <- CAIPADI_Vi_MODEL %>%
  mutate(age4 = age-age4 )



CAIPADI_Vi_MODEL <- CAIPADI_Vi_MODEL %>%
  select(bmi3,age4,hba1c3,glucose1,tg,hdl,age,waist,height,sex,sex_bin)

write.csv(CAIPADI_Vi_MODEL, paste(dir,"CAIPADI_",Vi,".csv",sep=""), row.names = FALSE)


####---- 1.3 Tomar las predicciones de la calculadora de omar (non insulin based model) y agregarlas a la tabla ----####
prediction_V1 <- read.csv(paste0(dir,"prediction1.csv")) %>%
  select("class","metsir","metsvf") %>%
  rename("class1" = "class","metsir1" = "metsir","metsvf1" = "metsvf" )

prediction_V4 <- read.csv(paste0(dir,"prediction4.csv") )%>%
  select("class","metsir","metsvf") %>%
  rename("class4" = "class","metsir4" = "metsir","metsvf4" = "metsvf" )

prediction_V5 <- read.csv(paste0(dir,"prediction5.csv")) %>%
  select("class","metsir","metsvf") %>%
  rename("class5" = "class","metsir5" = "metsir","metsvf5" = "metsvf" )

prediction_V6 <- read.csv(paste0(dir,"prediction6.csv")) %>%
  select("class","metsir","metsvf") %>%
  rename("class6" = "class","metsir6" = "metsir","metsvf6" = "metsvf" )

  
prediction_V7 <- read.csv(paste0(dir,"prediction7.csv") )%>%
      select("class","metsir","metsvf") %>%
  rename("class7" = "class","metsir7" = "metsir","metsvf7" = "metsvf" )


CAIPADI_MASTER <- bind_cols(CAIPADI_MASTER,prediction_V1,prediction_V4,prediction_V5,prediction_V6,prediction_V7)

####---- 1.4 Agregar la tasa de filtrado glomerular (TGFe) usando la formula ----####
CAIPADI_MASTER <- CAIPADI_MASTER %>%
  mutate(TGFe_V1 = ifelse(SEXO == 1, 142*pmin(ENF_V1_Creatinina/0.9,1)^(-0.302)*pmax(ENF_V1_Creatinina/0.9,1)^(-1.2)*0.9938^EDAD, 
                          142*pmin(ENF_V1_Creatinina/0.7,1)^(-0.241)*pmax(ENF_V1_Creatinina/0.7,1)^(-1.2)*0.9938^EDAD*1.012))
CAIPADI_MASTER <- CAIPADI_MASTER %>%
  mutate(TGFe_V5 = ifelse(SEXO == 1, 142*pmin(ENF_V5_Creatinina/0.9,1)^(-0.302)*pmax(ENF_V5_Creatinina/0.9,1)^(-1.2)*0.9938^EDAD, 
                          142*pmin(ENF_V5_Creatinina/0.7,1)^(-0.241)*pmax(ENF_V5_Creatinina/0.7,1)^(-1.2)*0.9938^EDAD*1.012))
CAIPADI_MASTER <- CAIPADI_MASTER %>%
  mutate(TGFe_V6 = ifelse(SEXO == 1, 142*pmin(ENF_V6_Creatinina/0.9,1)^(-0.302)*pmax(ENF_V6_Creatinina/0.9,1)^(-1.2)*0.9938^EDAD, 
                          142*pmin(ENF_V6_Creatinina/0.7,1)^(-0.241)*pmax(ENF_V6_Creatinina/0.7,1)^(-1.2)*0.9938^EDAD*1.012))
CAIPADI_MASTER <- CAIPADI_MASTER %>%
  mutate(TGFe_V7 = ifelse(SEXO == 1, 142*pmin(ENF_V7_Creatinina/0.9,1)^(-0.302)*pmax(ENF_V7_Creatinina/0.9,1)^(-1.2)*0.9938^EDAD, 
                          142*pmin(ENF_V7_Creatinina/0.7,1)^(-0.241)*pmax(ENF_V7_Creatinina/0.7,1)^(-1.2)*0.9938^EDAD*1.012))

####---- 2. Datos generales de los subgrupos de diabetes (V1) ----####
# Contar el número de personas por subgrupo
subgroup_counts <- table(CAIPADI_MASTER$CLASE_V1)

group_Hba1c <-grep(paste0("^", "ENF_V[0-9]+_HbA1c", "$"), names(CAIPADI_MASTER), value = TRUE)
group_colesterol_LDL <-grep(paste0("^", "ENF_V[0-9]+_Colesterol LDL", "$"), names(CAIPADI_MASTER), value = TRUE)
group_TA_sistelica <-grep(paste0("^", "ENF_V[0-9]+_TA sistê_lica", "$"), names(CAIPADI_MASTER), value = TRUE)
group_peso <-grep(paste0("^", "NUT_V[0-9]+_Peso", "$"), names(CAIPADI_MASTER), value = TRUE)
group_cintura <-grep(paste0("^", "NUT_V[0-9]+_Circunferencia de cintura", "$"), names(CAIPADI_MASTER), value = TRUE)
group_colesterol_total <-grep(paste0("^", "ENF_V[0-9]+_Colesterol total", "$"), names(CAIPADI_MASTER), value = TRUE)
group_colesterol_hdl <-grep(paste0("^", "ENF_V[0-9]+_Colesterol HDL", "$"), names(CAIPADI_MASTER), value = TRUE)
group_ast <-grep(paste0("^", "ENF_V[0-9]+_AST", "$"), names(CAIPADI_MASTER), value = TRUE)
group_alt <-grep(paste0("^", "ENF_V[0-9]+_ALT", "$"), names(CAIPADI_MASTER), value = TRUE)
group_ggt <-grep(paste0("^", "ENF_V[0-9]+_GGT", "$"), names(CAIPADI_MASTER), value = TRUE)
group_alb_cre <-grep(paste0("^", "ENF_V[0-9]+_Indice de albuminuria/creatinuria", "$"), names(CAIPADI_MASTER), value = TRUE)
group_creatinina <-grep(paste0("^", "ENF_V[0-9]+_Creatinina", "$"), names(CAIPADI_MASTER), value = TRUE)
group_TGFe <-grep(paste0("^", "TGFe_V[0-9]", "$"), names(CAIPADI_MASTER), value = TRUE)
group_imc <-grep(paste0("^", "NUT_V[0-9]+_IMC", "$"), names(CAIPADI_MASTER), value = TRUE)
group_perdida5 <-grep(paste0("^", "NUT_V[0-9]+_Perdida de peso del 5%", "$"), names(CAIPADI_MASTER), value = TRUE)
group_retino_izq <-grep(paste0("^", "OFT_V[0-9]+_ClasificacionRetinopatia ojo izquierdo", "$"), names(CAIPADI_MASTER), value = TRUE)
group_retino_der <-grep(paste0("^", "OFT_V[0-9]+_ClasificacionRetinopatia ojo derecho", "$"), names(CAIPADI_MASTER), value = TRUE)
group_fotoc_izq <-grep(paste0("^", "OFT_V[0-9]+_Fotocuagulacion ojo izquierdo", "$"), names(CAIPADI_MASTER), value = TRUE)
group_fotoc_der <-grep(paste0("^", "OFT_V[0-9]+_Fotocuagulacion ojo derecho", "$"), names(CAIPADI_MASTER), value = TRUE)
group_EM_izq <-grep(paste0("^", "OFT_V[0-9]+_ClasificacionEM ojo izquierdo", "$"), names(CAIPADI_MASTER), value = TRUE)
group_EM_der <-grep(paste0("^", "OFT_V[0-9]+_ClasificacionEM ojo derecho", "$"), names(CAIPADI_MASTER), value = TRUE)
group_glaucoma_izq <-grep(paste0("^", "OFT_V[0-9]+_Glaucoma ojo izquierdo", "$"), names(CAIPADI_MASTER), value = TRUE)
group_glaucoma_der <-grep(paste0("^", "OFT_V[0-9]+_Glaucoma ojo derecho", "$"), names(CAIPADI_MASTER), value = TRUE)
group_catarata_izq <-grep(paste0("^", "OFT_V[0-9]+_Catarata ojo izquierdo", "$"), names(CAIPADI_MASTER), value = TRUE)
group_catarata_der <-grep(paste0("^", "OFT_V[0-9]+_Catarata ojo derecho", "$"), names(CAIPADI_MASTER), value = TRUE)
group_metsir <-grep(paste0("^", "metsir[0-9]", "$"), names(CAIPADI_MASTER), value = TRUE)
group_metsvf <-grep(paste0("^", "metsvf[0-9]", "$"), names(CAIPADI_MASTER), value = TRUE)



all_columns = c("CLASE_V1","SEXO", "EDAD","ANIOS_DX_A_PRIMERA_VISITA",group_Hba1c,group_colesterol_LDL,group_TA_sistelica,
                group_peso,group_cintura,group_colesterol_total,group_colesterol_hdl,
                group_ast,group_alt,group_ggt,group_alb_cre,group_creatinina,
                group_TGFe,group_imc,group_perdida5,group_retino_izq,group_retino_der,
                group_fotoc_izq,group_fotoc_der,group_EM_izq,group_EM_der,
                group_glaucoma_izq,group_glaucoma_der, group_catarata_izq,group_catarata_der,
                group_metsir,group_metsvf)



# Genera tabla pero unicamente con las variables deseadas
individual_data <- CAIPADI_MASTER[all_columns]
individual_data <- individual_data %>%
  mutate(SEXO = ifelse(SEXO == 1, "Hombre", "Mujer"))


individual_data <- individual_data %>%
  rename("Edad" = "EDAD",
         "HbA1c Visita 1" = "ENF_V1_HbA1c",
         "HbA1c Visita 2" = "ENF_V2_HbA1c",
         "HbA1c Visita 3" = "ENF_V3_HbA1c",
         "HbA1c Visita 4" = "ENF_V4_HbA1c",
         "HbA1c Visita 5" = "ENF_V5_HbA1c",
         "HbA1c Visita 6" = "ENF_V6_HbA1c",
         "HbA1c Visita 7" = "ENF_V7_HbA1c",
         "Colesterol LDL Visita 1" = "ENF_V1_Colesterol LDL",
         "Colesterol LDL Visita 2" = "ENF_V2_Colesterol LDL",
         "Colesterol LDL Visita 3" = "ENF_V3_Colesterol LDL",
         "Colesterol LDL Visita 4" = "ENF_V4_Colesterol LDL",
         "Colesterol LDL Visita 5" = "ENF_V5_Colesterol LDL",
         "Colesterol LDL Visita 6" = "ENF_V6_Colesterol LDL",
         "Colesterol LDL Visita 7" = "ENF_V7_Colesterol LDL",
         "TA Sistólica Visita 1" = "ENF_V1_TA sistê_lica",
         "TA Sistólica Visita 2" = "ENF_V2_TA sistê_lica",
         "TA Sistólica Visita 3" = "ENF_V3_TA sistê_lica",
         "TA Sistólica Visita 4" = "ENF_V4_TA sistê_lica",
         "TA Sistólica Visita 5" = "ENF_V5_TA sistê_lica",
         "TA Sistólica Visita 6" = "ENF_V6_TA sistê_lica",
         "TA Sistólica Visita 7" = "ENF_V7_TA sistê_lica",
         "Peso Visita 1" = "NUT_V1_Peso",
         "Peso Visita 2" = "NUT_V2_Peso",
         "Peso Visita 3" = "NUT_V3_Peso",
         "Peso Visita 4" = "NUT_V4_Peso",
         "Peso Visita 5" = "NUT_V5_Peso",
         "Peso Visita 6" = "NUT_V6_Peso",
         "Peso Visita 7" = "NUT_V7_Peso",
         "Circunferencia de la cintura visita 1" = "NUT_V1_Circunferencia de cintura",
         "Circunferencia de la cintura visita 4" = "NUT_V4_Circunferencia de cintura",
         "Circunferencia de la cintura visita 5" = "NUT_V5_Circunferencia de cintura",
         "Circunferencia de la cintura visita 6" = "NUT_V6_Circunferencia de cintura",
         "Circunferencia de la cintura visita 7" = "NUT_V7_Circunferencia de cintura",
         "Colesterol Total Visita 1" = "ENF_V1_Colesterol total",
         "Colesterol Total Visita 2" = "ENF_V2_Colesterol total",
         "Colesterol Total Visita 3" = "ENF_V3_Colesterol total",
         "Colesterol Total Visita 4" = "ENF_V4_Colesterol total",
         "Colesterol Total Visita 5" = "ENF_V5_Colesterol total",
         "Colesterol Total Visita 6" = "ENF_V6_Colesterol total",
         "Colesterol Total Visita 7" = "ENF_V7_Colesterol total",
         "Colesterol HDL Visita 1" = "ENF_V1_Colesterol HDL",
         "Colesterol HDL Visita 2" = "ENF_V2_Colesterol HDL",
         "Colesterol HDL Visita 3" = "ENF_V3_Colesterol HDL",
         "Colesterol HDL Visita 4" = "ENF_V4_Colesterol HDL",
         "Colesterol HDL Visita 5" = "ENF_V5_Colesterol HDL",
         "Colesterol HDL Visita 6" = "ENF_V6_Colesterol HDL",
         "Colesterol HDL Visita 7" = "ENF_V7_Colesterol HDL",
         "AST Visita 1" = "ENF_V1_AST",
         "AST Visita 5" = "ENF_V5_AST",
         "AST Visita 6" = "ENF_V6_AST",
         "AST Visita 7" = "ENF_V7_AST",
         "ALT Visita 1" = "ENF_V1_ALT",
         "ALT Visita 5" = "ENF_V5_ALT",
         "ALT Visita 6" = "ENF_V6_ALT",
         "ALT Visita 7" = "ENF_V7_ALT",
         "GGT Visita 1" = "ENF_V1_GGT",
         "GGT Visita 5" = "ENF_V5_GGT",
         "GGT Visita 6" = "ENF_V6_GGT",
         "GGT Visita 7" = "ENF_V7_GGT",
         "Indice de albumuminuria-creatinuria Visita 1" =  "ENF_V1_Indice de albuminuria/creatinuria",
         "Indice de albumuminuria-creatinuria Visita 3" =  "ENF_V3_Indice de albuminuria/creatinuria",
         "Indice de albumuminuria-creatinuria Visita 4" =  "ENF_V4_Indice de albuminuria/creatinuria",
         "Indice de albumuminuria-creatinuria Visita 5" =  "ENF_V5_Indice de albuminuria/creatinuria",
         "Indice de albumuminuria-creatinuria Visita 6" =  "ENF_V6_Indice de albuminuria/creatinuria",
         "Indice de albumuminuria-creatinuria Visita 7" =  "ENF_V7_Indice de albuminuria/creatinuria",
         "Creatinina Visita 1" = "ENF_V1_Creatinina",
         "Creatinina Visita 5" = "ENF_V5_Creatinina",
         "Creatinina Visita 6" = "ENF_V6_Creatinina",
         "Creatinina Visita 7" = "ENF_V7_Creatinina",
         "TGFe Visita 1" = "TGFe_V1",
         "TGFe Visita 5" = "TGFe_V5",
         "TGFe Visita 6" = "TGFe_V6",
         "TGFe Visita 7" = "TGFe_V7",
         "IMC Visita 1" = "NUT_V1_IMC",
         "IMC Visita 2" = "NUT_V2_IMC",
         "IMC Visita 3" = "NUT_V3_IMC",
         "IMC Visita 4" = "NUT_V4_IMC",
         "IMC Visita 5" = "NUT_V5_IMC",
         "IMC Visita 6" = "NUT_V6_IMC",
         "IMC Visita 7" = "NUT_V7_IMC",
         "Pérdida de peso del cinco por ciento Visita 1" = "NUT_V1_Perdida de peso del 5%",
         "Pérdida de peso del cinco por ciento Visita 4" = "NUT_V4_Perdida de peso del 5%",
         "Pérdida de peso del cinco por ciento Visita 5" = "NUT_V5_Perdida de peso del 5%",
         "Pérdida de peso del cinco por ciento Visita 6" = "NUT_V6_Perdida de peso del 5%",
         "Pérdida de peso del cinco por ciento Visita 7" = "NUT_V7_Perdida de peso del 5%",
         "Clasificación Retinoparía ojo izquierdo Visita 2" = "OFT_V2_ClasificacionRetinopatia ojo izquierdo",
         "Clasificación Retinoparía ojo izquierdo Visita 5" = "OFT_V5_ClasificacionRetinopatia ojo izquierdo",
         "Clasificación Retinoparía ojo izquierdo Visita 6" = "OFT_V6_ClasificacionRetinopatia ojo izquierdo",
         "Clasificación Retinoparía ojo izquierdo Visita 7" = "OFT_V7_ClasificacionRetinopatia ojo izquierdo",
         "Clasificación Retinoparía ojo derecho Visita 2" = "OFT_V2_ClasificacionRetinopatia ojo derecho",
         "Clasificación Retinoparía ojo derecho Visita 5" = "OFT_V5_ClasificacionRetinopatia ojo derecho",
         "Clasificación Retinoparía ojo derecho Visita 6" = "OFT_V6_ClasificacionRetinopatia ojo derecho",
         "Clasificación Retinoparía ojo derecho Visita 7" = "OFT_V7_ClasificacionRetinopatia ojo derecho",
         "Fotocuagulación ojo izquierdo Visita 2" = "OFT_V2_Fotocuagulacion ojo izquierdo",
         "Fotocuagulación ojo izquierdo Visita 5" = "OFT_V5_Fotocuagulacion ojo izquierdo",
         "Fotocuagulación ojo izquierdo Visita 6" = "OFT_V6_Fotocuagulacion ojo izquierdo",
         "Fotocuagulación ojo izquierdo Visita 7" = "OFT_V7_Fotocuagulacion ojo izquierdo",
         "Fotocuagulación ojo derecho Visita 2" = "OFT_V2_Fotocuagulacion ojo derecho",
         "Fotocuagulación ojo derecho Visita 5" = "OFT_V5_Fotocuagulacion ojo derecho",
         "Fotocuagulación ojo derecho Visita 6" = "OFT_V6_Fotocuagulacion ojo derecho",
         "Fotocuagulación ojo derecho Visita 7" = "OFT_V7_Fotocuagulacion ojo derecho",
         "Clasificación NEM ojo izquierdo Visita 2" = "OFT_V2_ClasificacionEM ojo izquierdo",
         "Clasificación NEM ojo izquierdo Visita 4" = "OFT_V4_ClasificacionEM ojo izquierdo",
         "Clasificación NEM ojo izquierdo Visita 5" = "OFT_V5_ClasificacionEM ojo izquierdo",
         "Clasificación NEM ojo izquierdo Visita 6" = "OFT_V6_ClasificacionEM ojo izquierdo",
         "Clasificación NEM ojo izquierdo Visita 7" = "OFT_V7_ClasificacionEM ojo izquierdo",
         "Clasificación NEM ojo derecho Visita 2" = "OFT_V2_ClasificacionEM ojo derecho",
         "Clasificación NEM ojo derecho Visita 4" = "OFT_V4_ClasificacionEM ojo derecho",
         "Clasificación NEM ojo derecho Visita 5" = "OFT_V5_ClasificacionEM ojo derecho",
         "Clasificación NEM ojo derecho Visita 6" = "OFT_V6_ClasificacionEM ojo derecho",
         "Clasificación NEM ojo derecho Visita 7" = "OFT_V7_ClasificacionEM ojo derecho",
         "Glaucoma ojo izquierdo Visita 3" = "OFT_V3_Glaucoma ojo izquierdo",
         "Glaucoma ojo izquierdo Visita 5" = "OFT_V5_Glaucoma ojo izquierdo",
         "Glaucoma ojo izquierdo Visita 6" = "OFT_V6_Glaucoma ojo izquierdo",
         "Glaucoma ojo izquierdo Visita 7" = "OFT_V7_Glaucoma ojo izquierdo",
         "Glaucoma ojo derecho Visita 3" = "OFT_V3_Glaucoma ojo derecho",
         "Glaucoma ojo derecho Visita 5" = "OFT_V5_Glaucoma ojo derecho",
         "Glaucoma ojo derecho Visita 6" = "OFT_V6_Glaucoma ojo derecho",
         "Glaucoma ojo derecho Visita 7" = "OFT_V7_Glaucoma ojo derecho",
         "Catarata ojo izquierdo Visita 3" = "OFT_V3_Catarata ojo izquierdo",
         "Catarata ojo izquierdo Visita 5" = "OFT_V5_Catarata ojo izquierdo",
         "Catarata ojo izquierdo Visita 6" = "OFT_V6_Catarata ojo izquierdo",
         "Catarata ojo izquierdo Visita 7" = "OFT_V7_Catarata ojo izquierdo",
         "Catarata ojo derecho Visita 3" = "OFT_V3_Catarata ojo derecho",
         "Catarata ojo derecho Visita 5" = "OFT_V5_Catarata ojo derecho",
         "Catarata ojo derecho Visita 6" = "OFT_V6_Catarata ojo derecho",
         "Catarata ojo derecho Visita 7" = "OFT_V7_Catarata ojo derecho",
         "Metsir Visita 1" = "metsir1",
         "Metsir Visita 4" = "metsir4",
         "Metsir Visita 5" = "metsir5",
         "Metsir Visita 6" = "metsir6",
         "Metsir Visita 7" = "metsir7",
         "Metsvf Visita 1" = "metsvf1",
         "Metsvf Visita 4" = "metsvf4",
         "Metsvf Visita 5" = "metsvf5",
         "Metsvf Visita 6" = "metsvf6",
         "Metsvf Visita 7" = "metsvf7"
         )


#Generar una tabla descriptiva con promedio y SS de las varaibles anteriores 
Tabla_descriptiva_1 <- individual_data %>%
  group_by(CLASE_V1, SEXO) %>%
  summarise(across( .cols = everything(),
    .fns = list(promedio = ~mean(. , na.rm = TRUE), 
                               std = ~sd(. , na.rm = TRUE)),
                   .names = "{col}_{fn}"))



# LISTA de las posibles variables son las columnas de individual_data
all_variables <- data.frame(
  Variable = names(individual_data),
  Tipo = unlist(lapply(individual_data, function(x) {
    if (is.numeric(x)) {
      return("numérica")
    } else if (is.character(x)) {
      return("texto")
    } else if (inherits(x, "Date") || inherits(x, "POSIXct")) {
      return("fecha")
    } else if (is.logical(x)) {
      return("lógica")
    } else {
      return("otro")
    }
  }))
)



####---- 2.1 Distribuciones de subgrupos para cualquier variable----####


# Elegir una de las variables de entre "all_variables" (la lista de las variables a elegir)

#Elegir una y poner el nombre de la varibale 
variable <- "HbA1c Visita 7"

ggplot(individual_data, aes(x = CLASE_V1, y = !!sym(variable), fill = SEXO)) +
  # Añadir cajas y bigotes
  geom_boxplot(
    position = position_dodge(width = 0.8), 
               alpha = 0.7, 
               width = 0.7) +
  # Personalizar colores
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  # Etiquetas y título
  labs(
    title = paste(variable, "por Subgrupo de Diabetes y Sexo"),
    subtitle = "CAIPaDi 2013-2019",
    x = "Subgrupo de Diabetes",
    y =  variable,
    fill = "Género",
    color = "Género"
  ) +
  
  #DESCOMENTAR ESTA LINEA PARA CAMBIARL LOS LIMITES en Y
  # POR ejemplo, en RACU que el promedio es bajito pero algunos tienen màs de 4 mil de RACU
 
  ##scale_y_continuous(limits = c(0, 30))

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )



####---- 2.2 Distribuciones de subgrupos para TODAS las variables----####

library(rlang)

# Creando las columnas llamadas VARIABLES
variables <- setdiff(colnames(individual_data),c("CLASE_V1","SEXO"))

# El código de 2.1 pero generado a todas las variables y después guardarlos en una carpeta
plots <- lapply(variables, function(variable) {
  ggplot(individual_data, aes(x = CLASE_V1, y = !!sym(variable), fill = SEXO)) +
    geom_boxplot(position = position_dodge(width = 0.8), 
                 alpha = 0.7, 
                 width = 0.7) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
    labs(
      title = paste(variable, "por Subgrupo de Diabetes y Sexo"),
      subtitle = "CAIPaDi 2013-2019",
      x = "Subgrupo de Diabetes",
      y = variable,
      fill = "Género",
      color = "Género"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.background = element_rect(fill = "white", color = NA),  # QUE EL FONDO SEA BLANCO :) 
      panel.background = element_rect(fill = "white", color = NA)  
    )
})

# GUARDARLO
for (i in seq_along(variables)) {
  ggsave(filename = paste0(dir,"figuras_distribuciones/",variables[i], ".png"), plot = plots[[i]], width = 10, height = 6, dpi = 300)
}





####---- 3. Grafica de una variable como funcion del numero de visita ----####



# Cambiar para variable desdeada para graficar

library(dplyr)
variable <- "HbA1c"

variable_data <- individual_data %>%
  select(matches(variable))


variable_data_long <- variable_data %>%
  pivot_longer(cols = starts_with(paste0(variable, " Visita")), names_to = "Visita", values_to = variable) %>%
  mutate(Visita = as.numeric(str_extract(Visita, " \\d+")))

variable_data_mean <- variable_data_long %>%
  group_by(CLASE_V1, SEXO, Visita) %>%
  # Con cuartiles 
  summarise(mean = mean(!!sym(variable), na.rm=TRUE),
            lower_limit = quantile(!!sym(variable),0.25,na.rm = TRUE),
            upper_limit = quantile(!!sym(variable),0.75, na.rm=TRUE),
             .groups = "drop")

library(tidyverse)

custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

variable_data_mean <- variable_data_mean %>%
  mutate(group= interaction(CLASE_V1, SEXO, sep = "_")) %>%
  mutate(group = factor(group, levels = names(custom_colors) ))
         




ggplot(variable_data_mean, aes(x = Visita, y = mean, color = group, group = group)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +  # Thicker lines with transparency
  geom_point(size = 3, alpha = 0.8) +   # Bigger points with transparency
  geom_pointrange(aes(ymin = lower_limit, ymax = upper_limit), 
                  size = 1.3, alpha = 0.8,
                  linewidth=1.7) +  # Thicker point range lines
  scale_color_manual(values = custom_colors) + 
  scale_x_continuous(breaks = 1:7) +  
  labs(x = "Número de visita", y = paste0("Promedio de ",variable), title = paste0("Promedio de ", variable, " por visita"), color = "Clase_V1 y SEXO") +
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = "white", color = "white"),  # White panel background
    plot.background = element_rect(fill = "white", color = "white"),   # White plot background
    panel.grid.major = element_line(color = "grey80"),  # Light grid lines
    panel.grid.minor = element_line(color = "grey90")   # Lighter minor grid lines
  )



####---- 3.1 Hacer lo anterior para TODAS las variables ----####
variables <- c("HbA1c","Colesterol LDL","TA Sistólica","Peso","Circunferencia de la cintura",
             "Colesterol Total","Colesterol HDL","AST","ALT","GGT",
             "Indice de albumuminuria-creatinuria","Creatinina","TGFe","IMC",
             "Pérdida de peso del cinco por ciento","Clasificación Retinoparía ojo izquierdo",
             "Clasificación Retinoparía ojo derecho","Fotocuagulación ojo izquierdo",
             "Fotocuagulación ojo derecho","Clasificación NEM ojo izquierdo","Clasificación NEM ojo derecho",
             "Glaucoma ojo izquierdo","Glaucoma ojo derecho","Catarata ojo izquierdo","Catarata ojo derecho",
             "Metsir","Metsvf")

library(tidyverse)

custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

for (variable in variables) {
  
  
  
  variable_data <- individual_data %>%
    select(matches(variable))
  
  
  variable_data_long <- variable_data %>%
    pivot_longer(cols = starts_with(paste0(variable, " Visita")), names_to = "Visita", values_to = variable) %>%
    mutate(Visita = as.numeric(str_extract(Visita, " \\d+")))
  
  variable_data_mean <- variable_data_long %>%
    group_by(CLASE_V1, SEXO, Visita) %>%
    # Con cuartiles 
    summarise(mean = mean(!!sym(variable), na.rm=TRUE),
              lower_limit = quantile(!!sym(variable),0.25,na.rm = TRUE),
              upper_limit = quantile(!!sym(variable),0.75, na.rm=TRUE),
              .groups = "drop")
  
  
  variable_data_mean <- variable_data_mean %>%
    mutate(group= interaction(CLASE_V1, SEXO, sep = "_")) %>%
    mutate(group = factor(group, levels = names(custom_colors) ))
  
  
  
  pl<-ggplot(variable_data_mean, aes(x = Visita, y = mean, color = group, group = group)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +  # Thicker lines with transparency
    geom_point(size = 3, alpha = 0.8) +   # Bigger points with transparency
    geom_pointrange(aes(ymin = lower_limit, ymax = upper_limit), 
                    size = 1.3, alpha = 0.8,
                    linewidth=1.7) +  # Thicker point range lines
    scale_color_manual(values = custom_colors) + 
    scale_x_continuous(breaks = 1:7) +  
    labs(x = "Número de visita", y = paste0("Promedio de ",variable), title = paste0("Promedio de ", variable, " por visita"), color = "Clase_V1 y SEXO") +
    theme_minimal()+
    theme(
      panel.background = element_rect(fill = "white", color = "white"),  # White panel background
      plot.background = element_rect(fill = "white", color = "white"),   # White plot background
      panel.grid.major = element_line(color = "grey80"),  # Light grid lines
      panel.grid.minor = element_line(color = "grey90")   # Lighter minor grid lines
    )
  
  ggsave(filename = paste0(dir,"figuras_visita/",variable,".png"), plot = pl, width = 10, height = 6, dpi = 300)
  
  print(paste("Saved plot for", variable))
  
  
}








####---- 4.1 Tabla de Medicamentos ----####

# El V0 significa el tratamiento previo (antes en el excel viene como TX previo, pero estaba muy mal nombrado)

# Aquí se está leyendo únicmaente en las columnas donde pusieron la dosis tomando en cuenta la familia del medicamento pero no tal cual cada medicamento individual (por ejemplo, IECA como familia y no tal cual enalapril)

tabla_medicamentos <- CAIPADI_MASTER %>%
  select(
"END_V0_Dosis Inhibidores SGLT2 indicado",
"END_V1_Dosis Inhibidores SGLT2 indicado",
"END_V2_Dosis Inhibidores SGLT2 indicado",
"END_V3_Dosis Inhibidores SGLT2 indicado",
"END_V4_Dosis Inhibidores SGLT2 indicado",
"END_V5_Dosis Inhibidores SGLT2 indicado",
"END_V6_Dosis Inhibidores SGLT2 indicado",
"END_V7_Dosis Inhibidores SGLT2 indicado",
"END_V0_Dosis Sulfonilureas indicado",
"END_V1_Dosis Sulfonilureas indicado",
"END_V2_Dosis Sulfonilureas indicado",
"END_V3_Dosis Sulfonilureas indicado",
"END_V4_Dosis Sulfonilureas indicado",
"END_V5_Dosis Sulfonilureas indicado",
"END_V6_Dosis Sulfonilureas indicado",
"END_V7_Dosis Sulfonilureas indicado",
"END_V0_Dosis Anêçlogos de GLP-1 indicado",
"END_V1_Dosis Anêçlogos de GLP-1 indicado",
"END_V2_Dosis Anêçlogos de GLP-1 indicado",
"END_V3_Dosis Anêçlogos de GLP-1 indicado",
"END_V4_Dosis Anêçlogos de GLP-1 indicado",
"END_V5_Dosis Anêçlogos de GLP-1 indicado",
"END_V6_Dosis Anêçlogos de GLP-1 indicado",
"END_V6_Dosis Anêçlogos de GLP-1 indicado",
"END_V7_Dosis Anêçlogos de GLP-1 indicado",
"END_V0_Dosis Inhibidores de DPP-IV Indicado",
"END_V1_Dosis Inhibidores de DPP-IV Indicado",
"END_V2_Dosis Inhibidores de DPP-IV Indicado",
"END_V3_Dosis Inhibidores de DPP-IV Indicado",
"END_V4_Dosis Inhibidores de DPP-IV Indicado",
"END_V5_Dosis Inhibidores de DPP-IV Indicado",
"END_V6_Dosis Inhibidores de DPP-IV Indicado",
"END_V7_Dosis Inhibidores de DPP-IV indicado",

"END_V0_Dosis Insulina basal desayuno indicado",
"END_V0_Dosis Insulina basal comida indicado",
"END_V0_Dosis Insulina basal cena indicado",
"END_V1_Dosis Insulina basal desayuno indicado",
"END_V1_Dosis Insulina basal comida indicado",
"END_V1_Dosis Insulina basal cena indicado",
"END_V2_Dosis Insulina basal desayuno indicado",
"END_V2_Dosis Insulina basal comida indicado",
"END_V2_Dosis Insulina basal cena indicado",
"END_V3_Dosis Insulina basal desayuno indicado",
"END_V3_Dosis Insulina basal comida indicado",
"END_V3_Dosis Insulina basal cena indicado",
"END_V4_Dosis Insulina basal desayuno indicado",
"END_V4_Dosis Insulina basal comida indicado",
"END_V4_Dosis Insulina basal cena indicado",
"END_V5_Dosis Insulina basal desayuno indicado",
"END_V5_Dosis Insulina basal comida indicado",
"END_V5_Dosis Insulina basal cena indicado",
"END_V6_Dosis Insulina basal desayuno indicado",
"END_V6_Dosis Insulina basal comida indicado",
"END_V6_Dosis Insulina basal cena indicado",
"END_V7_Dosis Insulina basal desayuno indicado",
"END_V7_Dosis Insulina basal comida indicado",
"END_V7_Dosis Insulina basal cena indicado",

"END_V0_Dosis desayuno en bolo indicado",
"END_V0_Dosis comida en bolo indicado",
"END_V0_Dosis cena en bolo indicado",
"END_V1_Dosis desayuno en bolo indicado",
"END_V1_Dosis comida en bolo indicado",
"END_V1_Dosis cena en bolo indicado",
"END_V2_Dosis desayuno en bolo indicado",
"END_V2_Dosis comida en bolo indicado",
"END_V2_Dosis cena en bolo indicado",
"END_V3_Dosis desayuno en bolo indicado",
"END_V3_Dosis comida en bolo indicado",
"END_V3_Dosis cena en bolo indicado",
"END_V4_Dosis desayuno en bolo indicado",
"END_V4_Dosis comida en bolo indicado",
"END_V4_Dosis cena en bolo indicado",
"END_V5_Dosis desayuno en bolo indicado",
"END_V5_Dosis comida en bolo indicado",
"END_V5_Dosis cena en bolo indicado",
"END_V6_Dosis desayuno en bolo indicado",
"END_V6_Dosis comida en bolo indicado",
"END_V6_Dosis cena en bolo indicado",
"END_V7_Dosis desayuno en bolo indicado",
"END_V7_Dosis comida en bolo indicado",
"END_V7_Dosis cena en bolo indicado",

"END_V0_Dosis desayuno premezcla indicado",
"END_V0_Dosis comida premezcla indicado",
"END_V0_Dosis cena premezcla indicado",
"END_V1_Dosis desayuno premezcla indicado",
"END_V1_Dosis comida premezcla indicado",
"END_V1_Dosis cena premezcla indicado",
"END_V2_Dosis desayuno premezcla indicado",
"END_V2_Dosis comida premezcla indicado",
"END_V2_Dosis cena premezcla indicado",
"END_V3_Dosis desayuno premezcla indicado",
"END_V3_Dosis comida premezcla indicado",
"END_V3_Dosis cena premezcla indicado",
"END_V4_Dosis desayuno premezcla indicado",
"END_V4_Dosis comida premezcla indicado",
"END_V4_Dosis cena premezcla indicado",
"END_V5_Dosis desayuno premezcla indicado",
"END_V5_Dosis comida premezcla indicado",
"END_V5_Dosis cena premezcla indicado",
"END_V6_Dosis desayuno premezcla indicado",
"END_V6_Dosis comida premezcla indicado",
"END_V6_Dosis cena premezcla indicado",
"END_V7_Dosis desayuno premezcla indicado",
"END_V7_Dosis comida premezcla indicado",
"END_V7_Dosis cena premezcla indicado",

"END_V0_Dosis IECA indicado",
"END_V1_Dosis IECA indicado",
"END_V2_Dosis IECA indicado",
"END_V3_Dosis IECA indicado",
"END_V4_Dosis IECA indicado",
"END_V5_Dosis IECA indicado",
"END_V6_Dosis IECA indicado",
"END_V7_Dosis IECA indicado",
"END_V0_Dosis Calcio antagonistas indicado",
"END_V1_Dosis Calcio antagonistas indicado",
"END_V2_Dosis Calcio antagonistas indicado",
"END_V3_Dosis Calcio antagonistas indicado",
"END_V4_Dosis Calcio antagonistas indicado",
"END_V5_Dosis Calcio antagonistas indicado",
"END_V6_Dosis Calcio antagonistas indicado",
"END_V7_Dosis Calcio antagonistas indicado",
"END_V0_Dosis ê_ bloqueadores indicado",
"END_V1_Dosis ê_ bloqueadores indicado",
"END_V2_Dosis ê_ bloqueadores indicado",
"END_V3_Dosis ê_ bloqueadores indicado",
"END_V4_Dosis ê_ bloqueadores indicado",
"END_V5_Dosis ê_ bloqueadores indicado",
"END_V6_Dosis ê_ bloqueadores indicado",
"END_V7_Dosis ê_ bloqueadores indicado",
"END_V0_Dosis ARA2 indicado",
"END_V1_Dosis ARA2 indicado",
"END_V2_Dosis ARA2 indicado",
"END_V3_Dosis ARA2 indicado",
"END_V4_Dosis ARA2 indicado",
"END_V5_Dosis ARA2 indicado",
"END_V6_Dosis ARA2 indicado",
"END_V7_Dosis ARA2 indicado",
"END_V0_Dosis diuerê©tico indicado",
"END_V1_Dosis diuerê©tico indicado",
"END_V2_Dosis diuerê©tico indicado",
"END_V3_Dosis diuerê©tico indicado",
"END_V4_Dosis diuerê©tico indicado",
"END_V5_Dosis diuerê©tico indicado",
"END_V6_Dosis diuerê©tico indicado",
"END_V7_Dosis diuerê©tico indicado",
"END_V0_Dosis estatinas indicado",
"END_V1_Dosis estatinas indicado",
"END_V2_Dosis estatinas indicado",
"END_V3_Dosis estatinas indicado",
"END_V4_Dosis estatinas indicado",
"END_V5_Dosis estatinas indicado",
"END_V6_Dosis estatinas indicado",
"END_V7_Dosis estatinas indicado",
"END_V0_Dosis fibrato indicado",
"END_V1_Dosis fibrato indicado",
"END_V2_Dosis fibrato indicado",
"END_V3_Dosis fibrato indicado",
"END_V4_Dosis fibrato indicado",
"END_V5_Dosis fibrato indicado",
"END_V6_Dosis fibrato indicado",
"END_V7_Dosis fibrato indicado",
"END_V0_Aspirina indicado",
"END_V1_Aspirina indicado",
"END_V2_Aspirina indicado",
"END_V3_Aspirina indicado",
"END_V4_Aspirina indicado",
"END_V5_Aspirina indicado",
"END_V6_Aspirina indicado",
"END_V7_Aspirina indicado",
"END_V0_Ezetimibe indicado",
"END_V1_Ezetimibe indicado",
"END_V2_Ezetimibe indicado",
"END_V3_Ezetimibe indicado",
"END_V4_Ezetimibe indicado",
"END_V5_Ezetimibe indicado",
"END_V6_Ezetimibe indicado",
"END_V7_Ezetimibe indicado")

tabla_medicamentos <- tabla_medicamentos %>%
  mutate("Insulina Basal Visita 0" = `END_V0_Dosis Insulina basal desayuno indicado` + 
           `END_V0_Dosis Insulina basal comida indicado` + 
           `END_V0_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V0_Dosis Insulina basal desayuno indicado`, 
         -`END_V0_Dosis Insulina basal comida indicado`, 
         -`END_V0_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 1" = `END_V1_Dosis Insulina basal desayuno indicado` + 
           `END_V1_Dosis Insulina basal comida indicado` + 
           `END_V1_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V1_Dosis Insulina basal desayuno indicado`, 
         -`END_V1_Dosis Insulina basal comida indicado`, 
         -`END_V1_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 2" = `END_V2_Dosis Insulina basal desayuno indicado` + 
           `END_V2_Dosis Insulina basal comida indicado` + 
           `END_V2_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V2_Dosis Insulina basal desayuno indicado`, 
         -`END_V2_Dosis Insulina basal comida indicado`, 
         -`END_V2_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 3" = `END_V3_Dosis Insulina basal desayuno indicado` + 
           `END_V3_Dosis Insulina basal comida indicado` + 
           `END_V3_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V3_Dosis Insulina basal desayuno indicado`, 
         -`END_V3_Dosis Insulina basal comida indicado`, 
         -`END_V3_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 4" = `END_V4_Dosis Insulina basal desayuno indicado` + 
           `END_V4_Dosis Insulina basal comida indicado` + 
           `END_V4_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V4_Dosis Insulina basal desayuno indicado`, 
         -`END_V4_Dosis Insulina basal comida indicado`, 
         -`END_V4_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 5" = `END_V5_Dosis Insulina basal desayuno indicado` + 
           `END_V5_Dosis Insulina basal comida indicado` + 
           `END_V5_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V5_Dosis Insulina basal desayuno indicado`, 
         -`END_V5_Dosis Insulina basal comida indicado`, 
         -`END_V5_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 6" = `END_V6_Dosis Insulina basal desayuno indicado` + 
           `END_V6_Dosis Insulina basal comida indicado` + 
           `END_V6_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V6_Dosis Insulina basal desayuno indicado`, 
         -`END_V6_Dosis Insulina basal comida indicado`, 
         -`END_V6_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina Basal Visita 7" = `END_V7_Dosis Insulina basal desayuno indicado` + 
           `END_V7_Dosis Insulina basal comida indicado` + 
           `END_V7_Dosis Insulina basal cena indicado`) %>%
  select(-`END_V7_Dosis Insulina basal desayuno indicado`, 
         -`END_V7_Dosis Insulina basal comida indicado`, 
         -`END_V7_Dosis Insulina basal cena indicado`) %>%
  mutate("Insulina en Bolo Visita 0" = `END_V0_Dosis desayuno en bolo indicado` + 
           `END_V0_Dosis comida en bolo indicado` + 
           `END_V0_Dosis cena en bolo indicado`) %>%
  select(-`END_V0_Dosis desayuno en bolo indicado`, 
         -`END_V0_Dosis comida en bolo indicado`, 
         -`END_V0_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 1" = `END_V1_Dosis desayuno en bolo indicado` + 
           `END_V1_Dosis comida en bolo indicado` + 
           `END_V1_Dosis cena en bolo indicado`) %>%
  select(-`END_V1_Dosis desayuno en bolo indicado`, 
         -`END_V1_Dosis comida en bolo indicado`, 
         -`END_V1_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 2" = `END_V2_Dosis desayuno en bolo indicado` + 
           `END_V2_Dosis comida en bolo indicado` + 
           `END_V2_Dosis cena en bolo indicado`) %>%
  select(-`END_V2_Dosis desayuno en bolo indicado`, 
         -`END_V2_Dosis comida en bolo indicado`, 
         -`END_V2_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 3" = `END_V3_Dosis desayuno en bolo indicado` + 
           `END_V3_Dosis comida en bolo indicado` + 
           `END_V3_Dosis cena en bolo indicado`) %>%
  select(-`END_V3_Dosis desayuno en bolo indicado`, 
         -`END_V3_Dosis comida en bolo indicado`, 
         -`END_V3_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 4" = `END_V4_Dosis desayuno en bolo indicado` + 
           `END_V4_Dosis comida en bolo indicado` + 
           `END_V4_Dosis cena en bolo indicado`) %>%
  select(-`END_V4_Dosis desayuno en bolo indicado`, 
         -`END_V4_Dosis comida en bolo indicado`, 
         -`END_V4_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 5" = `END_V5_Dosis desayuno en bolo indicado` + 
           `END_V5_Dosis comida en bolo indicado` + 
           `END_V5_Dosis cena en bolo indicado`) %>%
  select(-`END_V5_Dosis desayuno en bolo indicado`, 
         -`END_V5_Dosis comida en bolo indicado`, 
         -`END_V5_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 6" = `END_V6_Dosis desayuno en bolo indicado` + 
           `END_V6_Dosis comida en bolo indicado` + 
           `END_V6_Dosis cena en bolo indicado`) %>%
  select(-`END_V6_Dosis desayuno en bolo indicado`, 
         -`END_V6_Dosis comida en bolo indicado`, 
         -`END_V6_Dosis cena en bolo indicado`) %>%
  mutate("Insulina en Bolo Visita 7" = `END_V7_Dosis desayuno en bolo indicado` + 
           `END_V7_Dosis comida en bolo indicado` + 
           `END_V7_Dosis cena en bolo indicado`) %>%
  select(-`END_V7_Dosis desayuno en bolo indicado`, 
         -`END_V7_Dosis comida en bolo indicado`, 
         -`END_V7_Dosis cena en bolo indicado`) %>%
  mutate("Insulina Premezcla Visita 0" = `END_V0_Dosis desayuno premezcla indicado` + 
           `END_V0_Dosis comida premezcla indicado` + 
           `END_V0_Dosis cena premezcla indicado`) %>%
  select(-`END_V0_Dosis desayuno premezcla indicado`, 
         -`END_V0_Dosis comida premezcla indicado`, 
         -`END_V0_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 1" = `END_V1_Dosis desayuno premezcla indicado` + 
           `END_V1_Dosis comida premezcla indicado` + 
           `END_V1_Dosis cena premezcla indicado`) %>%
  select(-`END_V1_Dosis desayuno premezcla indicado`, 
         -`END_V1_Dosis comida premezcla indicado`, 
         -`END_V1_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 2" = `END_V2_Dosis desayuno premezcla indicado` + 
           `END_V2_Dosis comida premezcla indicado` + 
           `END_V2_Dosis cena premezcla indicado`) %>%
  select(-`END_V2_Dosis desayuno premezcla indicado`, 
         -`END_V2_Dosis comida premezcla indicado`, 
         -`END_V2_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 3" = `END_V3_Dosis desayuno premezcla indicado` + 
           `END_V3_Dosis comida premezcla indicado` + 
           `END_V3_Dosis cena premezcla indicado`) %>%
  select(-`END_V3_Dosis desayuno premezcla indicado`, 
         -`END_V3_Dosis comida premezcla indicado`, 
         -`END_V3_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 4" = `END_V4_Dosis desayuno premezcla indicado` + 
           `END_V4_Dosis comida premezcla indicado` + 
           `END_V4_Dosis cena premezcla indicado`) %>%
  select(-`END_V4_Dosis desayuno premezcla indicado`, 
         -`END_V4_Dosis comida premezcla indicado`, 
         -`END_V4_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 5" = `END_V5_Dosis desayuno premezcla indicado` + 
           `END_V5_Dosis comida premezcla indicado` + 
           `END_V5_Dosis cena premezcla indicado`) %>%
  select(-`END_V5_Dosis desayuno premezcla indicado`, 
         -`END_V5_Dosis comida premezcla indicado`, 
         -`END_V5_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 6" = `END_V6_Dosis desayuno premezcla indicado` + 
           `END_V6_Dosis comida premezcla indicado` + 
           `END_V6_Dosis cena premezcla indicado`) %>%
  select(-`END_V6_Dosis desayuno premezcla indicado`, 
         -`END_V6_Dosis comida premezcla indicado`, 
         -`END_V6_Dosis cena premezcla indicado`) %>%
  mutate("Insulina Premezcla Visita 7" = `END_V7_Dosis desayuno premezcla indicado` + 
           `END_V7_Dosis comida premezcla indicado` + 
           `END_V7_Dosis cena premezcla indicado`) %>%
  select(-`END_V7_Dosis desayuno premezcla indicado`, 
         -`END_V7_Dosis comida premezcla indicado`, 
         -`END_V7_Dosis cena premezcla indicado`)

  
  
tabla_medicamentos <- tabla_medicamentos %>%
  rename("SGLT2 Visita 0" = "END_V0_Dosis Inhibidores SGLT2 indicado",
    "SGLT2 Visita 1" = "END_V1_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 2" = "END_V2_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 3" = "END_V3_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 4" = "END_V4_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 5" = "END_V5_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 6" = "END_V6_Dosis Inhibidores SGLT2 indicado",
         "SGLT2 Visita 7" = "END_V7_Dosis Inhibidores SGLT2 indicado",
    "Sulfonilureas Visita 0" = "END_V0_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 1" = "END_V1_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 2" = "END_V2_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 3" = "END_V3_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 4" = "END_V4_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 5" = "END_V5_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 6" = "END_V6_Dosis Sulfonilureas indicado",
        "Sulfonilureas Visita 7" = "END_V7_Dosis Sulfonilureas indicado",
    "Análogos de GLP-1 Visita 0" = "END_V0_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 1" = "END_V1_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 2" = "END_V2_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 3" = "END_V3_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 4" = "END_V4_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 5" = "END_V5_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 6" = "END_V6_Dosis Anêçlogos de GLP-1 indicado",
       "Análogos de GLP-1 Visita 7" = "END_V7_Dosis Anêçlogos de GLP-1 indicado",
    "iDPP4 Visita 0" = "END_V0_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 1" = "END_V1_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 2" = "END_V2_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 3" = "END_V3_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 4" = "END_V4_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 5" = "END_V5_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 6" = "END_V6_Dosis Inhibidores de DPP-IV Indicado",
      "iDPP4 Visita 7" = "END_V7_Dosis Inhibidores de DPP-IV indicado",
    "IECA Visita 0" = "END_V0_Dosis IECA indicado",
    "IECA Visita 1" = "END_V1_Dosis IECA indicado",
    "IECA Visita 2" = "END_V2_Dosis IECA indicado",
    "IECA Visita 3" = "END_V3_Dosis IECA indicado",
    "IECA Visita 4" = "END_V4_Dosis IECA indicado",
    "IECA Visita 5" = "END_V5_Dosis IECA indicado",
    "IECA Visita 6" = "END_V6_Dosis IECA indicado",
    "IECA Visita 7" = "END_V7_Dosis IECA indicado",
    "Calcio Antagonista Visita 0" = "END_V0_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 1" = "END_V1_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 2" = "END_V2_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 3" = "END_V3_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 4" = "END_V4_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 5" = "END_V5_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 6" = "END_V6_Dosis Calcio antagonistas indicado",
  "Calcio Antagonista Visita 7" = "END_V7_Dosis Calcio antagonistas indicado",
  "Beta Bloqueador Visita 0" = "END_V0_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 1" = "END_V1_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 2" = "END_V2_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 3" = "END_V3_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 4" = "END_V4_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 5" = "END_V5_Dosis ê_ bloqueadores indicado",
   "Beta Bloqueador Visita 6" = "END_V6_Dosis ê_ bloqueadores indicado",
  "Beta Bloqueador Visita 7" = "END_V7_Dosis ê_ bloqueadores indicado",
  "ARA2 Visita 0" =  "END_V0_Dosis ARA2 indicado",
  "ARA2 Visita 1" =  "END_V1_Dosis ARA2 indicado",
  "ARA2 Visita 2" =  "END_V2_Dosis ARA2 indicado",
  "ARA2 Visita 3" =  "END_V3_Dosis ARA2 indicado",
  "ARA2 Visita 4" =  "END_V4_Dosis ARA2 indicado",
  "ARA2 Visita 5" =  "END_V5_Dosis ARA2 indicado",
  "ARA2 Visita 6" =  "END_V6_Dosis ARA2 indicado",
  "ARA2 Visita 7" =  "END_V7_Dosis ARA2 indicado",
  "Diurético Visita 0" = "END_V0_Dosis diuerê©tico indicado",
    "Diurético Visita 1" = "END_V1_Dosis diuerê©tico indicado",
    "Diurético Visita 2" = "END_V2_Dosis diuerê©tico indicado",
    "Diurético Visita 3" = "END_V3_Dosis diuerê©tico indicado",
    "Diurético Visita 4" = "END_V4_Dosis diuerê©tico indicado",
    "Diurético Visita 5" = "END_V5_Dosis diuerê©tico indicado",
    "Diurético Visita 6" = "END_V6_Dosis diuerê©tico indicado",
  "Diurético Visita 7" = "END_V7_Dosis diuerê©tico indicado",
  "Estatinas Visita 0" = "END_V0_Dosis estatinas indicado",
    "Estatinas Visita 1" = "END_V1_Dosis estatinas indicado",
  "Estatinas Visita 2" = "END_V2_Dosis estatinas indicado",
  "Estatinas Visita 3" = "END_V3_Dosis estatinas indicado",
  "Estatinas Visita 4" = "END_V4_Dosis estatinas indicado",
  "Estatinas Visita 5" = "END_V5_Dosis estatinas indicado",
  "Estatinas Visita 6" = "END_V6_Dosis estatinas indicado",
  "Estatinas Visita 7" = "END_V7_Dosis estatinas indicado",
  "Fibrato Visita 0" = "END_V0_Dosis fibrato indicado",
    "Fibrato Visita 1" = "END_V1_Dosis fibrato indicado",
  "Fibrato Visita 2" = "END_V2_Dosis fibrato indicado",
  "Fibrato Visita 3" = "END_V3_Dosis fibrato indicado",
  "Fibrato Visita 4" = "END_V4_Dosis fibrato indicado",
  "Fibrato Visita 5" = "END_V5_Dosis fibrato indicado",
  "Fibrato Visita 6" = "END_V6_Dosis fibrato indicado",
  "Fibrato Visita 7" = "END_V7_Dosis fibrato indicado",
  "Aspirina Visita 0" =  "END_V0_Aspirina indicado",
    "Aspirina Visita 1" =  "END_V1_Aspirina indicado",
  "Aspirina Visita 2" =  "END_V2_Aspirina indicado",
  "Aspirina Visita 3" =  "END_V3_Aspirina indicado",
  "Aspirina Visita 4" =  "END_V4_Aspirina indicado",
  "Aspirina Visita 5" =  "END_V5_Aspirina indicado",
  "Aspirina Visita 6" =  "END_V6_Aspirina indicado",
  "Aspirina Visita 7" =  "END_V7_Aspirina indicado",
  "Ezetimibe Visita 0"=  "END_V0_Ezetimibe indicado",
    "Ezetimibe Visita 1"=  "END_V1_Ezetimibe indicado",
  "Ezetimibe Visita 2" = "END_V2_Ezetimibe indicado",
  "Ezetimibe Visita 3"=  "END_V3_Ezetimibe indicado",
  "Ezetimibe Visita 4" = "END_V4_Ezetimibe indicado",
  "Ezetimibe Visita 5"=  "END_V5_Ezetimibe indicado",
  "Ezetimibe Visita 6" =  "END_V6_Ezetimibe indicado",
  "Ezetimibe Visita 7" =  "END_V7_Ezetimibe indicado"
              )%>%
  mutate(SEXO = recode(SEXO,"0"="Mujer","1" = "Hombre")) 

####---- 4.2 Porcentajes y Dosis ----####

# Calcula los % de pacientes que toman el medicamento y dosis dividido por grupo (dosis promedio de los sujetos que toman el medicamento)

#Tallas toma la informacion de las tallas en cada visita, como forma de contar cuantos sujetos fueron a cada visita
tallas = read_xlsx(paste0(dir,"1. Base CAIPaDi 2013-2019 master.xlsx"),guess_max = 3200)%>% select("CLASE_V1","SEXO",contains("Talla m"))
tallas$`NUT_V0_Talla m` = tallas$`NUT_V1_Talla m`

#Conteos cuenta cuantos sujetos fueron a cada visita, separados por grupo y sexo.
conteos <- tallas %>%
       group_by(CLASE_V1,SEXO) %>%
       summarise(across(contains("Talla m"), ~sum(!is.na(.)))) %>%
       pivot_longer(cols = starts_with("NUT_V"), names_to = "variable", values_to = "conteo") %>%
       mutate(variable = gsub("NUT_V", "", variable)) %>%
  mutate(SEXO = recode(SEXO,"0"="Mujer","1" = "Hombre")) %>%
  rename(visita = variable) %>%
  mutate(visita = as.numeric(gsub("_Talla m", "", visita))) %>%
  mutate(conteo = as.numeric(conteo))



medicamentos_summary <- tabla_medicamentos %>%
  group_by(CLASE_V1, SEXO) %>%
  summarise(
    total_personas = n(),  # Calculate total_personas separately
    across(where(is.numeric) & !c(total_personas), list(
      porcentaje = ~sum(!is.na(.))/total_personas * 100,
      promedio = ~mean(., na.rm = TRUE)
    ), .names = "{.col}_{.fn}") )%>%
      mutate(across(everything(), ~replace_na(., 0)))
  



get_conteo <- function(clase, sexo, vis) {
  res <- conteos %>%
    filter(CLASE_V1 == clase, SEXO == sexo, visita == vis) %>%
    pull(conteo)
  
  return(res)
}

for (i in 1:nrow(medicamentos_summary)) {
  for (col in names(medicamentos_summary)) {
    if (str_detect(col, "_porcentaje$")) {  
      vis <- str_extract(col, "\\d+(?=_porcentaje$)")  
      
      conteo <- get_conteo(medicamentos_summary$CLASE_V1[i], medicamentos_summary$SEXO[i], vis)
      conteo_viejo <- get_conteo(medicamentos_summary$CLASE_V1[i], medicamentos_summary$SEXO[i], "0")
      # Apply multiplication to the column value
      medicamentos_summary[i, col] <- medicamentos_summary[i, col] / conteo * conteo_viejo
    }
  }
}





####---- 4.3 Crear tabla de los promedios de dosis para un dado numero de visita ----####

# Poner aqui el numero de visita (0 corresponde a tratamiento previo a la visita 1):

Vi <- "1"

library(stringr)

df_visita0 <- medicamentos_summary %>%
  select(CLASE_V1, SEXO, contains(paste0("Visita ",Vi,"_promedio")))

df_long <- df_visita0 %>%
  pivot_longer(cols = -c(CLASE_V1, SEXO), 
               names_to = "Medication", 
               values_to = "Valor") %>%
  mutate(Medication = str_remove(Medication, paste0(" Visita ",Vi,"_promedio")))  

library(ggplot2)

# Create a heatmap-like visualization
pl <- ggplot(df_long, aes(x = interaction(CLASE_V1, SEXO), y = Medication, fill = Valor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Valor, 1)), size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(x = "Grupo (CLASE_V1 y SEXO)", y = "Medicación", fill = "Promedio",
       title = paste0("Dosis promedio para Visita ", Vi))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Diagonal labels
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),  # White surrounding area
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
pl

ggsave(filename = paste0(dir,"tablas_dosis_promedio/visita",Vi, ".png"), plot = pl, width = 6, height = 10, dpi = 300)

## AUTOMATIZAR para guardar TODAS las gráficas de una vez 




for (Vi in seq(0,7,1)) {
  
  
  library(stringr)
  
  df_visita0 <- medicamentos_summary %>%
    select(CLASE_V1, SEXO, contains(paste0("Visita ",Vi,"_promedio")))
  
  df_long <- df_visita0 %>%
    pivot_longer(cols = -c(CLASE_V1, SEXO), 
                 names_to = "Medication", 
                 values_to = "Valor") %>%
    mutate(Medication = str_remove(Medication, paste0(" Visita ",Vi,"_promedio")))  
  
  library(ggplot2)
  
  # Create a heatmap-like visualization
  pl <- ggplot(df_long, aes(x = interaction(CLASE_V1, SEXO), y = Medication, fill = Valor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Valor, 1)), size = 4, color = "black") +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    labs(x = "Grupo (CLASE_V1 y SEXO)", y = "Medicación", fill = "Promedio",
         title = paste0("Dosis promedio para Visita ", Vi))+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Diagonal labels
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),  # White surrounding area
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )
  ggsave(filename = paste0(dir,"tablas_dosis_promedio/visita",Vi, ".png"), plot = pl, width = 6, height = 10, dpi = 300)
}



####---- 4.4 Crear tabla de los porcentajes de personas con un medicamento para un dado numero de visita ----####

# Poner aqui el numero de visita (0 corresponde a previo a la visita 1):
Vi <- "7"

library(stringr)

df_visita0 <- medicamentos_summary %>%
  select(CLASE_V1, SEXO, contains(paste0("Visita ",Vi,"_porcentaje")))

df_long <- df_visita0 %>%
  pivot_longer(cols = -c(CLASE_V1, SEXO), 
               names_to = "Medication", 
               values_to = "Valor") %>%
  mutate(Medication = str_remove(Medication, paste0(" Visita ",Vi,"_porcentaje")))  

library(ggplot2)

# Create a heatmap-like visualization
ggplot(df_long, aes(x = interaction(CLASE_V1, SEXO), y = Medication, fill = Valor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Valor, 1)), size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(x = "Grupo (CLASE_V1 y SEXO)", y = "Medicación", fill = "Porcentaje",
       title = paste0("Porcentaje de pacientes para Visita ", Vi))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Diagonal labels
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),  # White surrounding area
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )






## AUTOMATIZAR para guardar TODAS las gráficas de una vez 




for (Vi in seq(0,7,1)) {
  
  df_visita0 <- medicamentos_summary %>%
    select(CLASE_V1, SEXO, contains(paste0("Visita ",Vi,"_porcentaje")))
  
  df_long <- df_visita0 %>%
    pivot_longer(cols = -c(CLASE_V1, SEXO), 
                 names_to = "Medication", 
                 values_to = "Valor") %>%
    mutate(Medication = str_remove(Medication, paste0(" Visita ",Vi,"_porcentaje")))  
  
  # Create a heatmap-like visualization
  pl <- ggplot(df_long, aes(x = interaction(CLASE_V1, SEXO), y = Medication, fill = Valor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Valor, 1)), size = 4, color = "black") +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    labs(x = "Grupo (CLASE_V1 y SEXO)", y = "Medicación", fill = "Porcentaje",
         title = paste0("Porcentaje de pacientes para Visita ", Vi))+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Diagonal labels
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),  # White surrounding area
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )
  
  ggsave(filename = paste0(dir,"tablas_porcentaje/visita",Vi, ".png"), plot = pl, width = 6, height = 10, dpi = 300)
}

####---- 4.5 Ver el cambio de porcentaje de pacientes como fucion del numero de visitas ----####


# Permite ver el % de sujetos que toman un medicaemnto como función del número de visita
#PONER AQUI EL MEDICAMENTO:
variable <- "Fibrato"

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
# Custom colors
custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

# Transform data to long format
variable_cols <- colnames(medicamentos_summary) %>%
  .[str_detect(., paste0("^", variable, " Visita \\d+_porcentaje$"))]


variable_data_long <- medicamentos_summary %>%
  select(all_of(c("CLASE_V1", "SEXO", variable_cols))) %>%
  pivot_longer(cols = all_of(variable_cols), 
               names_to = "Visita", 
               values_to = "Porcentaje") %>%
  mutate(Visita = as.numeric(str_extract(Visita, "\\d+(?=[^\\d]*$)")),
         Group = paste(CLASE_V1, SEXO, sep = "_")) 

# Plot the data
ggplot(variable_data_long, aes(x = Visita, y = Porcentaje, group = Group, color = Group)) +
  geom_line(size = 1.2) +  # Make lines thicker
  geom_point(size = 3) +   # Make points bigger
  scale_color_manual(values = custom_colors) +  # Set custom colors
  scale_x_continuous(breaks = unique(variable_data_long$Visita)) +  # Ensure all visits have ticks
  labs(
    title = paste0("Porcentaje de pacientes: ",variable),
    x = "Número de Visita",
    y = "Porcentaje de pacientes",
    color = "Grupo (CLASE_V1 y SEXO)"
  ) +
  theme_classic()

####---- 4.6 Ver el cambio de dosis promedio como fucion del numero de visitas ----####


# Es lo mismo que el 4.5 pero para ver el CAMBIO DE LA DOSIS 
#PONER AQUI EL MEDICAMENTO
variable <- "Aspirina"

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
# Custom colors
custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

# Transform data to long format
variable_cols <- colnames(medicamentos_summary) %>%
  .[str_detect(., paste0("^", variable, " Visita \\d+_promedio$"))]


variable_data_long <- medicamentos_summary %>%
  select(all_of(c("CLASE_V1", "SEXO", variable_cols))) %>%
  pivot_longer(cols = all_of(variable_cols), 
               names_to = "Visita", 
               values_to = "Promedio") %>%
  mutate(Visita = as.numeric(str_extract(Visita, "\\d+(?=[^\\d]*$)")),
         Group = paste(CLASE_V1, SEXO, sep = "_")) 

# Plot the data
ggplot(variable_data_long, aes(x = Visita, y = Promedio, group = Group, color = Group)) +
  geom_line(size = 1.2) +  # Make lines thicker
  geom_point(size = 3) +   # Make points bigger
  scale_color_manual(values = custom_colors) +  # Set custom colors
  scale_x_continuous(breaks = unique(variable_data_long$Visita)) +  # Ensure all visits have ticks
  labs(
    title = paste0("Dosis Promedio: ",variable),
    x = "Número de Visita",
    y = "Dosis Promedio",
    color = "Grupo (Subgrupo de Diabetes y SEXO)"
  ) +
  theme_classic()

####---- 4.7 Crear todas las graficas de porcentaje de pacientes como funcion del numero de visita y guardarlas -----####
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


variables <- c("Sulfonilureas", "SGLT2", "Insulina Premezcla",
               "Insulina en Bolo","Insulina Basal",
               "IECA", "iDPP4","Fibrato","Ezetimibe","Estatinas",
               "Diurético","Calcio Antagonista","Beta Bloqueador",
               "Aspirina","ARA2","Análogos de GLP-1")  

# Custom colors
custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

# Iterate over all variables and create plots
for (variable in variables) {
  variable_cols <- colnames(medicamentos_summary) %>%
    .[str_detect(., paste0("^", variable, " Visita \\d+_porcentaje$"))]
  
  
  variable_data_long <- medicamentos_summary %>%
    select(all_of(c("CLASE_V1", "SEXO", variable_cols))) %>%
    pivot_longer(cols = all_of(variable_cols), 
                 names_to = "Visita", 
                 values_to = "Porcentaje") %>%
    mutate(Visita = as.numeric(str_extract(Visita, "\\d+(?=[^\\d]*$)")),
           Group = paste(CLASE_V1, SEXO, sep = "_")) 
  
  # Plot the data
  p<-ggplot(variable_data_long, aes(x = Visita, y = Porcentaje, group = Group, color = Group)) +
    geom_line(size = 1.2) +  # Make lines thicker
    geom_point(size = 3) +   # Make points bigger
    scale_color_manual(values = custom_colors) +  # Set custom colors
    scale_x_continuous(breaks = unique(variable_data_long$Visita)) +  # Ensure all visits have ticks
    labs(
      title = paste0("Porcentaje de Sujetos: ",variable),
      x = "Número de Visita",
      y = "Porcentaje de Sujetos",
      color = "Grupo (Subgrupo de Diabetes y SEXO)"
    ) +
    theme_classic()
  
  # Save the plot to a file
  ggsave(filename = paste0(dir,"porcentaje_de_pacientes/",variable,".png"), plot = p, width = 10, height = 6, dpi = 300)
  
  print(paste("Saved plot for", variable))
}



#### 4.8 Crear todas las graficas de dosis promedio como funcion del numero de visita y guardarlas
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# List of variables to iterate over
variables <- c("Sulfonilureas", "SGLT2", "Insulina Premezcla",
               "Insulina en Bolo","Insulina Basal",
               "IECA", "iDPP4","Fibrato","Ezetimibe","Estatinas",
               "Diurético","Calcio Antagonista","Beta Bloqueador",
               "Aspirina","ARA2","Análogos de GLP-1")  # Example variables

# Custom colors
custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)

# Iterate over all variables and create plots
for (variable in variables) {
  variable_cols <- colnames(medicamentos_summary) %>%
    .[str_detect(., paste0("^", variable, " Visita \\d+_promedio$"))]
  
  
  variable_data_long <- medicamentos_summary %>%
    select(all_of(c("CLASE_V1", "SEXO", variable_cols))) %>%
    pivot_longer(cols = all_of(variable_cols), 
                 names_to = "Visita", 
                 values_to = "Promedio") %>%
    mutate(Visita = as.numeric(str_extract(Visita, "\\d+(?=[^\\d]*$)")),
           Group = paste(CLASE_V1, SEXO, sep = "_")) 
  
  # Plot the data
  p<-ggplot(variable_data_long, aes(x = Visita, y = Promedio, group = Group, color = Group)) +
    geom_line(size = 1.2) +  # Make lines thicker
    geom_point(size = 3) +   # Make points bigger
    scale_color_manual(values = custom_colors) +  # Set custom colors
    scale_x_continuous(breaks = unique(variable_data_long$Visita)) +  # Ensure all visits have ticks
    labs(
      title = paste0("Dosis Promedio: ",variable),
      x = "Número de Visita",
      y = "Dosis Promedio",
      color = "Grupo (CLASE_V1 y SEXO)"
    ) +
    theme_classic()
  
  # Save the plot to a file
  ggsave(filename = paste0(dir,"dosis_promedio/",variable,".png"), plot = p, width = 10, height = 6, dpi = 300)
  
  print(paste("Saved plot for", variable))  # Print message to track progress
}










####---- 5. Ver porcentajes que cumplieron metas----####


variable <- "HbA1c"

individual_data_long <- individual_data %>%
  select("CLASE_V1","SEXO", starts_with(paste0(variable," Visita")))%>%
  pivot_longer(cols = starts_with(paste0(variable," Visita")),
               names_to = "Visita",
               values_to = variable) %>%
  mutate(Visita = as.numeric(str_extract(Visita, "\\d$")))
               
cumplidores <- individual_data_long %>%
  group_by(CLASE_V1, SEXO, Visita) %>%
  # Aqui abajo cambiar Hba1c por la condicion que se quiera
  summarise(Percentaje_cumple = mean(.data[[variable]] < 6.6, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(Group = paste(CLASE_V1, SEXO, sep = "_"))

custom_colors <- c(
  "mard_Hombre" = "#8B0000",  
  "mard_Mujer"  = "#FF6347",  
  "mord_Hombre" = "#00008B",  
  "mord_Mujer"  = "#4682B4",  
  "siid_Hombre" = "#006400",  
  "siid_Mujer"  = "#32CD32",  
  "sird_Hombre" = "#8B008B",  
  "sird_Mujer"  = "#DA70D6"   
)


ggplot(cumplidores, aes(x = Visita, y = Percentaje_cumple,  group = Group,
                       color = Group )) +
  geom_line(size = 1.5) +
  geom_point(size = 4, stroke = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = unique(cumplidores$Visita)) +
  labs(title = "Porcentaje de sujetos que cumplen la condición",
       x = "Visita", 
       y = "Porcentaje (%)",
       color = "Group (CLASE_V1 & SEXO)") +
  theme_minimal(base_size = 14) +  #
  theme(panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_blank()) 
