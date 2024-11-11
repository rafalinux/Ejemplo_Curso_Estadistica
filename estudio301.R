###
# Clase 2 - Curso de Estadística
# author: Rafa Carretero
# date: 23 Abril, 2024
###
#
#######################################
#                                     #
# CARGA DATOS Y LIBRERIAS             #
#                                     #
#######################################
library(dplyr)
library(gtsummary)
setwd("~/Python_Projects/Curso_Bioestadistica")
#setwd("~/Documentos/Curso_Estadística_Mostoles/Clase_2")
pacientes <- read.table("estudio301.csv", header=TRUE, sep="&")
#pacientes <- read.csv("estudio301_libreoffice.csv", header=TRUE, sep=",")

pacientes$tratamiento <- as.factor(pacientes$tratamiento)
pacientes$DiseaseLo <- as.factor(pacientes$DiseaseLo)
pacientes$muerto <- as.factor(pacientes$muerto)
pacientes$recaida <- as.factor(pacientes$recaida)

######################################
#                                    #
# AÑADIR ETIQUETAS                   #
#                                    #
######################################
library(Hmisc)
label(pacientes$EDAD) <- "Edad"
label(pacientes$DiseaseLo) <- "Localización de la enfermedad"
label(pacientes$BPI) <- "Escala de dolor BPI"
label(pacientes$C_ECOG) <- "Escala ECOG basal"
label(pacientes$TRATprevio) <- "N. de tratamiento previos"
label(pacientes$Visceral) <- "Enf. viceral"
label(pacientes$race) <- "Raza"
label(pacientes$muerto) <- "Muertes"
label(pacientes$recaida) <- "Recaída"
levels(pacientes$DiseaseLo)=c("Hueso","Node","Hígado")
levels(pacientes$tratamiento)=c("Placebo","Acetato de Abiraterona")
levels(pacientes$muerto)=c("Vivo","Fallecido")
levels(pacientes$recaida)=c("No recaída","Recaída")

#
#######################################
#                                     #
# TABLA CARACTERISTICAS GENERALES     #
#                                     #
#######################################
#
pacientes %>% tbl_summary()

pacientes %>%
  tbl_summary(by = muerto) %>%
  add_p()

pacientes %>%
  select(tratamiento, EDAD, DiseaseLo, BPI, C_ECOG, TRATprevio, Visceral, race) %>%
  
  tbl_summary(by = tratamiento,
              label = list(EDAD ~ "Edad",
                           DiseaseLo ~ "Localización de la enfermedad",
                           BPI ~ "Escala de dolor BPI",
                           C_ECOG ~ "Escala ECOG basal", 
                           TRATprevio ~ "N. de tratamiento previos", 
                           Visceral ~ "Enf. viceral", 
                           race ~ "Raza"),
              missing = "no",
              missing_text = NULL,
              percent="row") %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  add_overall() %>%
  modify_footnote(
    all_stat_cols() ~ "Medidas expresadas en Mediana (RIQ) o Frecuencia (%)"
  ) %>%
  #add_n() %>%
  modify_caption("**Tabla 1. Características de la cohorte**") %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Tipo de Tratamiento**") %>%
  bold_labels()

######################################
#                                    #
# TABLA con TABLEONE                 #
#                                    #
######################################
library(tableone)

CreateTableOne(data = pacientes)

## Get variables names
dput(names(pacientes))

myVars <- c("EDAD", 
            "DiseaseLo", "C_ECOG", "BPI", "TRATprevio", "Visceral", "race")
## Vector of categorical variables that need transformation
catVars <- c("DiseaseLo", "C_ECOG", "TRATprevio", "race", 
             "tratamiento")
## Create a TableOne object
tab2 <- CreateTableOne(vars = myVars, 
                       data = pacientes, 
                       factorVars = catVars,
                       strata = "tratamiento")
print(tab2, varLabels = TRUE) 

#######################################
#                                     #
# REGRESION LOGISTICA                 #
#                                     #
#######################################
#
# build logistic regression model
model1 <- glm(muerto ~ tratamiento, pacientes, family = binomial)
tbl_regression(model1, exponentiate = TRUE)

table(pacientes$muerto, pacientes$tratamiento)

######################################
#                                    #
# SUPERVIVENCIA                      #
#                                    #
######################################
library(survival)
library(ggsurvfit)
library(survminer)

#las variables "muerto" y "recaida" tienen que volver a ser
# números, porque si no, el análisis de supervivencia no sale
pacientes$muerto <- as.numeric(pacientes$muerto)
pacientes$recaida <- as.numeric(pacientes$recaida)

#######################################
#                                     #
# CURVA DE KAPLAN-MEIER (I)           #
#                                     #
#######################################
#
km_fit <- survfit(Surv(time_Global, muerto) ~ 1, data=pacientes)

ggsurvplot(km_fit, data = pacientes, risk.table = TRUE, conf.int = TRUE)


curva_1 <- ggsurvplot(km_fit, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = FALSE,
                      #pval=TRUE,
                      ggtheme = theme_bw(),
                      #xlim = c(0, 6400),
                      #xscale = "d_y",
                      legend = c("bottom"),
                      title = "Curva de Kaplan-Meier (supervivencia global)",
                      xlab = "Tiempo (meses)",
                      ylab = "Supervivencia (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      censor = FALSE
)
curva_1
# https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
curva_1$plot + scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  scale_colour_manual(values = c("steelblue","red","red","red")) 
#ggsave("kaplan_overall_alternativa.png", height = 6, width = 8)

#######################################
#                                     #
# CURVA DE KAPLAN-MEIER (II)          #
#                                     #
#######################################
#
km_fit <- survfit(Surv(time_f, recaida) ~ 1, data=pacientes)

ggsurvplot(km_fit, data = pacientes, risk.table = TRUE, conf.int = TRUE)


curva_1 <- ggsurvplot(km_fit, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = FALSE,
                      #pval=TRUE,
                      ggtheme = theme_bw(),
                      #xlim = c(0, 6400),
                      #xscale = "d_y",
                      legend = c("bottom"),
                      title = "Curva de Kaplan-Meier (supervivencia libre de enfermedad)",
                      xlab = "Tiempo (meses)",
                      ylab = "Supervivencia (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      censor = FALSE
)
curva_1
# https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
curva_1$plot + scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  scale_colour_manual(values = c("steelblue","red","red","red")) 
#ggsave("kaplan_overall_alternativa.png", height = 6, width = 8)

#######################################
#                                     #
# CURVA DE KAPLAN-MEIER (III)         #
#                                     #
#######################################
#
km_fit <- survfit(Surv(time_Global, muerto) ~ tratamiento, data=pacientes)

ggsurvplot(km_fit, data = pacientes, risk.table = TRUE, conf.int = TRUE)

curva_2 <- ggsurvplot(km_fit, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = FALSE,
                      #pval=TRUE,
                      ggtheme = theme_bw(),
                      #xlim = c(0, 6400),
                      #xscale = "d_y",
                      legend = c("bottom"),
                      title = "Curva de Kaplan-Meier (supervivencia global) \n (Abiterona vs Placebo)",
                      xlab = "Tiempo (meses)",
                      ylab = "Supervivencia (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      pval = TRUE,
                      censor = FALSE
)
curva_2
#ggsave("kaplan_overall.pdf", width = 8, height = 6)
#ggsave("kaplan_respondedor.png", height = 6, width = 8)
# https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
superv_1 <- curva_2$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  scale_colour_manual(values = c("#2E9FDF","#E7B800","steelblue","red","red","red")) +
  labs(title = "Curva de Kaplan-Meier",
       subtitle="Supervivencia global")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))
superv_1
#ggsave("kaplan_respondedor_1.png", height = 6, width = 8)

#######################################
#                                     #
# CURVA DE KAPLAN-MEIER (IV)          #
#                                     #
#######################################
#
km_fit <- survfit(Surv(time_f, recaida) ~ tratamiento, data=pacientes)

ggsurvplot(km_fit, data = pacientes, risk.table = TRUE, conf.int = TRUE)


curva_3 <- ggsurvplot(km_fit, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = FALSE,
                      #pval=TRUE,
                      ggtheme = theme_bw(),
                      #xlim = c(0, 6400),
                      #xscale = "d_y",
                      legend = c("bottom"),
                      title = "Curva de Kaplan-Meier (libre de enfermedad) \n (Abiterona vs placebo)",
                      xlab = "Tiempo (meses)",
                      ylab = "Supervivencia (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      pval = TRUE,
                      censor = FALSE
)
curva_3
#ggsave("kaplan_overall.pdf", width = 8, height = 6)
#ggsave("kaplan_respondedor.png", height = 6, width = 8)
# https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
superv_2 <- curva_3$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  scale_colour_manual(values = c("#2E9FDF","#E7B800","steelblue","red","red","red")) +
  labs(title = "Curva de Kaplan-Meier",
       subtitle="Supervivencia libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))
superv_2
#ggsave("kaplan_respondedor_2.png", height = 6, width = 8)

#######################################
#                                     #
# REGRESION DE COX                    #
#                                     #
#######################################
#
cox.clasico.1 <- coxph(Surv(time_Global, muerto)~ tratamiento, 
                     data = pacientes,na.action=na.exclude)
cox.clasico.2 <- coxph(Surv(time_f, recaida)~ tratamiento, 
                       data = pacientes,na.action=na.exclude)

t1 <- tbl_regression(cox.clasico.1, exponentiate = TRUE) %>%
  add_n()
t1

#https://www.danieldsjoberg.com/gtsummary/
t2 <-
  coxph(Surv(time_f, recaida) ~ tratamiento, pacientes) %>%
  tbl_regression(exponentiate = TRUE) %>%
  add_n()

# merge tables
tbl_merge_ex1 <-
  tbl_merge(
    tbls = list(t1, t2),
    tab_spanner = c("**Supervivencia global**", "**Tiempo libre de enfermedad**")
  ) %>%
  modify_caption("**Tabla 2. Hazard Ratios**") %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()
tbl_merge_ex1

cox.clasico.3 <- coxph(Surv(time_Global, muerto)~ tratamiento + C_ECOG + BPI, 
                       data = pacientes,na.action=na.exclude)

t3 <- tbl_regression(cox.clasico.3, exponentiate = TRUE) %>%
  add_n()
t3

#######################################
#                                     #
# PROBABILIDAD DE SUPERVIVENCIA       #
#                                     #
#######################################
#
pacientes$trt.libre <- pacientes$tratamiento
pacientes$trt.global <- pacientes$tratamiento
label(pacientes$trt.libre) <- "Libre de enfermedad"
label(pacientes$trt.global) <- "Global"

tbl_survfit_ex1 <- 
  list(
    survfit(Surv(time_Global, muerto) ~ trt.global, pacientes),
    survfit(Surv(time_f, recaida) ~ trt.libre, pacientes)
  ) %>%
  tbl_survfit(
  times = c(6, 12),
  label_header = "**Mes {time}**"
) %>% 
  #add_p() %>%
  bold_labels() %>%
  modify_header(label ~ "**Supervivencia**") %>%
  italicize_levels()
tbl_survfit_ex1

survfit(Surv(time_Global, muerto) ~ 1, pacientes) %>% 
  tbl_survfit(
    times = 365.25,
    label_header = "**Supervivencia (global) a 1 año (IC 95%)**"
  )

survfit(Surv(time_f, recaida) ~ 1, pacientes) %>% 
  tbl_survfit(
    times = 365.25,
    label_header = "**Supervivencia (libre de enfermedad) a 1 año (IC 95%)**"
  )

survfit(Surv(time_Global, muerto) ~ 1, pacientes) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Mediana de supervivencia (IC 95%)**"
  )
#######################################
#                                     #
# REGRESION LOGISTICA                 #
#                                     #
#######################################
#
pacientes$muerto[pacientes$muerto==1] <- 0
pacientes$muerto[pacientes$muerto==2] <- 1
pacientes$muerto_factor <- as.numeric(pacientes$muerto_factor)
# build logistic regression model
model1 <- glm(muerto ~ tratamiento, pacientes, family = binomial)
tbl_regression(model1, exponentiate = TRUE)

table(pacientes$muerto, pacientes$tratamiento)