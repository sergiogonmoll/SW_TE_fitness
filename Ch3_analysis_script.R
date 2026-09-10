library(lubridate)
library(coxme)
library(performance)
library(lme4)
library(tidyverse)
library(lmerTest)
library(readxl)
library(glmmTMB)
library(ggeffects)
library(survminer)
library(ggforestplot)
library(emmeans)
library(ggplot2)
library(DHARMa)
library(ggnewscale)

#### Survival analysis using coxme ####
survival_data<- read.csv("Survival_Data_Ch3.csv")

coxme_model <- coxme(Surv(PeriodEnd-PeriodStart, survival) ~
                                    logratD + logratS +lGCF + RainVar + Dam_RainVar + Sire_RainVar + FFL + FDL + FSL + dam_age + helpers_cat +
                                    age + Sex*lGCF+ Sex*(dam_age) + PeriodYear +
                                    (1|PeriodYear) + (1|BirdID)+ (1|dam) + (1|sire), data=survival_data)
extract_coxme_table <- function (mod){
  beta <- mod$coefficients
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

#### Extract coxme coefficient table ####

coxme_coef <- extract_coxme_table(mod = coxme_model)
coxme_coef$names <- c("Mother TE contribution", "Father TE contribution", "Offspring TE genomic content", "Rainfall variance","Maternal rainfall variance",
                      "Paternal rainfall variance", "Inbreeding coefficient",
                      "Mother inbreeding coefficient", "Father inbreeding coefficient",
                      "Mother's age", "Helper presence in natal territory",
                      "Age", "Sex",
                      "Period Year","GTEC:Sex(Male)","Mother's age:Sex(Male)")

#### Plot model results ####
coxme_hazards_plot<- ggforestplot::forestplot(coxme_coef, name = names, estimate = beta, se=se, pvalue = p, ci = 0.95, logodds = TRUE, xlab = "Hazard ratio")

#### Run sex-specific posthoc models ####
survival_data_F <- subset(survival_data, survival_data$Sex == 0)
survival_data_M <- subset(survival_data, survival_data$Sex == 1)

surv_F_m <- coxme(Surv(PeriodEnd-PeriodStart, survival) ~
                    logratD + logratS + lGCF +RainVar+ Dam_RainVar + Sire_RainVar+ FFL + FDL + FSL + dam_age + helpers_cat +
                    age + dam_age+ PeriodYear + logratD:RainVar +
                    (1|PeriodYear) + (1|BirdID)+ (1|dam) + (1|sire), data=survival_data_F)
surv_M_m <- coxme(Surv(PeriodEnd-PeriodStart, survival) ~
                    logratD + logratS + lGCF +RainVar+ Dam_RainVar + Sire_RainVar+ FFL + FDL + FSL + dam_age + helpers_cat +
                    age + dam_age + PeriodYear +
                    (1|PeriodYear) + (1|BirdID)+ (1|dam) + (1|sire), data=survival_data_M)

#### Plot sex-specific coxme model results ####
coxme_coef_M <- extract_coxme_table(mod = surv_M_m)
coxme_coef_M$names <- c("Mother TE contribution", "Father TE contribution", "Offspring TE genomic content", "Rainfall variance", "Maternal rainfall variance",
                        "Paternal rainfall variance", "Inbreeding coefficient",
                        "Mother inbreeding coefficient", "Father inbreeding coefficient",
                        "Mother's age", "Helper presence in natal territory",
                        "Age",
                        "Period Year")
coxme_hazards_plot_M<- ggforestplot::forestplot(coxme_coef_M, name = names, estimate = beta, se=se, pvalue = p, ci = 0.95, logodds = TRUE, xlab = "Hazard ratio")
coxme_coef_F <- extract_coxme_table(mod = surv_F_m)
coxme_coef_F$names <- c("Mother TE contribution", "Father TE contribution", "Offspring TE genomic content", "Rainfall variance","Maternal rainfall variance",
                        "Paternal rainfall variance", "Inbreeding coefficient",
                        "Mother inbreeding coefficient", "Father inbreeding coefficient",
                        "Mother's age", "Helper presence in natal territory",
                        "Age",
                        "Period Year", "Mother TE contribution:Rainfall variance")

coxme_hazards_plot_F<- ggforestplot::forestplot(coxme_coef_F, name = names, estimate = beta, se=se, pvalue = p, ci = 0.95, logodds = TRUE, xlab = "Hazard ratio")

#### Lifetime reproductive success analysis ####
LRS_data <- read.csv("LRS_Data_Ch3.csv")
LRS_data$Sex <- as.factor(LRS_data$Sex)

LRS_model <- glmmTMB::glmmTMB(offspring ~ logratD + logratS + FFL + FDL + FSL+
                                           Sex+ dam_age + LayYear+ RainVar +  Dam_RainVar + Sire_RainVar + helpers_cat+
                                           Sex*logratD + Sex*logratS  + (1|dam), ziformula = ~0, REML = T,
                                         family = 'poisson',control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")), data = LRS_data)
summary(LRS_model)
check_singularity(LRS_model)
check_collinearity(LRS_model)
check_overdispersion(LRS_model)
plot(simulateResiduals(LRS_model))

#### Obtain LRS model predictions ####
d.predLRS_D_dam <- ggpredict(LRS_model, terms = c("logratD [all]","Sex"), bias_correction = T)
d.predLRS_D_sire <- ggpredict(LRS_model, terms = c("logratS [all]","Sex"), bias_correction = T)
d.predLRS_FFL <- ggpredict(LRS_model, terms = "FFL [all]", bias_correction = T)

#### Plot LRS model predictions ####
plotDdam_Sex <- ggplot(d.predLRS_D_dam,aes(x=x,y=predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_smooth(aes(color = group),se = FALSE, linewidth=0.6, alpha = 0.6) +
  geom_point(data = LRS_data, width = 1, height = 0, size = 1.0, alpha = 0.6, aes(x=logratD, y=offspring, color = Sex)) +
  scale_x_continuous(limits = c(-2.5,3), breaks = seq(-3,3,1))+
  scale_y_continuous(limits = c(0,22), breaks = seq(0,30,10))+
  theme_minimal() +
  theme(
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.title = element_text(size = 14),
    axis.line = element_line(color="black", size = 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6),
    panel.background = element_rect(fill = "white", colour = "darkgrey")
  ) +
  xlab("Maternal TE contribution to offspring (SD)") +
  ylab("Offspring lifetime reproductive success") +
  annotate(geom = "text", x=0.2, y=20, label="**", size = 5)+
  scale_fill_viridis_d(name = "Sex", labels = c("Female", "Male")) +
  scale_color_viridis_d(name = "Sex", labels = c("Female", "Male"))

plotDsire_LRS <- ggplot(d.predLRS_D_sire,aes(x=x,y=predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_smooth(aes(color = group),se = FALSE, linewidth=0.6, alpha = 0.6) +
  geom_point(data = LRS_data, width = 1, height = 0, size = 1.0, alpha = 0.6, aes(x=logratD, y=offspring, color = Sex)) +
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1))+
  scale_y_continuous(limits = c(0,22), breaks = seq(0,30,10))+
  theme_minimal() +
  theme(
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.title = element_text(size = 14),
    axis.line = element_line(color="black", size = 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6),
    panel.background = element_rect(fill = "white", colour = "darkgrey")
  ) +
  xlab("Paternal TE contribution to offspring (SD)") +
  ylab("Offspring lifetime reproductive success") +
  annotate(geom = "text", x=0.2, y=21, label=".", size = 5)+
  scale_fill_viridis_d(name = "Sex", labels = c("Female", "Male")) +
  scale_color_viridis_d(name = "Sex", labels = c("Female", "Male"))

plotFrohLRs <- ggplot(d.predLRS_FFL,aes(x=x,y=predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_smooth(se = FALSE, linewidth=0.6, alpha = 0.6) +
  geom_point(data = LRS_data, width = 1, height = 0, size = 1.0, alpha = 0.6, aes(x=FFL, y=offspring)) +
  #scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1))+
  #scale_y_continuous(limits = c(0,20), breaks = seq(0,20,5))+
  theme_minimal() +
  theme(
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.title = element_text(size = 14),
    axis.line = element_line(color="black", size = 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6),
    panel.background = element_rect(fill = "white", colour = "darkgrey")
  ) +
  xlab("Inbreeding coefficient (SD)") +
  ylab("Offspring lifetime reproductive success") +
  annotate(geom = "text", x=-0.5, y=20, label="**", size = 5)

#### Reproductive senescence analysis #####
F_ARS_data <- read.csv("ARS_F_data_Ch3.csv")
M_ARS_data <- read.csv("ARS_M_data_Ch3.csv")

ARS_m_F<- glmmTMB(offspring ~ logratD + logratS + DeltaAgeParent +I(DeltaAgeParent^2) + MeanAgeParent + FFL +
                           FDL + FSL + HelpersF + GroupSize
                         + Focal_RainVarBlood + Dam_RainVar + Sire_RainVar 
                         + (1|BirdID),ziformula = ~0, REML = T,
                         family = 'genpois', data = F_ARS_data, control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

ARS_m_M <- glmmTMB(offspring ~ logratD+ logratS +DeltaAgeParent +I(DeltaAgeParent^2) + MeanAgeParent + FFL +
                            FDL*DeltaAgeParent + FSL + HelpersF + GroupSize
                          + Focal_RainVarBlood + Dam_RainVar + Sire_RainVar
                          + (1|BirdID),ziformula = ~0, REML = T,
                          family = 'genpois', data = M_ARS_data,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

summary(ARS_m_F)
check_singularity(ARS_m_F)
check_collinearity(ARS_m_F)
check_overdispersion(ARS_m_F)
plot(simulateResiduals(ARS_m_F))

summary(ARS_m_M)
check_singularity(ARS_m_M)
check_collinearity(ARS_m_M)
check_overdispersion(ARS_m_M)
plot(simulateResiduals(ARS_m_M))

#### Get reproductive senscence model predictions per sex ####
d.predARS_F <- ggpredict(ARS_m_F, terms = c("DeltaAgeParent [all]"), bias_correction = T)
d.predARS_F$Sex<- 0
d.predARS_M <- ggpredict(ARS_m_M, terms = c("DeltaAgeParent [all]"), bias_correction = T)
d.predARS_M$Sex<- 1
d.predARS_Sen <- rbind(d.predARS_F,d.predARS_M)
d.predARS_Sen$Sex<- as.factor(d.predARS_Sen$Sex)
M_ARS_data$Sex<- as.factor(M_ARS_data$Sex)
F_ARS_data$Sex<- as.factor(F_ARS_data$Sex)

d.predARS_FDL_M <- ggpredict(ARS_m_M, terms =c("DeltaAgeParent [all]","FDL"), bias_correction = T)

#### Plot reproductive senscence model predictions per sex ####
plotSexAgeARS <- ggplot(d.predARS_Sen,aes(x=x,y=predicted, color=Sex)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Sex), alpha = 0.25) +
  geom_smooth(aes(fill=Sex), se = FALSE, linewidth=0.6, alpha = 0.2) +
  geom_jitter(data = F_ARS_data, width = 1, height = 0, size = 1.0, alpha = 0.6, aes(x=DeltaAgeParent, y=offspring, color=Sex)) +
  geom_jitter(data = M_ARS_data, width = 1, height = 0, size = 1.0, alpha = 0.6, aes(x=DeltaAgeParent, y=offspring, color=Sex)) +
  scale_x_continuous(limits = c(-9,10), breaks = seq(-10,10,5))+
  scale_y_continuous(limits = c(0,7), breaks = seq(0,7,2))+
  theme_minimal() +
  theme(
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.title = element_text(size = 14),
    axis.line = element_line(color="black", size = 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6),
    panel.background = element_rect(fill = "white", colour = "darkgrey")
  ) +
  xlab("Age (mean-centered)") +
  ylab("Offspring produced yearly") +
  labs(title = "Maternal inbreeding effects on reproductive senescence") +
  scale_fill_viridis_d(name = "Sex", labels = c("Female","Male")) +
  scale_color_viridis_d(name = "Sex", labels = c("Female","Male"))+
  facet_wrap(~Sex, ncol = 1,labeller = labeller(Sex = c("0" = "Female","1" = "Male")))

plotFrohFARS <- ggplot(d.predARS_FDL_M,aes(x=x,y=predicted,)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill= group), alpha = 0.25) +
  geom_smooth(se = FALSE, linewidth=0.6, alpha = 0.2, aes(color=group)) +
  scale_fill_viridis_d(name = "Maternal Inbreeding", labels = c("Low","Medium","High")) +
  scale_color_viridis_d(name = "Maternal Inbreeding", labels = c("Low","Medium","High"))+
  new_scale(new_aes = "color")+
  geom_jitter(data = M_ARS_data, height = 0, width = 0.4, size = 1.0, alpha = 0.6, aes(x=DeltaAgeParent, y=offspring, color = Range_FD)) +
  scale_color_viridis_d(name = "Maternal Inbreeding", labels = c("Low","Medium","High"))+
  scale_y_continuous(limits = c(0,5), breaks = seq(0,7,1))+
  theme_minimal() +
  theme(
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.title = element_text(size = 14),
    axis.line = element_line(color="black", size = 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6),
    panel.background = element_rect(fill = "white", colour = "darkgrey")
  ) +
  xlab("Age (mean-centered)") +
  ylab("Offspring produced yearly") +
  labs(title = "Inbreeding depression on male annual reproductive success")

plot_coxme <- ggarrange(coxme_hazards_plot_F, coxme_hazards_plot_M, ncol = 2, labels = c("A. Female","B. Male"))
ggsave("Figure1_Survival.tiff", plot_coxme, dpi = 'retina', height = 15, width = 30, units = "cm")

LRS_top <- ggarrange(plotDdam_Sex, plotDsire_LRS, labels = c("A","B"),ncol=2, common.legend = T, legend = "right")
LRS_fig<-ggarrange(LRS_top,plotFrohLRs, labels = c("","C"), nrow = 2)
ggsave("Figure2_LRS.tiff", LRS_fig, dpi = 'retina', height = 26, width = 24, units = "cm")

ggsave("Figure3_Senescence_vs_maternal_inbreeding.tiff", plotFrohFARS, dpi = 'retina', height = 15, width = 21, units = "cm")

