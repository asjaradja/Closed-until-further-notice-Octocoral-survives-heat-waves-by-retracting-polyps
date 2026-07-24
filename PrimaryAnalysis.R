# 0 SET UP ENVIRONMENT ----
## 0.1 Load packages ----
library(cowplot)    # for figure organization
library(emmeans)    # for estimated marginal means
library(ggsignif)   # for brackets indicating statistically significant differences
library(multcomp)   # for compact letter displays
library(MuMIn)      # for linear mixed effects R-squared estimates
library(nlme)       # for fitting random-effects models
library(sjPlot)     # for checking model assumptions
library(tidyverse)  # for data wrangling 

## 0.2 Load datasets ----
# Dataset for heat-treated colonies
mastersheet=read_csv("mastersheet.csv",
                     col_types = 'cfffff')

# Dataset for extended and retracted polyps
CVOR=read_csv("CVOR.csv", col_types = 'ffn')    #Respiration
CVOP=read_csv("CVOP.csv", col_types = 'ffn')  #Photosynthesis

# 1 MORTALITY ----
## 1.1 Generate data table ----
tbl.mortality = table(mastersheet$Treatment, mastersheet$Mortality)

## 1.2 Run Fisher's Exact Test ----
fisher.test(tbl.mortality)

# 2 NECROSIS ----
## 2.1 Generate data table ----
tbl.necrosis = table(mastersheet$Treatment, mastersheet$Necrosis)

## 2.2 Run Fisher's Exact Test ----
fisher.test(tbl.necrosis)

# 3 BEHAVIOURAL ASSAY DATA ----
## 3.1 Create linear mixed effects model ----
modAvgT=lme(data=mastersheet,
             AvgT~Week*Treatment,
             random=~1|MC)

## 3.2 Check model assumptions ----
mod.plotsAvgT=plot_model(modAvgT,type="diag")
mod.plotsAvgT[[2]]
mod.plotsAvgT[[3]]

## 3.3 Generate model summary ----
anova.lme(modAvgT,            # ANOVA type summary
          type="marginal")

r.squaredGLMM(modAvgT)        # R-squared estimates

## 3.4 Generate plot ----
pBehavioral=ggplot(data=mastersheet,
                      aes(x=Week, y=AvgT, fill=Treatment)) + 
  geom_boxplot(color = "black",
               position = position_dodge(width = 0.8),
               width = 0.6,outlier.shape = NA,alpha = 0.6) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                dodge.width = 0.8,
                jitter.width = 0.2),
              shape = 21,size = 1.5,alpha = 0.9,
              show.legend = TRUE) +
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral"),
    labels = c("Ctrl" = "Control", "Exp" = "Simulated Heat Wave"),
    name = "Treatment") +
  scale_y_continuous(breaks = 1:6) +
  labs(x = "Week",
       y = "Average Extension Score") +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

# Add brackets
pBehavioral = pBehavioral +
  geom_signif(annotation=c("*","*","*","*"),
              xmin=c(0.8,1.8,2.8,3.8),
              xmax=c(1.2,2.2,3.2,4.2),
              y=3.1,
              tip_length = c(0.03,0.45,
                             0.03,0.55,
                             0.03,0.55,
                             0.03,0.7))

# Separate the legend
treatmentLegend=get_legend(pBehavioral+
                             theme(text = element_text(size=8)))
ggdraw(treatmentLegend)

ggsave("figures/treatmentLegend.pdf",
       width=3.25,
       height=2.6,
       units="cm")

#Separate the plot
pBehavioral+
  theme(legend.position="none")

ggsave("figures/pBehavioural.pdf",
       device="pdf",
       width = 17.75,
       height = 16,
       units = "cm")

# 4 RETRACTED-EXTENDED COMPARISON ----
## 4.1 Re-format data ----
### Net Photosynthesis
CVOP_wide <- CVOP %>%
  pivot_wider(names_from = S, values_from = PR)

### Dark Respiration
CVOR_wide <- CVOR %>%
  pivot_wider(names_from = S, values_from = RR)

## 4.2 T-tests ----
### Net Photosynthesis
t.test(CVOP_wide$C, CVOP_wide$O, paired = TRUE)
### Dark Respiration
t.test(CVOR_wide$C, CVOR_wide$O, paired = TRUE)


## 4.3 Summary statistics ----
### Net Photosynthesis
CVOP.sum = CVOP %>%
  group_by(S) %>%
  summarise(
    m = mean(PR),
    se = sd(PR) / sqrt(n()),
    .groups = "drop")

### Respiration
CVOR.sum= CVOR %>%
  group_by(S) %>%
  summarise(m=mean(RR),
            se=sd((RR)/sqrt(n())))

## 4.4 Generate plots ----
### Net Photosynthesis
pCVOP <- ggplot() +
  geom_jitter(data = CVOP,
              aes(x = S, y = PR, fill = S),
              position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.2),
              color = "black",shape = 21,size = 2,
              show.legend = FALSE) +
  scale_fill_manual(values = c("C" = "chartreuse3", "O" = "chartreuse3"),
                    guide = "none") +
  geom_errorbar(data = CVOP.sum,
                aes(x = factor(S, levels = c("C", "O")),
                    ymin = m - se,
                    ymax = m + se),
                color = "black",
                position = position_dodge(width = 0.8),
                width = 0.2,linewidth = 1,
                show.legend = FALSE) +
  geom_errorbar(data = CVOP.sum,
                aes(x = factor(S, levels = c("C", "O")),
                    ymin = m,
                    ymax = m),
                color = "black",
                width = 0.3,linewidth = 2,
                position = position_dodge(width = 0.8),
                show.legend = FALSE) +
  scale_x_discrete(limits = c("C", "O"),
                   labels = c("C" = "Retracted",
                              "O" = "Extended")) +
  labs(x = "Status",
       y = expression(italic(P)[net] ~ "(" * mu * mol ~ O[2] ~ cm^-2 ~ hr^-1 * ")")) +
  annotate("text",
           x = 1.5,
           y = max(CVOP$PR) * 1.05,
           label = "*",
           size = 10/.pt) +
  theme_classic()+
  theme(text=element_text(size=8))

pCVOP

ggsave("figures/CVOP.pdf",
  width = 11.3,
  height = 16.86,
  units="cm")

### Dark Respiration
pCVOR = ggplot() +
  geom_jitter(data = CVOR,
              aes(x = S, y = RR, fill = S),
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              shape = 21,size = 2,
              color = "black",
              show.legend = FALSE) +
  scale_fill_manual(values = c("C" = "cornflowerblue", "O" = "cornflowerblue"),
                    guide = "none") +
  geom_errorbar(data = CVOR.sum,
    aes(x = factor(S, levels = c("C", "O")),
        ymin = m - se, 
        ymax = m + se),
    color = "black",
    width = 0.2,linewidth = 1,
    position = position_dodge(width = 0.8),
    show.legend = FALSE) +
  geom_errorbar(data = CVOR.sum,
    aes(x = factor(S, levels = c("C", "O")),
        ymin = m,
        ymax = m),
    color = "black",
    width = 0.3,linewidth = 2,
    position = position_dodge(width = 0.8),
    show.legend = FALSE) +
  scale_x_discrete(labels = c("C" = "Retracted", "O" = "Extended")) +
  labs(x = "Status",
       y = expression(italic(R)[dark]~"(" * mu * mol~O[2]~cm^-2~hr^-1 * ")"))+
  annotate("text", 
           x = 1.5, 
           y = max(CVOR$RR) * 1.05, 
           label = "*", 
           size = 10/.pt) +
  theme_classic()+
  theme(text=element_text(size=8))

pCVOR

ggsave("figures/CVOR.pdf",
       width = 11.3,
       height = 16.86,
       units="cm")

# 5 NET PHOTOSYNTHESIS ----
## 5.1 Remove negative net photosynthesis value ----
mastersheetPn=mastersheet[-94, ]

## 5.2 Create linear mixed effects model ----
modPn=lme(data=mastersheetPn,
          log10(Pnet)~Week*Treatment,
          random=~1|Tank/MC)

## 5.3 Check model assumptions ----
mod.plotsPn=plot_model(modPn,type="diag")
mod.plotsPn[[2]]
mod.plotsPn[[3]]

## 5.4 Return model summary and R-squared values ----
anova(modPn, type="marginal")

r.squaredGLMM(modPn)

## 5.5 Pairwise comparisons ----
emmPn = emmeans(modPn, ~ Treatment * Week)

pairs(emmPn,simple="Treatment")
pairs(emmPn,simple="Week")
#heat-treated polyps have significantly lower 
#net photosynthesis at Weeks 2 and 4 than at week 1.

## 5.6 Summary statistics ----
Pn.sum = mastersheetPn %>%
  group_by(Treatment) %>%
  summarise(m=mean(Pnet,na.rm=TRUE),
            se=sd(Pnet,na.rm=TRUE)/sqrt(sum(!is.na(Pnet))))

## 5.7 Generate plot and save ----
pPnet = ggplot(data = mastersheetPn,
               aes(x = Week,
                   y = log(Pnet),
                   fill = Treatment)) + 
  geom_boxplot(color = "black",
               position = position_dodge(width = 0.8),
               width = 0.6,
               outlier.shape = NA,
               alpha = 0.6) +
  geom_jitter(color = "black",
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              shape = 21,size = 1.5,alpha = 0.9,
              show.legend = FALSE) +
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral")) +
  labs(x = "Week",
       y = expression(log[10]*italic(P)[net]*" ("*mu*"mol O"[2]*" "*cm^-2*" "*h^-1*")")) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none")+
  ylim(-7.5,0.8)+
  geom_signif(annotation = "*",
              xmin=c(1.2,1.2),
              xmax=c(2.2,4.2),
              y_position=c(-0.01,0.6),
              tip_length = c(0.05,0.23,
                             0.05,0.22))

pPnet

ggsave("figures/netPhotosynthesis.pdf",
       device="pdf",
       width=8.4,
       height=7.25,
       units="cm")

# 6 GROSS PHOTOSYNTHESIS ----
## 6.1 Remove rows with negative respiration or net photosynthesis ----
mastersheetPg=mastersheet[-c(81, 94, 110, 127), ]

## 6.2 Create linear mixed effects model ----
modPg= lme(data=mastersheetPg,
           Pgross~Week*Treatment,
           random=~1|Tank/MC)

## 6.3 Check model assumptions ----
mod.plotsPg=plot_model(modPg,type="diag")
mod.plotsPg[[2]]
mod.plotsPg[[3]]

## 6.4 Return model summary and R-squared values ----
anova(modPg, type="marginal")

r.squaredGLMM(modPg)

## 6.5 Pairwise comparisons ----
emmPg2 <- emmeans(modPg, ~ Treatment * Week)

pairs(emmPg2,simple="Treatment")
pairs(emmPg2,simple="Week")
#significant difference between weeks 1 and 2, 
#and 1 and 4 in heat treated polyps

## 6.6 Summary statistics ----
Pg.mean = mastersheetPg %>%
  group_by(Treatment) %>%
  summarise(m=mean(Pgross,na.rm=TRUE),
            se=sd(Pgross,na.rm=TRUE)/sqrt(sum(!is.na(Pgross))))

## 6.7 Generate plot and save ----
pPg = ggplot(data = mastersheetPg,
              aes(x = Week,y = Pgross,fill = Treatment)) + 
  geom_boxplot(color = "black",
               position = position_dodge(width = 0.8),
               width = 0.6,outlier.shape = NA,alpha = 0.6) +
  geom_jitter(color = "black",
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              shape = 21,size = 1.5,alpha = 0.9,
              show.legend = FALSE) +
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral"),
                    labels = c("Ctrl" = "Control", "Exp" = "Heat Wave")) +
  labs(x = "Week",
       y = expression(italic(P)[gross]*" ("*mu*"mol O"[2]*" "*cm^-2*" "*h^-1*")")) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.position = "none")+
  geom_signif(annotation = "*",
              xmin=c(1.2,1.2),
              xmax=c(2.2,4.2),
              y_position=c(1,1.1))+
  ylim(0,1.13)

pPg

ggsave("figures/grossPhotosynthesis.pdf",
       device="pdf",
       width=8.4,
       height=7.25,
       units="cm")

# 7 RESPIRATION ----
## 7.1 Remove negative respiration values  ----
mastersheetR=mastersheet[-c(81, 110, 127), ]

## 7.2 Create linear mixed effects model ----
modR=lme(data=mastersheetR,
         log10(R)~Week*Treatment,
         random=~1|Tank/MC)

## 7.3 Check model assumptions ----
mod.plotsR=plot_model(modR,type="diag")  # build assumptions plots
mod.plotsR[[2]]                          # check normality
mod.plotsR[[3]]                          # check equal variance

## 7.4 Return model summary and R-squared values ----
anova.lme(modR,               # ANOVA type summary
          type="marginal")

r.squaredGLMM(modR)           # R-squared estimates

## 7.5 Summary statistics ----
R.sum = mastersheetR %>%
  summarise(m=mean(R,na.rm=TRUE),
            se=sd(R,na.rm=TRUE)/sqrt(sum(!is.na(R))))

## 7.6 Generate plot and save ----
pR = ggplot(data = mastersheetR,
            aes(x = factor(Week, levels = c("1", "2", "3", "4")),
                y = log10(R),
                fill = Treatment)) + 
  geom_boxplot(color = "black",
               position = position_dodge(width = 0.8),
               width = 0.6,
               outlier.shape = NA,
               alpha = 0.6)+
  geom_jitter(color = "black",
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              shape = 21,
              size = 1.5,
              alpha = 0.9,
              show.legend = TRUE) +
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral")) +
  labs(x = "Week",
       y = expression(log[10]*italic(R)[dark]*" ("*mu*"mol O"[2]*" "*cm^-2*" "*h^-1*")")) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position= "none")

pR

ggsave("figures/respiration.pdf",
       device="pdf",
       width=8.4,
       height=7.25,
       units="cm")



# 8 P:R RATIO ----
## 8.1 Remove rows with negative respiration or net photosynthesis ----
mastersheetPR <- mastersheet[-c(81, 94, 110, 127), ]

## 8.2 Create linear mixed effects model ----
modPR= lme(data=mastersheetPR,
           log10(PRDAILY)~Week*Treatment,
           random=~1|Tank/MC)

## 8.3 Check model assumptions ----
mod.plotsPR=plot_model(modPR,type="diag")
mod.plotsPR[[2]]
mod.plotsPR[[3]]

## 8.4 Return model summary and R-squared values ----
anova.lme(modPR, type="marginal")

r.squaredGLMM(modPR)

## 8.5 Summary statistics ----
sumPR = mastersheetPR %>%
  summarize(m=mean(PRDAILY),
            se=sd(PRDAILY)/sqrt(sum(!is.na(PRDAILY))))

## 8.6 Generate plot and save ----
pPRD <- ggplot(data = mastersheetPR,
                aes(x = Week,
                    y = log(PRDAILY),
                    fill = Treatment,
                    color = Treatment)) + 
  geom_boxplot(color="black",
               position = position_dodge(width = 0.8),
               width = 0.6,outlier.shape = NA,alpha = 0.6) +
  geom_jitter(color="black",
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
    shape = 21,size = 1.5,alpha = 0.9,
    show.legend = FALSE) +
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral")) +   
  scale_y_continuous(breaks = 1:6) +
  geom_hline(yintercept = log(1), linetype = "dashed", color = "black", size = 1) +
  labs(x = "Week", 
    y = expression(log[10]*" "*italic(P)*":"*italic(R))) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = "none")

pPRD

ggsave("figures/prratio.pdf",
       device="pdf",
       width=8.4,
       height=7.25,
       units="cm")


# 9 SYMBIONT DENSITY ----
## 9.1 Create linear mixed effects model ----
modZoox= lme(data=mastersheet,
             cells~Week*Treatment,
             random=~1|Tank/MC)
## 9.2 Check model assumptions ----
mod.plotsZoox=plot_model(modZoox,type="diag")
mod.plotsZoox[[2]]
mod.plotsZoox[[3]]

## 9.3 Return model summary and R-squared values ----
anova(modZoox, type="marginal")

r.squaredGLMM(modZoox)

## 9.4 Pairwise comparisons ----
pairs(emmeans(modZoox,~Treatment))
cld(emmeans(modZoox,~Week))

# Heat treated are significantly lower density;
# density is also higher in weeks 1 and 2 than weeks 3 and 4

## 9.5 Summary statistics ----
Zoox.sum <- as.data.frame(summary(emmeans(modZoox,~Treatment*Week))) %>%
  dplyr::mutate(Week = as.integer(as.character(Week)),
    Treatment = factor(Treatment, levels = c("Ctrl", "Exp"))) %>%
  dplyr::select(Week, Treatment, emmean, SE, df, lower.CL, upper.CL) %>%
  mutate(alt_diff = emmean/lag(emmean, n = 1),
         alt_diff=replace(alt_diff,alt_diff>1,NA)) %>%  # calculate ratio between treatments
   summarise(m=mean(alt_diff,na.rm=T))                  # return the average ratio.

## 9.6 Generate plot and save ----
pZoox=ggplot(data = mastersheet,
                 aes(x = Week,
                     y = cells,
                     fill = Treatment)) + 
  geom_boxplot(color = "black",  # set outline color manually
               position = position_dodge(width = 0.8),
               width = 0.6,outlier.shape = NA,alpha = 0.6) +
  geom_jitter(color = "black",  # outline manually
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              shape = 21,size = 1,alpha = 0.9,
              show.legend = TRUE) +         
  scale_fill_manual(values = c("Ctrl" = "cornflowerblue", "Exp" = "coral"),
                    labels = c("Ctrl" = "Control", "Exp" = "Simulated Heat Wave"),
                    name = "Treatment") + 
  labs(x = "Week", 
       y = expression("Symbiodiniacaeae Density (cells "*mm^-2*")"),
       fill = "Treatment") +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none")+
  geom_signif(annotations="*",
              xmin=c(1.5,1.5),
              xmax=c(3,4),
              y_position=c(10000,10600),
              tip_length=c(0.05,0.4,
                           0.05,0.5),
              vjust=0.3)+
  annotate("segment",
  x=c(0.8,2.8,3.8),
  xend=c(2.2,3.2,4.2),
  y=c(9500,6200,5850))

pZoox

ggsave("figures/pZoox.pdf",
       width = 13.74,
       height = 9.88,
       units="cm")



