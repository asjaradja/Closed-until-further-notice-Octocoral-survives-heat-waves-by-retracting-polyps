# 0 Set-Up Environment ----
library(tidyverse)  # for data wrangling
library(vegan)      # for PERMANOVA

# 1 Relative Abundance Figures ----
## 1.1 Colours ----
fullBrown=c("#241710", "#301F15", "#3D271B", "#493021", "#563826",
            "#62402C", "#6F4932", "#7B5138", "#88593D", "#956243",
            "#A16A49", "#AE724F", "#BA7B55", "#C7835B", "#D38C61",
            "#E97451")
smallBrown=c("#241710","#62402C","#A16A49","#E97451")


## 1.2 Load data ----
ITS2rel = read_csv("ITS2DataFilteredRel.csv",
                   col_types="cfffc") %>%
  mutate(Time=fct_recode(Time,
                         "T0"="0",
                         "T1"="1"),
         Result=fct_recode(Result,
                           "Necrosed"="N",
                           "Healthy"="H"))

ITS2rel_long <- ITS2rel %>%
  pivot_longer(cols = where(is.numeric),
               names_to  = "Grouping",
               values_to = "Abundance") %>%
  mutate(Grouping = fct_recode(Grouping,
                               "Symbiodinium A2ba"="A2ba",
                               "Breviolum B8a-B1gb-B8f-B8i-B1ga-B1is"= "B8a-B1gb-B8f-B8i-B1ga-B1is",
                               "Breviolum B8a-B1ir-B1iq-B19bw-B8j-B8i-B8k" = "B8a-B1ir-B1iq-B19bw-B8j-B8i-B8k",
                               "Breviolum B8a-B1is-B8f-B8l-B8i" = "B8a-B1is-B8f-B8l-B8i"),
         Grouping = fct_relevel(Grouping,"Symbiodinium A2ba"))

relperm = read_csv("RelAbundanceITS2Only.csv",
                   col_types = "fffff") %>%
  mutate(Time=fct_recode(Time,
                         "T0"="0",
                         "T1"="1"),
         Result=fct_recode(Result,
                           "Necrosed"="N",
                           "Healthy"="H"))

relperm_long <- relperm %>%
  pivot_longer(cols = where(is.numeric),
               names_to  = "Grouping",
               values_to = "Abundance")

## Long format (create ITS2abs_long)


## 1.3 Relative Abundance ITS2 Figure Stacked ----
ITS2RelStack=ggplot(ITS2rel_long,
                    aes(x=Time,y=Abundance,fill=Grouping))+
  geom_col()+
  ggh4x::facet_nested(cols = vars(Result, Time),
                      scales    = "free_x",
                      space     = "free_x",
                      nest_line = TRUE,
                      strip     = ggh4x::strip_nested(),
                      labeller  = labeller(Time = c(T0 = "Week 1", T1 = "Week 4")))+
  scale_fill_manual(values=smallBrown,
                    breaks=c("Breviolum B8a-B1gb-B8f-B8i-B1ga-B1is",
                             "Breviolum B8a-B1ir-B1iq-B19bw-B8j-B8i-B8k",
                             "Breviolum B8a-B1is-B8f-B8l-B8i",
                             "Symbiodinium A2ba"),
                    labels = c(expression(italic("Breviolum")~"B8a-B1gb-B8f-B8i-B1ga-B1is"),
                               expression(italic("Breviolum")~"B8a-B1ir-B1iq-B19bw-B8j-B8i-B8k"),
                               expression(italic("Breviolum")~"B8a-B1is-B8f-B8l-B8i"),
                               expression(italic("Symbiodinium")~A2ba)))+
  ylab("ITS2 type profile relative abundance")+
  theme_classic(base_size = 12) +
  theme(panel.spacing.x      = unit(1.2, "lines"),
        axis.text.x          = element_blank(),
        axis.ticks.length.y  = unit(2, "pt"),
        axis.ticks.length.x  = unit(0,"pt"),
        strip.background     = element_blank(),
        strip.text           = element_text(face = "bold"),
        ggh4x.facet.nestline = element_line(),
        axis.title.x          = element_blank(),  # smaller X label
        axis.title.y          = element_text(size = 8),
        legend.position="bottom",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             title.position="top"))

## Display and save plot
ITS2RelStack
ggsave("figures/ITS2RelativeITS2Stack.pdf",
       width = 16,
       height = 10,
       units = "cm")

## 1.4 Relative Abundance Stacked Figure ----
relperm = read_csv("RelAbundanceITS2Only.csv",
                   col_types = "fffff") %>%
  mutate(Time=fct_recode(Time,
                         "T0"="0",
                         "T1"="1"),
         Result=fct_recode(Result,
                           "Necrosed"="N",
                           "Healthy"="H"))

relperm_long <- relperm %>%
  pivot_longer(cols = where(is.numeric),
               names_to  = "Grouping",
               values_to = "Abundance")

pRelStack=ggplot(relperm_long,
                 aes(x=Time,y=Abundance,fill=Grouping))+
  geom_col()+
  ggh4x::facet_nested(cols = vars(Result, Time),
                      scales    = "free_x",
                      space     = "free_x",
                      nest_line = TRUE,
                      strip     = ggh4x::strip_nested(),
                      labeller  = labeller(Time = c(T0 = "Week 1", T1 = "Week 4")))+
  scale_fill_manual(values=fullBrown,
                    breaks=c("B19bw","B1ga","B1gb","B1iq","B1ir",
                             "B1is","B8","B8a","B8d","B8f",
                             "B8h","B8i","B8j","B8k","B8l",
                             "A2ba"))+
  ylab("Relative Abundance")+
  theme_classic(base_size = 12) +
  theme(panel.spacing.x      = unit(1.2, "lines"),
        axis.text.x          = element_blank(),
        axis.ticks.length.y  = unit(2, "pt"),
        axis.ticks.length.x  = unit(0,"pt"),
        strip.background     = element_blank(),
        strip.text           = element_text(face = "bold"),
        ggh4x.facet.nestline = element_line(),
        axis.title.x          = element_blank(),  # smaller X label
        axis.title.y          = element_text(size = 8),
        legend.position="bottom",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             title.position="top"))

pRelStack
ggsave("figures/pRelativeStacked.pdf",
       width = 16,
       height = 10,
       units = "cm")         
# 2 Permanova of Bray-Curtis Distances for Relative Abundance between Week 1 and Week 4 & between Healthy and Necrosed nubbins ----

### Read in data
relperm = read_csv("RelAbundanceITS2Only.csv",
                   col_types = "fffff")

### Separate metadata and community matrix
metadatarel = relperm %>%
  select(where(is.factor))

cmatrel = relperm %>%
  select(where(is.numeric)) %>%
  filter(if_any(everything(), ~ .x != 0))  # Remove taxa with all zeros

### Square-root transform (common before Bray–Curtis)
cmatrel_sqrt <- sqrt(cmatrel)

### Compute Bray–Curtis dissimilarity
bcrel <- vegdist(cmatrel_sqrt, method = "bray")

dist(cmatrel_sqrt)

### PERMANOVA (adonis2)
set.seed(123)
permanova_rel <- adonis2(bcrel ~ Time+Result, data = metadatarel, permutations = 9999, by = "margin")

### Test for homogeneity of dispersion
disprel <- betadisper(bcrel, metadatarel$Time)
disp_testrel <- permutest(disprel, permutations = 9999)

### Output results
permanova_rel
disp_testrel