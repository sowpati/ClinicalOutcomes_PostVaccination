library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
library(paletteer)
library(ggrepel)
library(ggdist)
library(ggsignif)
library(gghalves)
library(showtext)
library(patchwork)

################################################################################
#       PREP DATA
################################################################################
remove(PV_data)
PV_data <- read.csv('./PostVaccData_Final.tsv', header = T, sep = "\t", stringsAsFactors = F)
PV_data$Date_Dose1 <- as.Date(PV_data$Date_Dose1, "%Y-%m-%d")
PV_data$Date_Dose2 <- as.Date(PV_data$Date_Dose2, "%Y-%m-%d")

PV_data <- PV_data %>% mutate(Vacc.Group =
                                case_when(Vaccine.Type == "Covaxin" & NumDoses == 1 ~ "Covaxin Dose 1", 
                                          Vaccine.Type == "Covaxin" & NumDoses == 2 ~ "Covaxin Dose 2",
                                          Vaccine.Type == "Covishield" & NumDoses == 1 ~ "Covishield Dose 1",
                                          Vaccine.Type == "Covishield" & NumDoses == 2 ~ "Covishield Dose 2",
                                          Vaccine.Type == "Unvaccinated" & NumDoses == 0 ~ "Unvaccinated"))
PV_data$Vacc.Group <- factor(PV_data$Vacc.Group, levels = c("Unvaccinated", "Covaxin Dose 1", "Covaxin Dose 2", "Covishield Dose 1", "Covishield Dose 2"))
PV_data$Vaccine.Type <- factor(PV_data$Vaccine.Type, levels = c("Unvaccinated", "Covaxin", "Covishield"))

PV_data$Mortality[PV_data$Mortality == "DISCHARGED"] <- "Recovered"
PV_data$Mortality[PV_data$Mortality == "DEATH"] <- "Deceased"
PV_data$Mortality <- factor(PV_data$Mortality, levels = c('Recovered', 'Deceased'))

head(PV_data)
colnames(PV_data)
PV_data$SpO2 <- as.numeric(PV_data$SpO2)
PV_data$Age <- as.numeric(PV_data$Age)
PV_data$Neutralizing.Antibodies <- as.numeric(PV_data$Neutralizing.Antibodies)
PV_data$CRP <- as.numeric(PV_data$CRP)
PV_data$Ferritin <- as.numeric(PV_data$Ferritin)
PV_data$D.Dimer <- as.numeric(PV_data$D.Dimer)
PV_data$Serum.LDH <- as.numeric(PV_data$Serum.LDH)
PV_data$WBC <- as.numeric(PV_data$WBC)
PV_data$Neutrophils <- as.numeric(PV_data$Neutrophils)
PV_data$Lymphocytes <- as.numeric(PV_data$Lymphocytes)
PV_data$Platelets <- as.numeric(PV_data$Platelets)

PV_data$Platelets <- PV_data$Platelets * 100 # Convert lakhs to thousands
PV_data$Lineage[PV_data$Lineage == "None"] <- "Others"


################################################################################
#       Basic Stats
################################################################################
PV_data %>%
  count(Vaccination.Status)

PV_data %>%
  count(Vaccine.Type)

PV_data %>%
  count(NumDoses)

PV_data %>%
  count(Vaccine.Type, NumDoses)


################################################################################
#       Table 1 - Vaccinated vs Unvaccinated
################################################################################

# Basic Stats - continuous
Table1_Basic_Continuous <- PV_data %>%
  melt(id.vars = "Vaccination.Status", measure.vars = c(
    "Age",
    "CRP",
    "Ferritin",
    "D.Dimer",
    "Serum.LDH",
    "WBC",
    "Neutrophils",
    "Lymphocytes",
    "Platelets",
    "Neutralizing.Antibodies"
  )) %>%
  group_by(Vaccination.Status, variable) %>%
  filter(!is.na(value)) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    count = n()
  ) %>%
  arrange(Vaccination.Status) %>%
  mutate(concat = paste(round(mean, 2), " ± ", round(sd, 2), " (n = ", count, ")", sep = ""))
write.table(Table1_Basic_Continuous, file="Table1_BasicStats_Continuous.tsv", quote = F, row.names = F, sep = "\t")

# Significance testing - Table 1 Continuous
t.test(Age ~ Vaccination.Status, data = PV_data)
t.test(CRP ~ Vaccination.Status, data = PV_data)
t.test(Ferritin ~ Vaccination.Status, data = PV_data)
t.test(D.Dimer ~ Vaccination.Status, data = PV_data)
t.test(Serum.LDH ~ Vaccination.Status, data = PV_data)
t.test(WBC ~ Vaccination.Status, data = PV_data)
t.test(Neutrophils ~ Vaccination.Status, data = PV_data)
t.test(Lymphocytes ~ Vaccination.Status, data = PV_data)
t.test(Platelets ~ Vaccination.Status, data = PV_data)
t.test(Neutralizing.Antibodies ~ Vaccination.Status, data = PV_data)

# Basic Stats - categorical
colnames(PV_data)
PV_data %>%
  group_by(Vaccination.Status) %>%
  count(Gender) %>%
  mutate(per = prop.table(n) * 100)

Table1_Basic_Categorical <- PV_data %>%
  melt(id.vars = "Vaccination.Status", measure.vars = c(
    22:28, # Comorbidities
    40:45 # Hospital parameters
  )) %>%
  group_by(Vaccination.Status, variable) %>%
  count(value) %>%
  mutate(per = prop.table(n) * 100) %>%
  filter(value == "Y") %>%
  mutate(concat = paste(n, " (", round(per, 2), "%)", sep = ""))
write.table(Table1_Basic_Categorical, file="Table1_BasicStats_Categorical.tsv", quote = F, row.names = F, sep = "\t")

PV_data %>%
  group_by(Vaccination.Status) %>%
  count(Severity == "Severe" | ICU.need.AT.Admission == "Y") %>%
  mutate(per = prop.table(n) * 100)

PV_data %>%
  group_by(Vaccination.Status, NumDoses) %>%
  count(Mortality) %>%
  mutate(per = prop.table(n) * 100)


# Significance testing - Table 1 Categorical

cont_table <- table(PV_data$Gender == "F", PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)  
  
cont_table <- table(PV_data$DM.HTN, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)  
  
cont_table <- table(PV_data$Hypothyroidism, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)
  
cont_table <- table(PV_data$KD, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$CLD, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$Respiratory, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$CD, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$Malignancy, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$Thrombotic.Complications, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$AKI, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$RRT, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$Mechanical.ventilation, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$Severity == "Severe" | PV_data$ICU.need.AT.Admission == "Y", PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

cont_table <- table(PV_data$ICU.admission.IN.HOSPITAL.STAY, PV_data$Vaccination.Status)
cont_table
fisher.test(cont_table)

# Mortality dose wise
filtered <- PV_data %>% 
  filter(NumDoses < 2) %>% 
  group_by(Mortality, Vaccination.Status)
cont_table <- table(filtered$Mortality, filtered$Vaccination.Status)
cont_table
fisher.test(cont_table)

filtered <- PV_data %>% 
  filter(NumDoses == 0 | NumDoses == 2) %>% 
  group_by(Mortality, Vaccination.Status)
cont_table <- table(filtered$Mortality, filtered$Vaccination.Status)
cont_table
fisher.test(cont_table)


################################################################################
#       Table S1 - Covaxin vs Covishield (2 doses)
################################################################################
CovaxinVsCovishield <- filter(PV_data, Vaccine.Type == "Covaxin" | Vaccine.Type == "Covishield") %>%
  filter(NumDoses == 2)
CovaxinVsCovishield$Vaccine.Type <- factor(CovaxinVsCovishield$Vaccine.Type, levels = c("Covaxin", "Covishield"))
count(CovaxinVsCovishield, Vaccine.Type)

# Basic Stats - continuous
TableS1_Basic_Continuous <- CovaxinVsCovishield %>%
  melt(id.vars = "Vaccine.Type", measure.vars = c(
    "CRP",
    "Ferritin",
    "D.Dimer",
    "Serum.LDH",
    "WBC",
    "Neutrophils",
    "Lymphocytes",
    "Platelets",
    "Neutralizing.Antibodies"
  )) %>%
  group_by(Vaccine.Type, variable) %>%
  filter(!is.na(value)) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    count = n()
  ) %>%
  arrange(Vaccine.Type) %>%
  mutate(concat = paste(round(mean, 2), " ± ", round(sd, 2), " (n = ", count, ")", sep = ""))
TableS1_Basic_Continuous
write.table(TableS1_Basic_Continuous, file="TableS1_BasicStats_Continuous.tsv", quote = F, row.names = F, sep = "\t")

# Significance testing - Table S1 Continuous
t.test(CRP ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Ferritin ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(D.Dimer ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Serum.LDH ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(WBC ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Neutrophils ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Lymphocytes ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Platelets ~ Vaccine.Type, data = CovaxinVsCovishield)
t.test(Neutralizing.Antibodies ~ Vaccine.Type, data = CovaxinVsCovishield)

# Basic Stats - categorical
TableS1_Basic_Categorical <- CovaxinVsCovishield %>%
  melt(id.vars = "Vaccine.Type", measure.vars = c(
    22:28, # Comorbidities
    40:45 # Hospital parameters
  )) %>%
  group_by(Vaccine.Type, variable) %>%
  count(value) %>%
  mutate(per = prop.table(n) * 100) %>%
  filter(value == "Y") %>%
  mutate(concat = paste(n, " (", round(per, 2), "%)", sep = ""))
write.table(TableS1_Basic_Categorical, file="TableS1_BasicStats_Categorical.tsv", quote = F, row.names = F, sep = "\t")

CovaxinVsCovishield %>%
  group_by(Vaccine.Type) %>%
  count(Severity) %>%
  mutate(per = prop.table(n) * 100)

CovaxinVsCovishield %>%
  group_by(Vaccine.Type) %>%
  count(Severity == "Severe" | ICU.need.AT.Admission == "Y") %>%
  mutate(per = prop.table(n) * 100)

CovaxinVsCovishield %>%
  group_by(Vaccine.Type) %>%
  count(Mortality) %>%
  mutate(per = prop.table(n) * 100)

# Significance testing - Table S1 Categorical
cont_table <- table(CovaxinVsCovishield$DM.HTN, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)  

cont_table <- table(CovaxinVsCovishield$Hypothyroidism, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$KD, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$CLD, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Respiratory, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$CD, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Malignancy, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Thrombotic.Complications, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$AKI, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$RRT, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Mechanical.ventilation, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Severity == "Severe", CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$ICU.need.AT.Admission, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Severity == "Severe" | CovaxinVsCovishield$ICU.need.AT.Admission == "Y", CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$ICU.admission.IN.HOSPITAL.STAY, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)

cont_table <- table(CovaxinVsCovishield$Mortality, CovaxinVsCovishield$Vaccine.Type)
cont_table
fisher.test(cont_table)


################################################################################
#       Pairwise Tests (with Holm correction)
################################################################################

pairwise.t.test(PV_data$CRP, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Ferritin, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$D.Dimer, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Serum.LDH, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$WBC, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Neutrophils, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Lymphocytes, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Platelets, PV_data$Vacc.Group, p.adjust.method = "holm")
pairwise.t.test(PV_data$Neutralizing.Antibodies, PV_data$Vacc.Group, p.adjust.method = "holm")


#####################################################################################
#      Configure theme of ggplot
#####################################################################################
font_add_google("Roboto", "roboto")
font_add_google("Noto Sans JP", "notosans")
font_add_google("Lato", "lato")

my_theme <- theme(
  text = element_text(family = "lato", size = 12),
  axis.text = element_text(size = 12),
  axis.title = element_text(size=14),
  legend.text = element_text(size=14),
  legend.title = element_text(size=16),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank()
)

#####################################################################################
#       Figure 1A - Lineages in vacc/unvacc
#####################################################################################
lineage_totals <-
  PV_data %>%
  filter(Sequenced == "Y") %>%
  group_by(Vaccination.Status) %>%
  summarise(total = n())
lineage_totals

Fig1A <- PV_data %>%
  filter(Sequenced == "Y") %>%
  ggplot(aes(x = Vaccination.Status, fill = Lineage)) +
  geom_bar(position = "fill") +
  geom_text(inherit.aes = F, aes(x = Vaccination.Status, y = 1, label = total), data = lineage_totals, vjust=-0.25, size = 5) +
  theme_minimal() +
  my_theme +
  # theme(axis.text.x = element_text(size = 10)) +
  labs(x='Vaccination Status', y='Proportion')
ggsave('./Fig1A_Lineages_StackedBar.png', plot = Fig1A, width=5, height=5, dpi=300, bg = "white")


#####################################################################################
#       Figure 1B - NAb values dose and type wise
#####################################################################################
Fig1B_signif <- data.frame(
  start = c(rep("Unvaccinated", 4)),
  end = c("Covaxin Dose 1", "Covaxin Dose 2", "Covishield Dose 1", "Covishield Dose 2"),
  label = c("n.s", "n.s", "***", "***"),
  y = rev(c(4.4, 4.1, 3.8, 3.5))
)

Fig1B <- PV_data %>%
  filter(!is.na(Vacc.Group)) %>%
  ggplot(aes(x = Vacc.Group, y = Neutralizing.Antibodies)) +
  geom_half_boxplot(aes(fill= Vacc.Group), outlier.shape = NA, nudge = 0.05) +
  geom_half_point(alpha = 0.75, aes(color = Vacc.Group)) +
  geom_signif(
    data = Fig1B_signif,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = T
  ) +
  scale_fill_paletteer_d("ggsci::default_jama", name='Category') +
  scale_color_paletteer_d("ggsci::default_jama", name='Category') +
  scale_y_log10() +
  theme_minimal() +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x='Vaccination Status', y='Neutralizing Ab levels\n(AU/ml, log scale)')
ggsave('./Fig1B_NAbLevels_DoseWise.png', plot = Fig1B, width=7, height=5, dpi=300, bg = "white")



#####################################################################################
#       Figure 1C - Ferritin levels
#####################################################################################
Fig1C_signif <- data.frame(
  start = c(rep("Unvaccinated", 4)),
  end = c("Covaxin Dose 1", "Covaxin Dose 2", "Covishield Dose 1", "Covishield Dose 2"),
  label = c("n.s", "**", "n.s", "***"),
  y = c(5000, 5300, 5600, 5900)
)

Fig1C <- PV_data %>%
  filter(!is.na(Vacc.Group)) %>%
  ggplot(aes(x = Vacc.Group, y = Ferritin)) +
  geom_half_boxplot(aes(fill= Vacc.Group), outlier.shape = NA, nudge = 0.05) +
  geom_half_point(alpha = 0.75, aes(color = Vacc.Group)) +
  geom_signif(
    data = Fig1C_signif,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = T
  ) +
  scale_fill_paletteer_d("ggsci::default_jama", name='Category') +
  scale_color_paletteer_d("ggsci::default_jama", name='Category') +
  theme_minimal() +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x='Vaccination Status', y='Ferritin (ng/ml)')
ggsave('./Fig1C_FerritinLevels_DoseWise.png', plot = Fig1C, width=7, height=5, dpi=300, bg = "white")




#####################################################################################
#       Figure 1D - LDH levels
#####################################################################################
Fig1D_signif <- data.frame(
  start = c(rep("Unvaccinated", 4)),
  end = c("Covaxin Dose 1", "Covaxin Dose 2", "Covishield Dose 1", "Covishield Dose 2"),
  label = c("n.s", "**", "n.s", "*"),
  y = c(1900,2100,2300,2500)
)

Fig1D <- PV_data %>%
  filter(!is.na(Vacc.Group)) %>%
  ggplot(aes(x = Vacc.Group, y = Serum.LDH)) +
  geom_half_boxplot(aes(fill= Vacc.Group), outlier.shape = NA, nudge = 0.05) +
  geom_half_point(alpha = 0.75, aes(color = Vacc.Group)) +
  geom_signif(
    data = Fig1D_signif,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = T
  ) +
  scale_fill_paletteer_d("ggsci::default_jama", name='Category') +
  scale_color_paletteer_d("ggsci::default_jama", name='Category') +
  scale_y_continuous(limits = c(0, 2500)) +
  theme_minimal() +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x='Vaccination Status', y='Serum LDH (Units/L)')
ggsave('./Fig1D_SerumLDHLevels_DoseWise.png', plot = Fig1D, width=7, height=5, dpi=300, bg = "white")


#####################################################################################
#       Figure 1E - NAb Levels in Vaccinated by Mortality
#####################################################################################

Fig1E <- PV_data %>%
  # filter(Vaccination.Status == "Vaccinated", NumDoses > 0) %>%
  filter(!is.na(Vacc.Group)) %>%
  ggplot(aes(x = as.character(NumDoses), y = Neutralizing.Antibodies)) +
  # geom_boxplot(aes(fill = Mortality)) +
  geom_half_boxplot(aes(fill = Mortality), nudge = 0.05, outlier.shape = NA) +
  geom_half_point(aes(color = Mortality)) +
  # geom_jitter(aes(color = as.character(NumDoses))) +
  scale_fill_paletteer_d("ggsci::default_aaas", name='Mortality') +
  scale_color_paletteer_d("ggsci::default_aaas", name='Mortality') +
  # scale_y_continuous(limits = c(0, 4500)) +
  scale_y_log10() +
  theme_minimal() +
  my_theme +
  labs(x='Number of doses', y='Neutralizing Ab levels\n(AU/ml, log scale)')
ggsave('./Fig1E_NAb_Mortality.png', plot = Fig1E, width=7, height=5, dpi=300, bg = "white")



#####################################################################################
#       Figure 1 Final with layout
#####################################################################################

Fig1_Layout <- '
AABBB
AABBB
AACCC
EECCC
EEDDD
EEDDD
'
wrap_plots(A = Fig1A, B = Fig1B, C = Fig1C, D = Fig1D, E = Fig1E, design = Fig1_Layout) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20))
ggsave('Fig1.png', width=18, height=12, dpi=300, bg = "white")



#####################################################################################
#       Figure S1 - Community prevalence of lineages
#####################################################################################
lineages.TG <- read.table("./lineages.TG.txt", header = T, sep = "\t")
lineages.TG %>%
  mutate(csum = rev(cumsum(rev(Count))), 
         pos = Count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Count/2, pos)) %>%
  ggplot(aes(x = "" , y = Count, fill = Lineage)) +
  geom_col(width =1 , color = c("#F8766D","#00BF7D","#00B0F6","#E76BF3")) +
  scale_fill_manual(values = c("#F8766D","#00BF7D","#00B0F6","#E76BF3")) + 
  coord_polar(theta = "y") +
  geom_label_repel(aes(x = 1.5,y = pos, 
                       label = paste(round(Count/792*100, 2), "%", sep = "")),
                   size = 4.5, nudge_x = 0.2, show.legend = F) +
  guides(fill = guide_legend(title = "Lineage",title.theme = element_text(size = 15),label.theme = element_text(size = 15))) +
  theme_void()
ggsave('FigS1.png', dpi=300, bg="white")



#####################################################################################
#       Figure S2 - NAb and Inflammatory markers - Vacc vs Unvacc
#####################################################################################
colnames(PV_data)
PV_melted_markers <- PV_data %>%
  filter(!is.na(Vaccination.Status)) %>%
  melt(id.vars = c('Vaccination.Status', 'Mortality'), measure.vars = c(5,31:48))
head(PV_melted_markers)

PV_melted_markers$value <- as.numeric(PV_melted_markers$value)
PV_melted_markers$value[PV_melted_markers$variable == "Neutralizing.Antibodies"] <- 
  log10(PV_melted_markers$value[PV_melted_markers$variable == "Neutralizing.Antibodies"]) # To plot in log scale in a facet wrap

PV_melted_markers %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = Vaccination.Status, y = value, fill=Vaccination.Status)) +
  geom_half_boxplot(outlier.shape = NA, nudge = 0) +
  # geom_violin() +
  geom_half_point(alpha=0.75, aes(color = Vaccination.Status)) +
  facet_wrap(~variable, nrow = 3, scales = "free", strip.position = "left", labeller = as_labeller(c(
    Neutralizing.Antibodies = "Log Neutralizing Ab \nlevels (AU/ml)",
    CRP = "CRP (mg/L)",
    Ferritin = "Ferritin (ng/ml)",
    Serum.LDH = "Serum LDH (Units/L)",
    D.Dimer = "D-Dimer (ng/ml)",
    WBC = "WBCs (cells/cmm)",
    Neutrophils = "Neutrophils (cells/cmm)",
    Lymphocytes = "Lymphocytes (cells/cmm)",
    Platelets = "Platelets\n(thousands/cmm)"
  ))) + 
  scale_fill_paletteer_d("ggsci::default_nejm", name='Category') +
  scale_color_paletteer_d("ggsci::default_nejm", name='Category') +
  # scale_y_log10() +
  theme_minimal() +
  my_theme +
  theme(legend.position = "top",
        strip.text = element_text(size = 18),
        strip.placement = "outside",
        axis.title = element_blank()
  )
ggsave('FigS2.png', width=12, height=10, dpi=300, bg="white")


#####################################################################################
#       Figure S3 - NAb levels Gender wise
#####################################################################################

PV_data %>%
  filter(!is.na(Neutralizing.Antibodies)) %>%
  ggplot(aes(x = Vaccination.Status, y = Neutralizing.Antibodies)) +
  geom_half_boxplot(aes(fill= Vaccination.Status), outlier.shape = NA, nudge = 0.05) +
  geom_half_point(alpha = 0.75, aes(color = Vaccination.Status)) +
  facet_wrap(~Gender, strip.position = "bottom", labeller = as_labeller(
    c(
      F = "Female",
      M = "Male"
    )
  )) +
  scale_fill_paletteer_d("ggsci::default_nejm", name='Category') +
  scale_color_paletteer_d("ggsci::default_nejm", name='Category') +
  scale_y_log10() +
  theme_minimal() +
  my_theme +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "top",
        strip.text = element_text(size = 18),
        strip.placement = "outside",
        axis.title.x = element_blank()
  ) +
  labs(y = 'Neutralizing Ab levels\n(AU/ml, log scale)')
ggsave('FigS3.png', width=7, height=5, dpi=300, bg="white")



#####################################################################################
#       Figure S4 - Inflammatory markers - Dose and Type wise
#####################################################################################

colnames(PV_data)
PV_melted_markers <- PV_data %>%
  filter(!is.na(Vacc.Group)) %>%
  melt(id.vars = c('Vacc.Group', 'Mortality'), measure.vars = c(31, 33, 35:38))
head(PV_melted_markers)

PV_melted_markers$value <- as.numeric(PV_melted_markers$value)

PV_melted_markers %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = Vacc.Group, y = value, fill=Vacc.Group)) +
  geom_half_boxplot(outlier.shape = NA, nudge = 0) +
  geom_half_point(alpha=0.75, aes(color = Vacc.Group)) +
  facet_wrap(~variable, nrow = 3, scales = "free", strip.position = "left", labeller = as_labeller(c(
    CRP = "CRP (mg/L)",
    D.Dimer = "D-Dimer (ng/ml)",
    WBC = "WBCs (cells/cmm)",
    Neutrophils = "Neutrophils (cells/cmm)",
    Lymphocytes = "Lymphocytes (cells/cmm)",
    Platelets = "Platelets\n(thousands/cmm)"
  ))) + 
  scale_fill_paletteer_d("ggsci::default_jama", name='Category') +
  scale_color_paletteer_d("ggsci::default_jama", name='Category') +
  theme_minimal() +
  my_theme +
  theme(legend.position = "top",
        strip.text = element_text(size = 18),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_blank()
  )
ggsave('FigS4.png', width=12, height=12, dpi=300, bg="white")

