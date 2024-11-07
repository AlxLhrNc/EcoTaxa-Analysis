##### Thesis 2 - Diversity
# 10/03/2023 - A.LHERIAU-NICE

# Setup
#####
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(lubridate)
library(openxlsx)
library(scales)
library(stringr)
library(tidyverse)
library(thematic)

rm(list = ls())
DIR_ALL <- dirname(rstudioapi::getSourceEditorContext()$path)
#####

# Functions
#####
library(RColorBrewer)
color.mkr <- function(pal, n){
  set.clr <- brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  return(set.clr[c(1:n)])
}

factor2num <- function(x){
  return(as.numeric(as.character(x)))
}

gamma_label <- function(.data){
  g_labels = as_labeller(
    setNames(
      .data %>% select(gamma.d, Date) %>% distinct() %>% mutate(label = paste0(Date, " (\u03b3 = ", gamma.d, ")")) %>% select(label) %>% unlist(),
      .data %>% select(gamma.d, Date) %>% distinct() %>% select(Date) %>% unlist()
      )
  )
  return(g_labels)
}

#####

# Preprocess
#####
# PHYTO
setwd(DIR_ALL)
TSV.FILE <- list.files(pattern = ".tsv")
tsv.raw <- read.csv(TSV.FILE, sep = "\t")

phyto.raw <- tsv.raw %>%
  filter(object_annotation_status == "validated") %>% 
  mutate(
    Date.UTC = ymd_hms(paste(object_date, object_time)),
    Date = factor(format(date(with_tz(Date.UTC, "Pacific/Auckland")), "%d %b %y"),
                  levels = c("04 Nov 22", "12 Jan 23", "13 Jan 23", "05 Apr 23", "06 Apr 23")),
    Sample = factor(str_sub(sample_id, start = -2)),

    #longitude = object_lon,
    class1 = str_match(object_annotation_hierarchy, "^(?:[A-z]+>){0,6}([A-z]*)")[,2], #"^(?:[A-z]+>){,6}([A-z]*)", 
    name = gsub("^.*>", "", object_annotation_hierarchy),
    name = case_when(name == "Protoperidinium(Oceanica)" ~ "Protoperidinium oceanica", 
                     name %in% c("part", "diatom centric") ~ class1,
                     .default = name),
    class2 = gsub(" .*$", "", name),
    ) %>%
  select(Date, Sample, name, class1, class2) %>%
  group_by(Date, Sample, name, class1, class2) %>%
  mutate(count = n()) %>%
  distinct()
#####

# MERGE %>% PURGE
#####
DATA <- phyto.raw %>%
  filter(!grepl("^[[:lower:]]", name)) %>%
  mutate(class2 = case_when(class2 == "Diatoms" ~ "Bacillariophyta", 
                            .default = class2),
         class1 = case_when(class1 %in% c("Bacillariophytina", "Diatoms", "Coscinodiscophytina") ~ "Bacillariophyta",
                            class1 %in% c("Protoperidinium") ~ "Dinophyceae",
                            .default = class1)
         ) %>%
  ungroup() %>%
  mutate(across(c(class1, class2, name), ~trimws(.x))) %>%
  arrange(class1, class2) %>%
  mutate(across(c(class1, class2, name), fct_inorder)) %>%
  filter(Date != "04 Nov 22")

levels(DATA$class1) <- c(levels(DATA$class1), "Others")
levels(DATA$class2) <- c(levels(DATA$class2), "Others")
levels(DATA$name) <- c(levels(DATA$name), "Others")
#####

# INDEXES (alpha, beta-sor, gamma-day, gamma)
#####

# BETA
tmp <- DATA %>%
  nest_by(Date) %>%
  mutate(data.wide = list(data %>%
                            select(Sample, name, count) %>%
                            pivot_wider(values_from = count) %>%
                            replace(is.na(.), 0) %>%
                            replace(.>1, 1)
                          )
         )
for(i in c(1:nrow(tmp))){
  beta.tbl <- data.frame("Sample"=NA, "beta.SIM"=NA,"beta.SNE"=NA,"beta.SOR"=NA)
  subdf <- tmp$data.wide[[i]]
  for(j in subdf$Sample){
    if(j == "01"){
      beta.res = betapart::beta.multi(subdf[c(1,2), -1])
    } else if(j == "16"){
      beta.res = betapart::beta.multi(subdf[c(15,16), -1])
    } else {
      j2 = factor2num(j)
      beta.res = betapart::beta.multi(subdf[c(j2-1,j2,j2+1), -1])
    }
    beta.tbl <- beta.tbl %>% rbind(data.frame(beta.res) %>% mutate(Sample = j))
  }
  tmp[i, "beta.tbl"] <- list(list(beta.tbl))
}

DATA <- tmp %>%
  mutate(full.data = list(data %>% left_join(beta.tbl, by = "Sample"))
         ) %>%
  select(-c(data, data.wide, beta.tbl)) %>%
  unnest(full.data)
rm(tmp)

# ALPHA + GAMMA.D + GAMMA
DATA <- DATA %>%
  filter(count != 0) %>%
  group_by(Date, Sample) %>%
  mutate(alpha = length(unique(name))
         ) %>%
  ungroup %>%
  group_by(Date) %>%
  mutate(gamma.d = length(unique(name))
  ) %>%
  mutate(gamma = length(unique(name))
  ) 
# Save the indexes
DATA %>% ungroup() %>%
  nest_by(Date) %>%
  mutate(data.wide = list(data %>%
                            select(-class1, -class2) %>%
                            pivot_wider(values_from = count) %>%
                            replace(is.na(.), 0) %>%
                            mutate(Sample = factor(Sample))
                          ),
         sheetName = paste(as.character(origin), as.character(Date))
         ) %>%
  ungroup() %>%
  select(sheetName, data.wide) %>%
  deframe() %>%
  write.xlsx(., file = paste(DIR_ALL, "diversity_results.xlsx", sep="/"))
  #purrrlyr::by_row(~write.xlsx(.$data.wide, file = paste(DIR_ALL, "diversity_results.xlsx", sep="/"), sheetName = .$sheetName, append = .$append))

#####

# SINGLE FIGURES
#####
# Visual thresholds
axfact <- 20
perctresh <- 0.01


# Important pre-shot
DATA <- DATA %>%
  group_by(Date, Sample) %>%
  mutate(perc = count / sum(count),
         class1 = case_when(perc<=perctresh ~ "Others", .default = class1),
         class2 = case_when(perc<=perctresh ~ "Others", .default = class2),
         name = case_when(perc<=perctresh ~ "Others", .default = name),
  ) %>%
  ungroup() %>%
  arrange(class1, class2) %>%
  mutate(across(c(class1, class2, name), fct_inorder)) %>%
  mutate(across(c(class1, class2, name), ~forcats::fct_relevel(., "Others", after = Inf)))


colours_Phyto <- setNames(color.mkr("Paired", DATA %>% distinct(class1) %>% nrow()),
                          DATA %>% arrange(class1) %>% distinct(class1) %>% unlist())

alpha_Phyto <- rep(1, DATA %>% distinct(class1, class2) %>% nrow())

#visuals
figs.bar <- DATA %>%
  ggplot(., aes(x = Sample, y = count, fill = class1)) +
                      #geom_bar(stat = "identity", position = "fill") +
                      geom_bar(stat = "identity", position = "fill")+
                      scale_fill_manual(values = colours_Phyto) +
                      scale_alpha_manual(values = alpha_Phyto) +
                      scale_y_continuous(labels = label_percent(accuracy = 1)) +
                      theme_minimal()+
                      theme(#plot.background = element_rect("white"),
                            legend.box = "horizontal") +
                      facet_grid(.~Date, labeller = gamma_label(DATA))+
                      labs(x="Station", y="Percentage", fill = "ID", alpha = "ID")
                    
figs.lin <- DATA %>%
  ggplot(., aes(x = Sample, group = Date), linewidth = 1.5) +
                      geom_line(aes(y = alpha, color = "\u03b1")) +
                      #geom_line(aes(y = beta.SOR * axfact, color = "\u03b2")) +
                      geom_line(aes(y = beta.SOR * axfact, color = "\u03b2", linetype = "Sorensen")) +
                      geom_line(aes(y = beta.SNE * axfact, color = "\u03b2", linetype = "Nestedness")) +
                      geom_line(aes(y = beta.SIM * axfact, color = "\u03b2", linetype = "Turnover")) +
                      scale_color_manual(values = c("forestgreen", "darkorchid3"))+
                      scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted"), c("Sorensen","Nestedness","Turnover")))+
                      theme_minimal()+
                      theme(#plot.background = element_rect("white"),
                            axis.title.y.left = element_text(colour = "forestgreen"),
                            axis.text.y.left = element_text(colour = "forestgreen"),
                            axis.title.y.right = element_text(colour = "darkorchid3"),
                            axis.text.y.right = element_text(colour = "darkorchid3"),
                      )+
                      scale_y_continuous(breaks = seq(0,25,5), name = "\u03b1 diversity",
                                         sec.axis = sec_axis(~ . / axfact,
                                                             name = "\u03b2 diversity",
                                                             breaks = seq(0,1.25,.25),
                                                             labels = seq(0,1.25,.25)
                                         )
                      )+
                      facet_grid(.~Date, scales = "free_y", labeller = gamma_label(DATA))+
                      labs(x="Station", y="Index value", color = "Diversity", linetype = "\u03b2 type")
    
  
#plots.table$figs.lin[[2]]

#Saves
#Barplot
ggsave(file = paste(DIR_ALL, "barplot.png", sep = "/"),
         plot = figs.bar +
           theme(plot.background = element_rect("white")),
         width = 16, height = 4, dpi = 300)

#Lineplot
ggsave(file = paste(DIR_ALL, "lineplot.png", sep = "/"),
         plot = figs.lin +
           theme(plot.background = element_rect("white")),
         width = 16, height = 4, dpi = 300)

#Both
ggarrange(figs.bar + theme(plot.margin = margin(r = 15, unit = "pt")),
          figs.lin + theme(plot.margin = margin(l = 15, unit = "pt")),
          labels = c("A", "B"),
          legend = "right", common.legend = FALSE,
          ncol = 1, nrow = 2) %>%
  ggsave(file = paste(DIR_ALL, "Full_results.png", sep = "/"), plot = .,
         width = 16, height = 8, dpi = 300)



#####


beep(0)

