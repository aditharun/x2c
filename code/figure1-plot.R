library(tidyverse)
library(readxl)

df <- "../data/fig-data/fig1.xlsx" %>% read_excel(sheet = 1) %>% slice(14:16)

df <- df %>% select(1:3)

fig1b <- df %>% magrittr::set_colnames(c("mutation", "EG", "NPC")) %>% pivot_longer(-mutation) %>% mutate(value = as.numeric(value), mutation = factor(mutation, levels = c("V600E", "V600K", "G12C"))) %>% ggplot(aes(x=mutation, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") + scale_fill_manual(values = c("maroon", "navy")) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("") + ylab("Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "top")

df <- "../data/fig-data/fig1.xlsx" %>% read_excel(sheet = 2) 

fig1c <- df %>% slice(6:8) %>% select(2:4) %>% magrittr::set_colnames(c("cancer", "Incidence", "Samples")) %>% pivot_longer(-cancer) %>% mutate(value = as.numeric(value), mutation = factor(cancer, levels = c("CRC", "LUAD", "PAAD"))) %>% ggplot(aes(x=cancer, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("") + ylab("Proportion (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "top")

fig1b %>% ggsave(., filename = "../results/fig1b.pdf", units = "in", device = cairo_pdf, height = 6, width = 4.5)

fig1c %>% ggsave(., filename = "../results/fig1c.pdf", units = "in", device = cairo_pdf, height = 6, width = 4.5)