library(tidyverse)
library(readxl)

df <- read_excel("../data/fig-data/fig2.xlsx")

df <- df %>% slice(1:8) %>% select(1:3) %>% magrittr::set_colnames(c("gene", "NPC", "EG")) 

fig2a <- df %>% pivot_longer(-gene) %>% mutate(value = as.numeric(value), gene = factor(gene, levels = df %>% mutate(EG = as.numeric(EG)) %>% arrange(desc(EG)) %>% pull(gene))) %>% ggplot(aes(x=gene, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") + scale_fill_manual(values = c("maroon", "navy")) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("") + ylab("Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "top")

fig2b <- df %>% type.convert(as.is = TRUE) %>% ggplot(aes(x=NPC, y=EG)) + geom_point(size = 2.5)+ theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("NPC Mutation Rate (%)") + ylab("EG Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "top") + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed", size = 0.8)

cowplot::plot_grid(fig2a, fig2b, nrow = 1, rel_widths = c(1.5, 1)) %>% ggsave(., filename = "../results/fig2.pdf", device = cairo_pdf, units = "in", height = 5, width = 11)

df <- df %>% type.convert(as.is = TRUE)
cor.test(df$NPC, df$EG)