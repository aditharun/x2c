library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../results"
if (!dir.exists(outdir)){
	dir.create(outdir)
}

tcga <- "../data/processed-data/Genomics_Output_Processed_specific_Cys.txt" %>% read_table()

#taken from supplementary material of stites 2021
denominator <- "../data/processed-data/Genomics_Output_Processed.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")

denom <- denominator %>% select(intersect(colnames(denominator), colnames(tcga))) %>% relocate(colnames(tcga))

#the total number of cases with at least 1 X2C mutation
total_cases <- tcga %>% filter(Hugo_Symbol == "Total")


#the total number of cases sequenced in all
tcga <- rbind(tcga %>% filter(Hugo_Symbol != "Total"), denom)

tcga_total_by_rosetta <- tcga %>% filter(Hugo_Symbol == "Total")
tcga_total_by_mut <- tcga %>% select(All)

tcga.df <- tcga %>% filter(Hugo_Symbol != "Total") %>% select(-All)

total_cases$Hugo_Symbol <- "cases with >= 1 acquired cysteine"
denom$Hugo_Symbol <- "total sequenced cases"
colnames(denom)[1] <- "status"
colnames(total_cases)[1] <- "status"
rbind(total_cases, denom) %>% as_tibble() %>% writexl::write_xlsx(., file.path(outdir, "cases-by-rosetta.xlsx"))


seer.df <- "../data/processed-data/SuppTable1_ROSETTA_Abundance_added12.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

total <- tcga_total_by_rosetta %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 

#renormalize weight vector to 1, adds to 0.93 now

total <- total %>% mutate(incidence = incidence / sum(incidence), incidence_frac = incidence_frac / sum(incidence_frac))

total <- total %>% magrittr::set_colnames(c("hugo", "all", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

tcga.muts <- tcga.df %>% pivot_longer(-Hugo_Symbol) %>% left_join(seer.df, by=c("name"="rosetta")) %>% distinct()

#create column with gene name

tcga.muts <- tcga.muts %>% mutate(gene = gsub("\\.(.*)", "", Hugo_Symbol)) %>% magrittr::set_colnames(c("hugo", "rosetta", "count", "cancer", "incidence", "incidence_frac", "gene"))

#proportion of cancer cases with each specific X2C mutation

tcga.gene <- tcga.muts %>% group_by(gene, cancer) %>% summarize(net_count=sum(count), incidence_frac = unique(incidence_frac), gene = unique(gene), rosetta = unique(rosetta)) %>% select(rosetta, incidence_frac, gene, net_count)


tcga.gene <- tcga.gene %>% mutate(weight.gene = incidence_frac * net_count)

total <- total %>% mutate(weight.tot = incidence_frac * count)

tcga.muts <- tcga.muts %>% mutate(weight.mut=incidence_frac * count) 

tcga.muts <- tcga.muts %>% mutate(idx = 1:n())

tcga.gene <- tcga.gene %>% ungroup() %>% mutate(idx = 1:n())


#perform confidence interval simulations

conf_int <- function(lambda, n=2000){

	values <- rpois(n=n, lambda=lambda)

	bounds <- quantile(values, c(0.025, 0.975)) %>% unname()

	return(list(lb = bounds[1], ub = bounds[2]))

}

tcga.muts.confint <- tcga.muts %>% filter(count > 0)  %>% mutate(result = purrr::map(count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)

tcga.muts <- tcga.muts %>% left_join(tcga.muts.confint %>% select(idx, lb, ub), by=c("idx"="idx"))

tcga.gene.confint <- tcga.gene %>% filter(net_count > 0)  %>% mutate(result = purrr::map(net_count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)

tcga.gene <- tcga.gene %>% left_join(tcga.gene.confint %>% select(idx, lb, ub), by=c("idx"="idx"))

tcga.gene <- tcga.gene %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)

tcga.muts <- tcga.muts %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)

tcga.gene <- tcga.gene %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))

tcga.muts <- tcga.muts %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))

# % of patients with each specific mutation across cancer types in US pop and TCGA

mut <- tcga.muts %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(hugo) %>% summarize(gene=unique(gene), pct_us = sum(weight.mut) / sum(weight.tot), pct_tcga = sum(count) / unique(all) , pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot)) %>% mutate(across(where(is.double), ~ . * 100))

# % of patients with mutations in a given gene across cancer types in US pop and TCGA

gene <- tcga.gene %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(gene) %>% summarize(pct_us = sum(weight.gene) / sum(weight.tot), pct_tcga = sum(net_count) / unique(all), pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot) ) %>% mutate(across(where(is.double), ~ . * 100))


all.muts <- mut %>% mutate(diff=(pct_us - pct_tcga)) %>% arrange(desc(abs(diff))) %>% mutate(gene = str_sub(gene, 1, -2)) %>% mutate(aa_change = gsub("(.*)\\.", "", hugo)) %>% mutate(mutation = paste0(gene, ".", aa_change)) %>% select(-hugo) %>% relocate(mutation, gene, aa_change, pct_tcga, pct_us, pct_lb, pct_ub, diff) %>% magrittr::set_colnames(c("mutation", "gene", "amino acid change", "Unweighted pan-cancer mutation rate (%)", "US population mutation rate (%)", "lower bound mutation rate in US pop (%)", "upper bound mutation rate in US pop (%)", "US - Unweighted pan-cancer Mutation Rate"))


all.genes <- gene %>% mutate(diff=(pct_us - pct_tcga)) %>% arrange(desc(abs(diff))) %>% mutate(gene = str_sub(gene, 1, -2)) %>% relocate(gene, pct_tcga, pct_us, pct_lb, pct_ub) %>% magrittr::set_colnames(c("gene", "Unweighted pan-cancer mutation rate (%)", "US population mutation rate (%)", "lower bound mutation rate in US pop (%)", "upper bound mutation rate in US pop (%)", "US - Unweighted pan-cancer Mutation Rate"))


x2c_points_mut <- mut %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% mutate(gene_fix = str_sub(gene, 1, -2), mut_fix = paste0(" (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")")) %>% mutate(label_text = paste0(gene_fix, mut_fix) ) %>% mutate(label = ifelse(us_rank <= 7, label_text, NA)) %>% select(-label_text) %>% filter(!is.na(label)) %>% mutate(position = ifelse(us_rank == 5, "right", "left")) %>% mutate(nx = ifelse(us_rank == 5, -0.02, 0.02)) %>% mutate(ny = ifelse(us_rank == 7, -0.01, 0.02)) %>% mutate(label = ifelse(us_rank == 5, paste0(gene_fix, "\n", mut_fix), label) ) 

x2c_mut_plot <- mut %>% ggplot(aes(x=pct_tcga, y=pct_us)) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + scale_y_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + ylab("Estimated Mutation Proportion\nfor U.S. Population (%)") + xlab("Proportion of Unweighted Pan-Cancer Samples (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=x2c_points_mut, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label, hjust = position), lineheight = 0.65)


#log-log plot mutation
x2c_points_mut_log_log <- x2c_points_mut %>% mutate(ny = ifelse(gene_fix == "FGFR3", ny + 0.02, ny))
x2c_mut_log_plot <- mut %>% ggplot(aes(x=log(pct_tcga), y=log(pct_us))) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + ylab("Log Estimated Mutation Proportion\nfor U.S. Population (log %)") + xlab("Log Proportion\nof Unweighted Pan-Cancer Samples (log %)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=x2c_points_mut_log_log, mapping=aes(x=log(pct_tcga+nx), y=log(pct_us+ny), label=label, hjust = position), lineheight = 0.65) + geom_vline(xintercept = c(log(1), log(0.5)), color = "grey60", linetype = "dotted", size = 0.65) +  geom_hline(yintercept = c(log(1), log(0.5)), color = "grey60", linetype = "dotted", size = 0.65)  + scale_x_continuous(limits = c(-5, 1.5), breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(limits = c(-16, 1), breaks = scales::pretty_breaks(n = 8)) + geom_text(data = data.frame(x=c(log(1), log(0.5), -5, -5), y = c(-16, -16, log(1), log(0.5)), label = c("1%", "0.5%", "1%", "0.5%")), mapping = aes(x = x, y = y, label = label))


x2c_points_gene <- gene %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% mutate(label = ifelse(us_rank <= 3, str_sub(gene, 1, -2), NA)) %>% filter(!is.na(label)) %>% mutate(nx = 0.05) %>% mutate(ny = 0.03)


x2c_gene_plot <- gene %>% ggplot(aes(x=pct_tcga, y=pct_us)) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,5), breaks=seq(0, 5, 1)) + scale_y_continuous(limits=c(0,5), breaks=seq(0, 5, 1)) + ylab("Estimated Mutation Proportion\nfor U.S. Population (%)") + xlab("Proportion of Unweighted Pan-Cancer Samples (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=x2c_points_gene, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label), hjust=0)

#log-log plot gene
x2c_points_gene_log_log <- x2c_points_gene %>% mutate(ny = ifelse(gene == "TTNp", ny + 1, ny)) %>% mutate(ny = ifelse(gene == "TP53p", ny - 0.5, ny)) %>% mutate(nx = ifelse(gene == "KRASp", nx + 0.1, nx))
x2c_gene_log_plot <- gene %>% ggplot(aes(x=log(pct_tcga), y=log(pct_us))) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + ylab("Log Estimated Mutation Proportion\nfor U.S. Population (log %)") + xlab("Log Proportion of Unweighted\nPan-Cancer Samples (log %)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=x2c_points_gene_log_log, mapping=aes(x=log(pct_tcga+nx), y=log(pct_us+ny), label=label),  hjust = 0) + geom_vline(xintercept = c(log(1), log(5), log(2)), color = "grey60", linetype = "dotted", size = 0.65) +  geom_hline(yintercept = c(log(1), log(5), log(2)), color = "grey60", linetype = "dotted", size = 0.65)  + scale_x_continuous(limits = c(-5, 2), breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(limits = c(-14, 2), breaks = scales::pretty_breaks(n = 8)) + geom_text(data = data.frame(x=c(log(1), log(5), -5, -5, -5, log(2)), y = c(-12, -12, log(1), log(5), log(2), -12), label = c("1%", "5%", "1%", "5%", "2%", "2%")), mapping = aes(x = x, y = y, label = label))

mut.bar <- mut %>% arrange(desc(pct_tcga)) %>% slice(1:50) 
mut.bar <- mut.bar %>% mutate(label_text = paste0(str_sub(gene, 1, -2), " (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")"))

mut.bar$label_text <- factor(mut.bar$label_text, levels = mut.bar$label_text[order(mut.bar$pct_us, decreasing = TRUE)])

x2c_mut_bar <- mut.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(0.5, 1, 1.5), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,1.8), breaks=seq(0,2,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = "Estimated Number of New Cases In 2023\nWith Acquired Cysteine Mutations In U.S.", breaks = seq(0, 30000, 30000/6))) + theme(axis.ticks.y = element_line(color="black")) + xlab("Top 50 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 14, angle = 90, hjust = 1)) 


gene.bar <- gene %>% arrange(desc(pct_tcga)) %>% slice(1:50) 
gene.bar <- gene.bar %>% mutate(label_text = str_sub(gene, 1, -2))
gene.bar$label_text <- factor(gene.bar$label_text, levels = gene.bar$label_text[order(gene.bar$pct_us, decreasing = TRUE)])

x2c_gene_bar <- gene.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(1,2,3,4), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,4.5), breaks=seq(0,4,1), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = "Estimated Number of New Cases In 2023\nWith X2C Mutation In U.S.", breaks = seq(0, 80000, 80000/8))) + theme(axis.ticks.y = element_line(color="black")) + xlab("Top 50 Mutated Genes") + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 14, angle = 90, hjust = 1)) + geom_hline(yintercept = 10000/(0.01 * 1958310), color="grey60", linetype = "dashed")



#mut.estimate.fig <- cowplot::plot_grid(x2c_mut_plot, x2c_mut_bar, nrow = 1, rel_widths = c(1, 2.85))
#ggsave(plot = mut.estimate.fig, filename = file.path(outdir, "mutation-pop-estimates.pdf"), device=cairo_pdf, units = "in", height = 7, width = 20)

#gene.estimate.fig <- cowplot::plot_grid(x2c_gene_plot, x2c_gene_bar, nrow = 1, rel_widths = c(1, 2.85))
#ggsave(plot = gene.estimate.fig, filename = file.path(outdir, "gene-pop-estimates.pdf"), device=cairo_pdf, units = "in", height = 6, width = 20)

ggsave(plot = x2c_mut_plot, filename = file.path(outdir, "mutation-scatterplot.pdf"), device = cairo_pdf, units = "in", height = 7, width = 7)
ggsave(plot = x2c_mut_bar, filename = file.path(outdir, "mutation-bar.pdf"), device = cairo_pdf, units = "in", height = 8, width = 15)
ggsave(plot = x2c_mut_log_plot, filename = file.path(outdir, "mutation-log-scatterplot.pdf"), units = "in", height = 7, width = 7)

ggsave(plot = x2c_gene_plot, filename = file.path(outdir, "gene-scatterplot.pdf"), device = cairo_pdf, units = "in", height = 7, width = 7)
ggsave(plot = x2c_gene_bar, filename = file.path(outdir, "gene-bar.pdf"), device = cairo_pdf, units = "in", height = 7, width = 15)
ggsave(plot = x2c_gene_log_plot, filename = file.path(outdir, "gene-log-scatterplot.pdf"), units = "in", height = 7, width = 7)


a_cor <- mut %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% arrange(desc(pct_tcga)) %>% mutate(tcga_rank = 1:n()) %>% filter(us_rank <= 25) %>% ggplot(aes(x=us_rank, y = tcga_rank)) + geom_point(size = 2.5) + theme_minimal() + theme(panel.background =  element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), axis.ticks = element_line(color = "grey60")) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) + xlab("Epidemio-Genomic Rank") + ylab("Unweighted Pan-Cancer Rank") + ggtitle("Top 25 Acquired Cysteine Mutations")  

b_cor <- mut %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% arrange(desc(pct_tcga)) %>% mutate(tcga_rank = 1:n()) %>% filter(us_rank <= 100) %>% ggplot(aes(x=us_rank, y = tcga_rank)) + geom_point(size = 2.5) + theme_minimal() + theme(panel.background =  element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), axis.ticks = element_line(color = "grey60")) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) + xlab("Epidemio-Genomic Rank") + ylab("Unweighted Pan-Cancer Rank") + ggtitle("Top 100 Acquired Cysteine Mutations")  


c_cor <- mut %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% arrange(desc(pct_tcga)) %>% mutate(tcga_rank = 1:n()) %>% filter(us_rank <= 1000) %>% ggplot(aes(x=us_rank, y = tcga_rank)) + geom_point(size = 2.5) + theme_minimal() + theme(panel.background =  element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), axis.ticks = element_line(color = "grey60")) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) + xlab("Epidemio-Genomic Rank") + ylab("Unweighted Pan-Cancer Rank") + ggtitle("Top 1000 Acquired Cysteine Mutations")  

cowplot::plot_grid(a_cor, b_cor, c_cor, nrow = 1, labels = c("A", "B", "C"), label_size = 25) %>% ggsave(plot = ., filename = file.path(outdir, "rank-correlation.pdf"), device = cairo_pdf, units = "in", height = 6, width = 18)
 
get_correlation <- function(df, n.cysteine){

	x <- df %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% arrange(desc(pct_tcga)) %>% mutate(tcga_rank = 1:n()) %>% filter(us_rank <= n.cysteine)

	ctest <- cor.test(x$us_rank, x$tcga_rank, method="pearson")

	list(n = n.cysteine, cor = ctest$estimate[["cor"]], p = ctest$p.value)
}
get_correlation(mut, 25)
get_correlation(mut, 100)
get_correlation(mut, 1000)

######### KRASG12C findings

mut %>% filter((grepl("RASp", gene)|grepl("GNASp", gene))) %>% select(hugo, gene, pct_us, pct_tcga) %>% mutate(total_pts = pct_us * 0.01 * 1958310) %>% arrange(desc(total_pts)) %>% magrittr::set_colnames(c("mutation", "gene", "pct_us", "pct_unweighted_pan_cancer", "total_pts")) %>% write_csv(., file = file.path(outdir, "RAS-mutations.csv"))


########## Tables 1 and 2

cancer_contr <- tcga.muts %>% select(hugo, rosetta, cancer, incidence, count) %>% distinct() %>% left_join(., total %>% select(cancer, count, rosetta) %>% magrittr::set_colnames(c("cancer", "count_cancer", "rosetta")), by=c("cancer"="cancer", "rosetta"="rosetta"))

ranked_pct_contr <- cancer_contr %>% group_by(hugo) %>% mutate(total = sum(count), frac = count / total) %>% arrange(desc(frac)) %>% slice(1:3)

t25_mut <- mut %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% filter(us_rank <= 25)

t25_mut_table <- left_join(ranked_pct_contr, t25_mut %>% select(pct_us, hugo, gene), by=c("hugo"="hugo")) %>% filter(!is.na(pct_us)) %>% select(hugo, cancer, frac, gene, pct_us) %>% mutate(cancer = str_replace_all(cancer, "_", " ")) %>% ungroup() %>% mutate(gene_fix = str_sub(gene, 1, -2), mut_fix = paste0(" (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")")) %>% mutate(label_text = paste0(gene_fix, mut_fix) ) %>% select(-c(gene_fix, gene, mut_fix, hugo)) %>% mutate(indiv_imp = paste0(cancer, " (", round(frac*100, 2), "%)")) %>% select(-c(cancer, frac)) %>% group_by(pct_us, label_text) %>% mutate(row_id = row_number()) %>% pivot_wider(names_from = row_id, values_from = indiv_imp, names_prefix = "V") %>% mutate(cancers = paste0(V1, ", ", V2, ", ", V3)) %>% select(-starts_with("V")) %>% arrange(desc(pct_us)) %>% relocate(label_text) %>% ungroup() %>% magrittr::set_colnames(c("Point Mutation", "Mutation Proportion (%)", "The most common cancers to harbor a mutation in this gene (ROSETTA classification, % of all instances of this mutation)"))

write_csv(t25_mut_table, file=file.path(outdir, "t25-mutations.csv"))


cancer_contr_gene <- tcga.gene %>% select(gene, net_count, rosetta) %>% distinct() %>% left_join(., total %>% select(cancer, count, rosetta) %>% magrittr::set_colnames(c("cancer", "count_cancer", "rosetta")), by=("rosetta"="rosetta"))

ranked_pct_contr_gene <- cancer_contr_gene %>% group_by(gene) %>% mutate(total = sum(net_count), frac = net_count / total) %>% arrange(desc(frac)) %>% slice(1:3)

t25_gene <- gene %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% filter(us_rank <= 25)

t25_gene_table <- left_join(ranked_pct_contr_gene, t25_gene %>% select(pct_us, gene), by=c("gene"="gene")) %>% filter(!is.na(pct_us)) %>% select(cancer, frac, gene, pct_us) %>% mutate(cancer = str_replace_all(cancer, "_", " ")) %>% ungroup() %>% mutate(gene = str_sub(gene, 1, -2)) %>% mutate(indiv_imp = paste0(cancer, " (", round(frac*100, 2), "%)")) %>% select(-c(cancer, frac)) %>% group_by(pct_us, gene) %>% mutate(row_id = row_number()) %>% pivot_wider(names_from = row_id, values_from = indiv_imp, names_prefix = "V") %>% mutate(cancers = paste0(V1, ", ", V2, ", ", V3)) %>% select(-starts_with("V")) %>% arrange(desc(pct_us)) %>% relocate(gene) %>% ungroup() %>% magrittr::set_colnames(c("Point Mutation", "Mutation Proportion (%)", "The most common cancers to harbor a mutation in this gene (ROSETTA classification, % of all instances of this mutation)"))

write_csv(t25_gene_table, file=file.path(outdir, "t25-genes.csv"))



##### TCGA by mutation analyses 


#rename rosetta cancer types
rename_rosetta <- function(df){
	df %>% mutate(cancer = str_replace_all(cancer, "_", " ") %>% str_to_title(.)) %>% mutate(cancer = ifelse(grepl("Nos", cancer), str_replace_all(cancer, "Nos", "NOS"), cancer)) %>% mutate(cancer = ifelse(cancer == "Acral Lentiginous Melanoma Malignant", "Malignant Acral Lentiginous Melanoma", cancer)) %>% mutate(cancer = ifelse(cancer == "Chronic Lymphocytic Leukemiasmall Lymphocytic Lymphoma", "Chronic Lymphocytic Leukemia", cancer)) %>% mutate(cancer = ifelse(cancer == "Renal Cell Carcinoma Chromophobe Type", "Chromophobe Renal Cell Carcinoma", cancer))
}


gene_by_cancer_tcga_rates <- cancer_contr_gene %>% group_by(gene) %>% mutate(total = sum(net_count), frac = net_count / total) %>% arrange(desc(frac)) %>% select(gene, rosetta, cancer, net_count, frac) %>% ungroup() %>% mutate(cancer = str_replace_all(cancer, "_", " ")) %>% mutate(cancer = str_to_title(cancer)) %>% select(-rosetta)

colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
              "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")

scramble <- function(x){
	sample(x, size = length(x), replace = FALSE)
}

hmap_tcga <- gene_by_cancer_tcga_rates %>% group_by(gene) %>% mutate(total_count = sum(net_count)) %>% ungroup() %>% rename_rosetta(.) %>% mutate(gene = str_sub(gene, 1, -2))


genelevels <- hmap_tcga %>% group_by(gene) %>% summarize(tc = unique(total_count)) %>% arrange(desc(tc)) %>% pull(gene)

cancerlevels <- hmap_tcga %>% group_by(cancer) %>% summarize(tc = sum(net_count)) %>% arrange(desc(tc)) %>% pull(cancer)

#top 50 genes by number of total mutations  and top 20 cancers by number of total mutations
hmap_tcga <- hmap_tcga %>% filter(gene %in% genelevels[1:50]) %>% filter(cancer %in% cancerlevels[1:20])

hmap <- hmap_tcga %>% mutate(gene = factor(gene, levels=genelevels), cancer = factor(cancer, levels=cancerlevels)) %>% ggplot(aes(x=gene, y=cancer, fill=net_count)) + geom_tile(color="black") + scale_fill_gradient2(name = "No. of Unweighted Pan-Cancer\nSamples with X2C Mutations", low = "black", mid = "grey80", high = colors[length(colors)], limits = c(0,120), guide = "colorbar", na.value = "grey80", breaks = seq(0, 120, 15)) + theme(panel.border = element_rect(color="transparent", fill='transparent'), panel.grid = element_blank()) + ylab("Top 20 Mutated Cancers") + xlab("Top 50 Mutated Genes") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) + theme(axis.title = element_text(size = 16))


by_cancer <- tcga.gene %>% select(rosetta, net_count, gene) %>% left_join(., total %>% select(rosetta, count, cancer), by=c("rosetta"="rosetta")) %>% select(-rosetta) %>% filter(net_count > 0) %>% group_by(cancer) %>% summarize(s = sum(net_count), cancer_count = unique(count)) 

#throw out less than 1% (14398*0.01)

muts_by_cancer <- by_cancer %>% rename_rosetta(.) %>% mutate(s_per = s / cancer_count) %>% mutate(s_per = round(s_per, 2)) %>% mutate(cancer = factor(cancer, levels = cancer[order(s_per, decreasing=FALSE)])) %>% ggplot(aes(y=cancer, x=s_per, label=s)) + geom_col(position = position_dodge(width = 7), color="grey60", alpha = 0.8, fill = "grey60") + theme_minimal() + theme(panel.border = element_rect(color="black", fill="transparent"), panel.grid = element_blank()) + xlab("# of X2C Mutations per Unweighted Pan-Cancer Sample") + ylab("") + theme(axis.ticks = element_line(color="black")) + geom_text(aes(x=s_per + 0.5, y = cancer), hjust = 0) + scale_x_continuous(limits = c(0, 62), breaks = seq(0, 60, 10)) + geom_vline(xintercept = seq(10,60,10), linetype="dashed", alpha = 0.85, color="grey70") + theme(axis.title = element_text(size = 16))


by.cancer <- by_cancer %>% rename_rosetta(.) %>% mutate(s_per = s / cancer_count) %>% mutate(s_per = round(s_per, 2)) %>% mutate(cancer = factor(cancer, levels = cancer[order(s_per, decreasing=FALSE)])) %>% magrittr::set_colnames(c("Cancer", "Total X2C Mutations", "# of Pan-Cancer Samples Sequenced", "X2C Mutations per Sample Sequenced"))

## amino acid analyses 

aa <- tcga.muts %>% select(hugo, rosetta, count, cancer, gene) %>% mutate(gene = str_sub(gene, 1, -2), aminoacid = gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, 3)) %>% filter(count > 0) %>% rename_rosetta(.)

aa_levels <- aa %>% group_by(aminoacid) %>% summarize(tc = sum(count)) %>% arrange(desc(tc)) %>% pull(aminoacid)

cancerlevels_aa <- aa %>% group_by(cancer) %>% summarize(tc = sum(count)) %>% arrange(desc(tc)) %>% pull(cancer)

hmap_aa <- aa %>% group_by(cancer, aminoacid) %>% summarize(tc = sum(count)) %>% ungroup() %>% mutate(aminoacid = factor(aminoacid, levels=aa_levels), cancer = factor(cancer, levels=cancerlevels_aa)) %>% ggplot(aes(x=aminoacid, y=cancer, fill=log(tc))) + geom_tile(color="black") + scale_fill_gradient2(name = "", low = "black", mid = "grey80", high = colors[length(colors)], limits = c(0,12), guide = "colorbar", na.value = "grey80", breaks = seq(0, 12, 2)) + theme(panel.border = element_rect(color="black", fill='transparent'), panel.grid = element_blank()) + ylab("") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))


aa_cancer <- aa %>% group_by(cancer, aminoacid) %>% summarize(tc = sum(count)) %>% ungroup() %>% mutate(aminoacid = factor(aminoacid, levels=aa_levels)) %>% ggplot(aes(x=aminoacid, y=log(tc))) + geom_boxplot(outlier.alpha = 0, alpha = 0.75, fill="grey80") + theme_minimal()  + theme(panel.border = element_rect(color="black", fill='transparent'), panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) + theme(axis.ticks = element_line(color = "black")) + geom_jitter(width = 0.15, alpha = 0.85, color = "grey60") + ylab("Log(number of X2C mutations with initial AA for each cancer type)") + xlab("") + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 16))



mut.burden.fig <- cowplot::plot_grid(muts_by_cancer, hmap, aa_cancer, nrow = 1, rel_widths = c(1, 1.25, 1), labels=c("A", "B", "C"), label_size = 22)

ggsave(filename = file.path(outdir, "mutation-burden-figure.pdf"), plot = mut.burden.fig, device = cairo_pdf, height = 12, width = 32, units = "in")


aa.transitions <- aa %>% group_by(cancer, aminoacid) %>% summarize(tc = sum(count)) %>% ungroup() %>% mutate(aminoacid = factor(aminoacid, levels=aa_levels)) %>% arrange(desc(tc)) %>% magrittr::set_colnames(c("Cancer", "Amino Acid Mutated to Cysteine", "Total Number of Mutations in Pan-Cancer Samples"))


writexl::write_xlsx(list(S1 = all.muts, S2 = all.genes, S3 = by.cancer, S4 = aa.transitions), file.path(outdir, "supplementary-tables.xlsx"))







