library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../results"
if (!dir.exists(outdir)){
	dir.create(outdir)
}

tcga <- "../data/braf-data/Genomics_Output_Processed_specific_Braf.txt" %>% read_table()

tcga <- tcga %>% filter(grepl("BRAF|Total", Hugo_Symbol))


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


seer.df <- "../data/processed-data/SuppTable1_ROSETTA_Abundance_added12.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

total <- tcga_total_by_rosetta %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 

#renormalize weight vector to 1, adds to 0.583 now

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

writexl::write_xlsx(
  x = list(V600_single = mut, V600_either = gene),
  path = "../results/braf-v600-mut-pct.xlsx"
)






