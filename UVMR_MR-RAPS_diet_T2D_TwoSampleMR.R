### Rscript for univariable MR analysis using MR-RAPS.
### Ref panel used for clumping: GRCh37, 1KG Phase3, 2013
### Dietary traits/preferences as exposures.
### T2D and related cardiometabolic traits as outcomes.

##Load packages for MR

library(devtools)
library(tidyverse)
library(TwoSampleMR)


##The genetic instruments (IVs) for exposure/dietary traits were selected via PLINK clumping procedure.

#Load the single summative MR info table consisting of exposure, outcome and potential mediator traits. 
tb = read.table("weak_instruments_mr_raps/final_single_MR_info_table_both_weak_strong_instruments_for_mr_raps_mediators_wo_pal.txt",sep="\t",header=T)
#This table contains no duplicate or palindromic variants for exposure traits.
#This table is made for MR-RAPS to test weak instruments of exposure traits of interest to determine whether causal association can be observed in presence of weak instrument bias. 


#Specify output dir
output_dir = "./MR_RAPS/"

#Specify exposure traits.
#The exposure traits of interest are: ALC, CHEESE, COF, COOKEDVEG, DRIEDFRU, MUESLI, PAL, RAWVEG, SAVCAL, SPREADS, STR, WHITEBRD.
exposure = c('ALC')

#Specify outcome (T2D) traits
#The outcome traits of interest are: ASAT, BMI, WHR, FI, FG, HDL, Liverfat, Liveriron, T2D, TG, ALP, ALT, GGT.
outcome = c('TG')

#Run this loop for UVMR analysis.
for (e in exposure){
  for(t in outcome){
    #Specify and select exposure and outcome traits of interest.
    #Choose to include both strong and weak instruments or only weak instruments.
    exp_data = tb %>% filter(Type == 'Exposure' & Trait == e)
    out_data = tb %>% filter(Type == 'Outcome' & Trait == t)
    #For weak instruments only:
       #exp_data = tb %>% filter(Type=='Exposure' & Trait==e & Strength == 'Weak')
       #out_data = tb %>% filter(Type=='Outcome' & Trait==t & Strength == 'Weak')
    exp <- format_data(exp_data, type = "exposure", log_pval = FALSE,snps = NULL, header=TRUE, phenotype_col = "Phenotype", chr_col = "CHR",pos_col = "POS", snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "EA",other_allele_col="NEA",eaf_col="EAF",pval_col="P",samplesize_col = "N")
    out <- format_data(out_data, type = "outcome", log_pval = FALSE,snps = NULL, header=TRUE, chr_col = "CHR",pos_col = "POS", snp_col = "SNP", beta_col = "BETA", se_col= "SE", effect_allele_col = "EA", other_allele_col= "NEA", eaf_col = "EAF",pval_col = "P", samplesize_col = "N")
    
    #Generate harmonized data for exposure-outcome trait pair.
    #action = 2: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    dat = harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
    
    #MR Steiger test for directionality
    st = directionality_test(dat)
    st$exposure = e
    st$outcome = t
    write.table(st, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_steiger_test.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Save the list of SNPs used for the analysis.
    dat2 = dat %>% mutate(Exposure = e, Outcome = t) 
    write.table(dat2,file = paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_SNP_data_table.txt"), sep="\t", row.names=FALSE, quote=FALSE)
   
    #Run MR analysis using MR-RAPS as a method.
    res <- mr(dat, method_list=c("mr_raps"))
    res2 = res %>% select(-exposure,-outcome)
    res2$exposure = e
    res2$outcome = t
    write.table(res2, file = paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    
    #Heterogeneity test.
    het = mr_heterogeneity(dat)
    het2 = het %>% select(-exposure,-outcome)
    het2$exposure = e
    het2$outcome = t
    write.table(het2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_heterogeneity.txt"), sep="\t", quote=F, row.names=F)
    
    #Pleiotropy test.
    plt = mr_pleiotropy_test(dat)
    plt2 = plt %>% select(-exposure,-outcome)
    plt2$exposure = e
    plt2$outcome = t
    write.table(plt2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_pleiotropy.txt"),sep="\t", quote=F, row.names=F)
    
    #MR analysis on each SNP individually.
    sin = mr_singlesnp(dat)
    sin2 = sin %>% select(-exposure,-outcome)
    sin2$exposure = e
    sin2$outcome = t
    write.table(sin2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_singleSNP_analysis.txt"),sep="\t", quote=F, row.names=F)
    
    #Add odds ratios (OR) to MR results. 
    or = generate_odds_ratios(res)
    or2 = or %>% select(-exposure,-outcome)
    or2$exposure = e
    or2$outcome = t
    write.table(or2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_OR_with_CI95.txt"),sep="\t", quote=F, row.names=F)
   
    #Save scatter and forest plots of MR results.
    p1 <- mr_scatter_plot(res2, dat)
    length(p1)
    ggsave(p1[[1]], file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_scatter_plot.pdf"), width=7, height=7)
    res_single <- mr_singlesnp(dat, all_method=c("mr_raps"))
    p2 <- mr_forest_plot(res_single)
    ggsave(p2[[1]], file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_forest_plot.pdf"), width=7, height=7)
    options(warn=-1)
  }
}


#Merge MR output files into single files.
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_results.txt   > Univariable_MR-RAPS_dietary_traits_results.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_OR_with_CI95.txt   > Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_heterogeneity.txt   > Univariable_MR-RAPS_dietary_traits_heterogeneity.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_pleiotropy.txt   > Univariable_MR-RAPS_dietary_traits_pleiotropy.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_steiger_test.txt   > Univariable_MR-RAPS_dietary_traits_steiger_test.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_SNP_data_table.txt   > Univariable_MR-RAPS_dietary_traits_SNP_data_table.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_singleSNP_analysis.txt   > Univariable_MR-RAPS_dietary_traits_singleSNP_analysis.txt"))

#Filter out the pairs that pass Bonferroni-adjusted pval.
d = read.table(file=paste0(output_dir,"Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"),sep='\t',header=T)

#Filter out exposure-outcome pairs based on the  pass the given Bonferroni-adjusted pval (P < 2.94e-3)
d2 = d %>% select(-id.exposure,-id.outcome) %>%
relocate(exposure,.before=method) %>%
relocate(outcome,.after=exposure) %>%
filter(pval < 2.94e-3) 

#Save significant associations from MR-RAPS.
write.table(d2, file=paste0(output_dir,"sig_Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"),sep='\t',row.names=F,quote=F)
