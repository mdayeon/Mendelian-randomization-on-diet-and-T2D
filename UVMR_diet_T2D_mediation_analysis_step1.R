### Rscript for univariable MR analysis using 3 standard methods (Inverse variance weighted, Weighted median, Egger regression.
### Ref panel used for clumping: GRCh37, 1KG Phase3, 2013
### Dietary traits/preferences as exposures.
### T2D and related cardiometabolic traits as outcomes.

##Load packages for MR

library(tidyverse)
library(TwoSampleMR)


##The genetic instruments (IVs) for exposure/dietary traits were selected via PLINK clumping procedure.

#Load the single summative MR info table consisting of exposure, outcome and potential mediator traits. 
tb = read.table("supp table 9 Table of genetic instruments for exposures, outcomes and potential mediators selected for UVMR and MVMR.txt",sep="\t",header=T)

#This table contains no duplicate or palindromic variants for exposure traits.
#It contains proxies for selected potential mediating traits to match the #SNPs used for MVMR analysis.  

#Specify output dir.
output_dir = "/project/voight_T2D_UM1_FGP/diane531/Diet/MR/MR_mediation/"

#Specify exposure traits. The below is the list of dietary exposure traits.
exposure = c("CHEESE")

#SKIMMLK has matching variants with T2D and 20 related traits, except AST, for harmonization. 
#LOWFAT and HEALTHY have zero IVs.
#FLORA and LOWFATMLK have less than 5 IVs.


#Specify mediators which will be tested as outcome traits of interest.
outcome = c("DBP")

#Run this loop for UVMR analysis.
for (e in exposure){
  for(t in outcome){
    #Specify and select exposure and outcome traits of interest.
    print(e)
    print(t)
    exp_data = tb %>% filter(Type=='Exposure' & Trait==e & Proxy == 'FALSE') %>% filter(!SNP %in% c('rs28569885','rs13107325'))
    out_data = tb %>% filter(Type=='Mediator' & Trait==t & Proxy == 'FALSE') 
    exp <- format_data(exp_data, type = "exposure", log_pval = FALSE,snps = NULL, header=TRUE, phenotype_col = "Phenotype", chr_col = "CHR",pos_col = "POS", snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "EA",other_allele_col="NEA",eaf_col="EAF",pval_col="P",samplesize_col = "N")
    out <- format_data(out_data, type = "outcome", log_pval = FALSE,snps = NULL, header=TRUE, chr_col = "CHR",pos_col = "POS", snp_col = "SNP", beta_col = "BETA", se_col= "SE", effect_allele_col = "EA", other_allele_col= "NEA", eaf_col = "EAF",pval_col = "P", samplesize_col = "N")
    #Generate harmonized data for exposure-outcome trait pair.
    #action = 2: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    dat = harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
    #MR Steiger test for directionality
    st = directionality_test(dat)
    st$exposure = e
    st$outcome = t
    write.table(st, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_steiger_test.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Save the list of SNPs used for the analysis.  
    dat2 = dat %>% mutate(Exposure = e, Outcome = t) 
    write.table(dat2,file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_SNP_data_table.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Run MR analysis using IVW, WM and Egger regression as methods.
    res <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    res2 = res %>% select(-exposure,-outcome)
    res2$exposure = e
    res2$outcome = t
    write.table(res2, file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Heterogeneity test.
    het = mr_heterogeneity(dat)
    het2 = het %>% select(-exposure,-outcome)
    het2$exposure = e
    het2$outcome = t
    write.table(het2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_heterogeneity.txt"), sep="\t", quote=F, row.names=F)
    #Pleiotropy test.
    plt = mr_pleiotropy_test(dat)
    plt2 = plt %>% select(-exposure,-outcome)
    plt2$exposure = e
    plt2$outcome = t
    write.table(plt2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_pleiotropy.txt"),sep="\t", quote=F, row.names=F)
    #MR analysis on each SNP individually.
    sin = mr_singlesnp(dat)
    sin2 = sin %>% select(-exposure,-outcome)
    sin2$exposure = e
    sin2$outcome = t
    write.table(sin2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_singleSNP_analysis.txt"),sep="\t", quote=F, row.names=F)
    #Add odds ratios (OR) to MR results. 
    or = generate_odds_ratios(res)
    or2 = or %>% select(-exposure,-outcome)
    or2$exposure = e
    or2$outcome = t
    write.table(or2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_OR_with_CI95.txt"),sep="\t", quote=F, row.names=F)
    #Leave one out sensitivity analysis using IVW and WV methods to determine whether single, particular SNP drives the causal effect.
    leave_ivw = mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
    leave_ivw2 = leave_ivw %>% select(-exposure,-outcome)
    leave_ivw2$exposure = e
    leave_ivw2$outcome = t
    write.table(leave_ivw2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_leaveoneout_ivw.txt"),sep="\t", quote=F, row.names=F)
    leave_wm =mr_leaveoneout(dat, parameters = default_parameters(), method = mr_weighted_median)
    leave_wm2 = leave_wm %>% select(-exposure,-outcome)
    leave_wm2$exposure = e
    leave_wm2$outcome = t
    write.table(leave_wm2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_leaveoneout_wm.txt"),sep="\t", quote=F, row.names=F)
    #Save scatter and forest plots of MR results.
    p1 <- mr_scatter_plot(res2, dat)
    length(p1)
    ggsave(p1[[1]], file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_scatter_plot.pdf"), width=7, height=7)
    res_single <- mr_singlesnp(dat, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    p2 <- mr_forest_plot(res_single)
    ggsave(p2[[1]], file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_forest_plot.pdf"), width=7, height=7)
    options(warn=-1)
  }
}


#Merge MR output files into single files.
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_results.txt   > Univariable_MR_dietary_traits_mediator_results.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_OR_with_CI95.txt   > Univariable_MR_dietary_traits_mediator_OR_with_CI95.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_heterogeneity.txt   > Univariable_MR_dietary_traits_mediator_heterogeneity.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_pleiotropy.txt   > Univariable_MR_dietary_traits_mediator_pleiotropy.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_ivw.txt   > Univariable_MR_dietary_traits_mediator_leaveoneout_ivw.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_wm.txt   > Univariable_MR_dietary_traits_mediator_leaveoneout_wm.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_steiger_test.txt   > Univariable_MR_dietary_traits_mediator_steiger_test.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_SNP_data_table.txt   > Univariable_MR_dietary_traits_mediator_SNP_data_table.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_singleSNP_analysis.txt   > Univariable_MR_dietary_traits_mediator_singleSNP_analysis.txt"))

#d = read.table(file=paste0(output_dir,"Univariable_MR_dietary_traits_mediator_OR_with_CI95.txt"),sep='\t',header=T)

#d2 = d %>% select(-id.exposure,-id.outcome) %>% relocate(exposure,.before=method) %>% relocate(outcome,.after=exposure) 

#Save significant associations from MR based on the 2 criteria above.
#write.table(d2, file=paste0(output_dir,"Univariable_MR_dietary_traits_mediator_OR_with_CI95.txt"),sep='\t',row.names=F,quote=F)


