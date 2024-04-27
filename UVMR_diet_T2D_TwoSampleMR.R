### Rscript for univariable MR analysis using 3 standard methods (Inverse variance weighted, Weighted median, Egger regression.
### Ref panel used for clumping: GRCh37, 1KG Phase3, 2013
### Dietary traits/preferences as exposures.
### T2D and related cardiometabolic traits as outcomes.

##Load packages for MR
library(devtools)
library(tidyverse)
library(TwoSampleMR)


##The genetic instruments (IVs) for exposure/dietary traits were selected via PLINK clumping procedure.

#Load the single summative MR info table consisting of exposure, outcome and potential mediator traits. 
#This table contains no duplicate or palindromic variants for exposure traits.
#It contains proxies for selected potential mediating traits to match the #SNPs used for MVMR analysis.
tb = read.table("final_single_MR_info_table_phase3_8mediators_FIunadj_proxies.txt",sep="\t",header=T)


 #Specify output dir.
output_dir = "./UVMR/"

#Specify exposure traits. The below is the list of dietary exposure traits.
exposure = c("ALCMEAL",
             "ALC",
             "BEEF",
             "BUTTER",
             "BUTMARG",
             "CARB",
             "CHAMPWH",
             "CHEESE",
             "COF",
             "COOKEDVEG",
             "CORNFLAK",
             "DRIEDFRU",
             "FAT",
             "FRESHFRU",
             "MUESLI",
             "NONOILYFSH",
             "PORK",
             "POULTRY",
             "PROTEIN",
             "RAWVEG",
             "REDWINE",
             "SPREADS",
             "SUGAR",
             "WHITEBRD",
             "WHOLEBRD",
             "WHOLEMLK",
             "ACQ",
             "CAFSWT",
             "COFALC",
             "DESS",
             "FATSALT",
             "LOWCAL",
             "PAL",
             "SAVCAL",
             "SAVOUR",
             "STR",
             "VEG")

#Note: SKIMMLK has matching variants only with T2D and 20 related traits, except AST, for harmonization. 
#LOWFAT and HEALTHY have zero IVs.
#FLORA and LOWFATMLK have less than 5 IVs.


#Specify outcome traits. Below contains T2D and 21 related cardiometabolic traits.
outcome = c("ALP",
            "ALT",
            "AST",
            "ASAT",
            "BMI",
            "FG",
            "FI",
            "GGT",
            "HbA1c",
            "HDL",
            "LDL",
            "Liverfat",
            "Liveriron",
            "Livervol",
            "Pancfat",
            "Panciron",
            "Pancvol",
            "TG",
            "T2D",
            "VAT",
            "WHR",
            "WHRadjBMI")

#Run this loop for UVMR analysis.
for (e in exposure){
  for(t in outcome){
    #Specify and select exposure and outcome traits of interest.
    exp_data = tb %>% filter(Type == 'Exposure' & Trait == e & Proxy== 'FALSE') 
    out_data = tb %>% filter(Type == 'Outcome' & Trait == t & Proxy == 'FALSE') 
    exp <- format_data(exp_data, 
                       type = "exposure",
                       log_pval = FALSE,
                       snps = NULL, 
                       header=TRUE, 
                       phenotype_col = "Phenotype", 
                       chr_col = "CHR",
                       pos_col = "POS", 
                       snp_col = "SNP",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "EA",
                       other_allele_col="NEA",
                       eaf_col="EAF",
                       pval_col="P",
                       samplesize_col = "N")
    
    out <- format_data(out_data,
                       type = "outcome",
                       log_pval = FALSE,
                       snps = NULL,
                       header=TRUE, 
                       chr_col = "CHR",
                       pos_col = "POS", 
                       snp_col = "SNP", 
                       beta_col = "BETA", 
                       se_col= "SE", 
                       effect_allele_col = "EA", 
                       other_allele_col= "NEA", 
                       eaf_col = "EAF",
                       pval_col = "P", 
                       samplesize_col = "N")
    
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

#Separate analysis on SKIMMLK because it has no matching variants with AST.
#The outcome includes only T2D and 20 cardiometabolic traits. 
exposure = 'SKIMMLK'
outcome = c("ALP",
              "ALT",
              "ASAT",
              "BMI",
              "FG",
              "FI",
              "GGT",
              "HbA1c",
              "HDL",
              "LDL",
              "Liverfat",
              "Liveriron",
              "Livervol",
              "Pancfat",
              "Panciron",
              "Pancvol",
              "TG",
              "T2D",
              "VAT",
              "WHR",
              "WHRadjBMI")

#The same code is used to run MR analysis on SKIMMLK.
for (e in exposure){
  for(t in outcome){ 
    exp_data = tb %>% filter(Type=='Exposure' & Trait == e) 
    out_data = tb %>% filter(Type=='Outcome' & Trait == t) 
    exp <- format_data(exp_data,
                       type = "exposure",
                       log_pval = FALSE,
                       snps = NULL, 
                       header=TRUE, 
                       phenotype_col = "Phenotype", 
                       chr_col = "CHR",
                       pos_col = "POS", 
                       snp_col = "SNP",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "EA",
                       other_allele_col="NEA",
                       eaf_col="EAF",
                       pval_col="P",
                       samplesize_col = "N")
    
    out <- format_data(out_data, 
                       type = "outcome", 
                       log_pval = FALSE,
                       snps = NULL, 
                       header=TRUE, 
                       chr_col = "CHR",
                       pos_col = "POS", 
                       snp_col = "SNP", 
                       beta_col = "BETA", 
                       se_col= "SE", 
                       effect_allele_col = "EA", 
                       other_allele_col= "NEA", 
                       eaf_col = "EAF",
                       pval_col = "P", 
                       samplesize_col = "N")
    
    #Generate harmonized data for exposure-outcome trait pairs
    #action = 2: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    dat = harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
    #MR Steiger test for the directionality of causal effect.
    st = directionality_test(dat)
    st$exposure = e
    st$outcome = t
    write.table(st, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_steiger_test.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    dat2 = dat %>% mutate(Exposure = e, Outcome = t) 
    write.table(dat2,file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_SNP_data_table.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    res <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    res2 = res %>% select(-exposure,-outcome)
    res2$exposure = e
    res2$outcome = t
    write.table(res2, file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    het = mr_heterogeneity(dat)
    het2 = het %>% select(-exposure,-outcome)
    het2$exposure = e
    het2$outcome = t
    write.table(het2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_heterogeneity.txt"), sep="\t", quote=F, row.names=F)
    plt = mr_pleiotropy_test(dat)
    plt2 = plt %>% select(-exposure,-outcome)
    plt2$exposure = e
    plt2$outcome = t
    write.table(plt2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_pleiotropy.txt"),sep="\t", quote=F, row.names=F)
    sin = mr_singlesnp(dat)
    sin2 = sin %>% select(-exposure,-outcome)
    sin2$exposure = e
    sin2$outcome = t
    write.table(sin2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_singleSNP_analysis.txt"),sep="\t", quote=F, row.names=F)
    or = generate_odds_ratios(res)
    or2 = or %>% select(-exposure,-outcome)
    or2$exposure = e
    or2$outcome = t
    write.table(or2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_OR_with_CI95.txt"),sep="\t", quote=F, row.names=F)
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
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_results.txt   > Univariable_MR_dietary_traits_results.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_OR_with_CI95.txt   > Univariable_MR_dietary_traits_OR_with_CI95.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_heterogeneity.txt   > Univariable_MR_dietary_traits_heterogeneity.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_pleiotropy.txt   > Univariable_MR_dietary_traits_pleiotropy.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_ivw.txt   > Univariable_MR_dietary_traits_leaveoneout_ivw.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_wm.txt   > Univariable_MR_dietary_traits_leaveoneout_wm.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_steiger_test.txt   > Univariable_MR_dietary_traits_steiger_test.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_SNP_data_table.txt   > Univariable_MR_dietary_traits_SNP_data_table.txt"))
system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_singleSNP_analysis.txt   > Univariable_MR_dietary_traits_singleSNP_analysis.txt"))

#Filter out exposure-outcome pairs based on 2 critera: (1) Pass the given Bonferroni-adjusted pval (5.99e-5) in at least 2 sensitivity analyses AND (2) more than 5 genetic instruments are used for MR.
d = read.table(file=paste0(output_dir,"Univariable_MR_dietary_traits_OR_with_CI95.txt"),sep='\t',header=T)

d2 = d %>% select(-id.exposure,-id.outcome) %>%
relocate(exposure,.before=method) %>%
relocate(outcome,.after=exposure) %>%
filter(pval <= 5.99e-5, nsnp >= 5) 
sig = d2 %>% group_by(exposure,outcome) %>% filter(n()>= 2)

#Save significant associations from MR based on the 2 criteria above.
write.table(sig, file=paste0(output_dir,"sig_Univariable_MR_dietary_traits_OR_with_CI95.txt"),sep='\t',row.names=F,quote=F)


