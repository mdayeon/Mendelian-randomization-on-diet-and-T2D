### Rscript for 2 multivariable MR (MVMR) using IVW: Single-mediator MVMR and Multiple-mediator MVMR.
### For further evaluation of significant exposure-outcome associations identified in UVMR.
### It first prepares and formats input data first to run via 2SampleMR package.
### Then,it uses the input data created to format the data and run multivariable MR via mvmr package.
### This approach allows harmonize exposure, mediator and outcome data and thereby compare multivariable MR results from between 2SampleMR and mvmr packages. 
### Reference: some of the codes are obtained and sourced from https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/
### The functions and steps used are also described and provided by https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#multivariable-mr

#Load package.
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)


#Load the single summative MR info table consisting of exposure, outcome and potential mediator (i.e., secondary exposure) traits, including proxies (for missing variants). 
tb = read.table("final_single_MR_info_table_phase3_8mediators_FIunadj_proxies.txt",sep="\t",header=T)

#Specify exposures(s) of interest.
exposure = c("MUESLI")

#Specify outcome(s) of interest.
outcome = c("T2D")

#The 8 potential mediators are: "BMI", "WHR", "FI", "EA", "PHY","SED","DBP","SBP".
#The table also includes "FIunadj" (unadjusted FI) for comparison because "FI" is adjusted for BMI.

#Specify potential mediator(s) of interest --> single mediator, non-BMI or all-inclusive.

#For single-mediator MR, specify 1 mediator trait at a time. 
	#mediator = c("BMI", "DBP")

#For multiple-mediator MR: non-BMI or all-inclusive
	#mediator = non_BMI = c("WHR","FI","EA","PHY","SED","DBP","SBP")
	mediator = all = c("BMI","WHR","FI","EA","PHY","SED","DBP","SBP")

#Specify path to output directories: Single-mediator MVMR, Non-BMI, All-inclusive MVMR.
output_dir = "./MVMR/"

##Create functions necessary to prepare input files and format output files.

#Function to create tidy outcome for results from TwoSampleMR.
tidy_pvals<-function(df){
  #The function rounds up output values and formats pvals in scientific notation.
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=4) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}


#Function to create tidy outcome for results from mvmr.
tidy_mvmr_output <- function(mvmr_res) {
  #The function is to tidy up the returned output.
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

#Function to create input file for MVMR. 
make_mvmr_input <- function(exposure_dat, outcome.data = outcome_dat){
  #The function provides exposure_dat created in the same way as for TwoSampleMR. 
  outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  #Harmonize datasets for exposure and outcome
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat) 
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  #Create variables fpr exposure and mediator (secondary exposures) for analysis.  
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  #Add "beta" and "se" columns for exposure.
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}


#Run this loop for multivariable MR analysis. In the code below, specify the filename based on the analysis selected.
for (e in exposure){
  for (t in outcome){
    for (m in mediator){
      #Load the data for exposure. 
      exp_dat = tb %>% filter(Type=='Exposure' & Trait==e)
      #Load the data for mediator, either single or multiple.
        #For single mediator/pairwise analysis.
           #med_dat = tb %>% filter(Type=='Mediator',Trait == m)
	#For multiple mediator analysis, either non-bmi or all-inclusive.
	   med_dat = tb %>% filter(Type == 'Mediator') %>% filter(Trait %in% c(mediator))			
      #Merge exposure and mediator (as a secodary exposure) data into one data table.
      exp_data = rbind(exp_dat,med_dat)
      #Format exposure/mediator data.
      exposure_dat <- format_data(exp_data, type = "exposure",
                                  log_pval = FALSE,
                                  snps = NULL,
                                  header=TRUE,
                                  phenotype_col = "Trait",
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
      #Load the data for outcome trait.
      out_data = tb %>% filter(Type=='Outcome' & Trait==t) 
      outcome_dat <- format_data(out_data, 
                                 type = "outcome",
                                 log_pval = FALSE,
                                 snps = NULL,
                                 header=TRUE,
                                 phenotype_col = "Trait",
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
      
      #Harmonise data to have all SNPs on the same reference allele. Also, harmonization removes SNPs that have incompatible alleles or are duplicates.
      mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
      #Save the SNPs harmonized and included for analysis.
         #For output filename, indicate if the analysis is single-mediator, non-bmi, or all-inclusive.
           #filename = paste0("Multivariable_MR_",e,"_",m,"_",t)
           #filename = paste0("Multivariable_MR_",e,"_nonbmi_",t)
           filename = paste0("Multivariable_MR_",e,"_all_",t)
      capture.output(df = data.frame(exposure.beta=mvdat$exposure_beta,
                                      exposure.pval=mvdat$exposure_pval,
                                      exposure.se=mvdat$exposure_se,
                                      outcome.beta = mvdat$outcome_beta,
                                      outcome.pval=mvdat$outcome_pval,
                                      outcome.se = mvdat$outcome_se),
                     file = paste0(output_dir,filename,"_SNP_data_table.txt")) 
      #Perform multivariable MR analysis.
      res <- mv_multiple(mvdat)
      
      #Create a tidy outcome.
      result_2smr <- res$result %>%
        split_outcome() %>%
        separate(outcome, "outcome", sep="[(]") %>% 
        mutate(outcome=stringr::str_trim(outcome))%>% 
        generate_odds_ratios() %>% 
        select(-id.exposure, -id.outcome) %>% 
        tidy_pvals()
      #Print out 2SampleMR multivariable MR results.
      result_2smr
      #Save MVMR results via TwoSampleMR package.
      write.table(result_2smr, file=paste0(output_dir,filename,"_2smr_ivw_results.txt"),sep='\t',row.names=F,quote=F)

      #Use "make_mvmr_input" function to create input data for analysis via mvmr package.
      mvmr_input <- make_mvmr_input(exposure_dat, outcome.data=outcome_dat)
      
      #Check the input data 
      glimpse(mvmr_input$XGs)
      glimpse(mvmr_input$Y)
      
      #Format data to be compatible for analysis via mvmr package.
      mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  # exposure betas
                              BYG = mvmr_input$YG$beta.outcome,                     # outcome beta
                              seBXGs = mvmr_input$XGs %>% select(contains("se")),  # exposure SEs
                              seBYG = mvmr_input$YG$se.outcome,                     # outcome SEs
                              RSID = mvmr_input$XGs$SNP)    
      
      head(mvmr_out)
      
      #Estimate causal effects using method via mvmr package.
      mvmr_res <- ivw_mvmr(r_input=mvmr_out)
      
      #Tidy up the output format
      #Specify the outcome trait
      result_mvmr <-
        mvmr_res %>% 
        tidy_mvmr_output() %>% 
        mutate(exposure = mvmr_input$exposures,
               outcome = t) %>% 
        select(exposure, outcome, everything()) %>% 
        tidy_pvals()
      
      #Print out multivariable MR results via mvmr package
      result_mvmr
      
      #Save MVMR results from mvmr package.
      write.table(result_mvmr, file=paste0(output_dir,filename,"_mvmr_ivw_results.txt"), sep='\t',row.names=F, quote=F)
    }
  }
}

##the example of 2SampleMR multivariable MR results 
#result_2smr
#exposure outcome nsnp       b     se    pval   lo_ci  up_ci     or or_lci95
#1      BMI     T2D    5  1.0515 0.2394 1.1e-05  0.5823 1.5206 2.8618   1.7902
#2   MUESLI     T2D   13 -0.1585 0.2078 4.5e-01 -0.5658 0.2488 0.8534   0.5679
#or_uci95
#1   4.5749
#2   1.2825






