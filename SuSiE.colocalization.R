
library(susieR)
library(TwoSampleMR)
library(biomaRt)



#load list of drug-target MR genes for colocalization (See Manuscript for details on selecting top genes)
longevity.top <- read.csv("path to/longevity.top.csv")


#get eQTLGen (Vosa et al. 2022, Nature Genetics: PMID:34475573)

ao <- available_outcomes()
ao.vosa <- subset(ao, author == "Vosa")

#load and format Timmers et al. multivariate longevity data (see Data availability for links)

longevity.format <- with(longeivty, data.frame(SNP=SNP, chr=CHR, pos=BP, effect_allele=A1, other_allele=A2, beta=BETA, se=SE, pval=P,samplesize=709709))
longevity.formatF <- format_data(longevity.format, type="outcome")




#get ensembl gene ID, start position, and end position 
ensembl=useMart("ensembl")
grch37 = useEnsembl(biomart="ensembl",GRCh=37)

esemblist <- as.data.frame(listDatasets(grch37))
ensembl = useDataset("hsapiens_gene_ensembl",mart= grch37)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)

timmers.forsusie <- merge(longevity.top, t2g, by= 'ensembl_gene_id')
timmers.forsusie$lower.100 <- timmers.forsusie$start_position - 100000
timmers.forsusie$upper.100 <- timmers.forsusie$end_position + 100000



#create exposure datasets eQTLGen eQTLs 
vosa.pheno.exp <- extract_instruments(outcomes=vosa.pheno$id, clump=F, p1=0.99)


vosa.pheno.exp$id <- sub(".*-", "",vosa.pheno.exp$exposure )


#write and save to folder 
sptDAT1 <- split(timmers.forsusie, timmers.forsusie$ensembl_gene_id) # for 1 outcome - where DAT is harmonized dataset of rosmap X on chronic pain

# Save 
setwd("~/path to instrument folder/")
lapply(names(sptDAT1), function(x){
  readr::write_csv(sptDAT1[[x]], path = paste(x, ".csv", sep = ""))
})




###creat cis-eQTL files to be harmonized with outcome data 
filenamesHD <- list.files(path = "~/Desktop/susie.longevity/",  
                          pattern = "csv",
                          full.names = TRUE)
#filenames <- filenamesHD[1:5]
filenames <- filenamesHD[1:length(filenamesHD)]
filenames

# eur - in calc.rho need pop="EAS" for eastern asianb
analyze3 <- function(filename) {
  dat <- read.csv(file = filename, header = TRUE)
  dat.exp <- subset(vosa.pheno.exp, vosa.pheno.exp$id == dat$ensembl_gene_id & vosa.pheno.exp$chr.exposure == dat$chromosome_name & (vosa.pheno.exp$pos.exposure >= dat$lower.100 & vosa.pheno.exp$pos.exposure <= dat$upper.100) )
    dat.exp$effect_allele.exposure <- ifelse(dat.exp$effect_allele.exposure == TRUE, "T", dat.exp$effect_allele.exposure)
  dat.exp$other_allele.exposure <- ifelse(dat.exp$other_allele.exposure == TRUE, "T", dat.exp$other_allele.exposure)
  if(nrow(dat.exp) > 0){
    dat.exp$id.exposure <- dat$ensembl_gene_id
    #print(dat)
    write.csv(dat.exp,file=paste0(filename, "instruments.csv"))
  }
}


for (f in filenames) {
  print(f)
  analyze3(f) 
  
}


####

#harmonize with longevity outcomes  
filenamesHD <- list.files(path = "~/path to instruments/",  
pattern = "csv",
full.names = TRUE)
#filenames <- filenamesHD[1:2]
filenames <- filenamesHD[1:length(filenamesHD)]
filenames

# eur - in calc.rho need pop="EAS" for eastern asianb
analyze3 <- function(filename) {
  dat <- read.csv(file = filename, header = TRUE)
  dat$exposure <- filename
  #dat <- subset(dat0, dat0$pval.exposure < 0.00000005)
  datout <- subset(longevity.formatF, SNP %in% dat$SNP)
   count <- nrow(datout)
  if(isTRUE(count > 0)) {
    dath <- TwoSampleMR::harmonise_data(dat, datout)
    dath$effect_allele.exposure <- ifelse(dath$effect_allele.exposure == TRUE, "T", dath$effect_allele.exposure)
    dath$other_allele.exposure <- ifelse(dath$other_allele.exposure == TRUE, "T", dath$other_allele.exposure)
    #print(dat)
    write.csv(dath,file=paste0(filename, "harmonized.csv"))
    
  }
  
}


for (f in filenames) {
  print(f)
  analyze3(f) 
  
}






###loop to perform SuSie colocaliztion over folder of harmonized files 
filenamesHD <- list.files(path = "~/path to harmonized/",  
                          pattern = ".csv",
                          full.names = TRUE)

filenames <- filenamesHD[1:length(filenamesHD)]
head(filenames) # to check to see if you pull in correct files


analyze6 <- function(filename) {
  dat <- read.csv(file = filename, header = TRUE)
  
  calc.rho <- ieugwasr::ld_matrix(
    variants = dat$SNP,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "/path/LDREF/EUR", with_alleles = FALSE)
  # NOT SAME ORDER AS DAT$SNP; ALSO DROPS SNPS NOT IN REF PANEL
  ## subset for only those found in ref file
  Msnp <- colnames(calc.rho)
  
  ##
  
  colnames(calc.rho)
  M <- calc.rho
  datrev <- subset(dat, SNP %in% Msnp) # revising harmonized data to exclude SNPs not found in Ref Panel (not in LD matrix)
  datrev$SNP
  
  # ordering snps in LD matrix (now called M) to align with datrev(ised)
  col.order <- datrev$SNP
  row.order <- datrev$SNP
  
  Mreord <- M[row.order, col.order]
  rownames(Mreord)
  
  dat <- datrev # changing name of dat revised back to dat (because I didn't want to rewrite subsequent code)
  head(dat)
  dat$SNP
  ### variance # make sure that dat has se.exposure and se.outcome (get rid of NAs)
  dat$var.exposure <- dat$se.exposure^2
  dat$var.outcome <- dat$se.outcome^2
  summary(dat)
  dat[,c("SNP","beta.exposure","var.exposure")]
  ### note that type, s and N should go in as numbers not vector from dat
  
  DS1=list(beta=dat$beta.exposure,varbeta=dat$var.exposure,snp=dat$SNP,position=dat$pos.exposure,pvalues=dat$pval.exposure,type="quant",N=samplesize.exposure,MAF=dat$eaf.exposure,LD=Mreord)
  DS2=list(beta=dat$beta.outcome,varbeta=dat$var.outcome,snp=dat$SNP,position=dat$pos.exposure,pvalues=dat$pval.outcome,type="quant",N=samplesize.outcome, MAF=dat$eaf.exposure,LD=Mreord)
  
  ### Mreord renamed LD (somewhere) - wait its in DS1 and DS2
  # COLOC ABF
  my.res <- coloc.abf(dataset1=DS1, dataset2=DS2)
  class(my.res)
  my.res
  sensitivity(my.res,"H4 > 0.9")
  # CHECKING DSs for LD matrix
  check_dataset(DS1,req="LD") # NULL means ok
  check_dataset(DS2,req="LD") # null this is with using eaf.exposure for both
  
  
  ### susie
  S1=runsusie(DS1,maxit=1000, estimate_prior_variance=FALSE)  # ok if converged = TRUE - can incresae maxit
  S2=runsusie(DS2,maxit=1000, estimate_prior_variance=FALSE)
  
  
  summary(S1)
  summary(S2)
  # well this is it
  susie1 <- if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(S1,S2)
    print(susie.res$summary)
    print(susie.res)
  }
  susie1
  res.summary <- susie1$summary
  res.summary$id <- filename
 write.csv(res.summary,file=paste0(filename, "susie.results.csv"))
  
  
}  

#, mc.cores =6)


for (f in filenames) {
  print(f)
  analyze6(f)
  
}

######################################################################
######################################################################
#######bindrows ... ########################################################
filenamesIVW <- list.files(path = "~/path to results/", 
                           pattern = "results.csv",
                           full.names = TRUE)
#
filenames <- filenamesIVW[1:length(filenamesIVW)]
filenames
datalist = list()

for(i in 1:length(filenamesIVW)) {
  read <- readr::read_csv(file = filenamesIVW[[i]])
  read$i <- i
  datalist[[i]] <- read
  
}
big_IVW = do.call(dplyr::bind_rows, datalist)
write.csv(big_IVW, file="~/path to results/RESULTS.susie.csv")

