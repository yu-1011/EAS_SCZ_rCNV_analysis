
## Load the R-packages
library("lme4qtl")
library("data.table")

## Read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1] 
input <- args[2] 
gene_file <- args[3]
phefile <- args[4] 
pi_hat_matrix <- args[5]
loci_file <- args[6]
output <- args[7]
perm <- 100

## Set working directory
setwd(wdir)

## Function to read and merge phenotype and fam files
read_and_merge_data <- function(input, phefile) {
  phe <- read.delim(phefile, header = TRUE, sep = " ")
  phe$matchID <- with(phe, paste(FID, ':', IID, sep = ''))
  phe$aff <- phe$AFF - 1
  
  fam <- read.table(paste0(input, ".fam"))
  colnames(fam) <- c("FID", "IID", "Father_ID", "Mother_ID", "Sex", "Affect")
  fam$matchID <- with(fam, paste(FID, ':', IID, sep = ''))
  
  merge(fam, phe, by = 'matchID')
}

## Function to print the start of the script
print_start <- function() {
  cat("Script started.
")
  cat(paste0("Working Directory: ", wdir, "
"))
  cat(paste0("Input File: ", input, "
"))
}

## Function to create ID files
create_id_files <- function(combined, data_name) {
  for (i in seq_along(data_name)) {
    ID <- combined[combined$CNV_platform == data_name[i], 2:3]
    write.table(ID, paste0(data_name[i], '.ID'), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '	')
  }
}

## Read and prepare data
print_start()
combined <- read_and_merge_data(input, phefile)
dat <- unique(sort(combined$CNV_platform))
data_name <- c('combined', as.character(dat))

## Create ID files for PLINK
create_id_files(combined, data_name)

## Load additional files and initialize matrices
GRM <- as.matrix(fread(pi_hat_matrix), rownames = 1)
idx <- intersect(rownames(GRM),combined$IID.x)
pi_h_mx <- GRM[idx,idx]

loci <- read.table(loci_file)

##set the start of datasets loop
startpoint=1
if(length(data_name)==2){
	startpoint = 2
}

## write out PLINK ID files 
for (i in 1:1) {
    ID <- combined[combined$CNV_platform==data_name[i],2:3]
    write.table(ID,paste(data_name[i],'.ID',sep=''),col=F,row=F,quo=F,sep='\t')
}

## set up loop parameter
type_name <- c('allCNV','del','dup')
region_name <- c(as.character(loci$V4))
size_name <- c('>20kb','>100kb','>200kb','>300kb','>400kb','>500kb','>600kb') 
freq_name <- c('Allfreq','singleton','2-5','6-10','11-20','21-40','41-80','81+')

data <- c('',paste(' --keep ', wdir,'/',data_name[2:length(data_name)],'.ID',sep=''))

type <- c('',
          '--cnv-del',
          '--cnv-dup')

region <- c(paste0('--cnv-intersect /home/yu/work/proj-cnv/burden_test/burden_known_loci/misc/hg19_loci/',loci$V4,'_region.txt'))

size <- c('',
          '--cnv-kb 100',
          '--cnv-kb 200',
          '--cnv-kb 300',
          '--cnv-kb 400',
          '--cnv-kb 500',
          '--cnv-kb 600') 

freq <- c('',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 1',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 5 --cnv-freq-exclude-below 2',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 10 --cnv-freq-exclude-below 6',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 20 --cnv-freq-exclude-below 11',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 40 --cnv-freq-exclude-below 21',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 80 --cnv-freq-exclude-below 41',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-below 81')

## PART 2: Variables: Descriptives and Statistics

## Note: we can get gene count numbers with --cnv-count

## Files:
##  - .cnv.indiv
##  - .cnv.grp.summary
##  - .cnv.summary.mperm
##  - .cnv with gene annotation
tmp_data <- list()
X <- 0 ## loop COUNTER
data_set <- NA; region_set <- NA; CNV_type <- NA; CNV_freq <- NA; CNV_size <- NA

## PLINK descriptives and stats:
FULL_SAMPLE_SIZE <- NA; FULL_cas_SAMPLE_SIZE <- NA; FULL_con_SAMPLE_SIZE <- NA; CNV_CARRIER_SAMPLE_SIZE <- NA; CNV_cas_CARRIER_SAMPLE_SIZE <- NA; CNV_con_CARRIER_SAMPLE_SIZE <- NA
COUNT_rate <- NA; COUNT_cas_rate <- NA; COUNT_con_rate <- NA; COUNT_cascon_ratio <- NA; COUNT_glm_OR <- NA; COUNT_glm_se <- NA; COUNT_glm_tval <- NA; COUNT_glm_pval <- NA; COUNT_perm_pval <- NA; COUNT_glm_lowerCI <- NA; COUNT_glm_upperCI <- NA
NGENE_rate <- NA; NGENE_cas_rate <- NA; NGENE_con_rate <- NA; NGENE_cascon_ratio <- NA; NGENE_glm_OR <- NA; NGENE_glm_se <- NA; NGENE_glm_tval <- NA; NGENE_glm_pval <- NA; NGENE_perm_pval <- NA; NGENE_glm_lowerCI <- NA; NGENE_glm_upperCI <- NA
KB_rate <- NA; KB_cas_rate <- NA; KB_con_rate <- NA; KB_cascon_ratio <- NA; KB_glm_OR <- NA; KB_glm_se <- NA; KB_glm_tval <- NA; KB_glm_pval <- NA; KB_perm_pval <- NA; KB_glm_lowerCI <- NA; KB_glm_upperCI <- NA
NSEG <- NA; NSEG_GENIC <- NA; NSEG_NONGENIC <- NA; CNV_GENIC_CARRIER_SAMPLE_SIZE <- NA; CNV_NONGENIC_CARRIER_SAMPLE_SIZE <- NA

## NSEG
NSEG_rate <- NA; NSEG_cas_rate <- NA; NSEG_con_rate <- NA; NSEG_cascon_ratio <- NA; NSEG_glm_OR <- NA; NSEG_glm_se <- NA; NSEG_glm_tval <- NA; NSEG_glm_pval <- NA; NSEG_perm_pval <- NA; NSEG_glm_lowerCI <- NA; NSEG_glm_upperCI <- NA 

## GENIC CNV
NSEG_GENIC_rate <- NA; NSEG_GENIC_cas_rate <- NA; NSEG_GENIC_con_rate <- NA; NSEG_GENIC_cascon_ratio <- NA; NSEG_GENIC_glm_OR <- NA; NSEG_GENIC_glm_se <- NA; NSEG_GENIC_glm_tval <- NA; NSEG_GENIC_glm_pval <- NA; NSEG_GENIC_perm_pval <- NA; NSEG_GENIC_glm_lowerCI <- NA; NSEG_GENIC_glm_upperCI <- NA
KB_GENIC_rate <- NA; KB_GENIC_cas_rate <- NA; KB_GENIC_con_rate <- NA; KB_GENIC_cascon_ratio <- NA; KB_GENIC_glm_OR <- NA; KB_GENIC_glm_se <- NA; KB_GENIC_glm_tval <- NA; KB_GENIC_glm_pval <- NA; KB_GENIC_perm_pval <- NA; KB_GENIC_glm_lowerCI <- NA; KB_GENIC_glm_upperCI <- NA

## FILTER sets:
for(a in 1:length(data_name)){
for(b in 1:length(type_name)){
for(c in 1:length(region_name)){
for(d in 1:length(freq_name)){
for(e in 1:length(size_name)){

## PART 3: CNV burden loop
## getting full count of burden analyses
X <- X+1
tot <- length(data_name)*length(type_name)*length(region_name)*length(freq_name)*length(size_name)

## PART 3.1: Setting up burden loop

## create temporary subdirectory within working directory (will not overwrite existing directory)
system(paste0("mkdir -p burden_loop"))


## assign specific loop parameters
data2 <- data[a]
## combining CNV type and region filtering
type2 <- type[b]
region2 <- region[c]
freq2 <- freq[d]
size2 <- size[e]

data_set[X] <- data_name[a]
CNV_type[X] <- type_name[b]
region_set[X] <- region_name[c]
CNV_freq[X] <- freq_name[d]
CNV_size[X] <- size_name[e]

## write out status report
cat('iteration =',X,'of',tot,'\n','data_set',a,'=',data_name[a],'\n','CNV_type',b,'=',type_name[b],'\n','region_set',c,'=',region_name[c],'\n','CNV_freq',d,'=',freq_name[d],'\n','CNV_size',e,'=',size_name[e],'\n',file=paste(output,".status",sep=''),append=F)

## PART 3.2: Applying PLINK filters


## =============== Filter by data, CNV type, and CNV regions first

  system(paste("plink --noweb --cfile ",input," ",data2," ",type2," ",region2," --cnv-write --out burden_loop/type_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/type_loop --cnv-make-map --out burden_loop/type_loop")


  ## --- Frequency pruning    
  system(paste("plink --noweb --cfile burden_loop/type_loop ",freq2," --cnv-write --out burden_loop/freq_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/freq_loop --cnv-make-map --out burden_loop/freq_loop")

  ## --- Size pruning
  system(paste("plink --noweb --cfile burden_loop/freq_loop ",size2," --cnv-write --out burden_loop/size_loop",sep=""))
  system("plink --noweb --cfile burden_loop/size_loop --cnv-make-map --out burden_loop/size_loop")
 
  ## Adding gene/exon count separately

  system(paste("plink --noweb --cfile burden_loop/size_loop --cnv-indiv-perm --mperm ",perm," --cnv-count ",gene_file," --out burden_loop/burden_loop",sep=""))
    
  system(paste("plink --noweb --cfile burden_loop/size_loop --cnv-freq-method2 0.5 --cnv-write-freq --cnv-track"," --cnv-count ",gene_file," --out burden_loop/burden_loop_", data_name[a],"_", type_name[b],"_", region_name[c],"_", freq_name[d],"_", size_name[e], sep=""))


## PART 3.3: Checking CNV burden coverage

## Check if any segments were matched

locus <- loci[c,4]
temp <- tryCatch(   
                     { cnv <- read.table('burden_loop/size_loop.cnv',h=T) 


## PART 3.4: Reading PLINK burden results

if (nrow(cnv) > 0) {
  
## == PLINK results
plink_sum <- read.table("burden_loop/burden_loop.cnv.grp.summary",h=T)
plink_emp <- read.table("burden_loop/burden_loop.cnv.summary.mperm",h=T)

## PART 3.5: CNV burden analysis in R

indiv <- read.table("burden_loop/burden_loop.cnv.indiv",h=T)
indiv$matchID <- paste(indiv$FID,':',indiv$IID,sep='')

cnv <- read.table("burden_loop/size_loop.cnv",h=T)
cnv$matchID <- paste(cnv$FID,':',cnv$IID,sep='')
cnv$bp <- cnv$BP2 - cnv$BP1

## split into genic and non-genic CNVs
genic <- subset(cnv,cnv$SCORE > 0)
nongenic <- subset(cnv,cnv$SCORE == 0)

# -- sum genes for each individual
gene.cnt <- tapply(cnv$SCORE,factor(cnv$matchID),sum)
indiv$NGENE <- 0
# Apply changes only to indexed subset
indx <- match(names(gene.cnt),indiv$matchID)
indiv$NGENE <- replace(indiv$NGENE,indx,gene.cnt)

## -- sum genic CNV count and KB for each individual
genic.tbl <- table(genic$matchID)
indiv$GENIC_CNV_COUNT <- 0
indx <- match(names(genic.tbl),indiv$matchID)
indiv$GENIC_CNV_COUNT <- replace(indiv$GENIC_CNV_COUNT,indx,genic.tbl)

genic.tbl <- tapply(genic$bp,genic$matchID,sum)/1000
indiv$GENIC_KB <- 0
indx <- match(names(genic.tbl),indiv$matchID)
indiv$GENIC_KB <- replace(indiv$GENIC_KB,indx,genic.tbl)


## -- merge with phenotype
comrg <- merge(indiv,phe,by='matchID')
comrg_cnv <- comrg
comrg_genic <- comrg
comrg_nongenic <- comrg
#comrg_cnv <- comrg[comrg$NSEG > 0,]
#comrg_genic <- comrg[comrg$GENIC_CNV_COUNT > 0,]
#comrg_nongenic <- comrg[comrg$NONGENIC_CNV_COUNT > 0,]
rownames(comrg) <- comrg$IID.x
idx <- intersect(rownames(GRM), rownames(comrg))
comrg <- comrg[idx,]
pi_h_mx <- GRM[idx,idx]
    
## ==== get SAMPLE_SIZES
FULL_SAMPLE_SIZE[X] <- nrow(comrg) 
FULL_cas_SAMPLE_SIZE[X] <- nrow(comrg[comrg$PHE==2,])
FULL_con_SAMPLE_SIZE[X] <- nrow(comrg[comrg$PHE==1,])

CNV_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_cnv)
CNV_cas_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_cnv[comrg_cnv$PHE==2,])
CNV_con_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_cnv[comrg_cnv$PHE==1,])

    
CNV_GENIC_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_genic)
CNV_NONGENIC_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_nongenic)

NSEG[X] <- sum(comrg$NSEG)
NSEG_GENIC[X] <- sum(comrg$GENIC_CNV_COUNT)
NSEG_NONGENIC[X] <- sum(comrg$NONGENIC_CNV_COUNT)

## ==== Genes covered
NGENE_rate[X] <- mean(comrg_cnv$NGENE)
NGENE_cas_rate[X] <- mean(comrg_cnv$NGENE[comrg_cnv$PHE==2])
NGENE_con_rate[X] <- mean(comrg_cnv$NGENE[comrg_cnv$PHE==1])
NGENE_cascon_ratio[X] <- NGENE_cas_rate[X]/NGENE_con_rate[X]

if(sum(comrg_cnv$aff==0)==0 | sum(comrg_cnv$aff==1)==0){
  NGENE_glm_OR[X] <- NA
  NGENE_glm_se[X] <- NA
  NGENE_glm_tval[X] <- NA
  NGENE_glm_pval[X] <- NA
  NGENE_glm_lowerCI[X] <- NA
  NGENE_glm_upperCI[X] <- NA }

if(sum(comrg_cnv$aff==0) > 0 & sum(comrg_cnv$aff==1) > 0){

if (data_set[X]=='combined') { 
    #NGENE.lm <- glm(aff ~ NGENE + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5,data=comrg_cnv,family='binomial') 
    NGENE.lm <- relmatGlmer(aff ~  NGENE + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)
    }
if (data_set[X]!='combined') { 
    #NGENE.lm <- glm(aff ~ NGENE + SEX + C1 + C2 + C3 + C4 + C5,data=comrg_cnv,family='binomial')
    NGENE.lm <- relmatGlmer(aff ~  NGENE + SEX + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)
}

NGENE.mod <- summary(NGENE.lm)
NGENE_glm_OR[X] <- exp(NGENE.mod$coefficients[2,1])    
NGENE_glm_se[X] <- NGENE.mod$coefficients[2,2]
NGENE_glm_tval[X] <- NGENE.mod$coefficients[2,3]
NGENE_glm_pval[X] <- NGENE.mod$coefficients[2,4]
NGENE_glm_lowerCI[X] <- exp(NGENE.mod$coefficients[2,1] - (1.96*NGENE_glm_se[X]))
NGENE_glm_upperCI[X] <- exp(NGENE.mod$coefficients[2,1] + (1.96*NGENE_glm_se[X])) }

## ==== CNV Count
NSEG_GENIC_rate[X] <- mean(comrg$GENIC_CNV_COUNT)
NSEG_GENIC_cas_rate[X] <- mean(comrg$GENIC_CNV_COUNT[comrg$PHE==2])
NSEG_GENIC_con_rate[X] <- mean(comrg$GENIC_CNV_COUNT[comrg$PHE==1])
NSEG_GENIC_cascon_ratio[X] <- NSEG_GENIC_cas_rate[X]/NSEG_GENIC_con_rate[X]

if (data_set[X]=='combined') { 
    #NSEG_GENIC.lm <- glm(aff ~ GENIC_CNV_COUNT + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5,data=comrg,family='binomial') 
    NSEG_GENIC.lm <- relmatGlmer(aff ~  GENIC_CNV_COUNT + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)}
if (data_set[X]!='combined') { 
    #NSEG_GENIC.lm <- glm(aff ~ GENIC_CNV_COUNT + SEX + C1 + C2 + C3 + C4 + C5,data=comrg,family='binomial') 
    NSEG_GENIC.lm <- relmatGlmer(aff ~  GENIC_CNV_COUNT + SEX + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)}


NSEG_GENIC.mod <- summary(NSEG_GENIC.lm)
if(nrow(NSEG_GENIC.mod$coefficients)==1){
  NSEG_GENIC_glm_OR[X] <- NA
  NSEG_GENIC_glm_se[X] <- NA
  NSEG_GENIC_glm_tval[X] <- NA
  NSEG_GENIC_glm_pval[X] <- NA
  NSEG_GENIC_glm_lowerCI[X] <- NA
  NSEG_GENIC_glm_upperCI[X] <- NA }
if(nrow(NSEG_GENIC.mod$coefficients) > 1){
NSEG_GENIC_glm_OR[X] <- exp(NSEG_GENIC.mod$coefficients[2,1])    
NSEG_GENIC_glm_se[X] <- NSEG_GENIC.mod$coefficients[2,2]
NSEG_GENIC_glm_tval[X] <- NSEG_GENIC.mod$coefficients[2,3]
NSEG_GENIC_glm_pval[X] <- NSEG_GENIC.mod$coefficients[2,4]
NSEG_GENIC_glm_lowerCI[X] <- exp(NSEG_GENIC.mod$coefficients[2,1] - (1.96*NSEG_GENIC_glm_se[X]))
NSEG_GENIC_glm_upperCI[X] <- exp(NSEG_GENIC.mod$coefficients[2,1] + (1.96*NSEG_GENIC_glm_se[X])) }

## ==== KB burden
comrg$GENIC_KB_100 <- (comrg$GENIC_KB)/100
KB_GENIC_rate[X] <- mean(comrg_genic$GENIC_KB)
KB_GENIC_cas_rate[X] <- mean(comrg_genic$GENIC_KB[comrg_genic$PHE==2])
KB_GENIC_con_rate[X] <- mean(comrg_genic$GENIC_KB[comrg_genic$PHE==1])
KB_GENIC_cascon_ratio[X] <- KB_GENIC_cas_rate[X]/KB_GENIC_con_rate[X]

if(sum(comrg_genic$aff==1)==0 | sum(comrg_genic$aff==1)==0){
  KB_GENIC_glm_OR[X] <- NA
  KB_GENIC_glm_se[X] <- NA
  KB_GENIC_glm_tval[X] <- NA
  KB_GENIC_glm_pval[X] <- NA
  KB_GENIC_glm_lowerCI[X] <- NA
  KB_GENIC_glm_upperCI[X] <- NA }

if(sum(comrg_genic$aff==1) > 0 & sum(comrg_genic$aff==1) > 0){

if (data_set[X]=='combined') { 
    #KB_GENIC.lm <- glm(aff ~ GENIC_KB + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5,data=comrg_genic,family='binomial')
    KB_GENIC.lm <- relmatGlmer(aff ~  GENIC_KB_100 + SEX + CNV_platform + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)}
if (data_set[X]!='combined') { 
    #KB_GENIC.lm <- glm(aff ~ GENIC_KB + SEX + C1 + C2 + C3 + C4 + C5,data=comrg_genic,family='binomial') 
    KB_GENIC.lm <- relmatGlmer(aff ~  GENIC_KB_100 + SEX + C1 + C2 + C3 + C4 + C5 + (1|IID.x), comrg, relmat = list(IID.x = pi_h_mx), family = binomial)}

KB_GENIC.mod <- summary(KB_GENIC.lm)
KB_GENIC_glm_OR[X] <- exp(KB_GENIC.mod$coefficients[2,1])    
KB_GENIC_glm_se[X] <- KB_GENIC.mod$coefficients[2,2] 
KB_GENIC_glm_tval[X] <- KB_GENIC.mod$coefficients[2,3]
KB_GENIC_glm_pval[X] <- KB_GENIC.mod$coefficients[2,4]
KB_GENIC_glm_lowerCI[X] <- exp(KB_GENIC.mod$coefficients[2,1] - (1.96*KB_GENIC_glm_se[X]))
KB_GENIC_glm_upperCI[X] <- exp(KB_GENIC.mod$coefficients[2,1] + (1.96*KB_GENIC_glm_se[X])) }


} ## End of mapped segments analysis 

tmp_data[[X]] <- cbind.data.frame(data_set=data_set[X],
                              region_set=region_set[X],
                              CNV_type=CNV_type[X],
                              CNV_freq=CNV_freq[X],
                              CNV_size=CNV_size[X],
                              FULL_SAMPLE_SIZE=FULL_SAMPLE_SIZE[X],
                              FULL_cas_SAMPLE_SIZE=FULL_cas_SAMPLE_SIZE[X],
                              FULL_con_SAMPLE_SIZE=FULL_con_SAMPLE_SIZE[X],
                              CNV_CARRIER_SAMPLE_SIZE=CNV_CARRIER_SAMPLE_SIZE[X],
                              CNV_cas_CARRIER_SAMPLE_SIZE=CNV_cas_CARRIER_SAMPLE_SIZE[X],
                              CNV_con_CARRIER_SAMPLE_SIZE=CNV_con_CARRIER_SAMPLE_SIZE[X],
                              CNV_GENIC_CARRIER_SAMPLE_SIZE=CNV_GENIC_CARRIER_SAMPLE_SIZE[X],
                              
                              NGENE_rate=NGENE_rate[X],    
                              NGENE_con_rate=NGENE_con_rate[X],
                              NGENE_cas_rate=NGENE_cas_rate[X],
                              NGENE_glm_OR=NGENE_glm_OR[X],
                              NGENE_glm_se=NGENE_glm_se[X],
                              NGENE_glm_tval=NGENE_glm_tval[X],
                              NGENE_glm_pval=NGENE_glm_pval[X],
                              NGENE_glm_lowerCI=NGENE_glm_lowerCI[X],
                              NGENE_glm_upperCI=NGENE_glm_upperCI[X],
                                  
                              KB_GENIC_rate=KB_GENIC_rate[X],
                              KB_GENIC_con_rate=KB_GENIC_con_rate[X],
                              KB_GENIC_cas_rate=KB_GENIC_cas_rate[X],
                              KB_GENIC_glm_OR=KB_GENIC_glm_OR[X],
                              KB_GENIC_glm_se=KB_GENIC_glm_se[X],
                              KB_GENIC_glm_tval=KB_GENIC_glm_tval[X],
                              KB_GENIC_glm_pval=KB_GENIC_glm_pval[X],
                              KB_GENIC_glm_lowerCI=KB_GENIC_glm_lowerCI[X],
                              KB_GENIC_glm_upperCI=KB_GENIC_glm_upperCI[X],
                              
                              NSEG_GENIC_rate=NSEG_GENIC_rate[X],    
                              NSEG_GENIC_con_rate=NSEG_GENIC_con_rate[X],
                              NSEG_GENIC_cas_rate=NSEG_GENIC_cas_rate[X],
                              NSEG_GENIC_glm_OR=NSEG_GENIC_glm_OR[X],
                              NSEG_GENIC_glm_se=NSEG_GENIC_glm_se[X],
                              NSEG_GENIC_glm_tval=NSEG_GENIC_glm_tval[X],
                              NSEG_GENIC_glm_pval=NSEG_GENIC_glm_pval[X],
                              NSEG_GENIC_glm_lowerCI=NSEG_GENIC_glm_lowerCI[X],
                              NSEG_GENIC_glm_upperCI=NSEG_GENIC_glm_upperCI[X])

## write to file
#system('rm -r burden_loop')  
                     },
                    # warning = function(w) { message('Waring @ ',locus) ; return(NA) },
                     error = function(e) { message('Error @ ',locus) ; return(NA) },
                     finally = { message('next...') }
                     )

}
}
}
}
}

## Wrap-up and save results
full_data <- do.call(rbind, tmp_data)
write.table(full_data, paste0(output, '.burden'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '	')
cat("Script completed.
")
