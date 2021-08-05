############ Mecp2 Fecal Microbiome Analysis #################

### loading required packages ###
library("limma")
library("phyloseq")
library("vegan")
library("edgeR")
library("ggplot2")
library("tidyr")

### setting working directory ###

setwd("/Users/karineier/Documents/Microbiome-Analyses/Output") 

### loading phyloseq data ###

load("/Users/karineier/Documents/Mecp2/Microbiome/phyloseq_nochim_silva.RData") 

ps = ps.silva.nochim

# remove Zymo/mock 
ps = subset_samples(ps, grepl("^JL.*Week", sample_names(ps)))
ps = prune_taxa(taxa_sums(ps)>0, ps)

head(tax_table(ps)) # checking taxa table


m = as(otu_table(ps), "matrix") # getting otu table and making it into a matrix - note: these are really ASVs
m = m + 1 # adding one to prevent log(0) issues
m = t(m) # transforming data so that columns = samples

m = m[,-c(50,53,65,223,242,682)] ## removing samples with very low ASV/read counts

taxonomy = tax_table(ps, errorIfNULL=FALSE)
if( !is.null(taxonomy)) {
  taxonomy = data.frame(as(taxonomy, "matrix"))
}

d = DGEList(counts=m, genes=taxonomy, remove.zeros=T) 
d = calcNormFactors(d, method="RLE") # normalizing ASV data using RLE method

if( !all(is.finite(d$samples$norm.factors))) {
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing 'method' argument")
}


### Preparing sample info ###

sample_info = read.csv("/Users/karineier/Documents/Mecp2/Microbiome/sample_info.csv") ## loading in sample data, separate from phyloseq data

pdata = sample_info[-c(50,53,65,223,242,682),] # removing samples with very low ASV/read counts
rownames(pdata) = pdata$Sample_ID

#### Setting up the model ####
### this model uses Genotype_sex_time which is a factor variable unique for each combo of Genotype, sex, and time and also includes a variable for the cross-foster dam to control for within litter effects ###

mm = model.matrix(~ 0 + Genotype_sex_time + CF_Dam, data = pdata)

# voom
y <- voom(d, mm, plot=T)

# estimate within-mouse correlations for longitudinal analysis #
cor <- duplicateCorrelation(y, mm, block=pdata$Mouse_ID)$consensus

# fitting the model

fit <- lmFit(y, mm, correlation = cor, block = pdata$Mouse_ID) 

### estimate mutant vs. wt at each age (cross-sectional analysis) ###
tidyr::crossing(age = levels(as.factor(pdata$Age_weeks)), 
                sex = c("M", "F")) %>%
  dplyr::slice(-c(16, 18, 20)) %>% #manually removing rows with ages that don't exist for males
  
  
  purrr::pwalk(function(age, sex) {
    
    contrast <- glue::glue("Genotype_sex_time{sex}_het_{age}-Genotype_sex_time{sex}_wt_{age}")
    contr <- do.call(makeContrasts, list(contrasts = contrast, levels = colnames(coef(fit))))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    tmp2 <- topTable(tmp, sort.by="p", n=Inf)
    tmp2$Taxa <- rownames(tmp2)
    tmp2 <- tmp2[,c("Taxa", "logFC", "AveExpr", "P.Value", "adj.P.Val")]
    cat("Number of differentially abundant taxa at age", age, "weeks in ", ifelse(sex=="F", "Females", "Males"), "is", length(which(tmp2$adj.P.Val<0.05)), "\n")
    sigs <- tmp2$Taxa[which(tmp2$adj.P.Val <0.05)]
    sigtab = cbind(as(tmp2, "data.frame"), as(tax_table(ps)[rownames(tmp2),], "matrix"))
    filename <- glue::glue("Genotype_Effect_Age_{age}_Weeks_{sex}.txt")
    write.table(file=filename, sigtab, col.names=T, row.names=F, sep="\t", quote=F)
    
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname="Set1", ...) {
      scale_fill_brewer(palette=palname, ...)
    }
    sigtabgen = subset(sigtab, !is.na(Genus))
    
    #Phylum order
    x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
    x = sort(x, TRUE)
    sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
    
    # Genus order
    x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
    x = sort(x, TRUE)
    sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
    p <- ggplot(sigtabgen, aes(x=Genus, y=logFC, color=Phylum)) + geom_point(size=2) +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size=7))
    p <- p + ggtitle(label = paste(age, "Weeks", ifelse(sex=="F", "Females", "Males"))) +
      theme(legend.position = "bottom")
    print(p)
    tiff(file=glue::glue("Fold_Change_Genotype_Effect_Age_{age}_Weeks_{sex}_by_Genus.tiff"), res=400, height=6, width=12, units="in")
    print(p)
    dev.off()
    return(sigs)
  })

###################################################################################################
##### Refit model for estimation of genotype and genotype x age as main effect longitudinally #####
###################################################################################################

# separate male and female data since males have fewer time points #

pdata$Age_weeks = as.factor(pdata$Age_weeks)

longitudinalModel = function(sex) {
  index = which(pdata$Sex == sex)
  pdata.0 = pdata[index,]
  pdata.0 = droplevels(pdata.0)
  m.0 = m[,index]
  
  d.0 = DGEList(counts=m.0, genes=taxonomy, remove.zeros=T)
  d.0 = calcNormFactors(d.0, method="RLE") # normalizing ASV data using RLE method
  
  if( !all(is.finite(d.0$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing 'method' argument")
  }
  
  mm.0 <- model.matrix( ~ Age_weeks*Genotype_cat + CF_Dam, data = pdata.0, contrasts=list(Age_weeks="contr.sum"))
  y.0 = voom(d.0, mm.0, plot=T)
  cor.0 = duplicateCorrelation(y.0, mm.0, block=pdata.0$Mouse_ID)$consensus
  
  fit.0 = lmFit(y.0, mm.0, correlation=cor.0, block=pdata.0$Mouse_ID)
  head(coef(fit.0))
  
  ### F-test of interaction coefficients ###
  
  pdata.0$CF_Dam = as.factor(pdata.0$CF_Dam)
  n.0 = 1 + nlevels(pdata.0$Age_weeks) + nlevels(pdata.0$CF_Dam)
  n.1 = ncol(coef(fit.0))
  tmp <- contrasts.fit(fit.0, coef=n.0:n.1)
  tmp <- eBayes(tmp)
  tmp2 <- topTable(tmp, sort.by="F", n=Inf)
  tmp2$Taxa <- rownames(tmp2)
  tmp2 <- tmp2[,c("Taxa", "AveExpr", "F", "P.Value", "adj.P.Val")]
  cat("Number of taxa with significant time-genotype interaction in", ifelse(sex=="F", "Females", "Males"), "is", length(which(tmp2$adj.P.Val<0.05)), "\n")
  
  #sigs <- c(sigs, tmp2$Taxa[which(tmp2$adj.P.Val <0.05)])
  sigtab = cbind(as(tmp2, "data.frame"), as(tax_table(ps)[rownames(tmp2), ], "matrix"))
  filename <- glue::glue("Time_Genotype_Interaction_{sex}.txt")
  write.table(file=filename, sigtab, col.names=T, row.names=F, sep="\t", quote=F)
  
  ### contrasts for genotype as main effect ###
  
  n.2 = nlevels(pdata.0$Age_weeks) + 1
  tmp.0 <- contrasts.fit(fit.0, coef = n.2)
  tmp.0 <- eBayes(tmp.0)
  tmp2.0 <- topTable(tmp.0, sort.by="P", n=Inf)
  tmp2.0$Taxa <- rownames(tmp2.0)
  tmp2.0 <- tmp2.0[,c("Taxa", "AveExpr", "B", "P.Value", "adj.P.Val")]
  cat("Number of taxa with significant main effect of genotype in", ifelse(sex=="F", "Females", "Males"), "is", length(which(tmp2.0$adj.P.Val<0.05)), "\n")
  
  sigtab.0 = cbind(as(tmp2.0, "data.frame"), as(tax_table(ps)[rownames(tmp2.0), ], "matrix"))
  filename.0 <- glue::glue("Main_Effect_Genotype_{sex}.txt")
  write.table(file=filename.0, sigtab.0, col.names=T, row.names=F, sep="\t", quote=F)
  
}

sexes = c("F", "M")

lapply(sexes, longitudinalModel)

