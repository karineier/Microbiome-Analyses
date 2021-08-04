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



