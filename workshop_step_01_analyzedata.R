# =============================================================================
#
# Serena Visit (2024)
#
# Microbiome workshop
# 
# =============================================================================
# load libs
# > tidy-R
library(tidyverse)
# > V-table summary tables
library(vtable)
# > GT summary table
library(gtsummary)
# > extra stuff for ggplot
library(ggpubr)
library(ggforce)
library(cowplot)
# > misc. statistics functions
library(coin)

# ==========================================
# ==========================================
#      MAIN
# ==========================================
# ==========================================

# set WD
setwd('~/UMCG/2024_Serena/Case_study_IBD/')
#  load extra codes
source('./helper_functions/R_Misc.R')
source('./helper_functions/R_Microbiome_scripts.R')
#   > modified summary table
source('./helper_functions/R_summarytable.R')

# =============================================================================
# =============================================================================
#   ANALYSIS: IBD vs HC microbiome data analysis
# =============================================================================

# =============================================================================
# LOAD AND PREP DATA
# =============================================================================
# > microbiome
inTaxa <- read.table('data_ready/mp406_ready.csv',sep=',',header = T) %>% column_to_rownames('ID')
# > metadata
inMeta <- read.table('data_ready/metadata.csv',sep=',',header = T) %>% column_to_rownames('ID')

# > humann: load & align to metadata
inHum <- read.table('data_ready/hum306_ready.csv',sep=',',header = T) %>% column_to_rownames('ID')

# > merge data layers
inDFm <- inMeta %>% merge(inTaxa,by='row.names') %>% column_to_rownames('Row.names')
inDFm <- inDFm %>% merge(inHum,by='row.names') %>% column_to_rownames('Row.names')

# --- small fixes and tweaks ---
inMeta$Diagnosis <- inMeta$Diagnosis %>% factor(levels = c('HC','UC','CD'))
inMeta$Diagnosis.2 <- NA
inMeta$Diagnosis.2[inMeta$Diagnosis %in% c("HC")] <- "HC"
inMeta$Diagnosis.2[inMeta$Diagnosis %in% c("UC","CD")] <- "IBD"

# =============================================================================
# 1. DATA QC & SUMMARY STATS [metadata by diagnosis]
# =============================================================================
# > clean & re-create output folder(s)
outF <- 'results_summaries'
unlink(outF,recursive = T)
if (!dir.exists(outF)) {dir.create(outF)}

# vtable (https://cran.r-project.org/web/packages/vtable/vignettes/sumtable.html)
sumstatsIBD <- inMeta %>% 
  mutate(Diagnosis.2 = as.factor(as.character(Diagnosis.2))) %>%
  #select(-any_of(c('Diagnosis'))) %>%
  rst(group = 'Diagnosis.2', group.test = TRUE, 
      summ = c("length(x)","mean(x)", "sd(x)"),
      summ.names = c("N","Mean","SD"),
      file = paste0(outF,'/summarystats_phenos_IBD_vtable.html'))

# gt summary table
sumstatsIBD <- inMeta %>%
  #select(-all_of(c('Diagnosis'))) %>%
  tbl_summary(by='Diagnosis.2',
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  add_p() %>%
  add_overall() %>% as_gt() %>% gt::gtsave(paste0(outF,'/summarystats_phenos_IBD_gt.html'))

# check NAs
for (cc in colnames(inMeta)) {
  nNAs <- sum(is.na(inMeta[[cc]]))
  nN <- nrow(inMeta)
  if (nNAs > 0) {
    print(paste0('>> ',cc,': ',class(inMeta[[cc]]),' : NAs: ',nNAs,'/',nN,' = ',round(nNAs/nN*100,2),'%'))
  }
}

# handle NAs: 
inMetaImp <- inMeta
# > drop data with NA for Diagnosis
inMetaImp <- inMetaImp %>% filter(!is.na(Diagnosis))
# > impute to the median for numeric variables
#  >> fix data types
inMetaImp$Fcal <- as.numeric(inMetaImp$Fcal)
inMetaImp$PFReads <- as.numeric(inMetaImp$PFReads)
#  >> replace NAs with median
inMetaImp <- inMetaImp %>% mutate(across(c(Age,BMI,Fcal,PFReads), 
                                         ~replace_na(., median(., na.rm=TRUE))))
#  >> replace NAs in medication use with N
inMetaImp <- inMetaImp %>% mutate(across(c(PPI,Laxatives,Antibiotics), 
                                         ~replace_na(., 'N')))
#  >> replace NA in Bristol.Chart with "hard"
inMetaImp <- inMetaImp %>% mutate(across(c(Bristol.Chart), 
                                         ~replace_na(., 'Hard')))
#  >> replace NA in Gender with F 
inMetaImp <- inMetaImp %>% mutate(across(c(Bristol.Chart), 
                                         ~replace_na(., 'F')))

# >>> ... to be adjusted appropriate to the study!

# re-check NAs
for (cc in colnames(inMetaImp)) {
  nNAs <- sum(is.na(inMetaImp[[cc]]))
  nN <- nrow(inMetaImp)
  if (nNAs > 0) {
    print(paste0('>> ',cc,': ',class(inMetaImp[[cc]]),' : NAs: ',nNAs,'/',nN,' = ',round(nNAs/nN*100,2),'%'))
  }
}

# re-make summaries
# ======================================
# vtable (https://cran.r-project.org/web/packages/vtable/vignettes/sumtable.html)
sumstatsIBD <- inMetaImp %>% 
  mutate(Diagnosis.2 = as.factor(as.character(Diagnosis.2))) %>%
  #select(-any_of(c('Diagnosis'))) %>%
  rst(group = 'Diagnosis.2', group.test = TRUE, 
      summ = c("length(x)","mean(x)", "sd(x)"),
      summ.names = c("N","Mean","SD"),
      file = paste0(outF,'/summarystats_phenos_imputed_IBD_vtable.html'))

# gt summary table
sumstatsIBD <- inMetaImp %>%
  #select(-all_of(c('Diagnosis'))) %>%
  tbl_summary(by='Diagnosis.2',
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  add_p() %>%
  add_overall() %>% as_gt() %>% gt::gtsave(paste0(outF,'/summarystats_phenos_imputed_IBD_gt.html'))

# Microbiome QC and summaries
# ======================================
# > BASED ON TAXA:
#   >> check UNCLASSIFIED data
ggplot(inTaxa) + aes(y=UNCLASSIFIED,x=1) + geom_violin() + geom_boxplot(width=0.2,alpha=0)
#   >> check Kingdoms
ggplot(inTaxa) + aes(y=k__Bacteria,x=1) + geom_violin() + geom_boxplot(width=0.2,alpha=0)
ggplot(inTaxa) + aes(y=k__Archaea,x=1) + geom_violin() + geom_boxplot(width=0.2,alpha=0)
#   >> check all-0 samples
sum(rowSums(inTaxa) == 0)
#   >> let's get rid of problematic samples
probIDs <- rownames(inTaxa)[rowSums(inTaxa) == 0 | inTaxa$UNCLASSIFIED > 0.50 | inTaxa$k__Bacteria < 0.5]
#   >> we can't do much with NA in outcomes, drop them!
probIDs <- c(probIDs,rownames(inMeta)[is.na(inMeta$Diagnosis)])

print(probIDs)

# REMOVE BAD SAMPLES
# ===========================================================================
inTaxaRdy <- inTaxa[!(row.names(inTaxa) %in% probIDs),]
inMetaRdy <- inMetaImp[!row.names(inMetaImp) %in% probIDs,]
inHumRdy <- inHum[!row.names(inHum) %in% probIDs,]
# RE-MERGE DATA
inDFm <- inMeta %>% merge(inTaxaRdy,by='row.names') %>% column_to_rownames('Row.names')
inDFm <- inDFm %>% merge(inHumRdy,by='row.names') %>% column_to_rownames('Row.names')

# check that everything is OK:
dim(inTaxaRdy)
dim(inMetaRdy)
dim(inHumRdy)
dim(inDFm)
# length(intersect(row.names(inMetaRdy),row.names(inTaxaRdy)))
# length(intersect(row.names(inMetaRdy),row.names(inHumRdy)))
# length(intersect(row.names(inTaxaRdy),row.names(inHumRdy)))

# RE-MAKE SUMMARIES OF METADATA
# ===========================================================================
# vtable (https://cran.r-project.org/web/packages/vtable/vignettes/sumtable.html)
sumstatsIBD <- inMetaRdy %>% 
  mutate(Diagnosis.2 = as.factor(as.character(Diagnosis.2))) %>%
  #select(-any_of(c('Diagnosis'))) %>%
  rst(group = 'Diagnosis.2', group.test = TRUE, 
      summ = c("length(x)","mean(x)", "sd(x)"),
      summ.names = c("N","Mean","SD"),
      file = paste0(outF,'/summarystats_phenos_imputed_ready_IBD_vtable.html'))

# gt summary table
sumstatsIBD <- inMetaRdy %>%
  #select(-all_of(c('Diagnosis'))) %>%
  tbl_summary(by='Diagnosis.2',
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  add_p() %>%
  add_overall() %>% as_gt() %>% gt::gtsave(paste0(outF,'/summarystats_phenos_imputed_ready_IBD_gt.html'))

# make TAXONOMY & PATHWAY summaries
# ======================================
# > TAXA
inDFsumm <- merge(inTaxaRdy,inMetaRdy[,c('Diagnosis.2'),drop=F],by='row.names') %>% column_to_rownames('Row.names')
mbsum <- makeMultiSummaryByResponse(inDF = inDFsumm,
                                    response = 'Diagnosis.2',
                                    ftrs = colnames(inDFsumm),
                                    doTestVsControl = F,
                                    includeTotals = F,
                                    doSort = F,
                                    doFDR = F) %>% 
  select(-any_of(c('niceOut','niceOutMSD','FDRfor','FDRsig','ResponseVar','Var.Type')))
write.table(mbsum,paste0(outF,'/summary_taxa.csv'),sep=',',row.names = F)

# > PATHWAYS
inDFsumm <- merge(inHumRdy,inMetaRdy[,c('Diagnosis'),drop=F],by='row.names') %>% column_to_rownames('Row.names')
mbsum <- makeMultiSummaryByResponse(inDF = inDFsumm,
                                    response = 'Diagnosis',
                                    ftrs = colnames(inHum),
                                    doTestVsControl = F,
                                    includeTotals = F,
                                    doSort = F,
                                    doFDR = F) %>% 
  select(-any_of(c('niceOut','niceOutMSD','FDRfor','FDRsig','ResponseVar','Var.Type')))
write.table(mbsum,paste0(outF,'/summary_PWYs.csv'),sep=',',row.names = F)

# =============================================================================
# 2. EXPLORATORY PLOTS
# =============================================================================
outF <- 'results_exploratory'
unlink(outF,recursive = T)
if (!dir.exists(outF)) {dir.create(outF)}

# > PLOT MICROBIOME SUMMARY
#  >> step-by-step: refer to previous workshop!

#  >> species
g <- plotMicrobiomeSummary(inDFm,toPlotLvl='S',groupBy='Diagnosis',topX=14)
print(g)
ggsave(plot=g,filename = paste0(outF,'/','microbiome_summary_lvl7_species.png'),width = 7,height = 5)
# > genera
g <- plotMicrobiomeSummary(inDFm,toPlotLvl='G',groupBy='Diagnosis',topX=9)
print(g)
ggsave(plot=g,filename = paste0(outF,'/','microbiome_summary_lvl6_genus.png'),width = 7,height = 5)
# > families
g <- plotMicrobiomeSummary(inDFm,toPlotLvl='F',groupBy='Diagnosis',topX=9)
print(g)
ggsave(plot=g,filename = paste0(outF,'/','microbiome_summary_lvl5_family.png'),width = 7,height = 5)
# > phyla
g <- plotMicrobiomeSummary(inDFm,toPlotLvl='P',groupBy='Diagnosis',topX=6)
print(g)
ggsave(plot=g,filename = paste0(outF,'/','microbiome_summary_lvl2_phyla.png'),width = 7,height = 5)

# > Firmicutes / Bacteroidetes ratio (F/B)
#  >> Calculate for Controls VS IBD
inDFfb <- inDFm
#   >>> calculation [we ln transform it for happiness and plotting]
inDFfb$F.B.Ratio <- log(inDFfb$k__Bacteria.p__Firmicutes/inDFfb$k__Bacteria.p__Bacteroidetes)
hist(inDFfb$F.B.Ratio)
# > digression: simple plot & statistics
ggplot(data = inDFfb,aes(x=Diagnosis.2,y=F.B.Ratio)) + geom_boxplot()
tst.FB <- wilcox.test(inDFfb$F.B.Ratio[inDFfb$Diagnosis.2=='HC'],
            inDFfb$F.B.Ratio[inDFfb$Diagnosis.2=='IBD'])
tst.FB

# > now let's automate it and make the plot nicer
g <- makePlotComparison(inDFfb,xVar = 'Diagnosis.2',yVar = 'F.B.Ratio',lbl = 'p.nice',returnData = T)
g[[1]] <- g[[1]] + ylab('Ln(Firmicutes / Bacteroides ratio)')
print(g)
ggsave(plot=g[[1]],filename = paste0(outF,'/Exploratory_FBratio.png'),height = 6,width = 4)
write.table(g[[2]],paste0(outF,'/FExploratory_FBratio.csv'),sep=',',row.names = F)
#  >> Calculate for Controls VS IBD types
g <- makePlotComparison(inDFfb,xVar = 'Diagnosis',yVar = 'F.B.Ratio',lbl = 'p.nice',returnData = T)
g[[1]] <- g[[1]] + ylab('Ln(Firmicutes / Bacteroides ratio)')
print(g)
ggsave(plot=g[[1]],filename = paste0(outF,'/Exploratory_FBratio_IBDtypes.png'),height = 6,width = 4)
write.table(g[[2]],paste0(outF,'/FExploratory_FBratio_IBDtypes.csv'),sep=',',row.names = F)

# =============================================================================
#   ================= DIVERSITY (ALPHA) ==============
# =============================================================================
# > NOTES:
#  - calculation of Alpha Diversity
#  - exploratory plots
# > folder for results
outF <- 'summaries_diversity'
unlink(outF,recursive = T)
if (!dir.exists(outF)) {dir.create(outF)}

# > digression calculation of diversity: 
#  >> subset microbiome to appropriate taxonomic level
#     (example here is Genus)
#   >>> get all "g__" columns == genera & lower
inDFm.g <- inDFm[,grep('g__',colnames(inDFm))]
#   >>> get rid of all "s__" columns == species & lower
inDFm.g <- inDFm.g[,grep('s__',colnames(inDFm.g),invert = T)]
#   >>> re-normalize (because we also excluded UNKNOWN)
inDFm.g <- inDFm.g / rowSums(inDFm.g)
#   >>> calculate diversity
inDFm.g.div.shannon <- diversity(inDFm.g,index = 'shannon')
#   check it:
head(inDFm.g.div.shannon)

# > Now let's automate this process!
inDFm2 <- inDFm
inDFm2$ID <- row.names(inDFm2)
divM <- calcDIVMetrics(inDFm2,
                       IDcol = 'ID',
                       DIVlvls=c("taxS","taxG"),
                       metrics = c("shannon","invsimpson","richness"))
inDFdiv <- merge(inDFm2,divM,by='ID')

# Exploratory plots: ALPHA diversity
# NOTES: 
#  - we do regression alongside with phenotypes later
# > make plots VS IBD subtype
for (dMet in c('DIV.S.shannon','DIV.S.invsimpson','DIV.S.richness',
               'DIV.G.shannon','DIV.G.invsimpson','DIV.G.richness')) {
  g <- makePlotComparison(inDFdiv,xVar = 'Diagnosis.2',yVar = dMet,returnData = F)
  print(g)
  ggsave(plot=g,filename = paste0(outF,'/Diversity_',dMet,'.png'),height = 6,width = 4)
}

# =============================================================================
#   =================== DIVERSITY (BETA) =======================
# =============================================================================
# NOTES:
#  - if using PCA, we use PCA (not PCoA) of CLR-transformed abundances
#    > PCA does euclidian distance by itself, so we should NEVER feed it distances!
# ALTERNATIVE:
#  - is PCoA of Aitchinson distances, results are very similar
# > folder for results
outF <- 'summaries_betadiversity'
unlink(outF,recursive = T)
if (!dir.exists(outF)) {dir.create(outF)}

# Step by step example for Genus-level beta-diversity Aitchinson PCoA plot
# > get Genera & Diagnosis!
#  >> subset microbiome to appropriate taxonomic level
#     (example here is Genus)
#   >>> get all "g__" columns == genera & lower
inDFm.t <- inDFm[,grep('s__',colnames(inDFm))]
#   >>> get rid of all "t__" columns == species & lower
inDFm.t <- inDFm.t[,grep('t__',colnames(inDFm.t),invert = T)]
#   >>>> CLR:
inDFmb.CLR <- decostand(inDFm.t,method = "clr",
                        pseudocount = min(inDFm.t[inDFm.t > 0])/2 )
#   >>>> Aitchinson distances (this is just euclidian distance mat. of CLR transformed values)
inDFmb.aitch <- vegdist(inDFmb.CLR,method='euclidian')
#   >>>> do PCoA transformation, keep first 10 components
pcoaMat <- cmdscale(inDFmb.aitch, eig = T, k=10)
#   >>>> calculate variance explained for first 2 PCs
variance <- head(eigenvals(pcoaMat)/sum(eigenvals(pcoaMat)))
var_pc1 <- as.integer(variance[1]*100)
var_pc2 <- as.integer(variance[2]*100)
#   >>>> now link it to phenotype(s)
pcoaDF <- as.data.frame(pcoaMat$points)
# annotate dataframe
colnames(pcoaDF)[1:2] <- c("PCo1","PCo2")
# merge with phenotypes for plotting
pcoaDF <- merge(pcoaDF,inDFm,by='row.names')
# and now plot it
g <- ggplot(pcoaDF,aes(x=PCo1,y=PCo2,col=Diagnosis.2)) + geom_point(size=3)
print(g)
# add ellipses to the plot
g <- g + geom_mark_ellipse(tol = 0.05,linetype='dashed',expand = unit(0, "mm"))
print(g)
# add centroids
ref <- reformulate("Diagnosis.2","cbind(PCo1,PCo2)")
centroids <- aggregate(ref,pcoaDF,mean)
g <- g + geom_point(data=centroids,shape=16,stroke=3,size=6,aes_string(x="PCo1",y="PCo2",col=ph),alpha=1) + 
  theme(legend.position="top")
print(g)

# Let's automate this process and make prettier plots!
# > variables to test
phenos <- c('Diagnosis','Diagnosis.2')
# > colors
colVec <- c('#5A98D6','#DEB045','#E05045')

# > make plot per phenotype
for (ph in phenos) {
  inDFm.noNA <- inDFm %>% drop_na(any_of(ph))
  inDFm.noNA$ID <- row.names(inDFm.noNA)
  # > prepare results
  inDFmb <- subsetMicrobiomeDF(inDFm.noNA,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,
                                getDivs = F,idToRowName = T,getID = T,idName = 'ID')
  inDFmeta <- subsetMicrobiomeDF(inDFm.noNA,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,
                                  getDivs = F,idToRowName = T,getID = T,idName = 'ID')
  inDFmb.spec <- filterMetaGenomeDF(inDFmb,presPerc = 0.05,minMRelAb = 1.0e-6,minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                     keepDomains = "All",keepLevels = c('S'),keepUnknown = F)
  # > CLR
  inDFmb.CLR <- decostand(inDFmb.spec,method = "clr",pseudocount = min(inDFmb.spec[inDFmb.spec > 0])/2 )

  # ======= SPECIES PCOA =============
  # ==================================
  # > calculate PCoA on Aitchinson distance [species]
  inDFmb.aitch <- vegdist(inDFmb.CLR,method='euclidian')
  pcoaMat <- cmdscale(inDFmb.aitch, eig = T, k=3)
  # calculate variance explained for first 3 PCs
  variance <- head(eigenvals(pcoaMat)/sum(eigenvals(pcoaMat)))
  var_pc1 <- as.integer(variance[1]*100)
  var_pc2 <- as.integer(variance[2]*100)
  var_pc3 <- as.integer(variance[3]*100)
  # convert to dataframe for happiness and plotting
  pcoaDF <- as.data.frame(pcoaMat$points)
  # annotate dataframe
  colnames(pcoaDF) <- c("PCo1","PCo2","PCo3")
  # merge with phenotypes for plotting
  pcoaDF <- merge(pcoaDF,inDFmeta,by='row.names')
  pcoaDFtoPlot <- pcoaDF
  # > plot
  print(paste0(' >> plotting results for ',ph))
  # plot PCo1 vs PCo2
  g <- ggplot(pcoaDFtoPlot,aes_string(x="PCo1",y="PCo2",col=ph)) + 
    geom_mark_ellipse(tol = 0.05,linetype='dashed',expand = unit(0, "mm")) + 
    geom_point(size=2,alpha=0.80) + 
    scale_color_manual(values=colVec) +
    theme_classic() + 
    theme(legend.position = 'top') +
    theme(text = element_text(size=16)) + 
    xlab(paste0("PCo1, ",var_pc1,"% variance")) +
    ylab(paste0("PCo2, ",var_pc2,"% variance")) 
  g <- g + ggtitle(paste0('Aitchison PCoA (Species) vs ',ph))
  # add centroids
  ref <- reformulate(ph,"cbind(PCo1,PCo2)")
  centroids <- aggregate(ref,pcoaDFtoPlot,mean)
  g <- g + geom_point(data=centroids,shape=16,stroke=3,size=6,aes_string(x="PCo1",y="PCo2",col=ph),alpha=1) + 
    theme(legend.position="top")
  # add density plots
  xdens <- 
    axis_canvas(g, axis = "x") + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo1', fill = ph, colour = ph), alpha = 0.5,size=0) + 
    scale_color_manual(values=colVec) + 
    scale_fill_manual(values=colVec)
  ydens <-
    axis_canvas(g, axis = "y", coord_flip = TRUE) + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo2', fill = ph, colour = ph), alpha = 0.5,size=0) +
    scale_color_manual(values=colVec) +
    scale_fill_manual(values=colVec) +
    coord_flip()
  # make prettier
  #g <- g + xlim(-60,70) + ylim(-30,60)
  # make even prettier
  g <- g %>% insert_xaxis_grob(xdens, grid::unit(1.25, "cm"), position = "top") %>% 
    insert_yaxis_grob(ydens, grid::unit(1.25, "cm"), position = "right")
  ggdraw(g)
  # save
  ggsave(plot=g,filename = paste0(outF,'/Aitchison_PCoA_spec_',ph,'_12.png'),height = 6,width = 6)
  
  # plot PCo2 vs PCo3
  g <- ggplot(pcoaDFtoPlot,aes_string(x="PCo2",y="PCo3",col=ph)) + 
    geom_mark_ellipse(tol = 0.05,linetype='dashed',expand = unit(0, "mm")) + 
    geom_point(size=3,alpha=0.65) + 
    scale_color_manual(values=colVec) +
    theme_classic() + 
    theme(legend.position = 'top') +
    theme(text = element_text(size=16)) + 
    xlab(paste0("PCo2, ",var_pc2,"% variance")) +
    ylab(paste0("PCo3, ",var_pc3,"% variance")) 
  g <- g + ggtitle(paste0('Aitchison PCoA (Species) vs ',ph))
  # add centroids
  ref <- reformulate(ph,"cbind(PCo2,PCo3)")
  centroids <- aggregate(ref,pcoaDFtoPlot,mean)
  g <- g + geom_point(data=centroids,shape=16,stroke=3,size=6,aes_string(x="PCo2",y="PCo3",col=ph),alpha=1) + 
    theme(legend.position="top")
  # add density plots
  xdens <- 
    axis_canvas(g, axis = "x") + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo2', fill = ph, colour = ph), alpha = 0.5,size=0) + 
    scale_color_manual(values=colVec) + 
    scale_fill_manual(values=colVec)
  ydens <-
    axis_canvas(g, axis = "y", coord_flip = TRUE) + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo3', fill = ph, colour = ph), alpha = 0.5,size=0) +
    scale_color_manual(values=colVec) +
    scale_fill_manual(values=colVec) +
    coord_flip()
  # make prettier
  g <- g %>% insert_xaxis_grob(xdens, grid::unit(1.25, "cm"), position = "top") %>% 
    insert_yaxis_grob(ydens, grid::unit(1.25, "cm"), position = "right")
  ggdraw(g)
  ggsave(plot=g,filename = paste0(outF,'/Aitchison_PCoA_spec_',ph,'_23.png'),height = 6,width = 6)
  # ===== PATHWAYS PCoA ======
  # ==================================
  inDFpwy <- subsetMicrobiomeDF(inDFm.noNA,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,
                               getDivs = F,idToRowName = T,getID = T,idName = 'ID')
  inDFmeta <- subsetMicrobiomeDF(inDFm.noNA,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,
                                 getDivs = F,idToRowName = T,getID = T,idName = 'ID')
  inDFpwy.CLR <- decostand(inDFpwy,method = "clr",
                           pseudocount = min(inDFpwy[inDFpwy > 0])/2 )
  # > calculate PCoA on Aitchinson distance [PWYs]
  inDFmb.aitch <- vegdist(inDFpwy.CLR,method='euclidian')
  pcoaMat <- cmdscale(inDFmb.aitch, eig = T, k=3)
  # calculate variance explained for first 3 PCs
  variance <- head(eigenvals(pcoaMat)/sum(eigenvals(pcoaMat)))
  var_pc1 <- as.integer(variance[1]*100)
  var_pc2 <- as.integer(variance[2]*100)
  var_pc3 <- as.integer(variance[3]*100)
  # convert to dataframe for happiness and plotting
  pcoaDF <- as.data.frame(pcoaMat$points)
  # annotate dataframe
  colnames(pcoaDF) <- c("PCo1","PCo2","PCo3")
  # merge with phenotypes for plotting
  pcoaDF <- merge(pcoaDF,inDFmeta,by='row.names')
  pcoaDFtoPlot <- pcoaDF
  # > plot
  print(paste0(' >> plotting PWY results for ',ph))
  # plot PCo1 vs PCo2
  g <- ggplot(pcoaDFtoPlot,aes_string(x="PCo1",y="PCo2",col=ph)) + 
    geom_mark_ellipse(tol = 0.05,linetype='dashed',expand = unit(0, "mm")) + 
    geom_point(size=2,alpha=0.80) + 
    scale_color_manual(values=colVec) +
    theme_classic() + 
    theme(legend.position = 'top') +
    theme(text = element_text(size=16)) + 
    xlab(paste0("PCo1, ",var_pc1,"% variance")) +
    ylab(paste0("PCo2, ",var_pc2,"% variance")) 
  g <- g + ggtitle(paste0('Aitchison PCoA (PWYs) vs ',ph))
  # add centroids
  ref <- reformulate(ph,"cbind(PCo1,PCo2)")
  centroids <- aggregate(ref,pcoaDFtoPlot,mean)
  g <- g + geom_point(data=centroids,shape=16,stroke=3,size=6,aes_string(x="PCo1",y="PCo2",col=ph),alpha=1) + 
    theme(legend.position="top")
  # add density plots
  xdens <- 
    axis_canvas(g, axis = "x") + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo1', fill = ph, colour = ph), alpha = 0.5,size=0) + 
    scale_color_manual(values=colVec) + 
    scale_fill_manual(values=colVec)
  ydens <-
    axis_canvas(g, axis = "y", coord_flip = TRUE) + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo2', fill = ph, colour = ph), alpha = 0.5,size=0) +
    scale_color_manual(values=colVec) +
    scale_fill_manual(values=colVec) +
    coord_flip()
  # make prettier
  #g <- g + xlim(-60,70) + ylim(-30,60)
  # make even prettier
  g <- g %>% insert_xaxis_grob(xdens, grid::unit(1.25, "cm"), position = "top") %>% 
    insert_yaxis_grob(ydens, grid::unit(1.25, "cm"), position = "right")
  ggdraw(g)
  # save
  ggsave(plot=g,filename = paste0(outF,'/Aitchison_PCoA_pwys_',ph,'_12.png'),height = 6,width = 6)
  
  # plot PCo2 vs PCo3
  g <- ggplot(pcoaDFtoPlot,aes_string(x="PCo2",y="PCo3",col=ph)) + 
    geom_mark_ellipse(tol = 0.05,linetype='dashed',expand = unit(0, "mm")) + 
    geom_point(size=3,alpha=0.65) + 
    scale_color_manual(values=colVec) +
    theme_classic() + 
    theme(legend.position = 'top') +
    theme(text = element_text(size=16)) + 
    xlab(paste0("PCo2, ",var_pc2,"% variance")) +
    ylab(paste0("PCo3, ",var_pc3,"% variance")) 
  g <- g + ggtitle(paste0('Aitchison PCoA (Species) vs ',ph))
  # add centroids
  ref <- reformulate(ph,"cbind(PCo2,PCo3)")
  centroids <- aggregate(ref,pcoaDFtoPlot,mean)
  g <- g + geom_point(data=centroids,shape=16,stroke=3,size=6,aes_string(x="PCo2",y="PCo3",col=ph),alpha=1) + 
    theme(legend.position="top")
  # add density plots
  xdens <- 
    axis_canvas(g, axis = "x") + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo2', fill = ph, colour = ph), alpha = 0.5,size=0) + 
    scale_color_manual(values=colVec) + 
    scale_fill_manual(values=colVec)
  ydens <-
    axis_canvas(g, axis = "y", coord_flip = TRUE) + 
    geom_density(data = pcoaDFtoPlot, 
                 aes_string(x = 'PCo3', fill = ph, colour = ph), alpha = 0.5,size=0) +
    scale_color_manual(values=colVec) +
    scale_fill_manual(values=colVec) +
    coord_flip()
  # make prettier
  g <- g %>% insert_xaxis_grob(xdens, grid::unit(1.25, "cm"), position = "top") %>% 
    insert_yaxis_grob(ydens, grid::unit(1.25, "cm"), position = "right")
  ggdraw(g)
  ggsave(plot=g,filename = paste0(outF,'/Aitchison_PCoA_pwys_',ph,'_23.png'),height = 6,width = 6)
}

# =============================================================================
#    =======   VARIANCE EXPLAINED BY PHENOTYPES (ADONIS) ======
# ===========================================================================
outF <- 'results_adonis'
if (!dir.exists(outF)) {
  dir.create(outF)
} 
#  >> prep data
inDFm$ID <- row.names(inDFm)
inDFhmb <- subsetMicrobiomeDF(inDFm,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,
                              getDivs = F,idToRowName = T,getID = T,idName = 'ID')
inDFhmeta <- subsetMicrobiomeDF(inDFm,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,
                                getDivs = F,idToRowName = T,getID = T,idName = 'ID')
inDFhmb.spec <- filterMetaGenomeDF(inDFhmb,presPerc = 0.05,minMRelAb = 1.0e-6,
                                   minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                   keepDomains = "All",keepLevels = c('S'),keepUnknown = F)
# >> prep phenotypes
phenosToTest <- colnames(inDFhmeta)
#  > covariates
covariates <- c('AGE','SEX','BMI','PFReads','PPI','Antibiotics')

# > RUN UNIVARIATE ADONIS (SPECIES, ALL SAMPLES)
aRes <- doAdonisTaxon(inPhenos = inDFhmeta,inTaxa = inDFhmb.spec,
                      adonisVarsTouse = phenosToTest, permNR = 1000)
write.table(aRes[[1]],paste0(outF,'/ADONIS_allsamples_species.csv'),
            sep=',',row.names = F)
g <- aRes[[2]] + theme(legend.position = "top") + theme(text = element_text(size = 21))
print(g)
ggsave(plot = g,filename = paste0(outF,'/ADONIS_allsamples_species.png'),width = 10, height = 12)

# > RUN UNIVARIATE ADONIS (PATHWAYS, ALL SAMPLES)
aRes <- doAdonisPWYs(inPhenos = inDFhmeta,inPWYs = inDFpwy,
                      adonisVarsTouse = phenosToTest, permNR = 1000)
write.table(aRes[[1]],paste0(outF,'/ADONIS_allsamples_PWYs.csv'),
            sep=',',row.names = F)
g <- aRes[[2]] + theme(legend.position = "top") + theme(text = element_text(size = 21))
print(g)
ggsave(plot = g,filename = paste0(outF,'/ADONIS_allsamples_PWYs.png'),width = 10, height = 12)

#TODO: add multivariable adonis!

# ================================================================================
# ================================================================================
# Differential abundance, diversity and pathways, Regression analysis [IBD vs HC]
# ================================================================================
# ================================================================================
# merge all data layers (diversity + taxa + pathways)
inDFm$ID <- row.names(inDFm)
# > prep Metadata
inDFmeta <- subsetMicrobiomeDF(inDFm,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,
                               getDivs = F,idToRowName = T,getID = T,idName = 'ID')
# > extract microbial taxa
inDFmb <- subsetMicrobiomeDF(inDFm,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,
                             getDivs = F,idToRowName = T,getID = T,idName = 'ID')
inDFmb.nonstrain <- filterMetaGenomeDF(inDFmb,presPerc = 0.1,minMRelAb = 1.0e-5,minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                       keepDomains = "All",keepLevels = c('P','C','O','F','G','S'),keepUnknown = F)
#   >> do CLR per taxonomy level
inDFmb.nonstrain.CLR <- NULL
for (taxlvl in c('P','C','O','F','G','S')) {
  print(paste0(' >> doing CLR for taxLvl = ',taxlvl))
  inDFmb.onelvl <- filterMetaGenomeDF(inDFmb.nonstrain,presPerc = 0.1,minMRelAb = 1.0e-5,minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                         keepDomains = "All",keepLevels = c(taxlvl),keepUnknown = F)
  inDFmb.onelvl.CLR <- decostand(inDFmb.onelvl,method = "clr",
                           pseudocount = min(inDFmb.onelvl[inDFmb.onelvl > 0])/2 )
  if (is.null(inDFmb.nonstrain.CLR)) {
    inDFmb.nonstrain.CLR <- inDFmb.onelvl.CLR
  } else {
    inDFmb.nonstrain.CLR <- merge(inDFmb.nonstrain.CLR,inDFmb.onelvl.CLR,by='row.names') %>% 
      column_to_rownames('Row.names')
  }
}

# > extract PWYs
inDFpwy <- filterHumannDF(inHum,presPerc = 0.1,minMRelAb = 1.0e-5,minMedRelAb = -1,rescale = T,
                             minSum = 1)
#   >> handle R-unfriendly pathway names
colnames(inDFpwy) <- make.names(colnames(inDFpwy))

#   >> do CLR
inDFpwy.CLR <- decostand(inDFpwy,method = "clr",
                        pseudocount = min(inDFpwy[inDFpwy > 0])/2 )
# diversities
inDFdiv <- calcDIVMetrics(inDFm,
                          DIVlvls=c("taxS","taxG"),
                          metrics = c("shannon","invsimpson","richness"))
# make master table
inDFmm <- inDFmeta %>% merge(inDFdiv,by='row.names') %>% column_to_rownames('Row.names') %>%
  merge(inDFmb.nonstrain,by='row.names') %>% column_to_rownames('Row.names') %>%
  merge(inDFpwy,by='row.names') %>% column_to_rownames('Row.names')

inDFmm.CLR <- inDFmeta %>% merge(inDFdiv,by='row.names') %>% column_to_rownames('Row.names') %>%
  merge(inDFmb.nonstrain.CLR,by='row.names') %>% column_to_rownames('Row.names') %>%
  merge(inDFpwy.CLR,by='row.names') %>% column_to_rownames('Row.names') 

# > save the data we used for regression analysis
saveRDS(inDFmm,'data_ready/data_rdyforregression.RDS')
saveRDS(inDFmm.CLR,'data_ready/data_rdyforregression.RDS')
write.table(inDFmm,'data_ready/data_rdyforregression.csv',sep=',')
write.table(inDFmm.CLR,'data_ready/CLR_data_rdyforregression.csv',sep=',')
write.table(inPwyNameTbl,'data_ready/pwy_names_conversion.csv',sep=',')
write.table(inDFdiv,'data_ready/data_diversities.csv',sep=',')

# > load the data for regression
# inDFmm <- readRDS('data_ready/LLD_1000IBD_data_rdyforregression.RDS')
# inDFmm.CLR <- readRDS('data_ready/LLD_1000IBD_CLR_data_rdyforregression.RDS')

# > prepare phenotypes to run regression on:
# =============================================
#  > phenos to test 
phenosToTest <- c('Diagnosis.2')
#  > covariates (to use)
mdlCovariates <- c('Age','Gender','BMI','PPI','Antibiotics','PFReads')
# > prepare features to test (and types of features)
varsToTest <- data.frame() %>% 
  bind_rows(data.frame(Feature=colnames(inDFdiv),Feature.Type='Diversity')) %>%
  bind_rows(data.frame(Feature=colnames(inDFmb.nonstrain.CLR),Feature.Type=getTaxLvlFromStringVec(colnames(inDFmb.nonstrain.CLR)))) %>%
  bind_rows(data.frame(Feature=colnames(inDFpwy),Feature.Type='Pathway'))

# output folder
outF <- "regression_models"
unlink(outF,recursive = T,force = T)
dir.create(outF)

# ==== RUN REGRESSION MODELS (IBD VS HC) ====
# ===========================================
# result accumulator
res <- NULL
#resEMMs <- NULL
for (ph in phenosToTest) {
  print('========================================================')
  print(paste0('>> phenotype = ',ph))
  print('========================================================')
  for (ftrN in c(1:nrow(varsToTest)) ) {
    t <- varsToTest$Feature[ftrN]
    # drop NAs and prep data for model building
    mdlData <- inDFmm.CLR[,c(ph,t,mdlCovariates)] %>% na.omit()
    # check for number of data points
    if (nrow(mdlData) <= 10) {
      print(paste0('WARNING, not enough data for model: ',t,' ~ ',ph))
    } else {
      # ============ FIT REGRESSION MODEL =======
      phClass <- class(mdlData[[ph]])
      frm <- reformulate(response=t,
                         termlabels = c(ph,mdlCovariates) )
      frmChar <- paste0(t,' ~ ',paste0(c(ph,mdlCovariates),collapse=' + ') )
      print(paste0(' > fitting model: ',frmChar))
      mdlFit <- lm(formula = frm,
                    data = mdlData)
      mdlFitS <- summary(mdlFit) %>% coef() %>% data.frame() %>% 
        rownames_to_column('Phenotype.Value') %>% filter(str_detect(Phenotype.Value,ph)) %>%
        mutate(Feature=t) %>% dplyr::rename(Std.Error = Std..Error) %>% dplyr::rename(P.value = Pr...t..) %>%
        mutate(Phenotype=ph,Model=frmChar,N=nrow(mdlData),Phenotype.class=phClass)
      res <- res %>% rbind(mdlFitS)
    }
  }
}
res <- res %>% mutate(FDR = p.adjust(P.value)) %>% arrange(P.value)
#resEMMs <- resEMMs %>% mutate(FDR = p.adjust(p)) %>% arrange(p) %>% mutate(p.adj = NULL)

write.table(res,paste0(outF,'/','results_regression.csv'),sep=',',row.names = F)
#write.table(resEMMs,paste0(outF,'/','1000IBD_IBD_vs_HC_results_regression_EMMs.csv'),sep=',',row.names = F)

# EXTRA: ESTIMATED MARGINAL MEANS (EMMs)
# ================================

if (class(mdlData[[ph]]) == 'factor' | class(mdlData[[ph]]) == 'character') {
  # calculate EMMs
  emmFrm <- reformulate(termlabels = ph)
  mdlEmm <- emmeans(mdlFit, emmFrm)
  #  > calculate contrasts
  ctr <- contrast(mdlEmm,"revpairwise",adjust='none')
  ctrsummary <- summary(ctr,type = "response", infer = c(TRUE, TRUE))
  #print(ctrsummary)
  #  > parse it
  emm_contrasts_p <- parseEMMContrasts(ctrsummary = ctrsummary,
                                       yname = ph,
                                       mtd = 'LM',
                                       p.adj = 'fdr')
  #emm_contrasts_sig <- emm_contrasts_p %>% filter(p < 0.05)
  emm_contrasts_p_tbl <- emm_contrasts_p %>%
    mutate(y.position = NULL, Feature=t, .y. = NULL, p.format = NULL,p.adj.f = NULL) %>%
    mutate(Phenotype = ph, Model = frmChar, N=nrow(mdlData),Phenotype.class=phClass) %>%
    relocate(Feature,Phenotype,group1,group2)
  resEMMs <- resEMMs %>% rbind(emm_contrasts_p_tbl)
}


# ==========================================================================
# ONE ALTERNATIVE: MAASLIN2
# ==========================================================================
# install: 
# > get Rtools
# > install BioConductor:
#install.packages("BiocManager")
# > Install Maaslin2
# BiocManager::install("Maaslin2")
# > load it
library(Maaslin2)

# prepare data:
# > taxonomy
inTaxa <- subsetMicrobiomeDF(inDFm,getTaxa = T,getPWYs = F,getVFs = F,
                             getCARDs = F,getPhenos = F,getDivs = F,idToRowName = F)
inTaxa <- filterMetaGenomeDF(inTaxa,presPerc = 0.1,minMRelAb = -1,
                             rescaleTaxa = T,keepLevels = c('G'))
# > metadata
inMeta <- subsetMicrobiomeDF(inDFm,getTaxa = F,getPWYs = F,getVFs = F,
                             getCARDs = F,getPhenos = T,getDivs = F,idToRowName = F)
# > run it
Maaslin2(input_data = inTaxa,
         input_metadata = inMeta,
         output = 'results_maaslin2',
         fixed_effects = c('Diagnosis'), #,'Age','Gender','BMI','PPI','Antibiotics','PFReads'),
         reference = c('Diagnosis,HC'),
         cores = 1,
         min_prevalence = 0.0,
         normalization = 'NONE',
         transform = 'NONE')
