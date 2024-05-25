# =========================================================
# by R.Gacesa (UMCG, 2019)
#
# MISC FUNCTIONS for plotting, microbiome stuff ...
# 
# NOTE:
# > also includes some functions  written by Dr. Alexander Kurilshikov (UMCG)
# =============================================================================
#library(pheatmap)

# annotates features types in a vector (usually microbiome results table column "Feature")
# returns data frame: feature, feature_type [e.g. Pathway, Taxon, Diversity], feature_type_narrow [e.g. Species]
# ========================================================================================
annotateMicrobiomeFeatures <- function(ftr) {
  ftrsWide <- c()
  ftrsNarrow <- c()
  ftrsShort <- c()
  for (f in ftr) {
    ftrTypeWide <- "Unknown"
    ftrTypeNarrow <- "Unknown"
    ftrShort <- NA
    if (grepl('t__',f)) {ftrTypeNarrow <- 'Strain'; ftrTypeWide <- 'Taxon'
    } else if (grepl('s__',f)) {ftrTypeNarrow <- 'Species'; ftrTypeWide <- 'Taxon'
    } else if (grepl('g__',f)) {ftrTypeNarrow <- 'Genus'; ftrTypeWide <- 'Taxon'
    } else if (grepl('f__',f)) {ftrTypeNarrow <- 'Family'; ftrTypeWide <- 'Taxon'
    } else if (grepl('o__',f)) {ftrTypeNarrow <- 'Order'; ftrTypeWide <- 'Taxon'
    } else if (grepl('c__',f)) {ftrTypeNarrow <- 'Class'; ftrTypeWide <- 'Taxon'
    } else if (grepl('p__',f)) {ftrTypeNarrow <- 'Phylum'; ftrTypeWide <- 'Taxon'
    } else if (grepl('k__',f)) {ftrTypeNarrow <- 'Kingdom'; ftrTypeWide <- 'Taxon'
    } else if (grepl('PWY',f,fixed = T)) {ftrTypeNarrow <- 'Pathway'; ftrTypeWide <- 'Pathway'
    } else if (grepl('DIV.G.',f,fixed = T)) {ftrTypeNarrow <- 'Genus'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.S.',f,fixed = T)) {ftrTypeNarrow <- 'Species'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.F.',f,fixed = T)) {ftrTypeNarrow <- 'Family'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.O.',f,fixed = T)) {ftrTypeNarrow <- 'Order'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.C.',f,fixed = T)) {ftrTypeNarrow <- 'Class'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.P.',f,fixed = T)) {ftrTypeNarrow <- 'Phylum'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.K.',f,fixed = T)) {ftrTypeNarrow <- 'Kingdom'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.T.',f,fixed = T)) {ftrTypeNarrow <- 'Strain'; ftrTypeWide <- 'Diversity'
    } else if (grepl('DIV.PWY.',f,fixed = T)) {ftrTypeNarrow <- 'Pathway'; ftrTypeWide <- 'Diversity'
    } 
    ftrsWide <- c(ftrsWide,ftrTypeWide)
    ftrsNarrow <- c(ftrsNarrow,ftrTypeNarrow)
    if (ftrTypeWide == 'Taxon') {ftrShort <- purgeMGNameOne(f)}
    ftrsShort <- c(ftrsShort,ftrShort)
  }
  res <- cbind(ftr,ftrsShort,ftrsWide,ftrsNarrow) %>% as.data.frame()
  #row.names(res) <- ftr
  colnames(res) <- c("Feature","Feature.Short","Feature.Type","Feature.Type.Specific")
  res
}

# get (lowest) Taxonomy level from metaphlan string
# ========================================
getTaxLvlFromString <- function(inStr) {
  ret <- NULL
  if (grepl('t__',inStr)) {ret <- 'Strain'} 
  else if (grepl('s__',inStr)) {ret <- 'Species'}
  else if (grepl('g__',inStr)) {ret <- 'Genus'}
  else if (grepl('f__',inStr)) {ret <- 'Family'}
  else if (grepl('o__',inStr)) {ret <- 'Order'}
  else if (grepl('c__',inStr)) {ret <- 'Class'}
  else if (grepl('p__',inStr)) {ret <- 'Phylum'}
  else if (grepl('k__',inStr)) {ret <- 'Kingdom'}
  return(ret)
}
getTaxLvlFromStringVec <- function(inStrV) {
  ret <- c()
  for (inStr in inStrV) {
    ret <- c(ret,getTaxLvlFromString(inStr))
  }
  return(ret)
}

# Helper functions
modelToFormula <- function(frmStr) {
  frmY <- gsub(' ','',unlist(strsplit(frmStr,'~',fixed = T))[[1]])
  frmX <- gsub(' ','',unlist(strsplit(frmStr,'~',fixed = T))[[2]])
  frmX <- strsplit(frmX,'+',fixed = T) %>% unlist() %>% unique()
  frm <- reformulate(response = frmY,termlabels = frmX)
  frm
}
modelToVector <- function(frmStr) {
  frmY <- gsub(' ','',unlist(strsplit(frmStr,'~',fixed = T))[[1]])
  frmX <- gsub(' ','',unlist(strsplit(frmStr,'~',fixed = T))[[2]])
  frmX <- strsplit(frmX,'+',fixed = T) %>% unlist() %>% unique()
  c(frmY,frmX)
}

## my happy plotter function
# ===========================================
makePlotComparison <- function(inDF,xVar='Diagnosis',yVar,yLog=F,autocolor=F,
                               lbl="p.signif",returnData=F,pointAlpha=0.2) {
  # >>> reads
  if (yLog) {
    inDF[[yVar]] <- log(inDF[[yVar]])
  }
  g <- ggplot(inDF,aes_string(x=xVar,y=yVar,col=xVar)) + 
    geom_jitter(width=0.2,height=0,alpha=pointAlpha) + 
    geom_violin(size=1.1,alpha=0.0) +
    geom_boxplot(outlier.alpha = 0,width=0.2,alpha=0.2,size=1.1) + 
    theme(legend.position = 'top') +
    theme_classic() + 
    theme(text = element_text(size=16)) + 
    theme(legend.position = 'none') 
  if (!autocolor) {
    g <- g + scale_color_manual(values=c('#5A98D6','#DEB045','#E05045')) 
  }
  if (yLog) {g <- g + ylab(paste0('ln(',yVar,')'))}
  
  frm <- reformulate(response=yVar,termlabels = xVar)
  tsts <- compare_means(frm,
                        data = inDF,
                        p.adjust.method = "fdr",
                        method = 'wilcox.test')
  tsts$p.nice <- paste0('p = ',tsts$p.format,' ',tsts$p.signif)
  tstsSig <- tsts[tsts$p.adj < 0.1,]
  my_comparisons <- getAllStatComparisons(tstsSig)
  g <- g + stat_compare_means(comparisons = my_comparisons,
                              label = lbl,
                              hide.ns = T)
  if (returnData) { list(g,tsts)
  } else {
    g
  }
}

# ============================================
## DO ADONIS ON PATHWAYS, UNIVARIATE
# ============================================
doAdonisPWYs <- function(inPhenos,inPWYs,adonisVarsTouse,permNR=1000,IDcol='row.names',
                          distmethod='aitchison') {
  adonisResults <- NULL
  spDF <- filterHumannDF(inPWYs,presPerc = -1,minMRelAb = -1)
  for (i in adonisVarsTouse) {
    print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
    print ('  >> collecting complete cases')
    inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"ID"),drop=F]
    allDF <- merge(x=inPhenosOneVarID,y=spDF,by=IDcol)
    rownames(allDF) <- allDF$Row.names
    allDF$Row.names <- NULL
    allDF$ID <- NULL
    allDF <- allDF[complete.cases(allDF),]
    av <- allDF[[i]]
    allDF[[i]] <- NULL
    print ('  >> calculating distance matrix')
    if (distmethod == 'aitchison') {
      inDmat <- vegdist(allDF,method = "aitchison",
                        pseudocount = min(allDF[allDF > 0])/2)
    } else {
      inDmat <- vegdist(allDF,method = "bray")
    }
    print ('  >> doing adonis')
    nrRows <- length(av)
    if (length(av) < 3 | length(unique(av)) < 2) {
      print(paste0(' >> WARNING: ',i,' has no useful data!!'))
    } else {
      # prepare model
      #frm <- reformulate(response = 'inDmat', termlabels = c(av,covariates) )
      #print(paste0(' NR NAs: ',sum(is.na(av))))
      ad <- adonis2(inDmat ~ av,permutations=permNR,parallel = 1)
      aov_table <- ad
      # accumulate results
      oneRow <- data.frame(Var=i,
                           NR_nonNA=nrRows,
                           DF=aov_table$Df[1],
                           SumsOfSqs=aov_table$SumOfSqs[1],
                           #MeanSqs=aov_table[1,3],
                           FModel=aov_table$F[1],
                           R2=aov_table$R2[1],
                           pval=aov_table$`Pr(>F)`[1],
                           #
                           FDR.BH=NA,
                           Significant=NA)
      print(oneRow)
      adonisResults <- rbind.data.frame(adonisResults,oneRow)
      print (paste0('--- ',i,' DONE! ---'))
    }
  }
  # calculate FDRs
  rownames(adonisResults) = adonisResults$Var
  adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
  adonisResults$Significance="No"
  adonisResults$Significance[adonisResults$pval<0.05]="Nominal"
  adonisResults$Significance[adonisResults$FDR.BH<0.05]="FDR"
  adonisResults$Significance <- factor(adonisResults$Significance,
                                       levels = c("FDR","Nominal","No"))
  # beautify for plotting and happiness
  adonisResults$Col='#FF8C00'
  adonisResults$Col[adonisResults$pval<0.05]="#3939F4"
  adonisResults$Col[adonisResults$FDR.BH<0.05]="#00B64F"
  adonisResults <- adonisResults[order(adonisResults$pval),]  
  colors <- distinct(adonisResults, Significance, Col)
  pal <- colors$Col
  names(pal) <- colors$Significance
  
  # make plot
  g <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
    ylab ("Explained variance (R^2) ") + xlab ("Factor")  + 
    theme(legend.position = 'bottom') +
    theme(text = element_text(size=16)) + 
    scale_fill_manual(values=pal)
  adonisResults$Col <- NULL
  # return results
  list(adonisResults,g)
}

# ============================================
## DO ADONIS ON ONE TAX LEVEL, UNIVARIATE
# ============================================
doAdonisTaxon <- function(inPhenos,inTaxa,adonisVarsTouse,permNR=1000,IDcol='row.names',
                          distmethod='aitchison',taxon='S') {
  adonisResults <- NULL
  spDF <- filterMetaGenomeDF(inTaxa,presPerc = -1,minMRelAb = -1,keepLevels = c(taxon))
  for (i in adonisVarsTouse) {
    print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
    print ('  >> collecting complete cases')
    inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"ID"),drop=F]
    allDF <- merge(x=inPhenosOneVarID,y=spDF,by=IDcol)
    rownames(allDF) <- allDF$Row.names
    allDF$Row.names <- NULL
    allDF$ID <- NULL
    allDF <- allDF[complete.cases(allDF),]
    av <- allDF[[i]]
    allDF[[i]] <- NULL
    print ('  >> calculating distance matrix')
    if (distmethod == 'aitchison') {
      inDmat <- vegdist(allDF,method = "aitchison",
                        pseudocount = min(allDF[allDF > 0])/2)
    } else {
      inDmat <- vegdist(allDF,method = "bray")
    }
    print ('  >> doing adonis')
    nrRows <- length(av)
    if (length(av) < 3 | length(unique(av)) < 2) {
      print(paste0(' >> WARNING: ',i,' has no useful data!!'))
    } else {
      # prepare model
      #frm <- reformulate(response = 'inDmat', termlabels = c(av,covariates) )
      #print(paste0(' NR NAs: ',sum(is.na(av))))
      ad <- adonis2(inDmat ~ av,permutations=permNR,parallel = 1)
      aov_table <- ad
      # accumulate results
      oneRow <- data.frame(Var=i,
                           NR_nonNA=nrRows,
                           DF=aov_table$Df[1],
                           SumsOfSqs=aov_table$SumOfSqs[1],
                           #MeanSqs=aov_table[1,3],
                           FModel=aov_table$F[1],
                           R2=aov_table$R2[1],
                           pval=aov_table$`Pr(>F)`[1],
                           #
                           FDR.BH=NA,
                           Significant=NA)
      print(oneRow)
      adonisResults <- rbind.data.frame(adonisResults,oneRow)
      print (paste0('--- ',i,' DONE! ---'))
    }
  }
  # calculate FDRs
  rownames(adonisResults) = adonisResults$Var
  adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
  adonisResults$Significance="No"
  adonisResults$Significance[adonisResults$pval<0.05]="Nominal"
  adonisResults$Significance[adonisResults$FDR.BH<0.05]="FDR"
  adonisResults$Significance <- factor(adonisResults$Significance,
                                      levels = c("FDR","Nominal","No"))
  # beautify for plotting and happiness
  adonisResults$Col='#FF8C00'
  adonisResults$Col[adonisResults$pval<0.05]="#3939F4"
  adonisResults$Col[adonisResults$FDR.BH<0.05]="#00B64F"
  adonisResults <- adonisResults[order(adonisResults$pval),]  
  colors <- distinct(adonisResults, Significance, Col)
  pal <- colors$Col
  names(pal) <- colors$Significance
  
  # make plot
  g <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
    ylab ("Explained variance (R^2) ") + xlab ("Factor")  + 
    theme(legend.position = 'bottom') +
    theme(text = element_text(size=16)) + 
    scale_fill_manual(values=pal)
  adonisResults$Col <- NULL
  # return results
  list(adonisResults,g)
}

# ============================================
## DO ADONIS ON ONE TAX LEVEL, MULTIVARIATE
# ============================================
doAdonisTaxonMultivariate <- function(inPhenos,inTaxa,adonisVarsTouse,permNR=1000,IDcol='row.names',
                                      distmethod='aitchison',taxon='S',covariates=c()) {
  adonisResults <- NULL
  spDF <- filterMetaGenomeDF(inTaxa,presPerc = -1,minMRelAb = -1,keepLevels = c(taxon))
  for (i in adonisVarsTouse) {
    print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
    if (i %in% covariates) {
      covariatesMdl <- covariates[covariates!=i] 
    } else {covariatesMdl <- covariates}
    print ('  >> collecting complete cases')
    inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"ID",covariates),drop=F]
    allDF <- merge(x=inPhenosOneVarID,y=spDF,by=IDcol)
    rownames(allDF) <- allDF$Row.names
    allDF$Row.names <- NULL
    allDF$ID <- NULL
    allDF <- allDF[complete.cases(allDF),]
    av <- allDF[,colnames(allDF) %in% c(i,"ID",covariatesMdl),drop=F]
    allDF <- allDF[,!(colnames(allDF) %in% c(i,"ID",covariatesMdl)),drop=F]
    print ('  >> calculating distance matrix')
    if (distmethod == 'aitchison') {
      inDmat <- vegdist(allDF,method = "aitchison",
                        pseudocount = min(allDF[allDF > 0])/2)
    } else {
      inDmat <- vegdist(allDF,method = "bray")
    }
    print ('  >> doing adonis')
    nrRows <- nrow(av)
    if (nrRows < 3 | length(unique(av[[i]])) < 2) {
      print(paste0(' >> WARNING: ',i,' has no useful data!!'))
    } else {
      # prepare model
      frm <- reformulate(response='inDmat',
                         termlabels = c(i,covariatesMdl) )
      frmChar <- paste0('DMAT ~ ',paste0(c(i,covariatesMdl),collapse=' + ') )
      #print(paste0(' NR NAs: ',sum(is.na(av))))
      ad <- adonis2(formula = frm,data = av,permutations=permNR,parallel = 4,by="margin")
      aov_table_row <- ad %>% data.frame() %>% rownames_to_column('VAR') %>% filter(VAR==i)
      # accumulate results
      oneRow <- data.frame(Var=i,
                           NR_nonNA=nrRows,
                           DF=aov_table_row$Df,
                           SumsOfSqs=aov_table_row$SumOfSqs,
                           FModel=aov_table_row$F,
                           R2=aov_table_row$R2,
                           pval=aov_table_row$Pr..F.,
                           Model=frmChar,
                           #
                           FDR.BH=NA,
                           Significant=NA)
      print(oneRow)
      adonisResults <- rbind.data.frame(adonisResults,oneRow)
      print (paste0('--- ',i,' DONE! ---'))
    }
  }
  # calculate FDRs
  rownames(adonisResults) = adonisResults$Var
  adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
  adonisResults$Significance="No"
  adonisResults$Significance[adonisResults$pval<0.05]="Nominal"
  adonisResults$Significance[adonisResults$FDR.BH<0.05]="FDR"
  adonisResults$Significance <- factor(adonisResults$Significance,
                                       levels = c("FDR","Nominal","No"))
  # beautify for plotting and happiness
  adonisResults$Col='#FF8C00'
  adonisResults$Col[adonisResults$pval<0.05]="#3939F4"
  adonisResults$Col[adonisResults$FDR.BH<0.05]="#00B64F"
  adonisResults <- adonisResults[order(adonisResults$pval),]  
  colors <- distinct(adonisResults, Significance, Col)
  pal <- colors$Col
  names(pal) <- colors$Significance
  
  # make plot
  g <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
    ylab ("Explained variance (R^2) ") + xlab ("Factor")  + 
    theme(legend.position = 'bottom') +
    theme(text = element_text(size=16)) + 
    scale_fill_manual(values=pal)
  adonisResults$Col <- NULL
  # return results
  list(adonisResults,g)
}

# MICROBIOME SUMMARY PLOT
# ============================================================================
# > inData = input dataframe (with taxa & phenotypes [if any])
# > toPlotLvl = which level to plot (P,O,C,F,G,S)
# > topX = how many taxa NOT to group (rest is grouped into OTHER)
# > groupBy = by which variable should data be grouped#
# ============================================================================
plotMicrobiomeSummary <- function(inData,toPlotLvl='P',groupBy=NULL,topX=9,minAb=1.0e-6,minPrev=0.00,getDataForPlot=F,
                                  prettierText=T) {
  inDataS <- filterMetaGenomeDF(inData,keepLevels = toPlotLvl,minMRelAb = minAb,presPerc = minPrev)
  inDataS <- purgeMGNamesRaw(inDataS)
  inDataS <- purgeMGNames(inDataS)
  inDataS$groupBy <- inDataS[[groupBy]]
  inDataS <- inDataS %>% select(c(contains('__'), groupBy))
  avgAb <- inDataS %>% select(contains('__')) %>% 
    apply(MARGIN = 2,FUN=function(x) {mean(x,na.rm = T)} ) %>% sort(decreasing = T)
  doGroup <- F
  if (topX < length(avgAb)) {
    toKeep <- avgAb[1:topX] %>% names()
    toGroup <- avgAb[(topX+1):length(avgAb)] %>% names()
    doGroup <- T
    toKeep <- c(toKeep,'OTHER')
  } else {
    toKeep <- avgAb %>% names()
    topX <- length(avgAb)
  }
  if (doGroup) {
    colOther <- inDataS %>% select(all_of(toGroup)) %>% rowSums()
    inDataS <- inDataS %>% mutate(OTHER=colOther) %>% select(!(all_of(toGroup)))
  }
  # summarize (get averages of taxa VS groupBy)
  if (!is.null(groupBy)) {
    inDataSS <- inDataS %>% group_by(groupBy) %>% 
      dplyr::summarise(across(c(toKeep), ~ mean(.x,na.rm = T) ))
  } else {
    inDataSS <- inDataS %>% dplyr::summarise(across(c(toKeep), ~mean(.x,na.rm = T) ))
  }
  # summarize (get total averages for sorting of factors)
  inDataSSa <- inDataS %>% dplyr::summarise(across(c(toKeep),~mean(.x,na.rm = T) )) %>% t() %>%
    as.data.frame() %>% rownames_to_column('Taxon') %>% arrange(-V1)
  
  # transform to long
  inDataSSL <- inDataSS %>% gather('Taxon','Abundance',c(toKeep))
  inDataSSL$Taxon <- fct_relevel(inDataSSL$Taxon,inDataSSa$Taxon)
  if (prettierText) {
    inDataSSL$Taxon2 <- as.character(inDataSSL$Taxon)
    inDataSSL$Taxon2 <- gsub('[spgfocpk]__','',inDataSSL$Taxon2)
    inDataSSL$Taxon2 <- gsub('_',' ',inDataSSL$Taxon2)
    inDataSSL$Taxon2 <- gsub('OTHER','Other',inDataSSL$Taxon2)
    
    inDataSSa$Taxon2 <- as.character(inDataSSa$Taxon)
    inDataSSa$Taxon2 <- gsub('[spgfocpk]__','',inDataSSa$Taxon2)
    inDataSSa$Taxon2 <- gsub('_',' ',inDataSSa$Taxon2)
    inDataSSa$Taxon2 <- gsub('OTHER','Other',inDataSSa$Taxon2)
    
    inDataSSL$Taxon <- as.factor(inDataSSL$Taxon2)
    inDataSSL$Taxon <- fct_relevel(inDataSSL$Taxon,inDataSSa$Taxon2)
  }
  
  if (!is.null(groupBy)) {
    # g <- ggplot(inDataSSL,aes_string(x=groupBy,y='Abundance',col='Taxon',fill='Taxon')) + geom_col() + 
    #   ylab('Relative Abundance') + theme_classic()
    g <- ggplot(inDataSSL,aes_string(x="groupBy",y='Abundance',fill='Taxon')) + geom_col(col='#222222') + 
      ylab('Relative Abundance') + theme_classic()
  } else {
    inDataSSL$Samples <- "All"
    g <- ggplot(inDataSSL,aes_string(x="Samples",y='Abundance',fill='Taxon')) + geom_col(col='#222222') + 
      ylab('Relative Abundance') + theme_classic() 
    if (!is(null(groupBy))) {
      g <- g + xlab(groupBy)
    }
    # g <- ggplot(inDataSSL,aes_string(x="Samples",y='Abundance',col='Taxon',fill='Taxon')) + geom_col() + 
    #   ylab('Relative Abundance') + theme_classic()
  }
  # prettify text a bit more
  if (prettierText) {
    g <- g + theme(legend.text = element_text(face = "italic"))
  }
  if (getDataForPlot) {  
    list(g,inDataSSL)
  } else {
    g
  }
}

# Parser for EMM contrasts summary
# ===========================================================================
# ===========================================================================
#  > turns the EMM summary table into data frame compatible with 
#   ggpubr function stat_pvalue_manual()

# Example: Toy code that demonstrates parser on cars dataset
# ==============================
# # > load data & libraries
# library(emmeans)
# library(ggpubr)
# library(tidyverse)
# data(cars)
# # > transforms cylinders into factor, turn values into string to avoid
# #   issues with EMM package value <-> name conversion#
# cars <- cars %>% mutate(Cyl.f = as.factor(paste0('Cyl.',Cylinder)))
# # > fit linear model (Price VS Cylinder category, adjust for Milage & NR. of doors)
# mdlFitLM <- lm(Price ~ Cyl.f + Mileage + Doors, cars)
# # > get EMMs
# mdlLmEmms <- emmeans(mdlFitLM, ~ Cyl.f)
# # > set up all pairs contrast, include FDR correction for multiple tests
# ctr <- contrast(mdlLmEmms,"revpairwise",adjust='fdr')
# # > summarize all pair-wise cylinder category tests
# ctrSummary <- summary(ctr,type = "response", infer = c(TRUE, TRUE))
# # > parse the data
# #  - Notes: yname is variable name we want to plot (Cyl.f in this case)
# #           ypos where brackets of p-values are plotted 
# ctr.parsed <- parseEMMContrasts(ctrSummary,
#                                 yname='Cyl.f',
#                                 ypos = max(as.data.frame(mdlLmEmms)$upper.CL) )
# # > plot EMMs with result significance
# #  - Notes: Cyl.f in X, marginal means on Y
# ggplot(as.data.frame(mdlLmEmms),
#        aes(x=Cyl.f,y=emmean,col=Cyl.f)) + geom_point() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   geom_linerange(aes(ymin = lower.CL, ymax = upper.CL),size=3,alpha=0.2) +
#   geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),width=0.2) +
#   ylab(paste0('EMM of Car price')) +
#   stat_pvalue_manual(ctr.parsed,
#                      label = "p.signif", step.increase = 0.035, tip.length = 0.0015)
# 
# # > plot EMMs with FDR correction:
# ggplot(as.data.frame(mdlLmEmms),
#        aes(x=Cyl.f,y=emmean,col=Cyl.f)) + geom_point() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   geom_linerange(aes(ymin = lower.CL, ymax = upper.CL),size=3,alpha=0.2) +
#   geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),width=0.2) +
#   ylab(paste0('EMM of Car price')) +
#   stat_pvalue_manual(ctr.parsed,
#                      label = "p.adj.f", step.increase = 0.055, tip.length = 0.005)
# ==========================================================================
parseEMMContrasts <- function(ctrsummary,yname,mtd='LM',p.adj='fdr',ypos=NA) {
  # .y.    group1   group2          p        p.adj p.format p.signif method   y.position
  # toTest [0, 10)  [60, 65) 3.25e- 4 0.0011       0.00032  ***      Wilcoxon         10
  ctrt <- as_tibble(ctrsummary)
  ctrt <- ctrt %>% 
    mutate(.y. = yname) %>% 
    separate_wider_delim(contrast," - ",names = c('group1','group2')) %>% 
    mutate(method = 'LM') %>%
    mutate(p.adj = p.adjust(p.value,method = p.adj)) %>%
    dplyr::rename(p = p.value) %>%
    mutate(p.signif = case_when(p.adj > 0.05 ~ '',
                                p.adj < 0.00005 ~ '****',
                                p.adj < 0.0005 ~ '***',
                                p.adj < 0.005 ~ '**',
                                p.adj <= 0.05 ~ '*')) %>% 
    mutate(p.adj.f = formatC(p.adj, format = "e", digits = 2))
  if (is.na(ypos)) {
    ctrt <- ctrt %>% 
      mutate(y.position = max(ctrsummary$estimate) ) 
  } else {
    ctrt <- ctrt %>% 
      mutate(y.position = ypos) }
  ctrt
}

# ==============================================================================
# pairwise proportion test with statistics
# ==============================================================================
pairwise_prop_test_wstats <- function (xtab, p.adjust.method = "fdr") 
{
  if (ncol(xtab) > 2 & nrow(xtab) == 2) {
    xtab <- t(xtab)
  }
  xtabdf <- as.data.frame(xtab)
  results <- stats::pairwise.prop.test(xtab, p.adjust.method = "none") %>%  
    tidy() %>% select(.data$group2, .data$group1, 
                      .data$p.value)
  colnames(results) <- c("group1", "group2", "p")
  results <- results %>% adjust_pvalue(method = p.adjust.method) %>% 
    add_significance("p.adj") %>% mutate(p = signif(.data$p, digits = 3), 
                                         p.adj = signif(.data$p.adj, digits = 3))
  N1 <- unique(xtabdf$Var2)[[1]]
  N2 <- unique(xtabdf$Var2)[[2]]
  results$V1 <- N1
  results$V2 <- N2
  results$group1V1 <- NA
  results$group1V2 <- NA
  results$group2V1 <- NA
  results$group2V2 <- NA
  for (g1 in results$group1) {
    for (g2 in results$group2) {
      results$group1V1[results$group1 == g1 & 
                         results$group2==g2] <- xtabdf$Freq[xtabdf$Var1==g1 & xtabdf$Var2==N1]
      results$group1V2[results$group1 == g1 & 
                         results$group2==g2] <- xtabdf$Freq[xtabdf$Var1==g1 & xtabdf$Var2==N2]
      results$group2V1[results$group1 == g1 & 
                         results$group2==g2] <- xtabdf$Freq[xtabdf$Var1==g2 & xtabdf$Var2==N1]
      results$group2V2[results$group1 == g1 & 
                         results$group2==g2] <- xtabdf$Freq[xtabdf$Var1==g2 & xtabdf$Var2==N2]
    }
  }
  
  results$group1V1.p <- results$group1V1/(results$group1V1+results$group1V2)
  results$group1V2.p <- results$group1V2/(results$group1V1+results$group1V2)
  results$group2V1.p <- results$group2V1/(results$group2V1+results$group2V2)
  results$group2V2.p <- results$group2V2/(results$group2V1+results$group2V2)
  results
}



# inverse rank transformation function by Alex
# ==============================================================================
invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
# ==============================================================================

# implementation of inverse rank transformation function by Alex to dataframe
# NOTE: it transforms only numeric and integer fields so it can eat mixed DF
# ==============================================================================
invrankDF <- function(inDF) {
  for (cn in colnames(inDF)) {
    if (class(inDF[[cn]]) == "numeric" | class(inDF[[cn]]) == "integer") {
      inDF[[cn]] <- invrank(inDF[[cn]] )
    }
  }
  inDF
}
# ==============================================================================


# ==================================================================
# KEGG BRITE MAPPER FOR HUMANN KO DATA
# INPUT
# - input: inBritePath = path to BRITE table
# - inKegg: cleaned KEGG humann2 output (obtained by running cleanHumannKegg function)
# -
# OUTPUT
# - list with:
#  [[1]] : BRITE lvl 3 grouped data
#  [[2]] : BRITE lvl 3 data summary
#  [[3]] : BRITE lvl 4 grouped data
#  [[4]] : BRITE lvl 4 data summary
#
# NOTES:
# - if verbose = T throws warnings if kos are missing (some missing ones are normal)

# ==================================================================
groupKEGG <- function(inBritePath, inKegg, verbose=T) {
  if (verbose) {
    print ('>>> STARTED KEGG-GROUPING FUNCTION <<<')
  }
  # load BRITE
  if (verbose) {
    print (paste0(' >> loading BRITE table [',inBritePath,']'))
  }
  inBrite <- read.table(inBritePath,sep=',',header=T,quote = '"',
                        stringsAsFactors = F, colClasses = "character")
  # load input data
  if (verbose) {
    print (paste0(' >> loading KEGG PATHWAYS table [',inKegg,']'))
  }
  inDF <- read.table(inKegg,sep='\t',header = T)
  
  # output collector
  outList = list()
  
  # merge keggs to lvl C (pathways) == level 3
  # ========================================
  # > build new DF, add C-levels one by one
  inDF$DUMMY <- NA
  rownames(inDF) <- inDF$ID
  inDFlvlC <- inDF[c("ID","DUMMY")]
  inDF$DUMMY <- NULL
  nrLvl3 <- 0
  for (lvlC in unique(inBrite$lvlC.ID)) {
    # get all KOs with this levelC ID
    kosToGet <- inBrite$lvlD.ID[inBrite$lvlC.ID==lvlC]
    # grep by it
    dfColsKO <- colnames(inDF)[colnames(inDF) %in% kosToGet]
    if (length(dfColsKO) == 0) {
      if (verbose) {
        print(paste0('   > WARNING: no KOs for lvl 3 [ko',lvlC,'] in data!'))
      }
    } else if (length(dfColsKO) == 1) {
      inDFlvlC[[paste0('ko',lvlC)]] <- inDF[[dfColsKO]]
      nrLvl3 <- nrLvl3 + 1
    } else {
      newCol <- rowSums(inDF[,dfColsKO])
      inDFlvlC[[paste0('ko',lvlC)]] <- newCol
      nrLvl3 <- nrLvl3 + 1
    }
  }
  if (verbose) {
    print (paste0(' >> merging of lvl 3s done, got total of ',nrLvl3,' pathways (at BRITE LVL 3)'))
  }
  # clean up gunk
  inDFlvlC$DUMMY <- NULL
  rownames(inDFlvlC) <- inDFlvlC$ID
  # re-normalize to relative values
  inDFlvlC$ID <- NULL
  inDFlvlCn <- inDFlvlC/rowSums(inDFlvlC)
  inDFlvlCn$ID <- rownames(inDFlvlCn)
  outList[[1]] = inDFlvlCn
  #write.table(inDFlvlCn,'kidney_kegg_filtered_lvl3.csv',sep=',',row.names = F)
  inDFlvlCn$ID <- NULL
  # do summaries
  summaries <- NULL
  for (cn in colnames(inDFlvlCn)) {
    oneRow <- makeOneSummary(inDFlvlCn[[cn]])
    oneRow$Var <- cn
    summaries <- rbind.data.frame(summaries,oneRow)
  }
  inBrite2 <- inBrite
  inBrite2$lvlC.ID <- paste0('ko',inBrite2$lvlC.ID)
  inBrite2 <- inBrite2[!duplicated(inBrite2$lvlC.ID),]
  summariesd <- merge(summaries, inBrite2, by.x="Var",by.y="lvlC.ID")
  summariesd$VarRaw <- gsub('ko','',summariesd$Var)
  summariesd <- summariesd[order(summariesd$Mean,decreasing = T),]
  
  #write.table(summariesd,'kidney_kegg_filtered_lvl3_summary.csv',sep=',',row.names = F)
  # make human readable summary
  summariesd$lvlB.ID <- NULL; summariesd$lvlA.ID <- NULL; summariesd$VarRaw <- NULL
  summariesd$lvlD.ID <- NULL; summariesd$lvlD.Name <- NULL; summariesd$dataType <- NULL
  summariesd$NrNAs <- NULL; summariesd$propNAs <- NULL
  colnames(summariesd)[1] <- "KEGG.ID"
  colnames(summariesd)[2] <- "SAMPLE.NR"
  colnames(summariesd)[12] <- "BRITE.LVLC"
  colnames(summariesd)[13] <- "BRITE.LVLB"
  colnames(summariesd)[14] <- "BRITE.LVLA"
  outList[[2]] = summariesd
  #write.table(summariesd,'kidney_kegg_filtered_lvl3_summary_hr.csv',sep=',',row.names = F)
  
  # =======> lvl 4 data
  inDF4 <- inDF
  inDF4$ID <- NULL
  summariesB <- NULL
  for (cn in colnames(inDF4)) {
    oneRow <- makeOneSummary(inDF4[[cn]])
    oneRow$Var <- cn
    summariesB <- rbind.data.frame(summariesB,oneRow)
  }
  inBrite1 <- inBrite#[,c("lvlD.ID","lvlD.Name")]
  inBrite1 <- inBrite1[!duplicated(inBrite1$lvlD.ID),]
  summariesBd <- merge(summariesB, inBrite1, by.x="Var",by.y="lvlD.ID")
  summariesBd <- summariesBd[order(summariesBd$Mean,decreasing = T),]
  outList[[3]] <- inDF4
  if (verbose) {
    print (paste0(' >> processing of lvl 4s done, got total of ',nrow(summariesBd),' kegg orthologs (at BRITE LVL 4)'))
  }
  
  #write.table(summariesBd,'dag3_kegg_filtered_lvl4_summary.csv',sep=',',row.names = F)
  # make human readable summary
  summariesBd$lvlB.ID <- NULL; summariesBd$lvlC.ID <- NULL; summariesBd$lvlA.ID <- NULL
  summariesBd$VarRaw <- NULL; summariesBd$lvlD.ID <- NULL; summariesBd$dataType <- NULL
  summariesBd$NrNAs <- NULL; summariesBd$propNAs <- NULL
  colnames(summariesBd)[1] <- "KEGG.ID"
  colnames(summariesBd)[2] <- "SAMPLE.NR"
  colnames(summariesBd)[12] <- "BRITE.LVLD"
  colnames(summariesBd)[13] <- "BRITE.LVLC"
  colnames(summariesBd)[14] <- "BRITE.LVLB"
  colnames(summariesBd)[15] <- "BRITE.LVLA"
  outList[[4]] = summariesBd
  #write.table(summariesBd,'dag3_kegg_filtered_lvl4_summary_hr.csv',sep=',',row.names = F)
  print ('>>> KEGG-GROUPING SUCCESSFUL, RETURNING RESULTS <<<')
  # return results
  outList
}


# recodes dataframe factor variables
# =========================================
recodeMulticlassFactorsDF <- function(inDF,minClasses=3,maxClasses=6,verbose=T,recodeInt=T,
                                      posClass='Y',negClass='N',idCol="") {
  
  if (verbose) {
    print(paste0('==================================================================================='))
    print(paste0('STARTING recodeMulticlassFactorsDF; minClasses = ',minClasses,', maxClasses = ',maxClasses))
    print(paste0(' recodeInt=',recodeInt))
    print(paste0('==================================================================================='))
  }
  # iterate over DF columns
  for (cn in colnames(inDF)) {
    tmpCol <- inDF[[cn]]
    # check if col is int, conditionally convert to factor
    if (class(tmpCol) == "integer" & recodeInt) {
      tmpCol <- factor(tmpCol)
    }
    # check if we want to recode it
    if (class(tmpCol) == "factor") {
      cclasses <- levels(tmpCol)
      if ( length(cclasses) >= minClasses & length(cclasses) <= maxClasses ) {
        if (verbose) {
          print(paste0(' > Column ',cn,' has ',length(cclasses),' classes, recoding!'))
        }
        newCols <- paste0(cn,'.',cclasses,'.',posClass,negClass)
        for (newColNR in c(1:length(cclasses)) ) {
          newCol <- inDF[[cn]]==cclasses[[newColNR]]
          newCol[newCol==TRUE] <- posClass
          newCol[newCol==FALSE] <- negClass
          newColName <- newCols[newColNR]
          inDF[[newColName]] <- newCol
        }
      }
    }
  } # end of column iteration loop
  # sort
  cnms <- colnames(inDF)
  cnmsWoIDs <- cnms[!cnms %in% c(idCol)]
  cnmsWoIDs <- cnmsWoIDs[order(cnmsWoIDs)]
  inDF <- inDF[,c("DAG3_sampleID",cnmsWoIDs) ]
  colnames(inDF) <- gsub('\\.\\.','\\.',colnames(inDF))
  inDF
}


# summary of one dataframe feature VS response
#TODO: extend to categorical variables
# ========================================================
makeOneSummaryByResponse <- function(inDF,ftr,response,doTestVsControl=T,controlClass="HC",
                                     formatSci=F,formatDig=2,formatRoundDig=3,dropNAs = T,verbose=T,
                                     includeTotals=T) {
  if (verbose) {
    print(paste0('==================================================================================='))
    print(paste0('STARTING makeOneSummaryByResponse; feature = ',ftr,', Response = ',response))
    print(paste0('==================================================================================='))
  }
  if (!response %in% colnames(inDF)) {
    stop(paste0('ERROR: ',response,' not in inDF columns! STOPPING!'))
  }
  if (!ftr %in% colnames(inDF)) {
    stop(paste0('ERROR: ',ftr,' not in inDF columns! STOPPING!'))
  }
  inDF <- inDF[,c(ftr,response)]
  if (dropNAs) {
    inDF <- inDF[complete.cases(inDF),]
  } else if (sum(is.na(inDF))) {
    print('WARNING: NAs in data and dropNAs is FALSE, CODE MIGHT CRASH!!!')
  }
  if (doTestVsControl) {
    if (!controlClass %in% inDF[[response]] & !controlClass == "OTHERS") {
      stop(paste0('ERROR: ',controlClass,' not in ',response, ' STOPPING! Use control = existing class or "OTHERS" to test VS rest of data'))
    }
  }
  ret = data.frame()
  retB = data.frame()
  retMG = data.frame()
  ftrsToTest <- unique(as.character(inDF[[response]]))
  if (includeTotals) {ftrsToTest <- c(ftrsToTest,"TOTAL")}
  for (rv in ftrsToTest) {
    if (verbose) {
      print(paste0(' >> making summary for ',ftr,' == ',rv,' VS ',response))
    }
    # INTEGER VARIABLE => transform to NUMERIC and treat as it
    if (class(inDF[[ftr]]) == "integer") {
      inDF[[ftr]] <- as.numeric(inDF[[ftr]])
    }
    # if (class(inDF[[ftr]]) == "character") {
    #   inDF[[ftr]] <- as.factor(inDF[[ftr]])
    # }
    # 
    # ==================================================
    # NUMERIC VARIABLE
    #  > use M/U/W test (wilcoxon implementation in R)
    #  > report mean/SD/min/Q1/Median/Q3/Max/NR/nonzero/prevalence (non na)
    # ==================================================
    if (class(inDF[[ftr]]) == "numeric") {
      if (verbose) {
        print(paste0('   -> ',ftr,' is numeric, making mean +/- sd and IQR range summary, M/U/W test for differences to controls!'))
      }
      if (rv == "TOTAL") {
        ftrData <- inDF[[ftr]]
      } else {
        ftrData <- inDF[inDF[[response]]==rv,][[ftr]]
      }
      mn <-  mean(ftrData)
      sd <-  sd(ftrData)
      md <-  median(ftrData)
      min <- min(ftrData)
      q1  <- quantile(ftrData)[[2]]
      q3  <- quantile(ftrData)[[4]]
      max <- max(ftrData)
      nr <- length(ftrData)
      nonzero <- sum(ftrData!=0)
      prev <- nonzero/nr
      if (is.na(prev)) {prev = 0.0}
      oneRow <- data.frame(Var=ftr,Var.Type=class(inDF[[ftr]]),Group=rv,Mean=mn,SD=sd,Min=min,Q1=q1,Median=md,Q3=q3,Max=max,
                           NR=nr,nonzero=nonzero,prevalence=prev,
                           niceOut=paste0(format(round(md,formatRoundDig),digits = formatDig,scientific = formatSci),
                                          ' [',format(round(q1,formatRoundDig),digits = formatDig,scientific = formatSci), '; ',
                                          format(round(q3,formatRoundDig),digits = formatDig,scientific = formatSci) ,']'),
                           niceOutMSD=paste0(format(round(mn,formatRoundDig),digits = formatDig,scientific = formatSci),
                                             ' +/- ',format(round(sd,formatRoundDig),digits = formatDig,scientific = formatSci))
      )
      print(doTestVsControl)
      if (!doTestVsControl) {
        print(paste0(' > ',response,' = ',rv,'; mean(',ftr,') = ',mn,'; sd = ',sd,'| median = ',md,' Q1 = ',q1,' Q3 = ',q3))
        ret <- rbind.data.frame(ret,oneRow)
      } else { 
        if (controlClass == "OTHERS") {
          if (rv != "TOTAL") {
            ctrData <- inDF[inDF[[response]] != rv,][[ftr]]
            tst <- wilcox.test(ctrData,ftrData)
            if (verbose) {
              print(paste0(" TESTING variable ",ftr,", Response variable(",response,") class " ,rv,' VS all other classes'))
            }
            sigStr = ""
            pV <- tst$p.value
            if (verbose) {
              print(paste0(' > P-value = ',pV))
            }
          }
        } else if (rv != controlClass) {
          ctrData <- inDF[inDF[[response]] == controlClass,][[ftr]]
          tst <- wilcox.test(ctrData,ftrData)
          if (verbose) {
            print(paste0(" TESTING variable ",ftr,", Response variable(",response,") class " ,rv,' VS ',controlClass))
          }
          sigStr = ""
          pV <- tst$p.value
          if (verbose) {
            print(paste0(' > P-value = ',pV))
          }
        } else {
          pV = NA; sig = ""; FDR = NA
        } 
        oneRow$pV <- pV
        ret <- rbind.data.frame(ret,oneRow)
      }
    } 
    # ===============================================
    # CATEGORICAL VARIABLE <BINARY>
    # > use chi/squared test
    # report: 
    #  > percent of group 1
    # TODO: add "others"
    # ===============================================
    if (class(inDF[[ftr]]) == "factor") {
      # BINARY VARIABLE
      if (length(unique(inDF[[ftr]])) == 2) {
        if (verbose) {
          print(paste0('   -> ',ftr,' is factor with 2 classes!!'))
        }
        g1 <- unique(inDF[[ftr]])[1]
        g2 <- unique(inDF[[ftr]])[2]
        if (rv == "TOTAL") {
          nr <- length(inDF[[response]])
          nonzero <- sum(inDF[[ftr]]!=0)
          nrg1 <- sum(inDF[[ftr]] == g1)
          nrg2 <- sum(inDF[[ftr]] == g2)
        } else {
          nr <- sum(inDF[[response]]==rv)
          nonzero <- sum(inDF[inDF[[response]]==rv,][[ftr]]!=0)
          nrg1 <- sum(inDF[[response]]==rv & inDF[[ftr]] == g1)
          nrg2 <- sum(inDF[[response]]==rv & inDF[[ftr]] == g2)
        }
        percg1 <- nrg1/nr
        percg2 <- nrg2/nr
        prev <- nonzero/nr
        if (is.na(prev)) {prev = 0.0}
        oneRowB <- data.frame(Var=ftr,Var.Type=class(inDF[[ftr]]),Group=rv,
                              NR=nr,nonzero=nonzero,prevalence=prev,
                              Lvl1=g1,nrLvl1=nrg1,percLvl1=percg1,
                              Lvl2=g2,nrLvl2=nrg2,percLvl2=percg2)
        if (!doTestVsControl) {
          retB <- rbind.data.frame(retB,oneRowB)
        } else if (controlClass == "OTHERS") {
          if (rv != "TOTAL") {
            toTestDF <- inDF
            toTestDF[[response]] <- as.character(toTestDF[[response]])
            toTestDF[[response]][toTestDF[[response]] != rv] <- "OTHERS"
            toTestDF[[response]] <- as.factor(toTestDF[[response]])
            
            toTestDF[[ftr]] <- as.factor(as.character(toTestDF[[ftr]] ))
            toTestDF[[response]] <- as.factor(as.character(toTestDF[[response]] ))
            tblTest <- table(toTestDF[[ftr]],toTestDF[[response]])
            tst <- chisq.test(tblTest)
            if (verbose) {
              print(tblTest)
              print(tst)
            }
            if (verbose) {
              print(paste0(" TESTING variable ",ftr,", Response variable(",response,") class " ,rv,' VS all other classes'))
            }
            sigStr = ""
            pV <- tst$p.value
            if (verbose) {
              print(paste0(' > P-value = ',pV))
            }
          }
        } else if (rv != controlClass & rv != "TOTAL") {
          toTestDF <- inDF[inDF[[response]] %in% c(rv,controlClass),]          
          toTestDF[[ftr]] <- as.factor(as.character(toTestDF[[ftr]] ))
          toTestDF[[response]] <- as.factor(as.character(toTestDF[[response]] ))
          tblTest <- table(toTestDF[[ftr]],toTestDF[[response]])
          tst <- chisq.test(tblTest)
          if (verbose) {
            print(tblTest)
            print(tst)
          }
          if (verbose) {
            print(paste0(" TESTING variable ",ftr,", Response variable(",response,") class " ,rv,' VS ',controlClass))
          }
          sigStr = ""
          pV <- tst$p.value
          if (verbose) {
            print(paste0(' > P-value = ',pV))
          }
        } else {
          pV = NA; sig = ""; FDR = NA
        } 
        oneRowB$pV <- pV
        retB <- rbind.data.frame(retB,oneRowB)
        ret <- retB
      }
    } else if (TRUE) {
      # ========================================
      # CATEGORICAL, > 2 groups
      # ========================================
    }
  }
  ret
}

# summary of one dataframe feature VS response
# > wraps makeOneSummaryByResponse, does testing for multiple features
# and does FDR and "nice" output
# > input is same as for
# makeOneSummaryByResponse, but ftrs is vector of variable names
# =====================================================================
makeMultiSummaryByResponse <- function(inDF,ftrs,response,doTestVsControl=T,controlClass="HC",
                                       formatSci=F,formatDig=3,formatRoundDig=3,doFDR=T,verbose=F,
                                       includeTotals=T,doSort=T,addResponseVarCol=T)
{
  ret <- data.frame()
  # debug / reality check
  if (sum(ftrs %in% colnames(inDF)) < length(ftrs)) {
    print(" WARNING: 1 or more feature variables missing from input dataframe (inDF): ")
    for (f in ftrs) {
      if (!(f %in% colnames(inDF))) {
        print(paste0(' > ',f,' is missing!'))
      }
    }
    print('  >> Continuing with remaining features ...')
  }
  for (f in ftrs) {
    ret <- rbind.data.frame(ret,makeOneSummaryByResponse(inDF,f,response,doTestVsControl = doTestVsControl,
                                                         controlClass,formatSci,
                                                         formatDig,formatRoundDig,verbose = verbose,
                                                         includeTotals=includeTotals))
  }
  if (doFDR) {ret$FDR <- p.adjust(ret$pV)}
  else {ret$FDR <- ret$pV}
  ret$FDRfor <- format(ret$FDR,digits = formatDig,scientific = T)
  ret$FDRsig <- as.character(ret$FDR <= 0.05)
  ret$FDRsig[ret$FDR <= 0.05] <- "*"
  ret$FDRsig[ret$FDR <= 1.0e-5] <- "**"
  #ret$FDRsig[ret$FDR <= 1.0e-10] <- "***"
  ret$FDRsig[ret$FDRsig=="FALSE"] <- ""
  ret$FDRsig[is.na(ret$FDRsig)] <- ""
  ret$niceOut <- paste0(ret$niceOut," ",ret$FDRsig)
  if (doSort) {
    ret <- ret[order(ret$pV,decreasing = F),]
  }
  if (addResponseVarCol) {
    ret$ResponseVar <- response
  }
  ret
}


# ============= makes summary of one vector ==============
# ========================================================
# NOTES:
# > fixer: if T, tries to fix some 'common' problems, such as spaces-only input being counted as valid data
#
makeOneSummary <- function(ftrV,varName='',fixer=T,infnullIsNA=T,nothingIsNA=T,nrClassesToList=5) {
  ret = data.frame()
  if (fixer) {
    backToFac <- F
    if (class(ftrV) == "factor") {
      ftrV <- as.character(ftrV)
      backToFac <- T
      ftrV <- trimws(ftrV)
      ftrV[ftrV == " "] <- ""
      ftrV[ftrV == "  "] <- ""
      ftrV[ftrV == "   "] <- ""
      ftrV[ftrV == "    "] <- ""
      ftrV[ftrV == "     "] <- ""
      ftrV[ftrV == "na"] <- NA
      ftrV[ftrV == "n/a"] <- NA
      ftrV[ftrV == "NA"] <- NA
      ftrV[ftrV == "N/A"] <- NA
      if (infnullIsNA) {
        ftrV[ftrV == "null"] <- NA
        ftrV[ftrV == "inf"] <- NA
        ftrV[ftrV == "NULL"] <- NA
        ftrV[ftrV == "INF"] <- NA
        ftrV[is.infinite(ftrV)] <- NA
        ftrV[is.null(ftrV)] <- NA
      }
      if (nothingIsNA) {
        ftrV[ftrV == ""] <- NA
      }
      ftrV <- as.factor(ftrV)
    }
    if (class(ftrV) == "integer") {
      ftrV <- as.numeric(ftrV)
    }
    if (class(ftrV) == "numeric") {
      if (infnullIsNA) {
        ftrV[is.infinite(ftrV)] <- NA
        ftrV[is.null(ftrV)] <- NA
      }
    }
  }
  if (class(ftrV) == "numeric") {
      mn <-  mean(ftrV,na.rm = T)
      sd <-  sd(ftrV,na.rm = T)
      md <-  median(ftrV,na.rm = T)
      min <- min(ftrV,na.rm = T)
      q1  <- quantile(ftrV,na.rm = T)[[2]]
      q3  <- quantile(ftrV,na.rm = T)[[4]]
      max <- max(ftrV,na.rm = T)
      nr <- length(ftrV)
      nonzero <- sum(ftrV!=0,na.rm = T)
      prev <- nonzero/nr
      nas <- sum(is.na(ftrV))
      naprev <- nas/nr
      if (is.na(prev)) {prev = 0.0}
      if (is.na(naprev)) {naprev = 0.0}
      ret <- data.frame(Var=varName,dataType=class(ftrV),NR=nr,Mean=mn,SD=sd,Min=min,Q1=q1,Median=md,Q3=q3,Max=max,
                        NrNonZero=nonzero,propNonZero=prev,NrNAs=nas,propNAs=naprev)
      
  } else if (class(ftrV) == "factor") {
    nr <- length(ftrV)
    nas <- sum(is.na(ftrV))
    naprev <- nas/nr
    something <- sum(!is.na(ftrV))
    somethingPrev <- something/nr
    nrClasses <- length(levels(ftrV))
    if (is.na(somethingPrev)) {somethingPrev = 0.0}
    if (is.na(naprev)) {naprev = 0.0}
    ret <- data.frame(Var=varName,dataType=class(ftrV),NR=nr,NrClasses=nrClasses,  
                      NrNonNA=something,propNonNA=somethingPrev,
                      NrNAs=nas,propNAs=naprev)
    if (nrClassesToList > 0) {
      ftrNRs <- as.data.frame(table(ftrV))
      ftrNRs <- ftrNRs[order(ftrNRs$Freq,decreasing = T),]
      for (nrC in c(1:nrClassesToList)) {
        if (nrC <= nrow(ftrNRs)) {
          toAddDF <- data.frame(A=ftrNRs$ftr[nrC],
                               B=ftrNRs$Freq[nrC],
                               C=ftrNRs$Freq[nrC]/nr)
          
        } else {
          toAddDF <- data.frame(A=NA,B=NA,C=NA)
        }
        colnames(toAddDF) <- c(paste0('C_',nrC,'_VALUE'),paste0('C_',nrC,'_NR'),paste0('C_',nrC,'_Prop'))
        ret <- cbind.data.frame(ret,toAddDF)
      }
    }
  }
  ret
}

# count taxa & pwys
# ===================================
countMicrobiomeFeatures <- function(inDF) {
  res <- data.frame()
  taxa <- purgeMGNames(subsetMicrobiomeDF(inDF = inDF,verbose = F,getTaxa = T,getPWYs = F,getVFs = F,getCARDs = F,getPhenos = F,getDivs = F))
  sumTaxa <- 0
  print ('---- TAXA -----------------------')
  for (tn in c("t__","s__","g__","f__","c__","o__","p__")) {
    nrTaxa <- sum(grepl(paste0("^",tn),colnames(taxa)))
    print(paste0(tn,' NR = ', nrTaxa))
    sumTaxa <- sumTaxa + nrTaxa
  }
  print (paste0(' TOTAL TAXA: ',sumTaxa))
  print ('---------------------------------')
  
  print ('---- PWYs -----------------------')
  pwys <- subsetMicrobiomeDF(inDF = inDF,verbose = F,getTaxa = F,getPWYs = T,getVFs = F,getCARDs = F,getPhenos = F,getDivs = F)
  print (paste0(' TOTAL PWYs: ',length(colnames(pwys))))
  
  taxa <- subsetMicrobiomeDF(inDF = inDF,verbose = F,getTaxa = T,getPWYs = F,getVFs = F,getCARDs = F,getPhenos = F,getDivs = F)
  taxa <- subsetMicrobiomeDF(inDF = inDF,verbose = F,getTaxa = T,getPWYs = F,getVFs = F,getCARDs = F,getPhenos = F,getDivs = F)
  taxa <- subsetMicrobiomeDF(inDF = inDF,verbose = F,getTaxa = T,getPWYs = F,getVFs = F,getCARDs = F,getPhenos = F,getDivs = F)
  
}


# ======================================================
# make ordination plot
# ======================================================
plotOrdination <- function(inDF,responseVar,doCentroids=F) {
  # ordinate on the distance matrix
  e.sir.pcoa <- cmdscale( data, eig = T )
  
  # variance explained 
  variance <- head(eigenvals(e.sir.pcoa)/sum(eigenvals(e.sir.pcoa)))
  x_variance <- as.integer(variance[1]*100)
  y_variance <- as.integer(variance[2]*100)
  
  # get scores for plotting
  e.sir.scores <- as.data.frame( e.sir.pcoa$points )
  
  # read in metadata file
  e.sir.meta <- read.delim( args_list$args[2], header = T, sep = "\t", row.names = 1 )
  # append to e.sir.scores
  e.sir.scores.meta <- merge( e.sir.scores, e.sir.meta, by = 'row.names' )
  # set first column as rownames and remove it
  rownames( e.sir.scores.meta ) <- e.sir.scores.meta[,1]
  e.sir.scores.meta[,1] <- NULL
  # change colnames
  colnames(e.sir.scores.meta) <- c( "PCo1", "PCo2", "SubjectID" )
  
  # plot ordination
  png( args_list$args[3], width = 750, height = 600, res = 150 )
  
  ggplot( e.sir.scores.meta, aes(PCo1, PCo2, color=SubjectID) ) + 
    geom_point(size = 4, alpha = 0.75) + theme_classic() + 
    theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.ticks = element_blank(), axis.text = element_blank()) + 
    xlab(paste("PCo1 (",x_variance,"% variance explained)")) + ylab(paste("PCo2 (",y_variance,"% variance explained)"))
  
  temp <- dev.off()
}


# cute cute function for PCA plotting regression plotting
# takes pcas table and responsevar (has to be in pcas dataset)
# ===========================================
plotPCAs <- function(pcas,nrPCs=3,responseVar,doCentroids=T) {
  c = 1
  res = list()
  # PC 1-2
  g <- ggplot(pcas ,aes_string(x = "PC1",y="PC2",col=responseVar )) + geom_point(size=2.25,alpha=0.8) 
  pcas$rV <- pcas[[responseVar]]
  if (doCentroids) {
    centroids <- aggregate(cbind(PC1,PC2)~rV,pcas,mean)
    colnames(centroids)[[1]] <- responseVar
    g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)
  }
  res[[c]] <- g
  if (nrPCs >= 3) {
    # PC 2-3
    g <- ggplot(pcas ,aes_string(x = "PC2",y="PC3",col=responseVar )) + geom_point(size=2.25,alpha=0.8) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC2,PC3)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3))
    }
    res[[c]] <- g
  }
  if (nrPCs >= 4) {
    # PC 3-4
    g <- ggplot(pcas ,aes_string(x = "PC3",y="PC4",col=responseVar )) + geom_point(size=2.25,alpha=0.4) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC3,PC4)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4))
    }
    res[[c]] <- g
  }
  if (nrPCs >= 5) {
    # PC 4-5
    g <- ggplot(pcas ,aes_string(x = "PC4",y="PC5",col=responseVar )) + geom_point(size=2.25,alpha=0.4) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC4,PC5)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar    
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5))
    }
    res[[c]] <- g
  }
  # return list of plots
  res
}

# ===========================================
# TSNE plotter
# ===========================================
ggplotTsne <- function(tSNEdf,feature) {
  g <- ggplot(tSNEdf, aes_string(x = "V1",y= "V2",col=feature )) + geom_point(size=2.5,alpha=0.9)
  ref <- reformulate(feature,"cbind(V1,V2)")
  centroids <- aggregate(ref,tSNEdf,mean)
  g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes_string(x="V1",y="V2",col=feature),alpha=1)
  g
}

# ================================================
# function for linear regression plotting
# takes linear model fit as argument
# ===========================================
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "blue") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       #"Intercept =",signif(fit$coef[[1]],2 ),
                       " Coef =",signif(fit$coef[[2]], 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))
}

# tester for linear fit
testsLinearFit <- function(inDF,responseVar,varsToTest,correctors=c(),nonZero=F) {
  tests <- data.frame()
  for (f in varsToTest) {
    inDFt <- inDF
    if (nonZero) {
      inDFt <- inDF[inDF[[f]] > 0,]
    }
    print(paste(' >> testing',f))
    frm <- reformulate(termlabels = c(responseVar,correctors), response = f)
    mdl <- lm(data = inDFt, formula =  frm)
    rs <- summary(mdl)$r.squared
    pv <- summary(mdl)$coefficients[responseVar,4]
    coeff <- summary(mdl)$coefficients[responseVar,1]
    tests <- rbind.data.frame(tests,data.frame(feature=f,rsquared=rs,coef=coeff,pvalue=pv,zeros=sum(inDFt[[f]] == 0),nonzeros=sum(inDFt[[f]] > 0)) )
  }
  tests <- tests[order(tests$pvalue),]
  tests$FDR <- p.adjust(tests$pvalue)
  tests$NonZero <- nonZero
  tests
}

# plotter
plotTestsLM <- function(inDF,tests,responseVarLM,responseVarX,FDRcutoff=0.2,pVcutoff=0.05,pointColVar=NA) {
  testor <- tests$feature[tests$FDR < FDRcutoff & tests$pvalue < pVcutoff]
  testor <- na.omit(testor)
  for (f in testor) {
    print(f)
    xannot = min(inDF[[responseVarLM]])
    yannot = max(inDF[[f]])
    tcoef <- format(tests$coef[tests$feature==f],digits=2)
    tpv <- format(tests$pvalue[tests$feature==f],digits = 2)
    trsq <- format(tests$rsquared[tests$feature==f],digits = 2)
    if (tests$NonZero[[1]]) {
      inDFt <- inDF[inDF[[f]] > 0,]
    } else {
      inDFt <- inDF
    }
    if (class(inDFt[[responseVarX]]) == "factor") {
      g <- ggplot(inDFt) +
        geom_boxplot(aes_string(x=responseVarX,y=f,col=responseVarX),outlier.alpha = 0) +
        geom_jitter(aes_string(x=responseVarX,y=f,col=responseVarX),width=0.2,alpha=0.25) +
        stat_smooth(aes_string(x=responseVarLM,y=f),method="lm",fullrange=F)
    } else if (is.na(pointColVar)) {
      g <- ggplot(inDFt) + geom_point(aes_string(x=responseVarX,y=f,col=responseVarX)) +
        stat_smooth(aes_string(x=responseVarLM,y=f),method="lm",fullrange=F) 
    } else {
      g <- ggplot(inDFt) + geom_point(aes_string(x=responseVarX,y=f,col=pointColVar)) +
        stat_smooth(aes_string(x=responseVarLM,y=f),method="lm",fullrange=F) 
      
    }
    g <- g + annotate(geom="text", y=yannot, hjust = 0, x = -Inf,
                      label=paste0("  LM fit: coef=",tcoef,"; R^2=",trsq,"; P=",tpv) )
    print(g)
  }
}

# color pallete maker
# ================================================
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

# save pheatmap
# ================================================
save_pheatmap_png <- function(x, filename, width=1200, height=1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_svg <- function(x, filename, width=5, height=20) {
  svg(filename, width=width, height=height)
  #grid::grid.newpage()
  #grid::grid.draw(x$gtable)
  print(x)
  dev.off()
}


# AGGREGATOR & PLOTTER for phenotype tests in DAG3
# =======================================================
# - agregates multiple tables and plots nice heatmap

# doFeatures can be: 
# - Abundances : abundances of taxa
# - Prevalences : prevalences of taxa
# - div :    diversity metrics
# - pwyAb :  abundances of PWYs

# NOTE:
# - currently works only on binary features (0 vs 1 OR n vs y)
# - metric should be : statsig, 

# NOTES:
# - returns list (heatmap,table of log transformed effect/p-values,non-transformed table)
#

aggregatePlotHeatmap <- function(tblFolder,doFeatures=c("Abundances"),nrToPlot=20,outFile="plot.png",pVcutoff=0.05,multiTest=F,
                                 clusterFeatures=T,clusterPhenos=F,subsetFeatures="",subsetPhenos="",
                                 subsetFeaturesGrep="",
                                 metric="statSig",enrichmentCutoff = 1.005,addXextra=0,addYextra=0,maxMax = 5,
                                 shortenTaxNames=T,invertYN=F,fixPhenos=T,cleanPhenoNamesGrep="") {
  if (!metric %in% c("statSig","meanEffectRatio","medEffectRatio","beta"))  {
    print (" >> metric MUST be one of [statSig, meanEffectRatio, medEffectRatio,beta]")
    break()
  }
  
  print ('>>> AGGREGATING RESULTS:')
  if (multiTest) {
    print (paste0(' doing multiple-testing correction (BH), using FDR = ',pVcutoff))
  } else {
    print (paste0(' using nominal p-value cutoff of ',pVcutoff))
  }

  hmDFsig <- NULL
  hmDFenrich <- NULL
  hmDFdirection <- NULL
  hmDFpv <- NULL
  mSS <- 0
  #maxMax <- 5
  # calculate NR of tests
  nrTests <- 0
  if (multiTest) {
    for (f in dir(tblFolder, pattern =".csv")) {
      doParse = F
      for (g in doFeatures) {if (grepl(g,f)) {doParse = T}}
      if (doParse) {
        #print(f)
        iT <- read.table(paste0(tblFolder,'/',f),sep=',',header = T,stringsAsFactors = F)
        #print(nrow(iT))
        if (subsetFeatures[[1]]!= "") {
          iT <- iT[iT$feature %in% subsetFeatures,]
        }
        if (subsetFeaturesGrep!="") {
          iT <- iT[grep(subsetFeaturesGrep,iT$feature),]
        }
        if (subsetPhenos[[1]]!= "")  {
          iT <- iT[iT$PHENOTYPE %in% subsetPhenos,]
        }
        nrTests <- nrTests+nrow(iT)
        #print (nrTests)
      }
    }
    print (paste0('  >> doing multiple-testing correction for ',nrTests,' tests!'))
  }
  
  print ('  >>> parsing tables ...')
  for (f in dir(tblFolder, pattern =".csv")) {
    doParse = F
    for (g in doFeatures) {if (grepl(g,f)) {doParse = T}}
    if (doParse) {
      print(f)
      iT <- read.table(paste0(tblFolder,'/',f),sep=',',header = T,stringsAsFactors = F)
      if (subsetFeatures[[1]]!= "") {
        iT <- iT[iT$feature %in% subsetFeatures,]
      }
      if (subsetFeaturesGrep!="") {
        iT <- iT[grep(subsetFeaturesGrep,iT$feature),]
      }
      if (subsetPhenos[[1]]!= "")  {
        iT <- iT[iT$PHENOTYPE %in% subsetPhenos,]
      }
      
      if (fixPhenos) {
        if (!"PHENOTYPE" %in% colnames(iT)) {
          iT$PHENOTYPE=gsub('_diversities','',gsub('phen_','',gsub('\\.csv','',f)))
        }
      }
      if (shortenTaxNames) {
        for (r in c(1:nrow(iT))) {
          if (grepl('__',iT$feature[r])) {
            iT$feature[r] <- purgeMGNameOne(iT$feature[r])
          }
        }
      }
      iT$V1 <- tolower(iT$V1)
      iT$V2 <- tolower(iT$V2)
      if (subsetFeatures[[1]] != "") {
        iT <- iT[subsetFeatures %in% iT$feature,]
      }
      if (subsetPhenos[[1]] != "")  {
        iT <- iT[iT$PHENOTYPE %in% subsetPhenos,]
      }
      if (multiTest) {
        if ("pCorr" %in% colnames(iT)) {
          iT$FDR <- p.adjust(iT$pCorr,n = nrTests)
        } else {
          iT$FDR <- p.adjust(iT$pValue,n = nrTests)
        }
      } else {
        if ("pCorr" %in% colnames(iT)) {
          iT$FDR <- iT$pCorr
        } else {
          iT$FDR <- iT$pValue
        }
      }
      if (!(iT$V2[1] %in% c("y","1","n","0"))) {
        print(paste0('WARNING: table ',f, ' / responseVar ',iT$V2[1],' : responseVar is not categorical 0/1 or n/y, skipping it ...'))
      } else {
        # add row(s) to dataframe
        # > pick what to plot
        if (metric == "statSig") {
          # > statistical significance 
          if (grepl("Prevalences", f)) {
            oneCol <- data.frame(feature=iT$feature,toplot=log(iT$FDR,base = 10)*-1*sign(1/iT$V1toV2-1)); colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
            oneColraw <- data.frame(feature=iT$feature,toplot=iT$FDR); colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          } else {
            oneCol <- data.frame(feature=iT$feature,toplot=log(iT$FDR,base = 10)*sign(iT$V1minusV2)); colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
            oneColraw <- data.frame(feature=iT$feature,toplot=iT$FDR); colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          }
        } else if (metric == "meanEffectRatio") {
          # > mean enrichment / depletion
          if (grepl("Prevalences", f)) {
            # > prevalences
            mS <- iT$V2prop/iT$V1prop
            mS <- mS[!is.na(mS)]; mS <- mS[!is.infinite(mS)]; mSS <- max(mS,mSS)
            eSize <- iT$V2prop/iT$V1prop
            eSize[is.na(eSize)] <- 0
            eSize[eSize == 1] <- 0
            eSize[eSize < 1] <- -1*log(1/eSize[eSize < 1])
            #print (eSize)
            eSize[eSize > 1] <- log(eSize[eSize > 1])
            eSize[is.infinite(eSize)] <- mSS*1.05
            oneCol <- data.frame(feature=iT$feature,toplot=(eSize)); 
            oneCol$toplot[iT$FDR > pVcutoff] <- 0.0
            #colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          } else {
            # > abundances
            mS <- iT$V2mean/iT$V1mean
            mS <- mS[!is.na(mS)]; mS <- mS[!is.infinite(mS)]; mSS <- max(mS,mSS)
            eSize <- iT$V2mean/iT$V1mean
            eSize[is.na(eSize)] <- 0
            eSize[eSize == 1] <- 0
            
            #eSize[eSize <= 0] <- abs(eSize[eSize <= 0])
            eSize[eSize <= 0] <- 0.0
            
            eSize[eSize < 1] <- -1*log(1/eSize[eSize < 1])
            eSize[eSize > 1] <- log(eSize[eSize > 1])
            eSize[is.infinite(eSize)] <- mSS*1.05
            oneCol <- data.frame(feature=iT$feature,toplot=(eSize))
            oneCol$toplot[iT$FDR > pVcutoff] <- 0.0
            #colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          }
        } else if (metric == "medEffectRatio") {
          # > median enrichment / depletion
          if (grepl("Prevalences", f)) {
            # > prevalences
            mS <- iT$V2prop/iT$V1prop
            mS <- mS[!is.na(mS)]; mS <- mS[!is.infinite(mS)]; mSS <- max(mS,mSS)
            eSize <- iT$V2prop/iT$V1prop
            eSize[is.na(eSize)] <- 0
            eSize[eSize == 1] <- 0
            eSize[eSize < 1] <- -1*log(1/eSize[eSize < 1])
            eSize[eSize > 1] <- log(eSize[eSize > 1])
            eSize[is.infinite(eSize)] <- mSS*1.05
            oneCol <- data.frame(feature=iT$feature,toplot=(eSize)); 
            oneCol$toplot[iT$FDR > pVcutoff] <- 0.0
            #colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          } else {
            # > abundances
            mS <- iT$V2median/iT$V1median
            mS <- mS[!is.na(mS)]; mS <- mS[!is.infinite(mS)]; mSS <- max(mS,mSS)
            eSize <- iT$V2median/iT$V1median
            eSize[eSize == 1] <- 0
            eSize[is.na(eSize)] <- 0
            #eSize[eSize <= 0] <- abs(eSize[eSize <= 0])
            eSize[eSize <= 0] <- 0.0
            eSize[eSize < 1] <- -1*log(1/eSize[eSize < 1])
            eSize[eSize > 1] <- log(eSize[eSize > 1])
            eSize[is.infinite(eSize)] <- mSS*1.05
            oneCol <- data.frame(feature=iT$feature,toplot=(eSize))
            oneCol$toplot[iT$FDR > pVcutoff] <- 0.0
            #colnames(oneCol) <- c("Feature",iT$PHENOTYPE[1])
          }
        }
        
        colnames(oneCol) <- c("Feature",gsub(cleanPhenoNamesGrep,"",iT$PHENOTYPE[1]))
        colnames(oneColraw) <- c("Feature",gsub(cleanPhenoNamesGrep,"",iT$PHENOTYPE[1]))
        
        if (is.null(hmDFsig)) {
          hmDFsig <- oneCol
          hmDFraw <- oneColraw
        } else {
          hmDFsig <- merge.data.frame(hmDFsig,oneCol,by = "Feature")
          hmDFraw <- merge.data.frame(hmDFraw,oneColraw,by = "Feature")
        }
      }
    }
  }
  hmDFsig <- hmDFsig[!duplicated(hmDFsig$Feature),]
  rownames(hmDFsig) <- hmDFsig$Feature
  hmDFsig$Feature <- NULL
  noSig = FALSE
  hmPlotDF <- hmDFsig
  if (metric == "statSig") {
    hmPlotDF[abs(hmPlotDF) < abs(log(pVcutoff,base=10)) ] <- 0.0
  } else if (metric == "meanEffectRatio") {
    hmPlotDF[abs(hmPlotDF) < log(enrichmentCutoff)] <- 0.0
  }
  hmPlotDF[is.na(hmPlotDF)] <- 0.0
  if (sum(hmPlotDF != 0) == 0) {
    hmPlotDF <- hmDFsig
    noSig = TRUE
    hmPlotDF <- hmPlotDF+rnorm(ncol(hmPlotDF)*nrow(hmPlotDF))/1000
    print(hmPlotDF)
  }
  hmPlotDF$sorter <- apply(hmPlotDF,MARGIN = 1,function(x) sum(abs(x)))
  
  hmPlotDF <- hmPlotDF[order(hmPlotDF$sorter,decreasing = T),]
  hmPlotDF$sorter <- NULL
  hmPlotDF <- hmPlotDF[1:min(nrToPlot,nrow(hmPlotDF)),]
  
  if (invertYN) {hmPlotDF = hmPlotDF * -1}
  hmDFdirection <- hmPlotDF 
  hmDFdirection[hmDFdirection > 0] <- 1
  hmDFdirection[hmDFdirection < 0] <- -1
  hmDFdirection[hmDFdirection == 0] <- 0
  hmDFdirection[hmDFdirection == 1] <- "+"
  hmDFdirection[hmDFdirection == -1] <- "-"
  hmDFdirection[hmDFdirection == 0] <- ""
  
  # set limits
  hmPlotDF[hmPlotDF >= maxMax] <- maxMax
  hmPlotDF[hmPlotDF <= -1*maxMax] <- -1*maxMax
  
  # set colors
  breaksList = seq(-max(abs(hmPlotDF)), +max(abs(hmPlotDF)),by=2 * max(abs(hmPlotDF))/55)
  cols <- colorRampPalette(c("Darkblue","White","Darkred"))(length(breaksList))
  cols[seq(length(cols)/2-1,length(cols)/2+1)] <- "White"

  if (noSig) {cols <- "white"; hmDFdirection[!is.na(hmDFdirection)] <- ""}
  # plot it!
  print ('  >>> Plotting')
  p <- pheatmap(hmPlotDF,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(hmPlotDF),legend = T,
                show_rownames = T, cluster_rows = clusterFeatures,cluster_cols = clusterPhenos,angle_col = 90,display_numbers = hmDFdirection,
                fontsize_number = 20,border_color = "#EEEEEE",na_col = "white",
                color = cols,treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",
                fontsize_col = 14,fontsize_row = 12,breaks = breaksList)
  
  
  
  xAdd = 0
  if ("pwyAbundance" %in% doFeatures) {xAdd=400+nrow(hmPlotDF)*4}
  print(addXextra)
  wtX = addXextra+xAdd+800+ncol(hmPlotDF)*35+nrow(hmPlotDF)*4
  print(wtX)
  save_pheatmap_png(p,outFile,width = wtX,
                    height = 800+nrow(hmPlotDF)*45+addYextra
                    ,res = 200)
  print (' >>> DONE')
  return( list(p,hmDFsig,hmDFraw) )
}

# ==================================================================
# ==================================================================
#
# AGGREGATOR AND PLOTTER FOR LINEAR MODELS
#
# ==================================================================
# ==================================================================

aggregatePlotHeatmapLM <- function(tblFolder,doFeatures=c("taxAbundances"),nrToPlot=20,outFile="plot.png",pVcutoff=0.05,multiTest=F,
                                 clusterFeatures=T,clusterPhenos=F,subsetFeatures="",subsetPhenos="",
                                 metric="statSig",enrichmentCutoff = 1.005,addXextra=0,addYextra=0,maxMax = 5,
                                 shortenTaxNames=F) {
  print ('>>> AGGREGATING RESULTS:')
  if (multiTest) {
    print (paste0(' doing multiple-testing correction (BH), using FDR = ',pVcutoff))
  } else {
    print (paste0(' using nominal p-value cutoff of ',pVcutoff))
  }
  if (!metric %in% c("statSig","beta","rsquared"))  {
    print (" >> metric MUST be one of [statSig, beta, rsquared]")
    break()
  }
  hmDFsig <- NULL
  hmDFenrich <- NULL
  hmDFdirection <- NULL
  mSS <- 0
  #maxMax <- 5
  nrTests <- 0
  if (multiTest) {
    for (f in dir(tblFolder, pattern =".csv")) {
      doParse = F
      for (g in doFeatures) {if (grepl(g,f)) {doParse = T}}
      if (doParse) {
        #print(f)
        iT <- read.table(paste0(tblFolder,'/',f),sep=',',header = T,stringsAsFactors = F)
        if (subsetFeatures != "") {
          iT <- iT[grep(subsetFeatures,iT$responseVar),]
        }
        if (subsetPhenos != "")  {
          iT <- iT[iT$Phenotype %in% subsetPhenos,]
        }
        nrTests <- nrTests+nrow(iT)
      }
    }
    print (paste0('  >> doing multiple-testing correction for ',nrTests,' tests!'))
  }
  print ('  >>> parsing tables ...')
  nrP <- 0
  for (f in dir(tblFolder, pattern =".csv")) {
    doParse = F
    for (g in doFeatures) {if (grepl(g,f)) {doParse = T}}
    if (doParse) {
      nrP <- nrP + 1
      #print(f)
      iT <- read.table(paste0(tblFolder,'/',f),sep=',',header = T,stringsAsFactors = F)
      if (subsetFeatures != "") {
        iT <- iT[grep(subsetFeatures,iT$responseVar),]
      }
      if (subsetPhenos != "")  {
        iT <- iT[iT$Phenotype %in% subsetPhenos,]
      }
      if (multiTest) {
        iT$FDR <- p.adjust(iT$pValue,n = nrTests)
      } else {
        iT$FDR <- iT$pValue
      }

      # add row(s) to dataframe
      # > pick what to plot
      if (metric == "statSig") {
        # > statistical significance 
        oneCol <- data.frame(feature=iT$responseVar,toplot=log(iT$FDR,base = 10)*-1*sign(iT$beta)); colnames(oneCol) <- c("Feature",iT$Phenotype[1])
      } else if (metric == "beta") {
        # > beta
        oneCol <- data.frame(feature=iT$responseVar,toplot=iT$beta ); colnames(oneCol) <- c("Feature",iT$Phenotype[1])
        # > test for significance
        oneCol[[2]][iT$FDR > pVcutoff] <- 0.0
      } else if (metric == "rsquared") {
        # > median enrichment / depletion
        iT$adjRsq[iT$adjRsq < 0] <- 0
        oneCol <- data.frame(feature=iT$responseVar,toplot=iT$adjRsq*1*sign(iT$beta) ); colnames(oneCol) <- c("Feature",iT$Phenotype[1])
        # > test for significance
        oneCol[[2]][iT$FDR > pVcutoff] <- 0.0
      }
      if (is.null(hmDFsig)) {
        hmDFsig <- oneCol
      } else {
        hmDFsig <- merge.data.frame(hmDFsig,oneCol,by = "Feature")
      }
    }
  }
  print (paste0('   >>> DONE, parsed, ',nrP,' tables'))
  
  hmDFsig <- hmDFsig[!duplicated(hmDFsig$Feature),]
  rownames(hmDFsig) <- hmDFsig$Feature
  hmDFsig$Feature <- NULL
  noSig = FALSE
  hmPlotDF <- hmDFsig
  if (metric == "statSig") {
    hmPlotDF[abs(hmPlotDF) < abs(log(pVcutoff,base=10)) ] <- 0.0
  } 
  if (sum(hmPlotDF != 0) == 0) {
    hmPlotDF <- hmDFsig
    noSig = TRUE
    hmPlotDF <- hmPlotDF+rnorm(ncol(hmPlotDF)*nrow(hmPlotDF))/1000
    #print(hmPlotDF)
  }
  hmPlotDF$sorter <- apply(hmPlotDF,MARGIN = 1,function(x) sum(abs(x)))
  
  hmPlotDF <- hmPlotDF[order(hmPlotDF$sorter,decreasing = T),]
  hmPlotDF$sorter <- NULL
  hmPlotDF <- hmPlotDF[1:min(nrToPlot,nrow(hmPlotDF)),]
  
  hmDFdirection <- hmPlotDF 
  hmDFdirection[hmDFdirection > 0] <- 1
  hmDFdirection[hmDFdirection < 0] <- -1
  hmDFdirection[hmDFdirection == 0] <- 0
  hmDFdirection[hmDFdirection == 1] <- "+"
  hmDFdirection[hmDFdirection == -1] <- "-"
  hmDFdirection[hmDFdirection == 0] <- ""
  
  # set limits
  hmPlotDF[hmPlotDF >= maxMax] <- maxMax
  hmPlotDF[hmPlotDF <= -1*maxMax] <- -1*maxMax
  
  # set colors
  breaksList = seq(-max(abs(hmPlotDF)), +max(abs(hmPlotDF)),by=2 * max(abs(hmPlotDF))/55)
  cols <- colorRampPalette(c("Darkblue","White","Darkred"))(length(breaksList))
  cols[seq(length(cols)/2-1,length(cols)/2+1)] <- "White"
  
  if (noSig) {cols <- "white"; hmDFdirection[!is.na(hmDFdirection)] <- ""}
  # plot it!
  print ('  >>> Plotting')
  p <- pheatmap(hmPlotDF,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(hmPlotDF),legend = T,
                show_rownames = T, cluster_rows = clusterFeatures,cluster_cols = clusterPhenos,angle_col = 90,display_numbers = hmDFdirection,
                fontsize_number = 20,border_color = "#EEEEEE",na_col = "white",
                color = cols,treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",
                fontsize_col = 14,fontsize_row = 12,breaks = breaksList)
  
  xAdd = 0
  if ("pwyAbundance" %in% doFeatures) {xAdd=400+nrow(hmPlotDF)*4}
  save_pheatmap_png(p,outFile,width = addXextra+xAdd+800+ncol(hmPlotDF)*35+nrow(hmPlotDF)*4,
                    height = 800+nrow(hmPlotDF)*45+addYextra, res = 200)
  
  print (' >>> DONE')
  return(p)
}

testFeatures <- function(inDFf,featureType='taxa',response='Diagnosis',plotsFolder='',dispPlot = 'Pn',xLab='Diagnosis',FDRcut=F,cutoff=0.05) {
  # find features with abundance changed [enriched or depleted]
  # ========================================================
  toTest <- colnames(inDFf[colnames(inDFf)!=response])
  resAbZ <- data.frame()
  # > run pair-wise mann-u-whitney tests for relative abundance
  for (t in toTest) {
    resAbZ <- rbind.data.frame(resAbZ,testOneFeature(inDFf,saveFolder = F,feature = t,responseVar = response,doPlots = F,discardZeros = F))
  }
  resAbZ <- na.omit(resAbZ)
  # correct for multiple testing
  resAbZ$FDR <- p.adjust(resAbZ$pValue)
  resAbZ <- resAbZ[order(resAbZ$pValue),]
  resAbZ$V1med_minus_V2med <- resAbZ$V1median - resAbZ$V2median
  # mark enrichment vs depletion
  resAbZ$Change <- "Depletion"
  resAbZ$Change[resAbZ$V1minusV2 < 0] <- "Enrichment"
  
  # > make pretty plots for depleted features only
  c = 0
  if (FDRcut) {
    for (t in resAbZ$feature[resAbZ$Change == "Depletion" & resAbZ$FDR < cutoff]) {
      c = c + 1
      testOneFeature(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_depleted_",c,"_",t),feature = t,
                     responseVar = response,doPlots = T,discardZeros = F,display = dispPlot,xLab = xLab)
    }
    c = 0
    # > make pretty plots for enriched features only
    for (t in resAbZ$feature[resAbZ$Change == "Enrichment" &  resAbZ$FDR < cutoff]) {
      c = c + 1
      testOneFeature(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_enriched_",c,"_",t),feature = t,
                     responseVar = response,doPlots = T,discardZeros = F,display = dispPlot,xLab = xLab)
    }
  } else {
    for (t in resAbZ$feature[resAbZ$Change == "Depletion" & resAbZ$pValue < cutoff]) {
      c = c + 1
      testOneFeature(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_depleted_",c,"_",t),feature = t,
                     responseVar = response,doPlots = T,discardZeros = F,display = dispPlot,xLab = xLab)
    }
    c = 0
    # > make pretty plots for enriched features only
    for (t in resAbZ$feature[resAbZ$Change == "Enrichment" &  resAbZ$pValue < cutoff]) {
      c = c + 1
      testOneFeature(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_enriched_",c,"_",t),feature = t,
                     responseVar = response,doPlots = T,discardZeros = F,display = dispPlot,xLab = xLab)
    }
  }
  resAbZ
}

testFeaturesPrev <- function(inDFf,featureType='taxa',response='Diagnosis',plotsFolder='',dispPlot = 'Pn',xLab='Diagnosis',FDRcut=F,cutoff=0.05) {
  # find features with abundance changed [enriched or depleted]
  # ========================================================
  toTest <- colnames(inDFf[colnames(inDFf)!=response])
  resAbZ <- data.frame()
  # > run pair-wise mann-u-whitney tests for relative abundance
  for (t in toTest) {
    print(paste0('   >>> testing ',t))
    resAbZ <- rbind.data.frame(resAbZ,testOneFeaturePrevalence(inDFf,saveFolder = F,feature = t,responseVar = response,doPlots = F))
  }
  resAbZ <- na.omit(resAbZ)
  # correct for multiple testing
  resAbZ$FDR <- p.adjust(resAbZ$pValue)
  resAbZ <- resAbZ[order(resAbZ$pValue),]
  #resAbZ$V1med_minus_V2med <- resAbZ$V1median - resAbZ$V2median
  # mark enrichment vs depletion
  resAbZ$Change <- "Depletion"
  resAbZ$Change[resAbZ$V1toV2 < 1] <- "Enrichment"
  
  # > make pretty plots for depleted features only
  c = 0
  if (FDRcut) {
    for (t in resAbZ$feature[resAbZ$Change == "Depletion" & resAbZ$FDR < cutoff]) {
      c = c + 1
      testOneFeaturePrevalence(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_depleted_",c,"_",t),feature = t,
                               responseVar = response,doPlots = T,display = dispPlot,xLab = xLab)
    }
    c = 0
    # > make pretty plots for enriched features only
    for (t in resAbZ$feature[resAbZ$Change == "Enrichment" &  resAbZ$FDR < cutoff]) {
      c = c + 1
      testOneFeaturePrevalence(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_enriched_",c,"_",t),feature = t,
                               responseVar = response,doPlots = T,display = dispPlot,xLab = xLab)
    }
  } else {
    for (t in resAbZ$feature[resAbZ$Change == "Depletion" & resAbZ$pValue < cutoff]) {
      c = c + 1
      testOneFeaturePrevalence(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_depleted_",c,"_",t),feature = t,
                               responseVar = response,doPlots = T,display = dispPlot,xLab = xLab)
    }
    c = 0
    # > make pretty plots for enriched features only
    for (t in resAbZ$feature[resAbZ$Change == "Enrichment" &  resAbZ$pValue < cutoff]) {
      c = c + 1
      testOneFeaturePrevalence(inDFf,saveFolder = paste0(plotsFolder,"/",featureType,"_enriched_",c,"_",t),feature = t,
                               responseVar = response,doPlots = T,display = dispPlot,xLab = xLab)
    }
  }
  resAbZ
}

# SUMMARY MAKER
# ===========================================================
# > makes summary statistics for dataframe
# > returns list[numeric summaries, categorical summaries]
# ===========================================================
makeMultiSummary <- function(inDF,ftrsToTest,nrClassesToList=2) {
  # > initialize variables
  summNum <- NULL
  summCat <- NULL
  for (ftr in ftrsToTest) {
    print(paste0(' > making summary for ',ftr))
    oneSumm <- makeOneSummary(ftr = inDF[[ftr]],
                              varName = ftr,
                              nrClassesToList = nrClassesToList)
    #print(oneSumm)
    oneSumm$Var <- ftr
    if (oneSumm$dataType=="factor") {
      summCat <- rbind.data.frame(summCat,oneSumm)
    } else if (oneSumm$dataType=="numeric") {
      summNum <- rbind.data.frame(summNum,oneSumm)
    }
  }
  list(S_NUM=summNum,S_CAT=summCat)
}

# =================================================================
# heatmap plotter for DAG3 phenotype-taxa associations
# =================================================================
# Parameters:
# 
#  phenosToPlot : which phenotypes to plot (must match phenotype column names, no default)
#  inData* : loaded dataframe with associations (produced by Alex)
#  statToPlot [def=Padj]: what to plot, other options: effect, R2, Chisq
#  featuresToPlot [def=G]: which features to consider, options: "Taxa" for all taxa, "S","G","F","O","P","K" for appropriate levels
#  clusterPhenos [def=T]: if True, clusters phenotypes
#  modelToUse [def=FullModel]: which model to use, options: FullModel, mu, sigma, nu, MuNu, 
#  sigStat [def=Padj]: which stat to consider for singificance, options: Padj, Pvalue
#  sigStatNominal [def=Pvalue]: which stat to consider for nominal significance, options: Padj, Pvalue
#  sigValue [def=0.05]: significance cutoff for sigStat
#  sigNominal [def=0.05]: significance cutoff for sigStatNominal
#  addText [def=Direction]: other options: "","Significance"
#  plotValuesSig [def=all]: what to plot & color: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue)
#  plotTextSig [def=all]:   what to put text on: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue) 
#  directionValue [def=effect.mu]: which value to use to determine direction (default: effect.mu), other options: effect.nu, effect.sigma
#  transformValues [def=T]: if T, does -1*log transform on p-values [Padj, Pvalue], no transformation on effect, R2, Chisq
#  logLimit [def=10]: limit log-values to entered number

plotTaxaAssociationsBEZIHeatmap <- function(phenosToPlot,
                                        inData,
                                        nrFeaturesToPlot = 20,
                                        featuresToPlot = "G",
                                        clusterPhenos = T,
                                        clusterFeatures = T,
                                        modelToUse = "FullModel",
                                        sigStat = "Padj",
                                        sigStatNominal = "Pvalue",
                                        sigValue = 0.05, 
                                        sigNominal = 0.05, 
                                        statToPlot = "Pvalue", 
                                        addText = "Direction",
                                        plotValuesSig = "all",
                                        plotTextSig = "padjSig", 
                                        directionValue = "effect.mu", 
                                        transformValues=T,
                                        logLimit = 10.0,
                                        sigAll = 0.9, 
                                        textSize=15, 
                                        colLow = "Orange", 
                                        colHigh = "Darkblue",
                                        colTextSize = 12,
                                        rowTextSize = 12,
                                        legendTextSize = 10,
                                        cellWidth=15, 
                                        cellHeight=15,
                                        nrColors = 13
)
{
  # other variables
  featuresToPlotGrep = ""
  # ================================================
  # ERROR CHECK:
  # =================================================
  print ('>> Function plotTaxaAssociationsHeatmap started ...')
  if (!modelToUse %in% c("FullModel","mu","sigma","nu","MuNu")) {
    stop("ERROR: modelToUse parameter MUST BE one of [FullModel, mu, sigma, nu, MuNu]")
  }
  if (!sigStat %in% c("Pvalue","Padj")) {
    stop("ERROR: sigStat parameter MUST BE one of [Pvalue, Padj]")
  }
  if (!statToPlot %in% c("Pvalue","Padj","effect","R2","Chisq")) {
    stop("ERROR: statToPlot parameter MUST BE one of [Pvalue, Padj, R2, Chisq, effect]")
  }
  if (!file.exists(dataPath)) {
    stop(paste0("ERROR: ",dataPath," does not exist!"))
  }
  if (!directionValue %in% c("effect.mu","effect.nu","effect.sigma")) {
    stop(paste0("ERROR: ",directionStat," MUST BE one of [effect.mu, effect.nu, effect.sigma]"))
  }
  if (!addText %in% c("","Direction","Significance")) {
    stop(paste0("ERROR: addText parameter MUST BE one of ['Direction', 'Significance', '']"))
  }
  if (!plotValuesSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotValuesSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  if (!plotTextSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotTextSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  
  if (featuresToPlot == "Taxa") {
    featuresToPlotGrep = "*"
  } else if (featuresToPlot == "S") {
    featuresToPlotGrep = "s__"
  } else if (featuresToPlot == "G") {
    featuresToPlotGrep = "g__"
  } else if (featuresToPlot == "F") {
    featuresToPlotGrep = "f__"
  } else if (featuresToPlot == "O") {
    featuresToPlotGrep = "o__"
  } else if (featuresToPlot == "P") {
    featuresToPlotGrep = 'p__'
  } else if (featuresToPlot == "K") {
    featuresToPlotGrep = 'k__'
  } else {
    stop(paste0("ERROR: featuresToPlot MUST BE one of [Taxa, S, G, F, O, P, K]"))
  }
  # ================================================
  
  # load data
  # ===========================
  print ('> loading data')
  inDF <- inData#read.table(dataPath,sep='\t',header=T,stringsAsFactors = F)
  # test for missing phenotypes
  if (sum(inDF$phenotype %in% phenosToPlot) == 0) {
    stop(paste0("ERROR: no phenotypes could be found, check phenosToPlot & dataPath parameters!"))
  }
  # fix column names
  inDFf <- inDF[inDF$phenotype %in% phenosToPlot,]
  colnames(inDFf) <- gsub('\\.full$','\\.FullModel',colnames(inDFf))
  toPlotValue = paste0(statToPlot,'.',modelToUse)
  nominalSigParameter = paste0(sigStatNominal,'.',modelToUse)
  sigParameter = paste0(sigStat,'.',modelToUse)
  # test for non-existing models (statistic to plot)
  if (!toPlotValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat "',statToPlot,'" + model "',modelToUse,'" does not exist in data, try different modelToUse or statToPlot!'))
  }
  # test for non-existing models (nominal significance)
  if (!nominalSigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat "',nominalSigParameter,'" + model "',modelToUse,'" does not exist in data, try different modelToUse or nominalSigParameter!'))
  }
  # test for non-existing direction values
  if (!directionValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: directionStat "',directionStat,'" + model "',modelToUse,'" does not exist in data, try directionStat = R2'))
  }
  # test for non-existing direction values
  if (!sigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: sigStat "',directionStat,'" + model "',modelToUse,'" does not exist in data, try different modelToUse or sigStat!'))
  }
  
  # collect data from results file
  inDFfs <- inDFf[,c("phenotype","taxon",toPlotValue,sigParameter,nominalSigParameter,directionValue)]
  # collect only values within nominal significance (as defined by parameters)
  inDFfs <- inDFfs[inDFfs[[nominalSigParameter]] < sigAll,]
  inDFfs <- inDFfs[complete.cases(inDFfs),]
  colnames(inDFfs) <- c("Phenotype","Taxon","ValueToPlot","Significance","NominalSig","Direction")
  if (nrow(inDFfs) == 0) {
    stop('WARNING: No significant results found, stopping!')
  }
  # shorten names
  for (r in c(1:nrow(inDFfs))) {
    #print(inDFfs$Taxon[r])
    inDFfs$Taxon[r] <- (purgeMGNameOne(inDFfs$Taxon[r]))
  }
  # make wide dataframe for heatmap plotting
  # ======================================
  inDFwide <- reshape(inDFfs, idvar = "Phenotype", timevar = "Taxon", direction = "wide")
  # subset taxa
  rownames(inDFwide) <- inDFwide$Phenotype
  inDFwide <- inDFwide[,grep(featuresToPlotGrep,colnames(inDFwide))]
  
  # extract individual dataframes
  # >>> DIRECTION DF
  inDirWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Direction',colnames(inDFwide)))]
  colnames(inDirWide) <- gsub('Direction\\.','',colnames(inDirWide))
  #  NAs for direction: values are 0 for purposes of plotting (no direction)
  for (cn in colnames(inDirWide)) {
    inDirWide[[cn]][is.na(inDirWide[[cn]])] <- 0.0
    inDirWide[[cn]] <- sign(as.numeric(inDirWide[[cn]]))
  }
  
  # >>> VALUE DF
  inValueWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^ValueToPlot',colnames(inDFwide)))]
  colnames(inValueWide) <- gsub('ValueToPlot\\.','',colnames(inValueWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  if (statToPlot %in% c("Pvalue","Padj")) {
    for (cn in colnames(inValueWide)) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 1.0
    }
  } else {
    for (cn in colnames(inValueWide)) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 0.0
    }
  }
  # >>> SIGNIFICANCE DATAFRAME (adjusted SIG)
  inSigAdjWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Significance',colnames(inDFwide)))]
  colnames(inSigAdjWide) <- gsub('Significance\\.','',colnames(inSigAdjWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigAdjWide)) {
    inSigAdjWide[[cn]][is.na(inSigAdjWide[[cn]])] <- 1.0
  }
  # >>> NOMINAL SIGNIFICANCE DATAFRAME (SIG)
  inSigNomWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^NominalSig',colnames(inDFwide)))]
  colnames(inSigNomWide) <- gsub('NominalSig\\.','',colnames(inSigNomWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigNomWide)) {
    inSigNomWide[[cn]][is.na(inSigNomWide[[cn]])] <- 1.0
  }
  # >>> PREP TEXT DATAFRAME (with +/- for direction or * for significance)
  doAddText = F
  if (addText == "Direction") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]][inTextWide[[cn]]=="1"] <- "+"
      inTextWide[[cn]][inTextWide[[cn]]=="-1"] <- "-"
      inTextWide[[cn]][inTextWide[[cn]]=="0"] <- ""
      # remove nonsignificant
      #inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  } else if (addText == "Significance") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]] <- "*"
      # remove nonsignificant
      inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  }
  # >>> filter dataframes by significance if needed
  #  > heatmap values
  #    > filter by nominal significance
  if (plotValuesSig == "nomSig") {
    if (statToPlot %in% c("Pvalue","Padj")) {
      inValueWide[inSigNomWide > sigNominal] <- 1.0
    } else {
      inValueWide[inSigNomWide > sigNominal] <- 0.0
    }
  } else if (plotValuesSig == "padjSig") {
    if (statToPlot %in% c("Pvalue","Padj")) {
      inValueWide[inSigAdjWide > sigValue] <- 1.0
    } else {
      inValueWide[inSigAdjWide > sigValue] <- 0.0
    }
  }
  # > heatmap text
  #    > filter by nominal significance or padj significance
  if (plotTextSig == "nomSig") {
    inTextWide[inSigNomWide > sigNominal] <- ""
  } else if (plotTextSig == "padjSig") {
    inTextWide[inSigAdjWide > sigValue] <- ""
  }
  
  # sort dataframes by number of associations
  dfNRAS <- as.data.frame(apply(inSigNomWide,MARGIN = 2,FUN = function(x) {sum(x <= sigNominal)}))
  dfNRAS$Var <- rownames(dfNRAS)
  colnames(dfNRAS) <- c("NR.Associations","Variable")
  dfNRAS <- dfNRAS[order(dfNRAS$NR.Associations,decreasing = T),]
  toKeep <- dfNRAS$Variable[1:nrFeaturesToPlot]
  
  # >>> transform values
  if (transformValues) {
    if (statToPlot %in% c("Pvalue","Padj")) {
      inValueWide <- log10(inValueWide)*-1
      for (cn in colnames(inValueWide)) {
        inValueWide[[cn]][inValueWide[[cn]] > logLimit] <- logLimit
      }
    } else {
      inValueWide <- abs(inValueWide)
    }
  }
  # >>> now multiply values by direction (for coloring of the plot)
  inValueMulWide <- inValueWide*inDirWide
  
  # >>> make heatmap
  # ==================================
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  
  # set colors
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  
  paletteLength <- nrColors
  myColor <- colorRampPalette(c(colLow, "white", colHigh))(paletteLength)
  myBreaks <- c(seq(min(plotDFValues), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(plotDFValues)/paletteLength, max(plotDFValues), length.out=floor(paletteLength/2)))
  
  p <- pheatmap(plotDFValues,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(plotDFValues),legend = T,
                show_rownames = T, cluster_rows = clusterPhenos,cluster_cols = clusterFeatures,angle_col = 90,
                fontsize_number = textSize,border_color = "#EEEEEE",na_col = "white",fontsize = legendTextSize,
                treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",color = myColor,
                fontsize_col = colTextSize,fontsize_row = rowTextSize,display_numbers = plotDFText,breaks = myBreaks,
                cellwidth = cellWidth,cellheight = cellHeight)
  p
}

# =================================================================
# heatmap plotter for DAG3 phenotype-taxa associations
# =================================================================
# Parameters:
# 
#  phenosToPlot : which phenotypes to plot (must match phenotype column names, no default)
#  inData* : loaded dataframe with associations (produced by Alex)
#  statToPlot [def=pval]: what to plot, other options: coef, qval
#  featuresToPlot [def=G]: which features to consider, options: "Taxa" for all taxa, "S","G","F","O","P","K" for appropriate levels
#  clusterPhenos [def=T]: if True, clusters phenotypes
#  sigStat [def=pval]: which stat to consider for singificance, options: qval
#  sigStatNominal [def=pval]: which stat to consider for nominal significance, options: qval
#  sigValue [def=0.05]: significance cutoff for sigStat
#  sigNominal [def=0.05]: significance cutoff for sigStatNominal
#  addText [def=Direction]: other options: "","Significance"
#  plotValuesSig [def=all]: what to plot & color: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue)
#  plotTextSig [def=all]:   what to put text on: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue) 
#  transformValues [def=T]: if T, does -1*log transform on p-values [qval, pval], no transformation on effect
#  logLimit [def=10]: limit log-values to entered number
plotTaxaAssociationsMaaslinHeatmap <- function(phenosToPlot,
                                            inData,
                                            nrFeaturesToPlot = 20,
                                            featuresToPlot = "G",
                                            clusterPhenos = T,
                                            clusterFeatures = T,
                                            sigStat = "qval",
                                            sigStatNominal = "pval",
                                            sigValue = 0.05, 
                                            sigNominal = 0.05, 
                                            statToPlot = "pval", 
                                            addText = "Direction",
                                            plotValuesSig = "all",
                                            plotTextSig = "all", 
                                            transformValues=T,
                                            logLimit = 10.0,
                                            sigAll = 0.9, 
                                            textSize=15, 
                                            colLow = "Orange", 
                                            colHigh = "Darkblue",
                                            colTextSize = 12,
                                            rowTextSize = 12,
                                            legendTextSize = 10,
                                            cellWidth=15, 
                                            cellHeight=15,
                                            nrColors = 13,
                                            flipCoords = T,
                                            shortenTaxNames = T,
                                            returnWideDFs = T,
                                            forcePhenosOrder=T
)
{
  # other variables
  featuresToPlotGrep = ""
  directionValue = "coef"
  # ================================================
  # ERROR CHECK:
  # =================================================
  print ('>> Function plotTaxaAssociationsHeatmap started ...')
  if (!sigStat %in% c("qval","pval")) {
    stop("ERROR: sigStat parameter MUST BE one of [pval, qval]")
  }
  if (!statToPlot %in% c("pval","qval","coef")) {
    stop("ERROR: statToPlot parameter MUST BE one of [pval, qval, coef]")
  }
  if (!addText %in% c("","Direction","Significance")) {
    stop(paste0("ERROR: addText parameter MUST BE one of ['Direction', 'Significance', '']"))
  }
  if (!plotValuesSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotValuesSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  if (!plotTextSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotTextSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  
  if (shortenTaxNames) {
    inData$feature <- as.character(inData$feature)
    for (cc in seq(1,nrow(inData))) {
      inData$feature[cc] <- purgeMGNameOne(inData$feature[cc])
    }
  }
  
  if (featuresToPlot == "Taxa" | featuresToPlot == "All") {
    featuresToPlotGrep = "*"
  } else if (featuresToPlot == "S") {
    featuresToPlotGrep = "s__"
  } else if (featuresToPlot == "G") {
    featuresToPlotGrep = "g__"
  } else if (featuresToPlot == "F") {
    featuresToPlotGrep = "f__"
  } else if (featuresToPlot == "O") {
    featuresToPlotGrep = "o__"
  } else if (featuresToPlot == "P") {
    featuresToPlotGrep = 'p__'
  } else if (featuresToPlot == "K") {
    featuresToPlotGrep = 'k__'
  } else {
    stop(paste0("ERROR: featuresToPlot MUST BE one of [Taxa | All, S, G, F, O, P, K]"))
  }
  # ================================================
  
  # load data
  # ===========================
  print ('> loading data')
  inDF <- inData
  # test for missing phenotypes
  if (sum(inDF$metadata %in% phenosToPlot) == 0) {
    stop(paste0("ERROR: no phenotypes could be found, check phenosToPlot & dataPath parameters!"))
  }
  # fix column names
  inDFf <- inDF[inDF$metadata %in% phenosToPlot,]
  #colnames(inDFf) <- gsub('\\.full$','\\.FullModel',colnames(inDFf))
  toPlotValue = statToPlot
  nominalSigParameter = sigStatNominal
  sigParameter = sigStat
  # test for non-existing models (statistic to plot)
  if (!toPlotValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat "',statToPlot,'" + model "',modelToUse,'" does not exist in data, try different modelToUse or statToPlot!'))
  }
  # test for non-existing models (nominal significance)
  if (!nominalSigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat ',nominalSigParameter,' does not exist in data, try different nominalSigParameter!'))
  }
  # test for non-existing direction values
  if (!directionValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: directionStat ',directionValue,' does not exist in data, try directionStat = R2'))
  }
  # test for non-existing direction values
  if (!sigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: sigStat "',sigParameter,' does not exist in data, try different sigStat!'))
  }
  
  # collect data from results file
  inDFfs <- inDFf[,c("metadata","feature",toPlotValue,sigParameter,nominalSigParameter,directionValue)]
  # collect only values within nominal significance (as defined by parameters)
  inDFfs <- inDFfs[inDFfs[[nominalSigParameter]] < sigAll,]
  inDFfs <- inDFfs[complete.cases(inDFfs),]
  colnames(inDFfs) <- c("Phenotype","Feature","ValueToPlot","Significance","NominalSig","Direction")
  if (nrow(inDFfs) == 0) {
    stop('WARNING: No significant results found, stopping!')
  }
  # shorten names
  # for (r in c(1:nrow(inDFfs))) {
  #   #print(inDFfs$Taxon[r])
  #   inDFfs$Feature[r] <- (purgeMGNameOne(inDFfs$Feature[r]))
  # }
  # make wide dataframe for heatmap plotting
  # ======================================
  inDFwide <- reshape(inDFfs, idvar = "Phenotype", timevar = "Feature", direction = "wide")
  # subset taxa
  rownames(inDFwide) <- inDFwide$Phenotype
  inDFwide <- inDFwide[,grep(featuresToPlotGrep,colnames(inDFwide))]
  
  # extract individual dataframes
  # >>> DIRECTION DF
  inDirWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Direction',colnames(inDFwide)))]
  colnames(inDirWide) <- gsub('Direction\\.','',colnames(inDirWide))
  #  NAs for direction: values are 0 for purposes of plotting (no direction)
  for (cn in colnames(inDirWide)[colnames(inDirWide)!="Phenotype"] ) {
    inDirWide[[cn]] <- as.character(inDirWide[[cn]])
    inDirWide[[cn]][is.na(inDirWide[[cn]])] <- "0.0"
    inDirWide[[cn]] <- as.numeric(inDirWide[[cn]])
    inDirWide[[cn]] <- sign(as.numeric(inDirWide[[cn]]))
  }
  rownames(inDirWide) <- inDirWide$Phenotype
  inDirWide$Phenotype <- NULL
  
  # >>> VALUE DF
  inValueWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^ValueToPlot',colnames(inDFwide)))]
  colnames(inValueWide) <- gsub('ValueToPlot\\.','',colnames(inValueWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  if (statToPlot %in% c("pval","qval")) {
    for (cn in colnames(inValueWide)[colnames(inValueWide)!="Phenotype"]) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 1.0
    }
  } else {
    for (cn in colnames(inValueWide)[colnames(inValueWide)!="Phenotype"]) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 0.0
    }
  }
  rownames(inValueWide) <- inValueWide$Phenotype
  inValueWide$Phenotype <- NULL
  
  # >>> SIGNIFICANCE DATAFRAME (adjusted SIG)
  inSigAdjWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Significance',colnames(inDFwide)))]
  colnames(inSigAdjWide) <- gsub('Significance\\.','',colnames(inSigAdjWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigAdjWide)[colnames(inSigAdjWide)!="Phenotype"]) {
    inSigAdjWide[[cn]][is.na(inSigAdjWide[[cn]])] <- 1.0
  }
  rownames(inSigAdjWide) <- inSigAdjWide$Phenotype
  inSigAdjWide$Phenotype <- NULL
  
  # >>> NOMINAL SIGNIFICANCE DATAFRAME (SIG)
  inSigNomWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^NominalSig',colnames(inDFwide)))]
  colnames(inSigNomWide) <- gsub('NominalSig\\.','',colnames(inSigNomWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigNomWide)[colnames(inSigNomWide)!="Phenotype"]) {
    inSigNomWide[[cn]][is.na(inSigNomWide[[cn]])] <- 1.0
  }
  rownames(inSigNomWide) <- inSigNomWide$Phenotype
  inSigNomWide$Phenotype <- NULL
  
  # >>> PREP TEXT DATAFRAME (with +/- for direction or * for significance)
  doAddText = F
  if (addText == "Direction") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)[colnames(inTextWide)!="Phenotype"]) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]][inTextWide[[cn]]=="1"] <- "+"
      inTextWide[[cn]][inTextWide[[cn]]=="-1"] <- "-"
      inTextWide[[cn]][inTextWide[[cn]]=="0"] <- ""
      # remove nonsignificant
      #inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  } else if (addText == "Significance") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)[colnames(inTextWide)!="Phenotype"]) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]] <- "*"
      # remove nonsignificant
      inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  }
  rownames(inTextWide) <- inTextWide$Phenotype
  inTextWide$Phenotype <- NULL
  inValueWideOrg <- inValueWide
  # >>> filter dataframes by significance if needed
  #  > heatmap values
  #    > filter by nominal significance
  if (plotValuesSig == "nomSig") {
    if (statToPlot %in% c("qval","pval")) {
      inValueWide[inSigNomWide > sigNominal] <- 1.0
    } else {
      inValueWide[inSigNomWide > sigNominal] <- 0.0
    }
  } else if (plotValuesSig == "padjSig") {
    if (statToPlot %in% c("qval","pval")) {
      inValueWide[inSigAdjWide > sigValue] <- 1.0
    } else {
      inValueWide[inSigAdjWide > sigValue] <- 0.0
    }
  }
  # > heatmap text
  #    > filter by nominal significance or padj significance
  if (plotTextSig == "nomSig") {
    inTextWide[inSigNomWide > sigNominal] <- ""
  } else if (plotTextSig == "padjSig") {
    inTextWide[inSigAdjWide > sigValue] <- ""
  }
  
  # sort dataframes by number of associations
  dfNRAS <- as.data.frame(apply(inSigAdjWide,MARGIN = 2,FUN = function(x) {sum(x <= sigValue)}))
  dfNRAS$Var <- rownames(dfNRAS)
  colnames(dfNRAS) <- c("NR.Associations","Variable")
  dfNRAS <- dfNRAS[order(dfNRAS$NR.Associations,decreasing = T),]
  toKeep <- dfNRAS$Variable[1:nrFeaturesToPlot]
  
  # >>> transform values
  if (transformValues) {
    if (statToPlot %in% c("pval","qval")) {
      inValueWide <- log10(inValueWide)*-1
      for (cn in colnames(inValueWide)) {
        inValueWide[[cn]][inValueWide[[cn]] > logLimit] <- logLimit
      }
    } else {
      inValueWide <- abs(inValueWide)
    }
  }
  # >>> now multiply values by direction (for coloring of the plot)
  inValueMulWide <- inValueWide*inDirWide
  
  # >>> make heatmap
  # ==================================
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  
  # set colors
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  rownames(plotDFText) <- rownames(plotDFValues)
  
  
  paletteLength <- nrColors
  myColor <- colorRampPalette(c(colLow, "white", colHigh))(paletteLength)
  myBreaks <- c(seq(min(plotDFValues), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(plotDFValues)/paletteLength, max(plotDFValues), length.out=floor(paletteLength/2)))
  
  if (forcePhenosOrder) {
    phenosToPlotF <- phenosToPlot[phenosToPlot %in% colnames(plotDFValues)]
    plotDFValues <- plotDFValues[levels(phenosToPlotF),]
    plotDFText <- plotDFText[levels(phenosToPlotF),]
  }
  
  if (flipCoords) {
    plotDFValues <- t.data.frame(plotDFValues)
    plotDFText <- t.data.frame(plotDFText)
  }
  p <- pheatmap(plotDFValues,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(plotDFValues),legend = T,
                show_rownames = T, cluster_rows = clusterPhenos,cluster_cols = clusterFeatures,angle_col = 90,
                fontsize_number = textSize,border_color = "#EEEEEE",na_col = "white",fontsize = legendTextSize,
                treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",color = myColor,
                fontsize_col = colTextSize,fontsize_row = rowTextSize,display_numbers = plotDFText,breaks = myBreaks,
                cellwidth = cellWidth,cellheight = cellHeight)
  if (!returnWideDFs) {
    p
  } else {
    list(Values.pv=inValueWideOrg,Values.qv=inSigAdjWide,Dir=inDirWide,Plot=p)
  }
}

# =================================================================
# heatmap plotter for DAG3 phenotype-feature associations (V2)
# =================================================================
# Notes:
#   V1 uses (now obsolete) BEZI models, this one uses new batch of Alex's results
# Parameters:
# 
#  phenosToPlot : which phenotypes to plot (must match phenotype column names, no default)
#  inData* : loaded dataframe with associations (produced by Alex)
#  statToPlot [def=pval]: what to plot, other options: coef, qval
#  featuresToPlot [def=G]: which features to consider, options: "Taxa" for all taxa, "S","G","F","O","P","K" for appropriate levels
#  clusterPhenos [def=T]: if True, clusters phenotypes
#  sigStat [def=pval]: which stat to consider for singificance, options: qval
#  sigStatNominal [def=pval]: which stat to consider for nominal significance, options: qval
#  sigValue [def=0.05]: significance cutoff for sigStat
#  sigNominal [def=0.05]: significance cutoff for sigStatNominal
#  addText [def=Direction]: other options: "","Significance"
#  plotValuesSig [def=all]: what to plot & color: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue)
#  plotTextSig [def=all]:   what to put text on: all (no significance filter), nomSig (nominal signifiance < sigNominal) or padjSig (padj significance < sigValue) 
#  transformValues [def=T]: if T, does -1*log transform on p-values [qval, pval], no transformation on effect
#  logLimit [def=10]: limit log-values to entered number
plotAssociationsDag3HeatmapV2 <- function(phenosToPlot,
                                          inData,
                                          nrFeaturesToPlot = 20,
                                          featuresToPlot = "G",
                                          dataType = "taxon",
                                          clusterPhenos = T,
                                          clusterFeatures = T,
                                          sigStat = "PadjBH",
                                          sigStatNominal = "Pvalue",
                                          sigValue = 0.05, 
                                          sigNominal = 0.05, 
                                          statToPlot = "Pvalue", 
                                          addText = "Direction",
                                          plotValuesSig = "nomSig",
                                          plotTextSig = "padjSig", 
                                          transformValues=T,
                                          logLimit = 10.0,
                                          sigAll = 0.9, 
                                          textSize=15, 
                                          colLow = "Orange", 
                                          colHigh = "Darkblue",
                                          colTextSize = 12,
                                          rowTextSize = 12,
                                          legendTextSize = 10,
                                          cellWidth=15, 
                                          cellHeight=15,
                                          nrColors = 13,
                                          flipCoords = T,
                                          doRefactor=F,
                                          sortPhenos = F,
                                          refactorIncludeAllClasses=F,
                                          trimNames = F,
                                          retData=F,
                                          keepOrder=F,
                                          directionValue = "effect.size"
)

{
  phenoOrder = phenosToPlot
  # other variables
  featuresToPlotGrep = ""
  # ================================================
  # ERROR CHECK:
  # =================================================
  print ('>> Function plotTaxaAssociationsHeatmap started ...')
  if (!sigStat %in% c("Pvalue",	"PadjBH","PadjBonferroni") ) {
    stop("ERROR: sigStat parameter MUST BE one of [Pvalue, PadjBH, PadjBonferroni]")
  }
  if (!statToPlot %in% c("Pvalue",	"PadjBH","PadjBonferroni","F.stat","R2","effect.size") ) {
    stop("ERROR: statToPlot parameter MUST BE one of [Pvalue, PadjBH, PadjBonferroni, R2, F.stat, effect.size]")
  }
  if (!addText %in% c("","Direction","Significance")) {
    stop(paste0("ERROR: addText parameter MUST BE one of ['Direction', 'Significance', '']"))
  }
  if (!plotValuesSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotValuesSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  if (!plotTextSig %in% c("all","nomSig","padjSig")) {
    stop(paste0("ERROR: plotTextSig parameter MUST BE one of ['all', 'nomSig', 'padjSig']"))
  }
  if (!dataType %in% c("taxon","pathway")) {
    stop(paste0("ERROR: dataType parameter MUST BE one of ['taxon', 'pathway']"))
  }
  
  if (featuresToPlot == "All") {
    featuresToPlotGrep = "*"
  } else if (featuresToPlot == "Taxa") {
    featuresToPlotGrep = "[sgfopk]__"
  } else if (featuresToPlot == "S") {
    featuresToPlotGrep = "s__"
  } else if (featuresToPlot == "G") {
    featuresToPlotGrep = "g__"
  } else if (featuresToPlot == "F") {
    featuresToPlotGrep = "f__"
  } else if (featuresToPlot == "O") {
    featuresToPlotGrep = "o__"
  } else if (featuresToPlot == "P") {
    featuresToPlotGrep = 'p__'
  } else if (featuresToPlot == "K") {
    featuresToPlotGrep = 'k__'
  } else {
    stop(paste0("ERROR: featuresToPlot MUST BE one of [Taxa | All, S, G, F, O, P, K]"))
  }
  # ================================================
  
  # load data
  # ===========================
  print ('> loading data')
  inDF <- inData
  # test for missing phenotypes
  if (sum(inDF$phenotype %in% phenosToPlot) == 0) {
    stop(paste0("ERROR: no phenotypes could be found, check phenosToPlot & dataPath parameters!"))
  }
  # fix column names
  inDFf <- inDF[inDF$phenotype %in% phenosToPlot,]
  #colnames(inDFf) <- gsub('\\.full$','\\.FullModel',colnames(inDFf))
  toPlotValue = statToPlot
  nominalSigParameter = sigStatNominal
  sigParameter = sigStat
  
  # test for non-existing models (statistic to plot)
  if (!toPlotValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat "',statToPlot,' does not exist in data, try different statToPlot!'))
  }
  # test for non-existing models (nominal significance)
  if (!nominalSigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: stat "',nominalSigParameter,' does not exist in data, try different nominalSigParameter!'))
  }
  # test for non-existing direction values
  if (!directionValue %in% colnames(inDFf)) {
    stop(paste0('ERROR: directionStat "',directionValue,' does not exist in data, try directionStat = R2'))
  }
  # test for non-existing direction values
  if (!sigParameter %in% colnames(inDFf)) {
    stop(paste0('ERROR: sigStat "',sigParameter,'" does not exist in data, try different sigStat!'))
  }
  

  # ================================================================
  # > refactor: keep 2nd and other factors, make new variables for each of there
  # NOTE: we only do this if we plot direction!
  # otherwise direction is average of all directions!
  # ================================================================
  if (doRefactor) {
    print ('reshaping factor results to long table')
    # refactor factor variables into multiple variables
    # > find factors
    phenosUnique <- unique(inDFf$phenotype)
    phenosFac <- c()
    for (ph in phenosUnique) {
      tt <- inDFf$levels[inDFf$phenotype == ph][1]
      #print(tt)
      if (grepl('\\:', tt)) {
        phenosFac <- c(phenosFac,ph)
      }
    }
    #inDFfbak <- inDFf
    print ('reshaping factor results to long table')
    cntr = 0
    newDFlist <- list() # accumulate in list for speed
    for (ph in phenosFac) {
      phDF <- inDFf[inDFf$phenotype==ph,]
      inDFf <- inDFf[inDFf$phenotype!=ph,]
      phLvls <- phDF$levels[phDF$phenotype == ph][1]
      print(paste0(' > refactoring ',ph,' <levels: ',phLvls,'>'))
      for (cc in c(1:nrow(phDF))) {
        oneRow <- phDF[cc,]
        splitLvls <- unlist(strsplit(oneRow$levels,':')) # levels split
        splitSS <- unlist(strsplit(oneRow$levels_SampleSize,':'))   # sample size split
        splitES <- unlist(strsplit(oneRow$effect.size,':'))  # effect size split
        if (refactorIncludeAllClasses) {startC = 1} else {startC = 2}
        if (dataType=="taxon") {dataTyp="taxon"
        } else if (dataType=="pathway") {dataTyp="pathway"}
        for (ln in c(startC:length(splitLvls))) {
          newRow <- data.frame(phenotype=paste0(oneRow$phenotype,'.C',ln-1,'_',splitLvls[ln]),
                               taxon=oneRow[[dataTyp]],
                               Nsamples=oneRow$Nsamples,
                               levels=splitLvls[ln],
                               levels_SampleSize=splitSS[ln],
                               effect.size=splitES[ln],
                               effect.size.asInteger=oneRow$effect.size.asInteger,
                               R2=oneRow$R2,
                               F.stat=oneRow$F.stat,
                               Pvalue=oneRow$Pvalue,
                               PadjBH=oneRow$PadjBH,
                               PadjBonferroni=oneRow$PadjBonferroni)
          cntr <- cntr+1
          newDFlist[[cntr]] <- newRow
        }
      }
    }
    print ('  > making dataframe')
    newDF <- do.call(rbind.data.frame, newDFlist)
    print ('  > sorting dataframe')
    newDF <- newDF[order(newDF$phenotype),]
    # band-aid
    if ("FDR" %in% colnames(inDFf) & !("FDR" %in% colnames(newDF))) {
      newDF$FDR <- newDF$PadjBH
    }
    print ('  > merging with original')
    if (dataType=="taxon") {colnames(newDF)[2] <- "taxon"
    } else if (dataType=="pathway") {colnames(newDF)[2] <- "pathway"}
    
    inDFf <- rbind.data.frame(inDFf,newDF)
    print ('reshaping done!')
  } else {
    # average directions: if 2 classes, take 2nd, if > 3 classes, take average
    for (cc in c(1:nrow(inDFf))) {
      oneRow <- inDFf[cc,]
      if (grepl('\\:',oneRow$effect.size)) {
        splitLvls <- unlist(strsplit(oneRow$levels,':')) # levels split
        #splitSS <- unlist(strsplit(oneRow$levels_SampleSize,':'))   # sample size split
        splitES <- unlist(strsplit(oneRow$effect.size,':'))  # effect size split
        if (length(splitLvls==2)) {
          oneRow$effect.size <- splitES[2]
          inDFf[cc,] <- oneRow
        } else {
          oneRow$effect.size <- oneRow$effect.size.asInteger
          inDFf[cc,] <- oneRow
        }
      }
    }
  }
  
  # collect data from results file
  if (!(dataType %in% colnames(inDFf))) {
    stop(paste0(' ERROR: ',dataType,' column not in input data, note: double-check the input data!'))
  }
  inDFfs <- inDFf[,c("phenotype",dataType,statToPlot,sigParameter,nominalSigParameter,directionValue)]
  # collect only values within nominal significance (as defined by parameters)
  inDFfs <- inDFfs[inDFfs[[nominalSigParameter]] < sigAll,]
  inDFfs <- inDFfs[complete.cases(inDFfs),]
  colnames(inDFfs) <- c("Phenotype","Feature","ValueToPlot","Significance","NominalSig","Direction")
  if (nrow(inDFfs) == 0) {
    stop('WARNING: No significant results found, stopping!')
  }
  # shorten names
  if (dataType=="taxon") {
    inDFfs$Feature <- as.character(inDFfs$Feature)
    for (r in c(1:nrow(inDFfs))) {
      #print(inDFfs$Taxon[r])
      inDFfs$Feature[r] <- (purgeMGNameOne(inDFfs$Feature[r]))
    }
  }
  
  # make wide dataframe for heatmap plotting
  # ======================================
  inDFwide <- reshape(inDFfs, idvar = "Phenotype", timevar = "Feature", direction = "wide")
  # kill class bug
  for (cc in colnames(inDFwide)) {
    inDFwide[[cc]] <- as.character(inDFwide[[cc]])
  }
  
  # subset taxa
  rownames(inDFwide) <- inDFwide$Phenotype
  inDFwide <- inDFwide[,grep(featuresToPlotGrep,colnames(inDFwide))]
  if ("Phenotype" %in% colnames(inDFwide)) {
    inDFwide$Phenotype <- NULL
  }
  
  # extract individual dataframes
  # >>> DIRECTION DF
  inDirWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Direction',colnames(inDFwide)))]
  colnames(inDirWide) <- gsub('Direction\\.','',colnames(inDirWide))
  #  NAs for direction: values are 0 for purposes of plotting (no direction)
  for (cn in colnames(inDirWide)[colnames(inDirWide)!="Phenotype"] ) {
    inDirWide[[cn]] <- as.character(inDirWide[[cn]])
    inDirWide[[cn]][is.na(inDirWide[[cn]])] <- "0.0"
    inDirWide[[cn]] <- as.numeric(inDirWide[[cn]])
    inDirWide[[cn]] <- sign(as.numeric(inDirWide[[cn]]))
  }
  # rownames(inDirWide) <- inDirWide$Phenotype
  # inDirWide$Phenotype <- NULL
  
  # >>> VALUE DF
  inValueWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^ValueToPlot',colnames(inDFwide)))]
  colnames(inValueWide) <- gsub('ValueToPlot\\.','',colnames(inValueWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  if (statToPlot %in% c("Pvalue",	"PadjBH","PadjBonferroni")) {
    for (cn in colnames(inValueWide)[colnames(inValueWide)!="Phenotype"]) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 1.0
    }
  } else {
    for (cn in colnames(inValueWide)[colnames(inValueWide)!="Phenotype"]) {
      inValueWide[[cn]] <- as.numeric(inValueWide[[cn]])
      inValueWide[[cn]][is.na(inValueWide[[cn]])] <- 0.0
    }
  }
  # rownames(inValueWide) <- inValueWide$Phenotype
  # inValueWide$Phenotype <- NULL
  
  # >>> SIGNIFICANCE DATAFRAME (adjusted SIG)
  inSigAdjWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^Significance',colnames(inDFwide)))]
  colnames(inSigAdjWide) <- gsub('Significance\\.','',colnames(inSigAdjWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigAdjWide)[colnames(inSigAdjWide)!="Phenotype"]) {
    inSigAdjWide[[cn]][is.na(inSigAdjWide[[cn]])] <- 1.0
    inSigAdjWide[[cn]] <- as.numeric(as.character(inSigAdjWide[[cn]]))
  }
  # rownames(inSigAdjWide) <- inSigAdjWide$Phenotype
  # inSigAdjWide$Phenotype <- NULL
  
  # >>> NOMINAL SIGNIFICANCE DATAFRAME (SIG)
  inSigNomWide <- inDFwide[,c(grep('^Phenotype',colnames(inDFwide)),grep('^NominalSig',colnames(inDFwide)))]
  colnames(inSigNomWide) <- gsub('NominalSig\\.','',colnames(inSigNomWide))
  #  NAs for values:
  #  if using p-values: values are 1 for purposes of plotting (p-value 1)
  #  otherwise values are 0 for purposes of plotting (effect size, chi-squared statistic, correlation)
  for (cn in colnames(inSigNomWide)[colnames(inSigNomWide)!="Phenotype"]) {
    inSigNomWide[[cn]][is.na(inSigNomWide[[cn]])] <- 1.0
    inSigNomWide[[cn]] <- as.numeric(as.character(inSigNomWide[[cn]]))
  }
  # rownames(inSigNomWide) <- inSigNomWide$Phenotype
  # inSigNomWide$Phenotype <- NULL
  
  # >>> PREP TEXT DATAFRAME (with +/- for direction or * for significance)
  doAddText = F
  if (addText == "Direction") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)[colnames(inTextWide)!="Phenotype"]) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]][inTextWide[[cn]]=="1"] <- "+"
      inTextWide[[cn]][inTextWide[[cn]]=="-1"] <- "-"
      inTextWide[[cn]][inTextWide[[cn]]=="0"] <- ""
      # remove nonsignificant
      #inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  } else if (addText == "Significance") {
    doAddText = T
    inTextWide <- as.data.frame(inDirWide)
    for (cn in colnames(inTextWide)[colnames(inTextWide)!="Phenotype"]) {
      inTextWide[[cn]] <- as.character(inTextWide[[cn]])
      inTextWide[[cn]] <- "*"
      # remove nonsignificant
      inTextWide[[cn]][inValueWide[[cn]] > sigValue] <- ""
    }
  }
  # rownames(inTextWide) <- inTextWide$Phenotype
  # inTextWide$Phenotype <- NULL
  # >>> filter dataframes by significance if needed
  #  > heatmap values
  #    > filter by nominal significance
  if (plotValuesSig == "nomSig") {
    if (statToPlot %in% c("Pvalue",	"PadjBH","PadjBonferroni")) {
      inValueWide[inSigNomWide > sigNominal] <- 1.0
    } else {
      inValueWide[inSigNomWide > sigNominal] <- 0.0
    }
  } else if (plotValuesSig == "padjSig") {
    if (statToPlot %in% c("Pvalue",	"PadjBH","PadjBonferroni")) {
      inValueWide[inSigAdjWide > sigValue] <- 1.0
    } else {
      inValueWide[inSigAdjWide > sigValue] <- 0.0
    }
  }
  # > heatmap text
  #    > filter by nominal significance or padj significance
  if (plotTextSig == "nomSig") {
    inTextWide[inSigNomWide > sigNominal] <- ""
  } else if (plotTextSig == "padjSig") {
    inTextWide[inSigAdjWide > sigValue] <- ""
  }
  
  # sort dataframes by number of associations
  dfNRAS <- as.data.frame(apply(inSigAdjWide,MARGIN = 2,FUN = function(x) {sum(x <= sigValue)}))
  dfNRAS$Var <- rownames(dfNRAS)
  colnames(dfNRAS) <- c("NR.Associations","Variable")
  dfNRAS <- dfNRAS[order(dfNRAS$NR.Associations,decreasing = T),]
  toKeep <- dfNRAS$Variable[1:nrFeaturesToPlot]
  
  # >>> transform values
  if (transformValues) {
    if (statToPlot %in% c("Pvalue",	"PadjBH","PadjBonferroni")) {
      inValueWide <- log10(inValueWide)*-1
      for (cn in colnames(inValueWide)) {
        inValueWide[[cn]][inValueWide[[cn]] > logLimit] <- logLimit
      }
    } else {
      inValueWide <- abs(inValueWide)
    }
  }
  # >>> now multiply values by direction (for coloring of the plot)
  inValueMulWide <- inValueWide*inDirWide
  
  # >>> make heatmap
  # ==================================
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  
  # set colors
  plotDFValues <- inValueMulWide[colnames(inValueMulWide) %in% toKeep]
  plotDFText <- inTextWide[colnames(inTextWide) %in% toKeep]
  
  paletteLength <- nrColors
  myColor <- colorRampPalette(c(colLow, "white", colHigh))(paletteLength)
  myBreaks <- c(seq(min(plotDFValues), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(plotDFValues)/paletteLength, max(plotDFValues), length.out=floor(paletteLength/2)))

  # sort (if sortPhenos = T)
  if (sortPhenos) {
    plotDFValues <- plotDFValues[order(rownames(plotDFValues)),]
    plotDFText <- plotDFText[order(rownames(plotDFText)),]
  }
  # keep order (if keepOrder = T)
  if (keepOrder) {
    plotDFValues <- plotDFValues[phenosToPlot,]
    plotDFText <- plotDFText[phenosToPlot,]
  }
  
  if (flipCoords) {
    plotDFValues <- t.data.frame(plotDFValues)
    plotDFText <- t.data.frame(plotDFText)
  }
  if (!flipCoords) {
    tt = clusterPhenos
    clusterPhenos = clusterFeatures
    clusterFeatures = tt
  }
  
  p <- pheatmap(plotDFValues,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(plotDFValues),legend = T,
                show_rownames = T, cluster_rows = clusterFeatures,cluster_cols = clusterPhenos,angle_col = 90,
                fontsize_number = textSize,border_color = "#EEEEEE",na_col = "white",fontsize = legendTextSize,
                treeheight_row = 0, treeheight_col = 0,legend_labels = "sig*log(p-value)",color = myColor,
                fontsize_col = colTextSize,fontsize_row = rowTextSize,display_numbers = plotDFText,breaks = myBreaks,
                cellwidth = cellWidth,cellheight = cellHeight)
  if (retData) {
    list(plotDFValues,plotDFText,p)
  } else {
    p
  }
}

# =========================================================
# do rarefaction
# =========================================================
#rarify (min abundance = 0.001)
doRarefaction <- function(inDF,replacements=F,doAll=F,
                          steps = c(10,20,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8229),
                          bootstraps=10,extrapolate=F,
                          doTaxa = T)
{
  if (doTaxa) {
    allCohTaxFt <- inDF[,grep('__',colnames(inDF))]
    allCohTaxFt$Cohort <- inDF$Cohort
  } else {
    isNum <- c()
    for (cc in colnames(inDF)) {
      isNum <- c(isNum,is.numeric(inDF[[cc]]))
    }
    allCohTaxFt <- inDF[,isNum]
    allCohTaxFt$Cohort <- inDF$Cohort
  }
  boots <- bootstraps
  rareDf2 <- data.frame()
  cohorts = unique(inDF$Cohort)
  if (doAll) {cohorts <- c("All",cohorts)}
  #for (n in seq(10,nrow(allCohTaxF),100)) {
  for (coh in cohorts) {
    print (paste(" >> Rarefying",coh))
    for (n in steps) {
      if (coh != "All" & n > sum(allCohTaxFt$Cohort==coh) & !extrapolate) {
        break()
      } else if (!extrapolate & n > length(allCohTaxFt$Cohort)) {
        break()
      }
      print(n)
      mns <- c()
      for (b in c(1:boots)) {
        if (coh != "All") {
          t <- allCohTaxFt[allCohTaxFt$Cohort==coh,]
          if (extrapolate & n > length(t$Cohort)) {
            smp <- t[sample(nrow(t), n,replace = T), ]
          } else {
            smp <- t[sample(nrow(t), n,replace = replacements), ]
          }
        } else {
          if (extrapolate & n > length(t$Cohort)) {
            smp <- t[sample(nrow(t), n,replace = T), ]
          }
          else {
            smp <- t[sample(nrow(t), n,replace = replacements), ]
          }
        }
        smp$Cohort <- NULL
        csums <- colSums(smp)
        nrSpec <- sum(csums>0)
        mns <- c(mns,nrSpec)
      }
      mn <- mean(mns)
      sdev <- sd(mns)
      rareDf2 <- rbind.data.frame(rareDf2,data.frame(nr=n,spec.nr.mn=mn,spec.nr.sd=sdev,cohort=coh))
    }
  }
  rareDf2
}

# GLM summary to p-value table
# ==============================================
glmToTbl <- function(fittedGLM) {
  gSmr <- summary(fittedGLM)
  gCoef <- as.data.frame(gSmr$coefficients)
  colnames(gCoef) <- c("Estimate","Std.Error","T.value","P.value")
  gCoef <- gCoef[row.names(gCoef) != "(Intercept)",]
  gCoef2 <- t(gCoef)
  colnames(gCoef2) <- paste0('pV.',colnames(gCoef2) )
  gCoef3 <- as_tibble(gCoef2)[4,]
}

# GLM summary to p-value table
# ==============================================
glmToTbl2 <- function(fittedGLM) {
  gSmr <- summary(fittedGLM)
  gCoef <- as.data.frame(gSmr$coefficients)
  colnames(gCoef) <- c("Estimate","Std.Error","T.value","P.value")
  gCoef <- gCoef[row.names(gCoef) != "(Intercept)",]
  gCoef2 <- t(gCoef)
  gCoef2 <- gCoef2[c("Estimate","P.value"),]
  gCoefPV <- as.tibble(gCoef2)[2,]
  colnames(gCoefPV) <- paste0(colnames(gCoefPV),'.pV' )
  gCoefEst <- as.tibble(gCoef2)[1,]
  colnames(gCoefEst) <- paste0(colnames(gCoefEst),'.Est' )
  gCoef3 <- cbind(gCoefPV,gCoefEst)
  gCoef3 <- gCoef3[,order(colnames(gCoef3))]
  gCoef3
}

## get all comparisons function
# =============================================================================
getAllStatComparisons <- function(.tbl) {
  .tbl %>%
    dplyr::select(.data$group1, .data$group2) %>%
    purrr::transpose() %>%
    purrr::modify_depth(1, unlist)
}
