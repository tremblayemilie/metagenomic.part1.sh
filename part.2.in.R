---
title: "MClaude_pollen.Rmd"
author: "E.T."
date: '2018-03-07'
output: html_document
---

```{r}

#Emilie's script for metagenomics by organism
#========================================================



```

tax.fill to remove uncertae sedis, unclassified...

```{r}
library("RAM")

pathDataAnMC <- "/isilon/cfia-ottawa-fallowfield/users/tremblaye/Debbie-Pollen-Feb_2018/RdataAnMclaude/"

dir.create(paste(pathDataAnMC, "taxFill/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)
taxFillPathMC <- paste(pathDataAnMC, "taxFill/", sep = "")


##############
#            #
#    ITS     #
#            #
############## works

temp <- read.table(paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  
 

##############
#            #
#    Oom     #
#            #
############## works
temp <- read.table(paste(pathDataAnMC, "DP_merged_OM.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  
##############
#            #
# Oom cured  #
#            #
############## works
temp <- read.table(paste(pathDataAnMC, "DP_merged_OM.curated.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill.old.version(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  
##############
#            #
#  Phy       #
#            #
############## works

temp <- read.table(paste(pathDataAnMC, "DP_merged_Phy.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill.old.version(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "Phy.table.taxFill.test.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
  
##############
#            #
#  Phy long  #
#            #
############## works

temp <- read.table(paste(pathDataAnMC, "DP_merged_Phy.table.mod.vignae2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill.old.version(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
 
##############
#            #
#  Phy NSP   #
#            #
############## works

temp <- read.table(paste(pathDataAnMC, "DP.Phy.NSP.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill.old.version(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "Phy.table.NSP.taxFill.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
       
  
##############
#            #
#    plants  #
#            #
############## works

temp <- read.table(paste(pathDataAnMC, "DP_merged_plants.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", 
                      comment.char = "", 
                      quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("taxonomy"="character"))

  temp$row.names <- NULL
  temp <- tax.fill.old.version(temp, downstream=TRUE)
  write.table(temp, file=paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote = FALSE)
  
```

Create diversity directory and path. List your datasets for diversity indexes. Adds a bunch of columns in a new diversity metadata file.
Indices obtained (in this specific order) are: Spec Number,  Simpson data,  Inv simpson data,	Shannon data,	Simpson eveness,	Shannon eveness,	Simpson true diversity,	shannon true diversity,	chao,	ACE.
```{r}
dir.create(paste(pathDataAnMC, "diversity/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)
diversityPathMC <- paste(pathDataAnMC, "diversity/", sep = "")


##############
#            #
#    ITS     #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "ITS.meta.div.2.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)

##############
#            #
#    Oom     #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "OM.meta.div.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)

####################################
#                                  #
#    Oom just phytophthora spp.    #
#                                  #
#################################### works

tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.phytophthora.only.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "OM.phytophthora.only.meta.div.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)

  
##############
#            #
#  Oom cured #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "OM.meta.div.cur.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)  
  

##############
#            #
# Phy        #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "Phy.meta.div.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)
  
  
##############
#            #
# Phy long   #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "Phy.meta.div.long.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)  
  

##############
#            #
# Phy NSP    #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.NSP.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "Phy.meta.div.NSP.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)  
    
 
##############
#            #
#    plants  #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

temp2 <- OTU.diversity(list(data=tabletemp), metatemp)

  write.table(temp2, file=paste(diversityPathMC, "DP.plants.meta.div.tsv", sep = ""),
              append    = FALSE,
              sep       = "\t",
              row.names = FALSE,
              quote=FALSE)
  
     
```

Whiskers plots (diversity)
```{r}

dir.create(paste(diversityPathMC, "diversityPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)
diversityPlotPathMC <- paste(diversityPathMC, "diversityPlot/", sep = "")

dir.create(paste(diversityPathMC, "EvennessPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)
EvennessPlotPathMC <- paste(diversityPathMC, "EvennessPlot/", sep = "")



##############
#            #
#    ITS     #
#            #
############## works


tabletemp <- read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "ITS.meta.div.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "ITS", ".Simpson.Evenness.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "ITS", ".Shannon.Evenness.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "ITS", ".Simpson.True.Div.2.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "ITS", ".Shannon.True.Div.2.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "ITS", ".Simpson.Div.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "ITS", ".Shannon.Div.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "ITS", ".Simpson.inv.Div.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "ITS", ".Species.nb.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


  
##############
#            #
#    Oom     #
#            #
############## works
  
tabletemp <- read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "OM.meta.div.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "OM", ".Simpson.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "OM", ".Shannon.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.True.Div.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Shannon.True.Div.2.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Shannon.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.inv.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Species.nb.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

####################################
#                                  #
#    Oom just phytophthora spp.    #
#                                  #
#################################### works
  
 tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.phytophthora.only.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
tabletemp$row.names <- NULL #if not using tXfilled tables you need it 
  
metatemp <-  read.table(paste(diversityPathMC, "OM.phytophthora.only.meta.div.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "OM", ".phytophthora.only.Shannon.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".phytophthora.only.Shannon.True.Div.2.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()  
  
  
#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".phytophthora.only.Shannon.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()
  
#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".phytophthora.only.Species.nb.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()  
  
  
##############
#            #
#  Oom cured #
#            #
############## works
  
tabletemp <- read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "OM.meta.div.cur.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "OM", ".Simpson.Evenness.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "OM", ".Shannon.Evenness.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.True.Div.cur.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Shannon.True.Div.cur.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.Div.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Shannon.Div.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Simpson.inv.Div.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "OM", ".Species.nb.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()
  
  
  
##############
#            #
#    Phy     #
#            #
############## works

  
tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "Phy.meta.div.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Simpson.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Shannon.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.True.Div.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.True.Div.2.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.inv.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Species.nb.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()
##############
#            #
#  Phy long  #
#            #
############## works

  
tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "Phy.meta.div.long.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Simpson.Evenness.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Shannon.Evenness.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.True.Div.long.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.True.Div.long.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.Div.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.Div.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.inv.Div.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Species.nb.long.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
##############
#            #
#  Phy NSP  #
#            #
############## works

  
tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.NSP.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "Phy.meta.div.NSP.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Simpson.Evenness.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "Phy", ".Shannon.Evenness.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.True.Div.NSP.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.True.Div.NSP.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.Div.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Shannon.Div.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Simpson.inv.Div.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "Phy", ".Species.nb.NSP.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()  

  
##############
#            #
#    plants  #
#            #
############## works


tabletemp <- read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "DP.plants.meta.div.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#Simpson Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "DP.plants", ".Simpson.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_even"), diversity.info=TRUE, compare="NULL", ylab="Simpson Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Evenness

# Open file-handle to get ready to make the png graph:
  png(filename=paste(EvennessPlotPathMC, "DP.plants", ".Shannon.Evenness.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_even"), diversity.info=TRUE, compare="NULL", ylab="Shannon Evenness", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

  
#Simpson True Diversity

# Open file-handle to get ready to make the png graph:

  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Simpson.True.Div.png", sep = ""),
      width = 1200, height = 600, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Simpson True Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon True Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Shannon.True.Div.png", sep = ""),
      width = 3000, height = 1200, units = "px", res = 300)

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan_trudiv"), diversity.info=TRUE, compare="NULL", ylab="Units of number of species", facet.x.cex = "10", facet.y.cex = "10", facet.y = TRUE)

# Close the png graph file handle:
  dev.off()

    
#Simpson Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Simpson.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("sim"), diversity.info=TRUE, compare="NULL", ylab="Simpson Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Shannon Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Shannon.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("shan"), diversity.info=TRUE, compare="NULL", ylab="Shannon Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()

#Simpson inverted Diversity

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Simpson.inv.Div.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("invsim"), diversity.info=TRUE, compare="NULL", ylab="Simpson inv Div", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()


#Species number

# Open file-handle to get ready to make the png graph:
  png(filename=paste(diversityPlotPathMC, "DP.plants", ".Species.nb.png", sep = ""),
      width = 1200, height = 400, units = "px")

myGroup.diversity(list(data=tabletemp), metatemp, factors = c("CollectionDate"), indices = c("spec"), diversity.info=TRUE, compare="NULL", ylab="Species nb", facet.x.cex = "10", facet.y.cex = "10")

# Close the png graph file handle:
  dev.off()
  
    
    
```


Rarefy your tables
```{r}

dir.create(paste(pathDataAnMC, "rarefaction/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

RarefactionPathMC <- paste(pathDataAnMC, "rarefaction/", sep = "")



##############
#            #
#    ITS     #
#            #
############## works

temp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character", "taxonomy"="character"))
  temp$row.names <- NULL 

rrf <- OTU.rarefy(list(data=temp), sample=NULL)
# Write out the dataframe:
  write.table(rrf, file=paste(RarefactionPathMC, "ITS.rarefy.2.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

  
temp2 <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
temp2$row.names <- NULL
temp2$taxonomy <- NULL

 


# Open file-handle to get ready to make the png graph:
  png(filename=paste(RarefactionPathMC, "ITS", ".RAREFplot.2.png", sep = ""),
      width = 1200, height = 400, units = "px")

rrfPLOT <- rarecurve(temp2, step =10, xlab = "Sequence nb per sample", ylab = "Species OTUs") 

# Close the png graph file handle:
  dev.off()
  
##############
#            #
#    Oom     #
#            #
############## works

temp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character", "taxonomy"="character"))
  temp$row.names <- NULL 

rrf <- OTU.rarefy(list(data=temp), sample=NULL)
# Write out the dataframe:
  write.table(rrf, file=paste(RarefactionPathMC, "OM.rarefy.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

  
temp2 <-  read.table(paste(pathDataAnMC, "DP_merged_OM.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
temp2$row.names <- NULL
temp2$taxonomy <- NULL

 


# Open file-handle to get ready to make the png graph:
  png(filename=paste(RarefactionPathMC, "OM", ".RAREFplot.png", sep = ""),
      width = 1200, height = 400, units = "px")

rrfPLOT <- rarecurve(temp2, step =10, xlab = "Sequence nb per sample", ylab = "Species OTUs") 

# Close the png graph file handle:
  dev.off()
  
##############
#            #
#    Oom  cur   #
#            #
############## works

temp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.curated.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character", "taxonomy"="character"))
  temp$row.names <- NULL 

rrf <- OTU.rarefy(list(data=temp), sample=NULL)
# Write out the dataframe:
  write.table(rrf, file=paste(RarefactionPathMC, "OM.rarefy.cur.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

  
temp2 <-  read.table(paste(pathDataAnMC, "DP_merged_OM.curated.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
temp2$row.names <- NULL
temp2$taxonomy <- NULL

 


# Open file-handle to get ready to make the png graph:
  png(filename=paste(RarefactionPathMC, "OM", ".RAREFplot.cur.png", sep = ""),
      width = 1200, height = 400, units = "px")

rrfPLOT <- rarecurve(temp2, step =10, xlab = "Sequence nb per sample", ylab = "Species OTUs") 

# Close the png graph file handle:
  dev.off()
    
  

##############
#            #
#    Phy     #
#            #
############## works

temp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character", "taxonomy"="character"))
  temp$row.names <- NULL 

rrf <- OTU.rarefy(list(data=temp), sample=NULL)
# Write out the dataframe:
  write.table(rrf, file=paste(RarefactionPathMC, "Phy.rarefy.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

  
temp2 <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
temp2$row.names <- NULL
temp2$taxonomy <- NULL

 


# Open file-handle to get ready to make the png graph:
  png(filename=paste(RarefactionPathMC, "Phy", ".RAREFplot.png", sep = ""),
      width = 1200, height = 400, units = "px")

rrfPLOT <- rarecurve(temp2, step =10, xlab = "Sequence nb per sample", ylab = "Species OTUs") 

# Close the png graph file handle:
  dev.off()


##############
#            #
#    plants  #
#            #
############## works

temp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character", "taxonomy"="character"))
  temp$row.names <- NULL 

rrf <- OTU.rarefy(list(data=temp), sample=NULL)
# Write out the dataframe:
  write.table(rrf, file=paste(RarefactionPathMC, "DP.plants.rarefy.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

  
temp2 <-  read.table(paste(pathDataAnMC, "DP_merged_plants.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
temp2$row.names <- NULL
temp2$taxonomy <- NULL

 


# Open file-handle to get ready to make the png graph:
  png(filename=paste(RarefactionPathMC, "DP.plants", ".RAREFplot.png", sep = ""),
      width = 1200, height = 400, units = "px")

rrfPLOT <- rarecurve(temp2, step =10, xlab = "Sequence nb per sample", ylab = "Species OTUs") 

# Close the png graph file handle:
  dev.off()  
  
```  
   ##TaxAbund
### group.abund.Taxa
# 1. generate the list of taxa desired with core taxa, then generate the plots with group.abund.Taxa
```{r}   
 

library("RAM")

# make new directory:

dir.create(paste(pathDataAnMC, "groupTaxaBarPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

taxaBarPlotPathMC <- paste(pathDataAnMC, "groupTaxaBarPlot/", sep = "")


##############
#            #
#    ITS     #
#            #
##############

tabletemp <-  read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".species.all.2.png", sep = ""),
      width = 10000, height = 10000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()

#genera_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".genera.all.bestres.png", sep = ""),
      width = 8000, height = 8000, units = "px", res = 300)

MyGroup.abund.taxa.version2(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 
  
  
#family_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "f", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".family.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "f",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Family_all")
                
# Close the png graph file handle:
  dev.off()   
  
 #order_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "o", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".order.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "o",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Order_all")
                
# Close the png graph file handle:
  dev.off()   
   
  
 #class_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "c", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".class.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "c",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Class_all")
                
# Close the png graph file handle:
  dev.off()     
  
 #phyla_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "p", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0, number=5)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "ITS", ".phyla.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "p",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "phyla_all")
                
# Close the png graph file handle:
  dev.off()     
              
##############
#            #
#    Oom     #
#            #
############## works
  
tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".species.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()

#genera_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".genera.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 

  
#phyla
  

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "p", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".phyla.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 
  
####################################
#                                  #
#    Oom just phytophthora spp.    #
#                                  #
#################################### works
  
 tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.phytophthora.only.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
tabletemp$row.names <- NULL #if not using tXfilled tables you need it 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".phytophthora.only.species.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()


##############
#            #
#    Oom cur    #
#            #
############## works
  
tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
#tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".species.all.cur.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()
  


#genera_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "OM", ".genera.all.cur.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 
  
  
    
##############
#            #
#    Phy     #
#            #
############## works  
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "Phy", ".species.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()

#genera_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "Phy", ".genera.all.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 
  
##############
#            #
# Phy long   #
#            #
############## works  
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses = c("CollectionDate"= "Date"))


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "Phy", ".species.all.long.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()

##############
#            #
# Phy NSP    #
#            #
############## works  
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.NSP.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
#tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses = c("CollectionDate"= "Date"))


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "Phy", ".species.all.NSP.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()
    
##############
#            #
#    plants  #
#            #
##############

tabletemp <-  read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("row.names"="character","taxonomy"="character"))
  
tabletemp$row.names <- NULL 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".species.all.2.png", sep = ""),
      width = 10000, height = 10000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "s",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Species_all")
                
# Close the png graph file handle:
  dev.off()

#genera_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".genera.all.bestres.png", sep = ""),
      width = 8000, height = 8000, units = "px", res = 300)

MyGroup.abund.taxa.version2(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "g",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Genera_all")
                
# Close the png graph file handle:
  dev.off() 
  
  
#family_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "f", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".family.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "f",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Family_all")
                
# Close the png graph file handle:
  dev.off()   
  
 #order_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "o", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".order.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "o",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Order_all")
                
# Close the png graph file handle:
  dev.off()   
   
  
 #class_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "c", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".class.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "c",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "Class_all")
                
# Close the png graph file handle:
  dev.off()     
  
 #phyla_all

coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "p", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

#taxaDate <- coreTaxaITSdate$data$CollectionDate$taxa


# Open file-handle to get ready to make the png graph:
  png(filename=paste(taxaBarPlotPathMC, "DP.plants", ".phyla.all.2.png", sep = ""),
      width = 3000, height = 2000, units = "px", res = 100)

myGroup.abund.Taxa(data = list(data = tabletemp),
                    is.OTU = TRUE,
                    rank = "p",    
                    drop.unclassified = FALSE,
                    meta = metatemp,
                    meta.factor = "CollectionDate",
                    taxa = coreTaxaITSdate,
                    main = "phyla_all")
                
# Close the png graph file handle:
  dev.off()       
  
```
  # GroupTaxabar

```{r}

dir.create(paste(pathDataAnMC, "groupAbundnbBarPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

abundNbBarPlotPathMC <- paste(pathDataAnMC, "groupAbundnbBarPlot/", sep = "")

##############
#            #
#    ITS     #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.species.top30.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()


  #genus
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.genera.top15.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="g", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the genus level ", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  #family
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.family.top15.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="f", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the family level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
   

   #order
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.order.top10.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="o", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the order level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  

  
 #class
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.class.top10.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="c", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the class level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  

  #phyla
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "ITS", ".gr.abund.phyla.top5.2.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="p", top = 5, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 50 OTU abundance (count) at the phyla level ", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
      
##############
#            #
#    Oom     #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.species.top30.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()


  #genus
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.genera.top15.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="g", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the genus level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  
  
  #order

# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.order.top10.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="o", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the order level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  
  #family

# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.family.top15.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="f", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the family level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off() 

  
 ####################################
#                                  #
#    Oom just phytophthora spp.    #
#                                  #
#################################### works
  
 tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.phytophthora.only.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
tabletemp$row.names <- NULL #if not using tXfilled tables you need it 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work
 
#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".phytophthora.only.gr.abund.species.top30.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 50, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 50 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
    
##############
#            #
#    Oom cur     #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.species.top30.cur.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()


  #genus
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.genera.top15.cur.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="g", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the genus level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  
  
  #order

# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.order.top10.cur.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="o", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the order level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  
  #family

# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".gr.abund.family.top15.cur.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="f", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the family level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off() 


####################################
#                                  #
#    Oom cur just phytophthora spp.    #
#                                  #
#################################### works
  
 tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.curated.phytophthora.only.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
tabletemp$row.names <- NULL #if not using tXfilled tables you need it 
  
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work
 
#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "OM", ".phytophthora.only.gr.abund.species.top30.cur.png", sep = ""),
      width = 3500, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 50, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 50 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()  
  
    
    
##############
#            #
#    Phy     #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "Phy", ".gr.abund.species.top30.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

#remove unidentified
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.taxFill.no.unidentified.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work
  
  
  #species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "Phy", ".gr.abund.species.top30.minus.unidentified.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = TRUE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  
########  
  #genus
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "Phy", ".gr.abund.genera.top10.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="g", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the genus level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  
##############
#            #
#    Phy  long   #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "Phy", ".gr.abund.species.top30.long.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

#remove unidentified
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.no.unidentified.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work
  
  
  #species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "Phy", ".gr.abund.species.top30.minus.unidentified.long.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = TRUE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  
##############
#            #
#    plants  #
#            #
############## works

tabletemp <-  read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work



#species
# Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.species.top30.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="s", top = 30, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 30 counts of taxonomic groups at the species level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()


  #genus
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.genera.top15.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="g", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the genus level ", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()

  #family
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.family.top15.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="f", top = 15, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 15 counts of taxonomic groups at the family level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
   

   #order
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.order.top10.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="o", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the order level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  

  
 #class
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.class.top10.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="c", top = 10, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 10 counts of taxonomic groups at the class level", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  

  #phyla
  # Open file-handle to get ready to make the plot:
  png(filename=paste(abundNbBarPlotPathMC, "DP.plants", ".gr.abund.phyla.top5.png", sep = ""),
      width = 3000, height = 3000, units = "px", pointsize = 12, res = 300)

group.abundance.meta(data=(list(data=tabletemp)), rank="p", top = 5, count =  TRUE, drop.unclassified = FALSE, cex.x = 10, main = "Top 50 OTU abundance (count) at the phyla level ", meta = metatemp, meta.factor = c("CollectionDate"))

# Close the png graph file handle:
  dev.off()
  
  
```



```{r}
dir.create(paste(pathDataAnMC, "PCoAPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

PcoAPlotPathMC <- paste(pathDataAnMC, "PCoAPlot/", sep = "")

##############
#            #
#    ITS     #
#            #
##############

tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, row.names = 1)
tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)


tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 


metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

myvec <- c("CollectionDate")

names(myvec) <- "CollectionDate"

#either euclidean and normalize or hellinger and bray or jaccard (pick the one that represent more % of your data)
#stand.method	One of "total", "max", "freq", "normalize", "range", "standardize", "pa", "chi.square", "hellinger" or "log"
#dist.method	one of "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower","morisita", "horn", "mountford", "raup", "binomial", "chao", or "cao"

# Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "ITS", ".PCOA.genera.top.100.dist.euc.standnorm.2.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "normalize", dist.method = "euclidean", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
 
  # Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "ITS", ".PCOA.genera.top.100.distHell.standBray.2.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "hellinger", dist.method = "bray", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
  
##############
#            #
#    Oom     #
#            #
############## works
  
tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_OM.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, row.names = 1)
tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)


tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 


metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

myvec <- c("CollectionDate")

names(myvec) <- "CollectionDate"

#either euclidean and normalize or hellinger and bray or jaccard (pick the one that represent more % of your data)
#stand.method	One of "total", "max", "freq", "normalize", "range", "standardize", "pa", "chi.square", "hellinger" or "log"
#dist.method	one of "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower","morisita", "horn", "mountford", "raup", "binomial", "chao", or "cao"

# Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "OM", ".PCOA.genera.top.100.dist.euc.standnorm.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "normalize", dist.method = "euclidean", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
 
  # Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "OM", ".PCOA.genera.top.100.distHell.standBray.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "hellinger", dist.method = "bray", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
  

##############
#            #
#    Phy     #
#            #
############## 

tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_Phy.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, row.names = 1)
tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)


tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 


metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

myvec <- c("CollectionDate")

names(myvec) <- "CollectionDate"


# Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "Phy", ".PCOA.genera.top.10.dist.raup.standnorm.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "normalize", dist.method = "raup", sample.labels = FALSE, top = 10, ellipse = FALSE, main = "PCOA for top 10 Phytophthora genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
 
  # Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "Phy", ".PCOA.genera.top.10.dist.raup.standhell.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "hellinger", dist.method = "raup", sample.labels = FALSE, top = 10, ellipse = FALSE, main = "PCOA for top 10 Phytophthora genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
   

##############
#            #
#    plants  #
#            #
##############

tabletemp <- read.table(paste(pathDataAnMC, "DP_merged_plants.table.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, row.names = 1)
tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)


tabletemp$row.names <- NULL #if not using tXfilled tables you need it
 


metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

myvec <- c("CollectionDate")

names(myvec) <- "CollectionDate"

#either euclidean and normalize or hellinger and bray or jaccard (pick the one that represent more % of your data)
#stand.method	One of "total", "max", "freq", "normalize", "range", "standardize", "pa", "chi.square", "hellinger" or "log"
#dist.method	one of "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower","morisita", "horn", "mountford", "raup", "binomial", "chao", or "cao"

# Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "DP.plants", ".PCOA.genera.top.100.dist.euc.standnorm.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "normalize", dist.method = "euclidean", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
 
  # Open file-handle to get ready to make the png graph:
png(filename=paste(PcoAPlotPathMC, "DP.plants", ".PCOA.genera.top.100.distHell.standBray.png", sep = ""),
      width = 1200, height = 1000, units = "px", res = 100)

pcoa.plot(data = tabletemp, is.OTU = TRUE, meta = metatemp, factors = (myvec), rank = "g", stand.method = "hellinger", dist.method = "bray", sample.labels = FALSE, top = 100, ellipse = FALSE, main = "PCOA for top 100 Fungal genera", ggplot2 = TRUE, bw = FALSE)

# Close the png graph file handle:
  dev.off() 
    
    
```


#plot Venn Venn Plot

```{r}
dir.create(paste(pathDataAnMC, "VennPlot/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

VennPathMC <- paste(pathDataAnMC, "VennPlot/", sep = "")

##############
#            #
#    ITS     #
#            #
############## 

tabletemp <-  read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.species.months.2.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#genera
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.genera.months.2.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#black and white a.k.a grey
  
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.genera.months.bw.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("grey10", "grey40", "grey70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()  
    
######get partition of venn plots 

  library("VennDiagram")
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])



#see rows 1 to 3
partitions$..values..[1:3]
  
  
  #####
  
#family
  
 coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "f", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.family.months.2.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
 
#order
  
 coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "o", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.order.months.2.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  
#phylum
  
   coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "p", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "ITS", ".VennPlot.phylum.months.2.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  
 
#############
#           #
#  oom      #
#           #
#############
  
tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.species.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#genera
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.genera.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  

#grey black and white
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.genera.months.bw.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("grey10", "grey40", "grey70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
    
    
    library("VennDiagram")
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])



#see rows 1 to 3
partitions$..values..[1:3]
  
  
#############
#           #
#  oom  cur    #
#           #
#############
  
tabletemp <-  read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.species.months.cur.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#genera
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.genera.months.cur.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  

#grey black and white
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "OM", ".VennPlot.genera.months.bw.cur.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("grey10", "grey40", "grey70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
    
    
    library("VennDiagram")
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])



#see rows 1 to 3
partitions$..values..[1:3]
    
  
#############
#           #
#  Phy      #
#           #
#############
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa


taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=Ventaxamay2)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "Phy", ".VennPlot.species.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  

#gray
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "Phy", ".VennPlot.species.months.bw.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("gray10", "gray40", "gray70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
    
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])
  
  
  

#genera
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa


taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=Ventaxamay2)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "Phy", ".VennPlot.genera.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()  
   
  
#############
#           #
#  Phy long     #
#           #
#############
  
tabletemp <-  read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("CollectionDate"="Date"))


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa


taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=Ventaxamay2)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "Phy", ".VennPlot.species.months.long.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  

#gray
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "Phy", ".VennPlot.species.months.bw.long.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("gray10", "gray40", "gray70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
    
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])
  
##############
#            #
#    plants  #
#            #
############## 

tabletemp <-  read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))

tabletemp$row.names <- NULL 
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE)


metatemp$row.names <- NULL 
rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work  

#species
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "s", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.species.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#genera
coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "g", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)


#colours
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.genera.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()

#black and white a.k.a grey
  
# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.genera.months.bw.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("grey10", "grey40", "grey70"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()  
    
######get partition of venn plots 

  library("VennDiagram")
    
partitions <- get.venn.partitions(taxVectorMC)
View(partitions)

#to see all rows
partitions$..values..

#to see row 1
partitions$..values..[1]

#rename the list as you want 
May_June_July <- as.list(partitions$..values..[1])
May_July <- as.list(partitions$..values..[2])
May_June <- as.list(partitions$..values..[3])
May <- as.list(partitions$..values..[4])
June_July <- as.list(partitions$..values..[5])
July <- as.list(partitions$..values..[6])
June <- as.list(partitions$..values..[7])



#see rows 1 to 3
partitions$..values..[1:3]
  
  
  #####
  
#family
  
 coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "f", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.family.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
 
#order
  
 coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "o", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.order.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
  
#phylum
  
   coreTaxaITSdate <- core.Taxa(data = list(data = tabletemp), is.OTU = TRUE, meta = metatemp, rank = "p", drop.unclassified = FALSE, meta.factor = "CollectionDate", percent = 0)

Ventaxamay1 <- coreTaxaITSdate$data$'2017-05-27'$taxa
Ventaxamay2 <- coreTaxaITSdate$data$'2017-05-17'$taxa
Ventaxajun1 <- coreTaxaITSdate$data$'2017-06-02'$taxa
Ventaxajun2 <- coreTaxaITSdate$data$'2017-06-13'$taxa
Ventaxajun3 <- coreTaxaITSdate$data$'2017-06-19'$taxa
Ventaxajun4 <- coreTaxaITSdate$data$'2017-06-28'$taxa
Ventaxajui1 <- coreTaxaITSdate$data$'2017-07-12'$taxa

taxamay <- c(unique(coreTaxaITSdate$data$'2017-05-27'$taxa, coreTaxaITSdate$data$'2017-05-17'$taxa))
taxajuneAll <- c(unique(fromLast = TRUE, coreTaxaITSdate$data$'2017-06-02'$taxa, coreTaxaITSdate$data$'2017-06-13'$taxa, coreTaxaITSdate$data$'2017-06-19'$taxa, coreTaxaITSdate$data$'2017-06-28'$taxa))
taxajuly <- coreTaxaITSdate$data$'2017-07-12'$taxa
             
taxVectorMC <- list(June_2017=taxajuneAll, July_2017=taxajuly, May_2017=taxamay)



# Open file-handle to get ready to make the png graph:
  png(filename=paste(VennPathMC, "DP.plants", ".VennPlot.phylum.months.png", sep = ""),
      width = 1500, height = 1500, units = "px", res = 100)

group.venn(vectors = taxVectorMC, cat.cex =3, cex = 3, label = FALSE, lab.cex = 3, fill = c("darkorchid1", "aquamarine1", "navajowhite4"), lab.col = "black")
                 
# Close the png graph file handle:
  dev.off()
      
  
```

#Tax abund tables generated at each taxonomic rank
```{r}

dir.create(paste(pathDataAnMC, "abundTables/", sep = ""),
           showWarnings = TRUE,
           recursive    = FALSE)

abundTablePath <- paste(pathDataAnMC, "abundTables/", sep = "")

##############
#            #
#    ITS     #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = NULL, mode = "percent")
write.table(tmp, file=paste(abundTablePath, "ITS.taxab.table.percent.2.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = 50, mode = "number")
write.table(tmp, file=paste(abundTablePath, "ITS.taxab.table.top50.number.2.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

##############
#            #
#    OOM     #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = NULL, mode = "percent")
write.table(tmp, file=paste(abundTablePath, "OM.taxab.table.percent.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = 50, mode = "number")
write.table(tmp, file=paste(abundTablePath, "OM.taxab.table.top50.number.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

##############
#            #
#    OOM cur    #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = NULL, mode = "percent")
write.table(tmp, file=paste(abundTablePath, "OM.taxab.table.percent.cur.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = 50, mode = "number")
write.table(tmp, file=paste(abundTablePath, "OM.taxab.table.top50.number.cur.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

##############
#            #
#    PHY     #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = NULL, mode = "percent")
write.table(tmp, file=paste(abundTablePath, "Phy.taxab.table.percent.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = 10, mode = "number")
write.table(tmp, file=paste(abundTablePath, "Phy.taxab.table.top10.number.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

##############
#            #
#    plants  #
#            #
############## works

tabletemp <- read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

#metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = NULL, mode = "percent")
write.table(tmp, file=paste(abundTablePath, "DP.plants.taxab.table.percent.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)

tmp <- tax.abund(tabletemp, otu2 = NULL, rank = NULL, drop.unclassified =FALSE, top = 50, mode = "number")
write.table(tmp, file=paste(abundTablePath, "DP.plants.taxab.table.top50.number.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = FALSE,
                               quote=FALSE)


                                 
```  
  #group temporal
  
```{r}  
  
temporalPathMC <- paste(pathDataAnMC, "temporal/", sep = "")

#https://stackoverflow.com/questions/41592862/scale-x-date-remove-extra-month-on-axis-ggplot2-2-2-0-9
##############
#            #
#    ITS     #
#            #
##############

tabletemp <- read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_ITS.temporal.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_genus)
View(revamptemp$data_genus)


write.table(revamptemp$data_genus, file=paste(revampPathMC, "DP18_datarevamp.ITS.genus.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 25
is.numeric(colourCount) #must be true				


#genus
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "ITS", ".group.temporal.genus.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("Aspergillus", "Itersonilia", "belongs.to.k.Fungi", "Alternaria", "belongs.to.p.Ascomycota", "Filobasidium", "Epicoccum", "Vishniacozyma", "belongs.to.No.blast.hit", "Mycosphaerella", "Verticillium", "Penicillium", "belongs.to.o.Hypocreales", "Lapidomyces", "belongs.to.o.Pleosporales", "Sporobolomyces", "belongs.to.c.Dothideomycetes", "Arthrocatena", "Lemonniera", "Blumeria", "Taphrina", "Septoriella", "Mycocentrospora", "Rhodotorula", "Fusarium"), group = c("Aspergillus", "Itersonilia", "belongs.to.k.Fungi", "Alternaria", "belongs.to.p.Ascomycota", "Filobasidium", "Epicoccum", "Vishniacozyma", "belongs.to.No.blast.hit", "Mycosphaerella", "Verticillium", "Penicillium", "belongs.to.o.Hypocreales", "Lapidomyces", "belongs.to.o.Pleosporales", "Sporobolomyces", "belongs.to.c.Dothideomycetes", "Arthrocatena", "Lemonniera", "Blumeria", "Taphrina", "Septoriella", "Mycocentrospora", "Rhodotorula", "Fusarium"), rank = "g")


# Close the png graph file handle:
  dev.off() 

  
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#colourCount <- 5
#is.numeric(colourCount) #must be true
  
  
##############
#            #
#    Oom     #
#            #
##############  

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.temporal.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_genus)
View(revamptemp$data_genus)

revampPathMC <- paste(pathDataAnMC, "revamp/", sep = "")

write.table(revamptemp$data_genus, file=paste(revampPathMC, "DP18_datarevamp.Om.genus.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 18
is.numeric(colourCount) #must be true												

#genus
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "OM", ".group.temporal.genus.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("Peronospora", "Hyaloperonospora", "belongs.to.c.uncultured.Oomycetes", "Plasmopara", "Plasmoverna", "belongs.to.No.blast.hit", "Bremia", "Pythium", "Basidiophora", "Phytophthora", "Saprolegnia", "Pseudoperonospora", "uncultured.Phytophthora", "Lagenidium", "Halophytophthora", "Albugo", "Pythiogeton", "Phytopythium"), group = c("Peronospora", "Hyaloperonospora", "belongs.to.c.uncultured.Oomycetes", "Plasmopara", "Plasmoverna", "belongs.to.No.blast.hit", "Bremia", "Pythium", "Basidiophora", "Phytophthora", "Saprolegnia", "Pseudoperonospora", "uncultured.Phytophthora", "Lagenidium", "Halophytophthora", "Albugo", "Pythiogeton", "Phytopythium"), rank = "g")


# Close the png graph file handle:
  dev.off() 

  
  
##############
#            #
#    Oom cur     #
#            #
##############  

tabletemp <- read.table(paste(taxFillPathMC, "OM.table.cur.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_OM.temporal.cur.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_genus)
View(revamptemp$data_genus)

revampPathMC <- paste(pathDataAnMC, "revamp/", sep = "")

write.table(revamptemp$data_genus, file=paste(revampPathMC, "DP18_datarevamp.Om.cur.genus.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 17
is.numeric(colourCount) #must be true												

#genus
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "OM", ".group.temporal.genus.cur.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

  
  my.group.temporal.2(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("Albugo",	"Basidiophora", "belongs_to_c_uncultured.Oomycetes","belongs_to_No.blast.hit",	"Bremia",	"Halophytophthora",	"Hyaloperonospora",	"Lagenidium",	"Peronospora",	"Phytophthora",	"Plasmopara",	"Plasmoverna",	"Pseudoperonospora",	"Pythiogeton",	"Pythium",	"Saprolegnia",	"uncultured.Phytophthora"), group = c("Albugo",	"Basidiophora",	"belongs_to_c_uncultured.Oomycetes",	"belongs_to_No.blast.hit",	"Bremia",	"Halophytophthora",	"Hyaloperonospora",	"Lagenidium",	"Peronospora",	"Phytophthora",	"Plasmopara",	"Plasmoverna",	"Pseudoperonospora",	"Pythiogeton",	"Pythium",	"Saprolegnia",	"uncultured.Phytophthora"), rank = "g")


#Close the png graph file handle:
  dev.off()   

  
  
##############
#            #
#    Phy     #
#            #
##############  

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.temporal.meta.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_genus)
View(revamptemp$data_genus)

revampPathMC <- paste(pathDataAnMC, "revamp/", sep = "")

write.table(revamptemp$data_species, file=paste(revampPathMC, "DP18_datarevamp.Phy.species.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 32
is.numeric(colourCount) #must be true												

																															

#genus
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "Phy", ".group.temporal.species.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("belongs.to.k.No.blast.hit", "VIGNAE", "MEGAKARYA", "IRANICA", "SULAWESIENSIS", "HIBERNALIS", "BAHAMENSIS", "SANSOMEA", "TRIFOLII", "FALLAX", "RASPTASM", "SALIXSOIL", "FOLIORUM", "COLOCASIAE", "RICHARDIAEMACROCHLAMYDOSPORA", "ALNISUBSPMULTIFORMIS", "KERNOVIAE", "PINIFOLIA", "MEADII", "PALMIVORA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "INUNDATA", "TAXRIVERSOIL", "OREGONENSIS", "CLANDESTINA", "ERWINII", "GONAPODYIDES", "MEGASPERMA", "PGCHLAMYDO", "QUERCINA", "TENTACULATA"), group = c("belongs.to.k.No.blast.hit", "VIGNAE", "MEGAKARYA", "IRANICA", "SULAWESIENSIS", "HIBERNALIS", "BAHAMENSIS", "SANSOMEA", "TRIFOLII", "FALLAX", "RASPTASM", "SALIXSOIL", "FOLIORUM", "COLOCASIAE", "RICHARDIAEMACROCHLAMYDOSPORA", "ALNISUBSPMULTIFORMIS", "KERNOVIAE", "PINIFOLIA", "MEADII", "PALMIVORA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "INUNDATA", "TAXRIVERSOIL", "OREGONENSIS", "CLANDESTINA", "ERWINII", "GONAPODYIDES", "MEGASPERMA", "PGCHLAMYDO", "QUERCINA", "TENTACULATA"), rank = "s")


# Close the png graph file handle:
  dev.off() 
  

  
##############
#            #
#    Phy long    #
#            #
##############  

tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.temporal.long.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_genus)

revampPathMC <- paste(pathDataAnMC, "revamp/", sep = "")
write.table(revamptemp$data_species, file=paste(revampPathMC, "DP18_datarevamp.Phy.species.long.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  
  
  
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 30
is.numeric(colourCount) #must be true												

																															
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "Phy", ".group.temporal.species.long.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal.2(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("No_blast_hit", "ALNISUBSPMULTIFORMIS", "BAHAMENSIS", "CLANDESTINA", "COLOCASIAE", "ERWINII", "FALLAX", "FOLIORUM", "GONAPODYIDES", "HIBERNALIS", "INUNDATA", "KERNOVIAE", "MEADII", "MEGAKARYA", "MEGASPERMA", "OREGONENSIS", "PALMIVORA","PGCHLAMYDO"  , "PINIFOLIA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "QUERCINA","RASPTASM", "RICHARDIAEMACROCHLAMYDOSPORA", "SALIXSOIL", "SANSOMEA", "SULAWESIENSIS", "TAXRIVERSOIL", "TENTACULATA", "TRIFOLII"), group = c("No_blast_hit", "ALNISUBSPMULTIFORMIS", "BAHAMENSIS", "CLANDESTINA", "COLOCASIAE", "ERWINII", "FALLAX", "FOLIORUM", "GONAPODYIDES", "HIBERNALIS", "INUNDATA", "KERNOVIAE", "MEADII", "MEGAKARYA", "MEGASPERMA", "OREGONENSIS", "PALMIVORA","PGCHLAMYDO"  , "PINIFOLIA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "QUERCINA","RASPTASM", "RICHARDIAEMACROCHLAMYDOSPORA", "SALIXSOIL", "SANSOMEA", "SULAWESIENSIS", "TAXRIVERSOIL", "TENTACULATA", "TRIFOLII"), rank = "s")


# Close the png graph file handle:
  dev.off() 
    

  
#no unidentified
  
  
  tabletemp <- read.table(paste(taxFillPathMC, "Phy.table.long.taxFill.no.unidentified.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_Phy.meta.temporal.long.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 29
is.numeric(colourCount) #must be true												

																															
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "Phy", ".group.temporal.species.long.no.unidentified.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal.2(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("ALNISUBSPMULTIFORMIS", "BAHAMENSIS", "CLANDESTINA", "COLOCASIAE", "ERWINII", "FALLAX", "FOLIORUM", "GONAPODYIDES", "HIBERNALIS", "INUNDATA", "KERNOVIAE", "MEADII", "MEGASPERMA", "OREGONENSIS", "PALMIVORA","PGCHLAMYDO"  , "PINIFOLIA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "QUERCINA","RASPTASM", "RICHARDIAEMACROCHLAMYDOSPORA", "LACUSTRIS", "SANSOMEA", "TAXRIVERSOIL", "TENTACULATA", "TRIFOLII"), group = c("ALNISUBSPMULTIFORMIS", "BAHAMENSIS", "CLANDESTINA", "COLOCASIAE", "ERWINII", "FALLAX", "FOLIORUM", "GONAPODYIDES", "HIBERNALIS", "INUNDATA", "KERNOVIAE", "MEADII", "MEGASPERMA", "OREGONENSIS", "PALMIVORA","PGCHLAMYDO"  , "PINIFOLIA", "PSEUDOSYRINGAE", "PSYCHROPHILA", "QUERCINA","RASPTASM", "RICHARDIAEMACROCHLAMYDOSPORA", "LACUSTRIS", "SANSOMEA", "TAXRIVERSOIL", "TENTACULATA", "TRIFOLII"), rank = "s")



# Close the png graph file handle:
  dev.off() 

  
##############
#            #
#    plants  #
#            #
##############

tabletemp <- read.table(paste(taxFillPathMC, "DP.plants.table.taxFill.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(pathDataAnMC, "DP_merged_plants.meta.temporal.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, colClasses=c("CollectionDate"="Date"))

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work

######perform revamp before to obtain taxa
revamptemp <- data.revamp(data=list(data=tabletemp), is.OTU=TRUE, ranks = ranks, stand.method = "normalize", top = NULL, mode = "number")

revamptemp
#as.data.frame(revamptemp)
names(revamptemp)
names(revamptemp$data_family)
View(revamptemp$data_genus)


write.table(revamptemp$data_genus, file=paste(revampPathMC, "DP18_datarevamp.DP.plants.genus.tsv", sep = ""),
                               append    = FALSE,
                               sep       = "\t",
                               row.names = TRUE,
                               col.names = NA,
                               quote=FALSE)
  
library("RColorBrewer")
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount <- 25
is.numeric(colourCount) #must be true				


#genus
# Open file-handle to get ready to make the png graph:
  png(filename=paste(temporalPathMC, "DP.plants", ".group.temporal.species.png", sep = ""),
      width = 2200, height = 1800, units = "px", res = 300)

my.group.temporal.2(data=tabletemp, meta=metatemp, date.col= "CollectionDate", factors = c("Anemone",	"belongs_to_p_Streptophyta",	"Betula",	"Brassica",	"Centaurea",	"Crataegus",	"Echium",	"Malus",	"Melilotus",	"Parthenocissus",	"Potentilla",	"Prunus",	"Pyrus",	"Quercus",	"Rhamnus",	"Rhus",	"Sambucus",	"Sorbus",	"Symphytum",	"Syringa",	"Taraxacum",	"Toxicodendron",	"Trifolium",	"Viburnum",	"Zanthoxylum"), group = c("Anemone",	"belongs_to_p_Streptophyta",	"Betula",	"Brassica",	"Centaurea",	"Crataegus",	"Echium",	"Malus",	"Melilotus",	"Parthenocissus",	"Potentilla",	"Prunus",	"Pyrus",	"Quercus",	"Rhamnus",	"Rhus",	"Sambucus",	"Sorbus",	"Symphytum",	"Syringa",	"Taraxacum",	"Toxicodendron",	"Trifolium",	"Viburnum",	"Zanthoxylum"), rank = "g")


# Close the png graph file handle:
  dev.off() 

  
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#colourCount <- 5
#is.numeric(colourCount) #must be true
    
  
    
      
```  
  
 #heatmap
```{r}   
 
heatPathMC <- paste(pathDataAnMC, "heatmap/", sep = "")

##############
#            #
#    ITS     #
#            #
############## works


tabletemp <- read.table(paste(taxFillPathMC, "ITS.table.taxFill.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE, check.names = FALSE, colClasses=c("taxonomy"="character"))
#tabletemp$row.names <- NULL if not using tXfilled tables you need it
 
metatemp <-  read.table(paste(diversityPathMC, "ITS.meta.div.2.tsv", sep = ""),
                      sep = "\t", header = TRUE, dec = ".", comment.char = "", quote = "", stringsAsFactors = TRUE,
                      as.is = TRUE)

metatemp$row.names <- NULL

rownames(metatemp) <- colnames(tabletemp)[-ncol(tabletemp)] #seems to work


#group.heatmap.simple

# Open file-handle to get ready to make the png graph:
  png(filename=paste(heatPathMC, "ITS", ".heatmapgroup..png", sep = ""),
      width = 1200, height = 400, units = "px")

group.heatmap.simple(data=tabletemp, is.OTU=TRUE, meta= metatemp, rank = NULL, row.factor = c("CollectionDate"), top = NULL, drop.unclassified = FALSE, dendro = "row", count = TRUE)

# Close the png graph file handle:
  dev.off()


 
 
``` 
 
 plot with concentric circles
 
```{r}   
 #mapping with concentric circles 
#trying for kinburn ontario with 5km circle  


d <- data.frame(lat = c(45.38319), lon = c(-76.201510))


library(ggplot2)
library(ggmap)

lon <- d$lon
lat <- d$lat

#plot without circle
mapkinb <- get_map(location = c(lon = mean(d$lon), lat = mean(d$lat)), zoom = 12, 
                     maptype = "terrain", color = c("bw"))

# Open file-handle to get ready to make the png graph:
png(filename=paste(pathDataAnMC, "Bees", ".geo.map.png", sep = ""),
    width = 3000, height = 3500, units = "px", pointsize = 6)

ggmap(mapkinb) +
  geom_point(data = d, aes(x = lon, y = lat, fill = "black", alpha = 2), size = 9, shape = 19) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)

# Close the png graph file handle:
dev.off()  


####add cicrle 

library(ggplot2)
library(ggmap)
data = data.frame(
  ID = as.numeric(c(1)),
  longitude = as.numeric(c(-76.201510)),
  latitude = as.numeric(c(45.383139))
)


#################################################################################
# create circles data frame from the centers data frame
make_circles <- function(centers, radius, nPoints = 100){
  # centers: the data frame of centers with ID
  # radius: radius measured in kilometer
  #
  meanLat <- mean(centers$latitude)
  # length per longitude changes with lattitude, so need correction
  radiusLon <- radius /111 / cos(meanLat/57.3) 
  radiusLat <- radius / 111
  circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
  angle <- seq(0,2*pi,length.out = nPoints)
  
  circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radiusLon * cos(angle)))
  circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radiusLat * sin(angle)))
  return(circleDF)
}

# here is the data frame for all circles
myCircles <- make_circles(data, 5)



##################################################################################
# https://stackoverflow.com/questions/34183049/plot-circle-with-a-certain-radius-around-point-on-a-map-in-ggplot2
#https://gis.stackexchange.com/questions/119736/ggmap-create-circle-symbol-where-radius-represents-distance-miles-or-km/119858
#https://www.rdocumentation.org/packages/ggforce/versions/0.1.1/topics/geom_circle
#http://ggplot2.tidyverse.org/reference/geom_point.html

#data of the dude
island = get_map(location = c(lon = -63.247593, lat = 17.631598), zoom = 13, maptype = "satellite")
islandMap = ggmap(island, extent = "panel", legend = "bottomright")
RL = geom_point(aes(x = longitude, y = latitude), data = data, color = "#ff0000")
islandMap + RL + 
  scale_x_continuous(limits = c(-63.280, -63.21), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(17.605, 17.66), expand = c(0, 0)) +
  ########### add circles
  geom_polygon(data = myCircles, aes(lon, lat, group = ID), color = "red", alpha = 0)

#################
#my data

#get map
mapkinb <- get_map(location = c(lon = d$lon, lat = d$lat), zoom = 11, 
                   maptype = "terrain", color = c("bw"))

#plot and have circle
# Open file-handle to get ready to make the png graph:
png(filename=paste(pathDataAnMC, "Bees", ".geo.map.circled.png", sep = ""),
    width = 3000, height = 3500, units = "px", pointsize = 6)

ggmap(mapkinb) +
  geom_point(data = d, aes(x = lon, y = lat, fill = "black", alpha = 2), size = 15, shape = 19) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE) + geom_polygon(data = myCircles, aes(lon, lat, group = ID), color = "grey10", alpha = 0.1, size = 3)

# Close the png graph file handle:
dev.off()  

  
``` 
 species query with metaresultextractor
 
 
```{r}  

---
  title: "species.query.marie-claude"
author: "E.T."
date: "5/6/2018"


#######################


#########################
#                       #
#    ITS by species     #
#                       #
#########################

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "genusQuery.txt", sep = "")
otu_in <- paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = "")
meta_in <- paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = "")
fasta_in <- paste(pathDataAnMC, "ref.seqs.ITS.DP.fasta", sep = "")
table_taxon_out <- paste(WordMatchDir, "ITS_species_query_by_taxon.MC.pollen.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "ITS_species_query_by_sample.MC.pollen.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/ITS_species_query_MC_pollen", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)

#then, use the cluster script
#########################
#                       #
#    ITS by genuses     #
#                       #
#########################

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "genusQuery.txt", sep = "")
otu_in <- paste(pathDataAn, "ITS.table.tsv", sep = "")
meta_in <- paste(pathDataAn, "ITS.meta.tsv", sep = "")
fasta_in <- paste(pathDataAn, "new_refseqs_ITS.fna", sep = "")
table_taxon_out <- paste(WordMatchDir, "ITS_genus_query_by_taxon2.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "ITS_genus_query_by_sample2.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/ITS_genus_query", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)

####
#more species/genera to the list

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "genus.and.species.quesry.pollen.txt", sep = "")
otu_in <- paste(pathDataAnMC, "DP_merged_ITS.table.tsv", sep = "")
meta_in <- paste(pathDataAnMC, "DP_merged_ITS.meta.tsv", sep = "")
fasta_in <- paste(pathDataAnMC, "ref.seqs.ITS.DP.fasta", sep = "")
table_taxon_out <- paste(WordMatchDir, "ITS_species_ag.query_by_taxon.MC.pollen.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "ITS_species_ag.query_by_sample.MC.pollen.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/ITS_species_query_MC_pollen", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)


#oomycetes
#more species/genera to the list

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "oomycete.query.pollen.txt", sep = "")
otu_in <- paste(pathDataAnMC, "DP_merged_OM.table.tsv", sep = "")
meta_in <- paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = "")
fasta_in <- paste(pathDataAnMC, "OM.ref.seqs.DP.fasta", sep = "")
table_taxon_out <- paste(WordMatchDir, "omycete_query_by_taxon.MC.pollen.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "oomycetes_query_by_sample.MC.pollen.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/ITS_species_query_MC_pollen", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)

###oom curated


#oomycetes
#more species/genera to the list

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "oomycete.query.pollen.txt", sep = "")
otu_in <- paste(pathDataAnMC, "DP_merged_OM.curated.tsv", sep = "")
meta_in <- paste(pathDataAnMC, "DP_merged_OM.meta.tsv", sep = "")
fasta_in <- paste(pathDataAnMC, "OM.ref.seqs.DP.cur.fasta", sep = "")
table_taxon_out <- paste(WordMatchDir, "omycete_query_by_taxon.MC.pollen.cur.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "oomycetes_query_by_sample.MC.pollen.cur.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/ITS_species_query_MC_pollen", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)


#plants
#oomycetes
#more species/genera to the list

#Define your 6 arguments to run the perl script
query_in <- paste(WordMatchDir, "plants.query.pollen.txt", sep = "")
otu_in <- paste(pathDataAnMC, "DP_merged_plants.table.tsv", sep = "")
meta_in <- paste(pathDataAnMC, "DP_merged_plants.meta.tsv", sep = "")
fasta_in <- paste(pathDataAnMC, "Plants.ref.seqs.DP.fasta", sep = "")
table_taxon_out <- paste(WordMatchDir, "plants_query_by_taxon.MC.pollen.tsv", sep = "")
table_sample_out <- paste(WordMatchDir, "plants_query_by_sample.MC.pollen.tsv", sep = "")
fasta_out_folder <- paste(WordMatchDir, "/FastaOut/plants_species_query_MC_pollen", sep = "")

cmd <- paste("perl", paste(ScriptsPath, "metaResultExtractorGOOD.pl", sep = ""),
             query_in,
             otu_in,
             meta_in,
             fasta_in,
             table_taxon_out,
             table_sample_out,
             fasta_out_folder,
             sep = " ")
system(cmd)

