## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(nMyo)

## ----loading------------------------------------------------------------------
data(nMyo_Data)
data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)

# listed components
names(data)

# the CorrCounts
data$Counts[1:5, 1:5]

# the design_table
data$Design[1:5, ]
            
# the DE
data$DEstats[1:5, ]

# the data stamp (for file storage)
data$Date_Stamp

## ----query1-------------------------------------------------------------------
data_Postn <- MarkerQuery(Data = data, marker = "Postn")

## ----query2-------------------------------------------------------------------
data_Postn_alt <- MarkerQuery(Data = data, marker = "ENSMUSG00000027750")

## ----query3-------------------------------------------------------------------
data_Postn_alt2 <- MarkerQuery(Data = data, marker = "ENSMUSG00000027750:Postn")

## ----query4-------------------------------------------------------------------
data_Postn_err <- MarkerQuery(Data = data, marker = "Post")

## ----query5-------------------------------------------------------------------
data_err <- MarkerQuery(Data = data, marker = "Tp53")

## ----tsne---------------------------------------------------------------------
scatter_Postn <- doScatterDR(Data = data_Postn,
                             highlight.by = "CellType",
                             show.plot = TRUE,
                             save.plot = FALSE)

## ----violin-------------------------------------------------------------------
violin_Postn <- doViolin(Data = data_Postn,
                         grouping.by = "CellType",
                         grouping2.by = "Condition",
                         show.plot = TRUE,
                         save.plot = FALSE) 

## ----violin2------------------------------------------------------------------
violin_Postn2 <- doViolin(Data = data_Postn,
                         grouping.by = "CellType",
                         grouping2.by = "",
                         show.plot = TRUE,
                         save.plot = FALSE) 

## ----volcano------------------------------------------------------------------
volcano_Postn <- doVolcano(Data = data_Postn,
                           filters = "Comparison=='Fibroblast: Sham-MI'",
                           logFC = 1,
                           FDR = 0.01,
                           show.plot = TRUE,
                           save.plot = FALSE)

## ----pairwise-----------------------------------------------------------------
# all pairwise comparisons (to be used in filter parameter)
names(table(data$DEstats$Comparison))

## ----ma-----------------------------------------------------------------------
ma_Postn <- doMA(Data = data_Postn,
                 filters = "Comparison=='Fibroblast: Sham-MI'",
                 logFC = 1,
                 FDR = 0.01,
                 show.plot = TRUE,
                 save.plot = FALSE)

## ----run_app, eval=FALSE------------------------------------------------------
#  run_nMyo()

## ----pressure, echo=FALSE, fig.cap="nMyo welcome screen", out.width = '100%'----
knitr::include_graphics("screen1.png")

## ----pressure2, echo=FALSE, fig.cap="nMyo welcome screen", out.width = '100%'----
knitr::include_graphics("screen2.png")

