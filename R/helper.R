
#' Color specifications
#'
#' color_hue
#' @param n numeric. The number of colors to be used from the pallete.
#' @keywords gg_color_hue
#' @return A vector of colors
#'
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Stop the analysis if the files are not provided or the filenames are wrong
#'
#' Helper function
#' @param Counts_file character. The location of the normalised marker expression file.
#' @param Design_file character. The location of the experimental design file.
#' @keywords firstStop
#' @return A stop or nothing
#'
firstStop<-function(Counts_file,Design_file){
  # if the necessary files are missing...
  if(is.null(Counts_file)){
    stop("
          ***** A Counts file is required in the application. *****
          ")
  }
  if(is.null(Design_file)){
    stop("
          ***** A Design file is required in the application. *****
          ")
  }
  if(!is.null(Counts_file)){
    if(!file.exists(Counts_file)){
      stop("
              ***** The Counts file is not found in the system. Check the provided name. *****
              ")
    }
  }
  if(!is.null(Design_file)){
    if(!file.exists(Design_file)){
      stop("
              ***** The Design file is not found in the system. Check the provided name. *****
              ")
    }
  }
}

#' Read the required data files
#'
#' Helper function
#' @param data_file character. The data file containing the normalised and corrected counts (CorrCounts), the design table (DesignTable) and the differential expressin analysis (DEs).
#' @keywords RequiredDataReader
#' @return A list of data
#' @import data.table stringr
#'
RequiredDataReader<-function(data_file){
  counts<-data_file$CorrCounts
  rownames(counts)<-counts[,1]
  cc<-unique(apply(matrix(as.character(rownames(counts)),ncol=1),1,str_count,":"))
  if(length(cc)>1 | cc[1]>1){
    stop("
               ***** EachThe marker ID of the Counts_file should contain only one : character. *****
               ")
  }
  counts<-counts[,-1]

  desi<-data_file$DesignTable
  rownames(desi)<-desi[,1]
  mm<-match(colnames(counts),rownames(desi),nomatch=0)
  if(length(mm[mm==0])>0){
    stop("
               ***** The SampleIDs of the Design file (rows) do not match exactly the SampleIDs of the Counts file (columns). *****
               ")
  }
  desi<-desi[mm,]

  annot<-matrix(rownames(counts),ncol=1)
  annot<-data.frame(cbind(annot,GeneNames(annot,what="Other"),GeneNames(annot,what="Name")))
  rownames(annot)<-annot[,1]
  colnames(annot)<-c("ID","ID1","ID2")

  # remove unnecessary characters from names
  colnames(counts)<-str_replace(colnames(counts),"_","")
  rownames(desi)<-str_replace(rownames(desi),"_","")
  desi$Condition<-str_replace(str_replace(desi$Condition,"-","_"),"&","_")
  desi$CellType<-str_replace(str_replace(desi$CellType,"-","_"),"&","_")

  de<-data_file$DEs
  data<-list(Counts=counts,Design=desi,DEstats=de,Annotation=annot)


  return(data)
}


#' Read the DE_file and tests the data validity
#'
#' Helper function
#' @param data list. The outcome of the RequiredDataReader() internal function.
#' @keywords testDataValidity
#' @return A list of data
#' @import data.table stringr
#'
testDataValidity<-function(data){

      data<-validity(data)

      uu<-unique(as.character(data$DEstats$CellType))
      if(length(grep("-",uu))>0){
        data$DEstats$CellType<-str_replace(data$DEstats$CellType,"[-|_@;/^:]","&")
      }
      dec<-unique(as.character(data$DEstats$Comparison))
      strc<-apply(matrix(dec,ncol=1),1,str_count,"-")
      if(length(which(strc==2))>0 | length(which(strc==0))>0){
        stop("
                       ***** Each entry of DE_file's Comparison column should contain only one - character. *****
                       ")
      }
      dec<-matrix(t(apply(t(matrix(unlist(strsplit(dec,"-")),nrow=2)),1,sort)),ncol=2)
      if(nrow(unique(dec))!=nrow(dec)){
        stop("
             ***** Reversed comparisons are not allowed in the DE file: <level1> - <level2> and <level2> - <level1> should be joined! *****")
      }


 return(data)
}


#' Test the data validity
#'
#' A helper (internal) function that reports an error if the data do not have the required format
#' @param data list. It must include four components: component 'Counts' includes the normalised count data, component 'Design' includes the experimental design information, component 'Annotation' include the marker annotation (e.g. gene names and Ensembl IDs) and component 'DEstats' includes the differential expression analysis results.
#' @keywords validity
#'
validity<-function(data){

  if(!identical(colnames(data$Counts),rownames(data$Design))){
    stop("
        ***** The sample names of the raw count and the design matrices differ! *****
        ")
  }

  g<-match(paste("Dim",1:3,sep=""),colnames(data$Design),nomatch=0)
  if(length(g[g>0])<=1){
    stop("
        ***** There are no dimensionality reduction components in the data. *****
        ***** Add the coordinates of a tSNE or UMAP or PCA with column names Dim1, Dim2 and Dim3 (if exists)! *****
        ")
  }

  if(!is.null(data$DEstats)){
    g1<-match(c("Marker","logFC","logCPM","PValue","FDR","Comparison","CellType"),colnames(data$DEstats),nomatch=0)
    if(length(g1[g1==0])>0){
      stop("
            ***** The DE statistics table is incomplete or the column names are not well-defined! *****
            ")
    }

    # appropriateness of DE comparison and cell type
    comps<-unique(as.character(data$DEstats$Comparison))
    g2<-grep("-",comps)
    if(length(g2)<length(comps)){
      stop("
            ***** Modify the data$DEstats$Comparison levels as <condition1> - <condition2> (omit the space). *****
            ")
    }

    cts<-unique(as.character(data$DEstats$CellType))
    g3<-grep("-",cts)
    if(length(g3)>0){
      stop("
            ***** The data$DEstats$CellType should not contain the - symbol. *****
            ***** Use a single or two words separated by & matching the Design Condition or CellType entries *****
            ")
    }

    # commons between Design and DE
    desiCellType<-unique(as.character(data$Design$CellType))
    desiCondition<-unique(as.character(data$Design$Condition))
    deCellType<-unique(as.character(data$DEstats$CellType))
    deCondition<-unique(as.character(data$DEstats$Comparison))

    mmctype1<-matchDEentries(deEntry=deCellType,the_other=desiCellType,the_character="&")
    mmctype2<-matchDEentries(deEntry=deCellType,the_other=desiCondition,the_character="&")
    mmctype<-unique(c(mmctype1,mmctype2))

    mmcond1<-matchDEentries(deEntry=deCondition,the_other=desiCondition,the_character="-")
    mmcond2<-matchDEentries(deEntry=deCondition,the_other=desiCellType,the_character="-")
    mmcond<-unique(c(mmcond1,mmcond2))

    if(length(mmctype)==0){
      stop("
            ***** There are no common cell types among those of the Design file and the DE file. *****
            ***** Check the data or consider removing the DE file from the analysis. *****
            ")
    }
    if(length(mmcond)==0){
      stop("
            ***** There are no common conditions among those of the Design file and the DE file (in Comparison column). *****
            ***** Check the data or consider removing the DE file from the analysis. *****
            ")
    }

  }

  g2<-match(c("SampleID","Condition","CellType"),colnames(data$Design),nomatch=0)
  if(length(g2[g2==0])>0){
    stop("
        ***** The information of the design table is incomplete! *****
        ***** At least one of the required column names (SampleID, Condition and CellType) is missing. *****
        ")
  }

  return(data)
}


#' Read the Markers file
#'
#' Helper function
#' @param Markers_file character. The location of the markers file (it can be NULL).
#' @keywords ReadMarkers
#' @return A vector of gene names or NULL
#'
ReadMarkers<-function(Markers_file){
  # if the Genes file is missing, set the data to an empty vector
  if(!is.null(Markers_file)){
    #status<-"Marker file"

    if(file.exists(Markers_file)){
      genes<-scan(Markers_file,what="character")
    } else {
      stop("
              ***** The Markers file is not found in the system. Check the provided name. *****
              ")
    }
  } else {
    #status<-"DE file"
    genes<-NULL
  }
 return(genes)
}



#' Read the DE_file and tests the data validity
#'
#' Helper function
#' @param data list. The outcome of the RequiredDataReader() internal function updated with other components (is.Exact, logFC_cutoff, FDR_cutoff, status).
#' @param genes character. A list of genes.
#' @keywords finaliseMarkerList
#' @return A list with the significant DE genes, the final list of genes and the missing genes
#'
finaliseMarkerList<-function(data,genes){

  missing<-NULL
  sigDE<-NULL
  if(is.null(genes)){
    if(is.null(data$DEstats)){
      data$Status<-"Manual input"
      genes<-rownames(data$Counts)
    } else {

      sigDE<-subset(data$DEstats,abs(logFC)>=data$logFC & FDR<=data$FDR)

      if(nrow(sigDE)>0){
        genes<-unique(as.character(sigDE[,1]))
        mm<-findGeneType(templ=data$Annotation,genes=genes,multiple=FALSE)
        genes<-rownames(data$Annotation)[mm]
        print("***** The data have been successfully loaded! *****")
      } else {
        stop("
                ***** The analysis returned 0 differentially expressed markers at the specified cutoffs. Relax the logFC and FDR criteria. *****
                ")
      }
    }
  } else {

    if(!data$Exact_Marker_Match){
      newgenes<-NULL
      gg<-as.list(rep(0,length(genes)))
      names(gg)<-genes
      for(i in 1:length(genes)){
        gg[[i]]<-grep(genes[i],rownames(data$Annotation))
        newgenes<-c(newgenes,rownames(data$Annotation)[gg[[i]]])
      }
      newgenes<-unique(newgenes)
      missing<-NULL
      ll<-lapply(gg,length)
      if(length(ll[ll==0])>0){
        missing<-names(gg)[ll==0]
      }
      genes<-newgenes

    } else {
      mm<-findGeneType(templ=data$Annotation,genes=genes,multiple=FALSE)
      newgenes<-rownames(data$Annotation)[mm]
      missing<-NULL
      if(length(mm[mm==0])>0){
        missing<-genes[mm==0]
      }

      genes<-newgenes
    }

    # get the DE genes when Genes file exists
    if(!is.null(data$DEstats)){
      sigDE<-subset(data$DEstats,abs(logFC)>=data$logFC & FDR<=data$FDR)
      x<-as.character(sigDE[,1])
      x<-cbind(x,GeneNames(x,what="Other"),GeneNames(x,what="Name"))
      mm<-findGeneType(templ=x,genes=genes,multiple=TRUE)
      sigDE<-sigDE[mm>0,]

      if(nrow(sigDE)==0){
        print("***** None of the markers was found to be differentially expressed in the DE_file *****")
      }
    }
  }
  return(list(significantDE=sigDE,Markers=genes,MissingMarkers=missing))
}


#' It creates the filteredData slot
#'
#' Helper function
#' @param data list. The outcome of the RequiredDataReader() internal function updated with other components (is.Exact, logFC_cutoff, FDR_cutoff, status,Output_folder and Dimensions).
#' @param sigDE.genes.status list. The outcome of finaliseMarkerList().
#' @keywords fixFilteredData
#' @return A list with the filteredData
#'
fixFilteredData<-function(data,sigDE.genes.status){

  redData<-data[which(names(data)=="Counts" | names(data)=="Design" | names(data)=="Annotation" | names(data)=="DEstats")]

  if(!is.null(sigDE.genes.status$Markers)){
    genes_index<-match(as.character(sigDE.genes.status$Markers),rownames(data$Counts))
    redData$Counts<-redData$Counts[genes_index,]
    redData$Annotation<-redData$Annotation[genes_index,]
  } else {
    redData$Counts<-data$Counts
    redData$Annotation<-data$Annotation
  }

  if(!is.null(data$DEstats)){
    redData$DEstats<-sigDE.genes.status$significantDE
  } else {
    redData$DEstats<-NULL
  }
  redData$Marker_List<-sigDE.genes.status$Markers
 return(redData)
}

#' GeneNames
#'
#' Extract gene IDs
#'
#' @param data character vector. A vector of IDs. Each entry may contain two ID separated by a symbol e.g. EnsemblID:MarkerName.
#' @param symbol character. The separator used. Default is the symbol ':'.
#' @param what character. Which of the above two IDs (if two) to keep. Options are 'Name' (default) to keep the marker name or 'Other' to keep the other ID. If each entry consists of one ID then this is the one kept.
#' @keywords GeneNames
#' @return A vector of IDs
#' @export
#'
GeneNames<-function(data,symbol=":",what="Name"){
    ss<-strsplit(as.character(data), symbol)
    k<-ifelse(what=="Name",2,1)
    res<-rep(0,length(ss))
    for(i in 1:length(ss)){
        res[i]<-ss[[i]][min(k,length(ss[[i]]))]
    }
    return(res)
}

#' Find the marker of interest in the data
#'
#' Helper function
#' @param templ matrix. It contains the marker annotations.
#' @param genes character. The IDs of the markers of interest.
#' @param correction logical. Whether the gene names will be automatically edited to a format that allows a case-insensitive marker search. Default is FALSE
#' @param multiple logical. Whether the temp contains repeats of the same gene name. Default is FALSE.
#' @keywords findGeneType
#' @return A vector of indices
#'
findGeneType<-function(templ,genes,correction=FALSE,multiple=FALSE){

    if(length(grep("ENSMUSG",toupper(genes)))==0){
        if(correction){
            s <- strsplit(tolower(genes), " ")[[1]]
            genes<-paste(toupper(substring(s, 1,1)), substring(s, 2),sep="", collapse=" ")
        }
    } else {
        if(correction){
            genes<-unlist(strsplit(toupper(genes),":"))[1]
        }
    }
    if(!multiple){
        mm<-matrix(0,4,length(genes))
        for(i in 1:3){
            mm[i,]<-match(genes,as.character(templ[,i]),nomatch=0)
        }
        for(i in 1:ncol(mm)){
            w<-which(mm[,i]>0)
            if(length(w)>0){
                mm[4,i]<-mm[w[1],i]
            } else {
                mm[4,i]<-0
            }
        }
    } else {
       mm<-matrix(0,4,nrow(templ))
       for(i in 1:3){
           mm[i,]<-match(as.character(templ[,i]),genes,nomatch=0)
       }
       for(i in 1:ncol(mm)){
           w<-which(mm[,i]>0)
           if(length(w)>0){
               mm[4,i]<-mm[w[1],i]
           } else {
               mm[4,i]<-0
           }
       }
    }
  return(mm[4,])
}


#' Modify the data if the subsequent steps of DIMER have been used
#'
#' Helper function
#' @param Data list. The outcome of readData().
#' @param func character. Which function is used.
#' @keywords Sample.DataTester
#' @return The Data list after the modification (removal of some slots)
#'
DataTester<-function(Data,func){

    if(func=="SampleQuery_types"){
        w<-which(names(Data)== "selectedData" |
            names(Data)=="Date_Stamp" |
            names(Data)=="DRScatterplot" |
            names(Data)=="Violin" |
            names(Data)=="MA" |
            names(Data)=="Volcano" |
            names(Data)=="Heatmap" |
            names(Data)=="TCMA" |
            names(Data)=="TCVolcano" |
            names(Data)=="Summary_Stats" |
            names(Data)=="Marker_Data" |
            names(Data)=="DE_Stats")

        Data$filteredData<-Data$temp
    }

    if(func=="MarkerQuery_types" | func=="DataIntegration_types"){
        w<-which(names(Data)== "selectedData" |
            names(Data)=="DRScatterplot" |
            names(Data)=="Violin" |
            names(Data)=="MA" |
            names(Data)=="Volcano" |
            names(Data)=="Heatmap" |
            names(Data)=="TCMA" |
            names(Data)=="TCVolcano" |
            names(Data)=="Summary_Stats" |
            names(Data)=="Marker_Data" |
            names(Data)=="DE_Stats")
    }

    if(func!="MarkerQuery_types" &
        func!="DataIntegration_types" &
        func!="SampleQuery_types"){

        w<-which(names(Data)==func |
                 names(Data)=="Summary_Stats" |
                 names(Data)=="Marker_Data" |
                 names(Data)=="DE_Stats")
    }


  if(length(w)>0){
    Data<-Data[-w]
  }

 return(Data)
}


#' Filter the Samples of the Design file based on the filters
#'
#' Helper function
#' @param designData data.frame. The data of the Design_file.
#' @param filters expression. The filters.
#' @keywords filterSampleData
#' @return The filtered data of the Design_file
#' @import dplyr
#'
filterSampleData<-function(designData,filters){
  res<-tryCatch({
    designData %>% filter(eval(filters[1]))},
    error=function(x){
      return(NULL)
    })
  if(is.null(res)){
    stop("*****
        One or more filter variables does not exist in the Design file. *****",sep="")
  }

  if(length(filters)>1){
    for(i in 2:length(filters)){

      res1<-tryCatch({
        designData %>% filter(eval(filters[i]))},
        error=function(x){
          return(NULL)
        })
      if(is.null(res1)){
        stop("*****
          One or more filter variables does not exist in the Design file. *****",sep="")
      }

      res<-rbind(res,res1)
    }
  }
 return(res)
}


#' Filter the Samples of the DE file based on the filters
#'
#' Helper function
#' @param Data list. The outcome of the readData().
#' @keywords filterSampleDEstatsData
#' @return The filtered data of the DE_file
#'
filterSampleDEstatsData<-function(Data){

  desiCellType<-unique(as.character(Data$filteredData$Design$CellType))
  desiCondition<-unique(as.character(Data$filteredData$Design$Condition))

  # get the cell types of the DEstats
  mm1<-matchDEentries(deEntry=as.character(Data$filteredData$DEstats$CellType),the_other=desiCellType,the_character="&")
  mm2<-matchDEentries(deEntry=as.character(Data$filteredData$DEstats$CellType),the_other=desiCondition,the_character="&")
  mm<-unique(c(mm1,mm2))
  if(length(mm)==0){
    print(paste("***** The DE file does not contain any of the kept cell types and it will be removed from the analysis. *****",sep=""))
    Data$filteredData$DEstats<-NULL
  } else {
    mmm<-match(as.character(Data$filteredData$DEstats$CellType),unique(as.character(Data$filteredData$DEstats$CellType))[mm],nomatch=0)
    Data$filteredData$DEstats<-Data$filteredData$DEstats[mmm>0,]
  }
  # get the conditions of the DEstats
  mm1<-matchDEentries(deEntry=as.character(Data$filteredData$DEstats$Comparison),the_other=desiCondition,the_character="-")
  mm2<-matchDEentries(deEntry=as.character(Data$filteredData$DEstats$Comparison),the_other=desiCellType,the_character="-")
  mm<-unique(c(mm1,mm2))
  if(length(mm)==0){
    print(paste("***** The DE file does not contain any of the kept conditions and it will be removed from the analysis. *****",sep=""))
    Data$filteredData$DEstats<-NULL
  } else {
    mmm<-match(as.character(Data$filteredData$DEstats$Comparison),unique(as.character(Data$filteredData$DEstats$Comparison))[mm],nomatch=0)
    Data$filteredData$DEstats<-Data$filteredData$DEstats[mmm>0,]
  }

 return(Data)
}

#' Create a date stamp
#'
#' Helper function
#' @param date character. The date.
#' @keywords DateStamp
#' @return A date to use in the filename
#'
DateStamp<-function(date){
  setDate<-unlist(strsplit(date," ",fixed=T))
  setDate<-setDate[c(1,2,3,5,4)]
  setDate[5]<-chartr(old = ":", new = ".", setDate[5])
  setDate<-paste(setDate,collapse="_")
 return(setDate)
}


#' Check the filtering variables
#'
#' Helper function
#' @param filters character. The filters expression.
#' @param dataMatrix data.frame. The data.frame whose column names the filters refer to.
#' @param mode character. The function that the filter is applied to.
#' @keywords testFilters
#' @return Filters and possible mistakes in the filter names
#'
testFilters<-function(filters,dataMatrix,mode){
  f1<-testFiltersFactor(filters=filters,dataMatrix=dataMatrix,mode=mode)
  f2<-testFiltersNumeric(filters=filters,dataMatrix=dataMatrix,mode=mode)
  f<-list(CellType=f1[[1]],Comparison=f1[[2]],FactorFilters=f1[[3]],NumericFilters=f2)
  if(length(unlist(lapply(f,is.null))==FALSE)==0){
    stop("
         ***** The filtering conditions for marker sorting are not well-defined. *****
         ")
  }
  return(f)
}


#' Check the factor filtering variables
#'
#' Helper function
#' @param filters character. The filters expression.
#' @param dataMatrix data.frame. The data.frame whose column names the filters refer to.
#' @param mode character. The function that the filter is applied to.
#' @keywords testFiltersFactor
#' @return Filters and possible mistakes in the filter names
#'
testFiltersFactor<-function(filters,dataMatrix,mode){
  f<-unlist(strsplit(unlist(strsplit(as.character(filters),"&")),"|",fixed=T))
  chars<-"=="
  mistakes<-c()
  nonmistakes<-c()
  g<-grep(chars,f)
  if(length(g)>0){
    f1<-f[g]
    f1<-t(matrix(unlist(strsplit(as.character(f1),chars)),nrow=2))
    f1[,1]<-gsub(" ", "", f1[,1], fixed = TRUE)
    f1[,2]<-gsub("\"", "", f1[,2], fixed = TRUE)
    f1[,2]<-gsub(" ", "", f1[,2], fixed = TRUE)

    mm<-match(f1[,1],colnames(dataMatrix),nomatch=0)
    if(length(mm[mm==0])>0){
      stop("
          ***** Some factor filters are not well-defined. Check the filters parameter! *****
          ")
    } else {
      w1<-which(f1[,1]=="CellType")
      ct<-unique(f1[w1,2])
      w2<-which(f1[,1]=="Comparison")
      co<-unique(f1[w2,2])
    }
  } else {
    f1<-c()
    ct<-c()
    co<-c()
    if(mode=="SampleSorting"){
      stop("
             ***** The factor filtering conditions for marker sorting are not well-defined. *****
             ")
    }
  }
  return(list(ct,co,f1))
}


#' Check the numeric filtering variables
#'
#' Helper function
#' @param filters character. The filters expression.
#' @param dataMatrix data.frame. The data.frame whose column names the filters refer to.
#' @param mode character. The function that the filter is applied to.
#' @keywords testFiltersNumeric
#' @return Filters and possible mistakes in the filter names
#'
testFiltersNumeric<-function(filters,dataMatrix,mode){
  f<-unlist(strsplit(unlist(strsplit(as.character(filters),"&")),"|",fixed=T))
  chars<-c(">=","<=",">","<")
  mistakes<-c()
  nonmistakes<-c()
  res<-matrix(0,1,2)
  for(i in 1:length(chars)){
    g<-grep(chars[i],f)
    if(length(g)>0){
      f1<-f[g]
      f1<-t(matrix(unlist(strsplit(as.character(f1),chars[i])),nrow=2))
      f1[,1]<-gsub(" ", "", f1[,1], fixed = TRUE)
      f1[,2]<-gsub("\"", "", f1[,2], fixed = TRUE)
      f1[,2]<-gsub(" ", "", f1[,2], fixed = TRUE)
      isN<-suppressWarnings(apply(matrix(as.numeric(as.character(f1[,2])),ncol=1),1,is.na))
      f1<-matrix(f1[which(!isN),],ncol=2)
      mm<-match(f1[,1],colnames(dataMatrix),nomatch=0)
      if(length(mm[mm==0])>0){
        stop("
              ***** Some numeric filters are not well-defined. Check the filters parameter! *****
              ")
      }
      res<-rbind(res,f1)
    }

  }
  res<-matrix(unique(res[-1,]),ncol=2)

  w1<-which(res[,1]=="diffFC")
  w2<-which(res[,1]=="diffCPM")
  w3<-which(res[,1]=="diffPV")

  if(mode=="SampleTime" & length(w1)==0 & length(w2)==0 & length(w3)==0){
    stop("
         ***** The numeric filtering conditions for marker timecourse are not well-defined. *****")
  }

  return(res)
}


#' Find the marker of interest in the differential expression results
#'
#' Helper function
#' @param deEntry character. A marker characteristics from the DE file.
#' @param the_other character. A marker characteristic from the design file.
#' @param the_character character. A character to split the marker characteristics from the DE file.
#' @keywords matchDEentries
#' @return A vector of indices
#'
matchDEentries<-function(deEntry,the_other,the_character){
    s<-strsplit(unique(as.character(deEntry)),the_character)
    if(length(s[[1]])==1){
        res<-matrix(unlist(s),ncol=1)
    }
    if(length(s[[1]])>1){
        res<-t(matrix(unlist(s),nrow=length(s[[1]])))
    }

    mm<-res
    for(i in 1:ncol(res)){
        mm[,i]<-match(res[,i],the_other,nomatch=0)
    }
    mm<-apply(mm,1,function(x) length(x[x>0]))
    mm<-which(mm>0)

  return(mm)
}

#' Find the marker of interest in the differential expression results
#'
#' Helper function
#' @param deData matrix. The differential expression analysis results.
#' @param filters character. The filtering variables.
#' @keywords DEplot_entries
#' @return A vector of indices
#' @import dplyr
#'
DEplot_entries<-function(deData,filters){

    if(is.null(filters)){
      stop("
           ***** Provide the filtering variable to select the Volcano / MA plotting data. ***** ")
    }

    if(length(filters)==0){
      stop("
             ***** Provide the filtering variable to select the Volcano / MA plotting data. ***** ")
    }

    if(length(filters)>1){
      stop("
               ***** The filtering variable of the Volcano / MA plot should be a single expression vector *****
               ***** of the form CellType==<level1> & Comparison==<level2>. *****
           ")
    }
    t<-testFilters(filters=filters,dataMatrix=deData,mode="DE")
    #if(length(t$CellType)==0 | length(t$Comparison)==0){
    #  stop("
    #       ***** The filtering variable of the Volcano / MA plot must specify the CellType and Comparison levels to be plotted. *****
    #       ")
    #}


    res<-tryCatch({
      deData<-deData %>% filter(eval(filters[1]))},
      error=function(x){
        return(NULL)
      })
    if(is.null(deData)){
      stop("*****
               One or more filter variables does not exist in the file. *****",sep="")
    }

    if(nrow(deData)==0){
        stop("
        ***** The DE filtering variable did not return any data for the Volcano / MA plot. *****
        ***** Make sure that the de.filters are correctly specified and match to SampleQuery() selection! *****
        ")
    }

 return(list(deData,t))
}

#' Join two strings
#'
#' Helper function
#' @param x string. The first string.
#' @param y string. The second string.
#' @keywords joinstring
#' @return A joined string
#'
'%&%' <- function(x, y)paste0(x,y)


#' Assign colors to the scatterplot
#'
#' Helper function
#' @param Data list. The outcome of the MarkerQuery().
#' @param plotData data.frame. The data to use for plotting.
#' @param feature character. The feature to be plotted.
#' @keywords makecolors_scatterplot
#' @return A list
#' @import RColorBrewer viridis viridisLite
#'
makecolors_scatterplot<-function(Data,plotData,feature){

    type<-c()
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    filteredData<-Data$selectedData$Design

    if(length(which(colnames(filteredData)==feature))==0){
        stop("
        ***** This feature name does not in the Design file! *****
        ")
    }

    if(which(colnames(filteredData)==feature)==ncol(filteredData)){
        my_palette <- colorRampPalette(c("#00aedb","#ff4e50"))(n = 50)
        column<-as.numeric(as.character(Data$Counts[match(feature,rownames(Data$Counts)),match(rownames(plotData),colnames(Data$Counts))]))
        type<-"gene_numeric"
        feature<-"various gene levels"
    } else {
        w1<-which(colnames(filteredData)==feature)
        if(length(which(paste("Colors_",feature,sep="")==colnames(filteredData)))>0){
            if(is.numeric(filteredData[1,w1])){
                column<-as.numeric(as.character(plotData[,w1]))
                type<-"cell_exist_numeric"
            } else {
                column<-factor(plotData[,w1])
                type<-"cell_exist_factor"
            }
            w2<-which(colnames(filteredData)==paste("Colors_",feature,sep=""))
            uc<-unique(as.matrix(plotData[,c(w1,w2)]))
            uc<-uc[sort.list(uc[,1]),]
            my_palette <- as.character(uc[,2])
            names(my_palette)<-uc[,1]
            feature<-as.character(uc[,1])
        } else {
            if(is.numeric(filteredData[1,w1])){
                column<-as.numeric(as.character(plotData[,w1]))
                my_palette<-viridis(50)
                type<-"cell_nonexist_numeric"
                feature<-"various other levels"
            } else {
                column<-factor(plotData[,w1])
                d1<-unique(as.matrix(plotData[,w1]))
                d1<-d1[sort.list(d1)]
                if(length(table(d1))>74){
                    stop("
                    ***** The feature contains more than 74 levels. *****
                    ***** Consider modifying it into a continuous variable in the Design file. *****
                    ")
                }
                d2<-col_vector[1:length(d1)]
                uc<-cbind(d1,d2)
                uc<-uc[sort.list(uc[,1]),]
                my_palette <- as.character(uc[,2])
                names(my_palette)<-uc[,1]
                type<-"cell_nonexist_factor"
                feature<-as.character(uc[,1])
            }
        }
    }
    return(list(column=column,features=feature,palette=my_palette,type=type))
}

#' Assign dot shapes and sizes to the scatterplot
#'
#' Helper function
#' @param plotData data.frame. The data to use for plotting.
#' @param my_colors list. The outcome of makecolors_scatterplot().
#' @param dot.size numeric. The size of the dots to be plotted.
#' @param feature character. The feature to be plotted.
#' @param grouping.by character. The grouping variable to be highlighted.
#' @keywords ShapeSize
#' @return A data.frame
#'
ShapeSize<-function(plotData,my_colors,dot.size,feature,grouping.by){

    # fix the shapes (genes)
    shap<-cbind(c(16,15,17:18,2:14),rep(dot.size,17))
    shap<-rbind(shap,cbind(c(16,15,17:18,2:14),rep(dot.size,17)))

    plotData.filt<-plotData[plotData$SelectedSamples=="Selected",]


    if(!is.null(grouping.by)){

        if(my_colors$type=="gene_numeric" | my_colors$type=="cell_nonexist_numeric" | my_colors$type=="cell_exist_numeric"){
            w<-which(colnames(plotData)==grouping.by)
            if(length(w)==0){
                stop("
                ***** The grouping.by variable does not exist in the Design file. *****
                ")
            }
            tt<-table(plotData[,w],plotData$SelectedSamples)
            tt<-tt[,colnames(tt)=="Selected"]>0
            tt<-names(tt[tt==TRUE])
            if(length(tt)>34){
                stop("
                ***** Using a grouping.by variable with more than 34 levels does not support this type of visualisation. *****
                ")
            }
            for(i in 1:length(tt)){
                plotData.filt$SelectedShapes[as.character(plotData.filt[,w])==tt[i]]<-shap[i,1]
                plotData.filt$SelectedSizes[as.character(plotData.filt[,w])==tt[i]]<-shap[i,2]
            }
            mm<-match(rownames(plotData.filt),rownames(plotData))
            plotData$SelectedShapes[mm]<-plotData.filt$SelectedShapes
            plotData$SelectedSizes[mm]<-plotData.filt$SelectedSizes

        } else {

            if(feature!=grouping.by){
                w1<-which(colnames(plotData)==feature)
                if(length(w1)==0){
                    stop("
                    ***** The feature variable does not exist in the Design file. *****
                    ")
                }
                w2<-which(colnames(plotData)==grouping.by)
                if(length(w2)==0){
                    stop("
                    ***** The grouping.by variable does not exist in the Design file. *****
                    ")
                }
                tt<-table(plotData[,w1],plotData$SelectedSamples,plotData.filt[,w2])
                tt<-tt[,colnames(tt)=="Selected",]>0
                wh<-which(tt==TRUE,arr.ind=TRUE)
                rr<-rownames(tt)[unique(wh[,1])]
                if((length(rr)*length(cc))>34){
                    stop("
                    ***** This visualisation module supports less than 35 groups (feature x grouping.by levels). *****
                    ***** Select a smaller number of feature / grouping.by levels from SampleQuery(). *****")
                }

                index<-0
                for(i in 1:length(rr)){
                    for(j in 1:length(cc)){
                        index<-index+1

                        plotData.filt$SelectedShapes[as.character(plotData.filt[,w1])==rr[i] & as.character(plotData.filt[,w2])==cc[j]]<-shap[index,1]
                        plotData.filt$SelectedSizes[as.character(plotData.filt[,w1])==rr[i] & as.character(plotData.filt[,w2])==cc[j]]<-shap[index,2]
                    }
                }
                mm<-match(rownames(plotData.filt),rownames(plotData))
                plotData$SelectedShapes[mm]<-plotData.filt$SelectedShapes
                plotData$SelectedSizes[mm]<-plotData.filt$SelectedSizes

            }

            if(feature==grouping.by){
                w<-which(colnames(plotData)==feature)
                if(length(w)==0){
                    stop("
                    ***** The feature variable does not exist in the Design file. *****
                    ")
                }
                tt<-table(plotData[,w],plotData$SelectedSamples)
                tt<-tt[,colnames(tt)=="Selected"]>0
                tt<-names(tt[tt==TRUE])
                if(length(tt)>34){
                    stop("
                    ***** Using a feature variable with more than 34 levels does not support this type of visualisation. *****
                    ")
                }
                for(i in 1:length(tt)){
                    plotData.filt$SelectedShapes[as.character(plotData.filt[,w])==tt[i]]<-shap2[i,1]
                    plotData.filt$SelectedSizes[as.character(plotData.filt[,w])==tt[i]]<-shap2[i,2]
                }
                mm<-match(rownames(plotData.filt),rownames(plotData))
                plotData$SelectedShapes[mm]<-plotData.filt$SelectedShapes
                plotData$SelectedSizes[mm]<-plotData.filt$SelectedSizes
            }
        }
    } else {
       plotData$SelectedShapes[plotData$SelectedSamples=="Selected"]<-shap[1,1]
       plotData$SelectedSizes[plotData$SelectedSamples=="Selected"]<-shap[1,2]
    }

 return(plotData)
}

#' Assign text annotation to the ggplot scatterplot
#'
#' Helper function
#' @param plotData data.frame. The data to use for plotting.
#' @param vars character. The variables whose levels will be depicted.
#' @param dims character. The dimensions to be plotted.
#' @keywords addggText
#' @return A vector of names
#'
addggText<-function(plotData,vars,dims){
    ll<-apply(plotData[,match(vars,colnames(plotData))],1,paste,collapse=":")
    x<-cbind(as.numeric(as.character(plotData[,dims[1]])),as.numeric(as.character(plotData[,dims[2]])),ll)
    a1<-aggregate(as.numeric(as.character(x[,1])),list(x[,3]),median)
    a2<-aggregate(as.numeric(as.character(x[,2])),list(x[,3]),median)
    res<-cbind(a1,a2[,2])

    mm<-match(x[,3],res[,1],nomatch=0)
    extra<-res[mm,]
    extra<-cbind(abs(as.numeric(as.character(x[,1])))-abs(as.numeric(as.character(extra[,2]))),
        abs(as.numeric(as.character(x[,2])))-abs(as.numeric(as.character(extra[,3]))),as.character(extra[,1]))
    extra<-cbind(abs(apply(matrix(abs(as.numeric(as.character(extra[,1:2]))),ncol=2),1,sum)),extra[,3])
    extra<-aggregate(as.numeric(as.character(extra[,1])),list(extra[,2]),which.min)

    for(i in 1:nrow(res)){
        res[i,2:3]<-x[which(x[,3]==res[i,1])[extra[i,2]],1:2]
    }

 return(res)
}



#' Generate the scatterplot
#'
#' The function that generates the scatterplot.
#' @param Data data.frame. The data to use for plotting.
#' @param feature character. The feature to be plotted.
#' @param grouping.by character.  The grouping variable to be highlighted automatically.
#' @param highlight.by character.  The second grouping variable to be highlighted on click.
#' @param showDims character. The dimensions to be plotted.
#' @param dot.size numeric. The size of the dots to be plotted.
#' @param interactive logical. If TRUE, an interactive (plotly) plot is generated.
#' @keywords Scatter
#' @return Data for plotting
#' @import plotly ggplot2 viridis ggpubr
#'
Scatter<-function(Data,feature,grouping.by,highlight.by,showDims,dot.size,interactive){

    # extract the data of interest
    plotData<-Data$Design
    filteredData<-Data$selectedData$Design
    mm<-match(rownames(filteredData),rownames(plotData),nomatch=0)
    sel<-rep("Non-Selected",nrow(plotData))
    sel[mm]<-"Selected"
    plotData<-cbind(plotData,SelectedSamples=sel,SelectedShapes=rep(1,nrow(plotData)),SelectedSizes=rep(1,nrow(plotData)))

    # record the dimensions
    dims<-Data$Dimensions


    # colors for gene (for cells features it is the viridis)
    my_colors<-makecolors_scatterplot(Data=Data,plotData=plotData,feature=feature)

    # fix the shapes and sizes
    plotData<-ShapeSize(plotData=plotData,my_colors=my_colors,dot.size=dot.size,feature=feature,grouping.by=grouping.by)

    # check for the dims
    if(length(showDims)<2){
        stop("
        ***** The scatterplot is requires at least two dimensions. Check the showDims parameter. *****
        ")
    }
    cc<-colnames(plotData)[as.numeric(dims[,2])]
    mm<-match(showDims,cc,nomatch=0)
    if(length(mm[mm==0])>0){
        stop("
        ***** One of the requested dimensions does not exist in the data. Consider revising the showDims parameter. *****
        ")
    }

    # check for the existence of the highlight variable
    if(!is.null(highlight.by)){
        ww<-which(colnames(plotData)==highlight.by)
        if(length(ww)==0){

            spl<-unlist(strsplit(highlight.by,":",fixed=TRUE))
            if(length(spl)>1){
                hb<-matrix(0,nrow(plotData),length(spl))
                for(i in 1:length(spl)){
                    hb[,i]<-plotData[,which(colnames(plotData)==spl[i])]
                }
                hb<-apply(hb,1,paste,collapse="_")
                plotData<-cbind(plotData,Interaction=hb)
                highlight.by<-"Interaction"
            } else {
                print(paste("***** The highlight.by variable does not exist on the data. Setting it to NULL! *****",sep=""))
                highlight.by<-1
            }
        }
    } else {
        highlight.by<-1
    }


    if(my_colors$type=="gene_numeric"){

            dot_label<-paste(rownames(plotData),": ",
                                plotData$Condition," ",
                                plotData$CellType," Expr=",
                                my_colors$column,sep="")

    } else {

            dot_label<-paste(rownames(plotData),": ",
                            plotData$Condition," ",
                            plotData$CellType,sep="")

    }



    # start the plot
    ppp1<-c()
    plot.int<-plot.noint<-as.list(rep(0,3))

    ti<- paste("tSNE plot of <b> Gene ID: <a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
                GeneNames(feature,what="Name"),"'>",
                GeneNames(feature,what="Name"),"</a> </b>",sep="")



    if(interactive){
        index<-0
        for(i in 1:(length(mm)-1)){
            for(j in (i+1):length(mm)){

                index<-index+1

                plot.int[[index]] <- plotData %>% highlight_key(formula(paste("~",highlight.by))) %>% plot_ly(
                x=as.numeric(as.character(plotData[,as.numeric(dims[mm[i],2])])),
                y=as.numeric(as.character(plotData[,as.numeric(dims[mm[j],2])])),
                color=my_colors$column,
                colors=my_colors$palette,
                type="scatter",
                marker=list(size=plotData$SelectedSizes,line=list(width=3)),
                text = dot_label,
                symbol=plotData$SelectedShapes,
                symbols=sort(as.numeric(as.character(unique(plotData$SelectedShapes))))
                ) %>%

                layout(title=list(text=ti,font=list(color="black",size=22),x=0.03,y=0.997),
                    xaxis=list(title=paste("<b>tSNE ",colnames(plotData)[as.numeric(dims[mm[i],2])],"</b> <br> <br> (<i>The color gradient bar shows the normalised expression range of the selected gene</i>) ",sep=""),
                    range = range(as.numeric(as.character(Data$Design[,as.numeric(dims[mm[i],2])])))+c(-0.5,0.5)),
                    yaxis=list(title=paste("<b>tSNE ",colnames(plotData)[as.numeric(dims[mm[j],2])],"</b>",sep=""),range = range(as.numeric(as.character(Data$Design[,as.numeric(dims[mm[j],2])])))+c(-0.5,0.5)),showlegend=F)



            }
        }

        if(length(mm)>2){
            ppp1<-subplot(plot.int[[1]],plot.int[[2]],plot.int[[3]],nrows=2,shareX=T,shareY=T)
        } else {
            ppp1<-plot.int[[1]]
        }


    } else {

        # the non-interactive Scatterplot genes is done here
        index<-0
        for(i in 1:(length(mm)-1)){
            for(j in (i+1):length(mm)){

                index<-index+1

                plot.noint[[index]]<-ggplot(plotData,aes_string(
                    x=as.numeric(as.character(plotData[,as.numeric(dims[mm[i],2])])),
                    y=as.numeric(as.character(plotData[,as.numeric(dims[mm[j],2])])),
                    color=my_colors$column)) +
                    #geom_point(size=(dot.size/10)) +
                    geom_point(size=(plotData$SelectedSizes/4),shape=plotData$SelectedShapes,stroke=0.8) +
                    xlim(range(as.numeric(as.character(Data$Design[,as.numeric(dims[mm[i],2])])))) +
                    ylim(range(as.numeric(as.character(Data$Design[,as.numeric(dims[mm[j],2])]))))
                if(my_colors$type=="gene_numeric"){
                    plot.noint[[index]] <- plot.noint[[index]] +
                        scale_colour_gradient(low = my_colors$palette[1], high = my_colors$palette[50])
                }
                if(my_colors$type=="cell_nonexist_numeric"){
                        plot.noint[[index]] <- plot.noint[[index]] +
                            scale_color_viridis()
                }
                if(my_colors$type=="cell_nonexist_factor"){
                        plot.noint[[index]] <- plot.noint[[index]] +
                            scale_color_manual(breaks=my_colors$features,values=my_colors$palette,name=feature)
                }
                if(my_colors$type=="cell_exist_numeric" | my_colors$type=="cell_exist_factor"){
                    plot.noint[[index]] <- plot.noint[[index]] +
                        scale_color_manual(breaks=my_colors$features,values=my_colors$palette,name=feature)
                }
                plot.noint[[index]] <- plot.noint[[index]] +
                    labs(x=colnames(plotData)[as.numeric(dims[mm[i],2])],y=colnames(plotData)[as.numeric(dims[mm[j],2])],
                        color=plotData[,which(colnames(plotData)==feature)]) +
                        labs(color=GeneNames(feature,what="Name")) +
                    theme_light()

                if(!is.null(grouping.by)){
                    if(my_colors$type=="gene_numeric" | my_colors$type=="cell_nonexist_numeric" | my_colors$type=="cell_exist_numeric"){
                        pp<-addggText(plotData=plotData,vars=c(grouping.by,"SelectedSamples"),dims=c(as.numeric(dims[mm[i],2]),as.numeric(dims[mm[j],2])))
                        g<-grep("Non-Selected",pp[,1])
                        if(length(g)>0){
                            pp<-as.matrix(pp[-g,])
                        }
                        pp[,1]<-t(matrix(unlist(strsplit(pp[,1],":")),nrow=2))[,1]
                    } else {
                        pp<-addggText(plotData=plotData,vars=c(feature,grouping.by,"SelectedSamples"),dims=c(as.numeric(dims[mm[i],2]),as.numeric(dims[mm[j],2])))
                        g<-grep("Non-Selected",pp[,1])
                        if(length(g)>0){
                            pp<-as.matrix(pp[-g,])
                        }
                        pp[,1]<-apply(t(matrix(unlist(strsplit(pp[,1],":")),nrow=3))[,1:2],1,paste,collapse=":")
                    }

                } else {

                    if(my_colors$type=="gene_numeric" | my_colors$type=="cell_nonexist_numeric" | my_colors$type=="cell_exist_numeric"){
                        pp<-NULL
                        print(paste(
                        "***** The scatterplot will not be annotated unless the grouping.by variable is specified. *****",sep="")
                        )
                    } else {
                        pp<-addggText(plotData=plotData,vars=c(feature,"SelectedSamples"),dims=c(as.numeric(dims[mm[i],2]),as.numeric(dims[mm[j],2])))
                        g<-grep("Non-Selected",pp[,1])
                        if(length(g)>0){
                            pp<-as.matrix(pp[-g,])
                        }
                        pp[,1]<-t(matrix(unlist(strsplit(pp[,1],":")),nrow=2))[,1]
                    }

                }

                if(!is.null(pp)){
                    plot.noint[[index]] <- plot.noint[[index]] + annotate("text",x=as.numeric(as.character(pp[,2])),y=as.numeric(as.character(pp[,3])),label=as.character(pp[,1]),size=2)
                }
            }
        }

        if(length(mm)>2){
            ppp1<-ggarrange(plot.noint[[1]],ggarrange(plot.noint[[2]], plot.noint[[3]], ncol = 2, labels = c(" ", " ")),nrow = 2,labels = " ")
        } else {
            ppp1<-plot.noint[[1]]
        }
    }

 return(list(Plot=ppp1,feature=feature))
}




#' Assign colors to the violin
#'
#' Helper function
#' @param plotData list. The outcome of the MarkerQuery().
#' @param feature character. The feature to be plotted.
#' @param grouping.by character.  The grouping variable to be highlighted in the x-axis.
#' @param split.by character. The variable to split each violin with.
#' @keywords makecolors_violin
#' @return A list
#' @import RColorBrewer
#'
makecolors_violin<-function(plotData,feature,grouping.by,split.by){

    type<-c()
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    wgrouping<-which(colnames(plotData)==grouping.by)
    wsplit<-which(colnames(plotData)==split.by)


    # various stops and warning
    if(length(which(colnames(plotData)==feature))==0){
        stop("
        ***** This feature name does not exist in the Design file! *****
        ")
    } else {
        if(suppressWarnings(is.na(as.numeric(plotData[1,which(colnames(plotData)==feature)])))){
            stop("
            ***** Character / factor features are not supported in the violin plot! *****
            ")
        }
    }

    if(is.null(grouping.by)){
        stop("
        ***** Define the grouping.by variable of the violin plot. *****
        ")
    }

    if(length(wgrouping)==0){
        stop("
        ***** The grouping.by variable does not exist in the Design file! *****
        ")
    }

    if(length(table(as.character(plotData[,wgrouping])))>74){
        stop("
        ***** The violin plot of more than 74 groups is not supported! *****
        ")
    }

    if(!is.null(split.by)){
        if(length(wsplit)==0){
            split.by<-NULL
            print(paste(
            "***** The grouping2.by variable does not exist in the Design file. Setting grouping2.by = NULL! *****",sep="")
            )
        }
        if(length(table(as.character(plotData[,wsplit])))==1){
            split.by<-NULL
            print(paste(
            "***** The grouping2.by variable contains only one level. Setting grouping2.by = NULL! *****",sep="")
            )
        }
        if(length(table(as.character(plotData[,wsplit])))>2){
            stop("
            ***** The split-type violin plot for more than 2 groups is not supported! *****
            ")
        }
        if(!is.null(split.by)){
          if(wsplit==wgrouping){
              stop("
              ***** The violin plot cannot be generated because the grouping.by variable and the grouping2.by variable for each violin are the same. *****
              ")
          }
        }
    }


    # start the plotting
    if(!is.null(split.by)){
        if(length(which(paste("Colors_",split.by,sep="")==colnames(plotData)))>0){
            colors<-plotData[,which(paste("Colors_",split.by,sep="")==colnames(plotData))]
        } else {
            uu<-sort(as.character(unique(plotData[,wsplit])))
            cc<-col_vector[1:length(unique(plotData[,wsplit]))]
            colors<-rep(0,nrow(plotData))
            for(i in 1:length(uu)){
                colors[as.character(plotData[,wsplit])==uu[i]]<-cc[i]
            }
        }
    } else {
        if(length(which(paste("Colors_",grouping.by,sep="")==colnames(plotData)))>0){
            colors<-plotData[,which(paste("Colors_",grouping.by,sep="")==colnames(plotData))]
        } else {
            uu<-sort(as.character(unique(plotData[,wgrouping])))
            cc<-col_vector[1:length(unique(plotData[,wgrouping]))]
            colors<-rep(0,nrow(plotData))
            for(i in 1:length(uu)){
                colors[as.character(plotData[,wgrouping])==uu[i]]<-cc[i]
            }
        }
    }
    return(list(colors=colors,wgrouping=wgrouping,wsplit=wsplit))
}


#' Generates the Violin plot
#'
#' The function that generates the violin plot
#' @param Data list. The outcome of the MarkerQuery().
#' @param feature character. The feature to be plotted.
#' @param grouping.by character.  The grouping variable to be highlighted in the x-axis.
#' @param split.by character. The variable to split each violin with.
#' @param dot.size numeric. The size of the dots to be plotted.
#' @param interactive logical. If TRUE, an interactive (plotly) plot is generated.
#' @keywords Violin
#' @return Data for plotting
#' @import plotly ggplot2
#'
Violin<-function(Data,feature,grouping.by,split.by,dot.size,interactive){

    # extract the data of interest
    plotData<-data.frame(as.matrix(Data$selectedData$Design),check.names=FALSE)
    dims<-Data$Dimensions

    my_colors<-makecolors_violin(plotData=plotData,feature=feature,grouping.by=grouping.by,split.by=split.by)
    if(length(my_colors$wsplit)==0){
      split.by<-NULL
    }


    ppp2<-c()


    if(interactive){



            if(is.null(split.by)){

                ti<- paste("Violin plot of <b> Gene ID: <a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
                GeneNames(feature,what="Name"),"'>",
                GeneNames(feature,what="Name"),"</a> </b>",sep="")

                ppp2<-plot_ly(plotData,
                x=~as.character(plotData[,my_colors$wgrouping]),
                y=~as.numeric(as.character(plotData[,ncol(plotData)])),
                split=~as.character(plotData[,my_colors$wgrouping]),
                type="violin",
                spanmode="hard",
                text=paste(rownames(plotData),": ",plotData$Condition," ",plotData$CellType,sep=""),
                box = list(visible = F),meanline = list(visible = T),
                points = 'all',
                bandwidth=0.999,
                marker=list(size=dot.size),
                color=I(as.character(plotData$Colors_CellType))) %>%
                    layout(title=list(text=ti,font=list(color="black",size=22),x=0.03,y=0.997),
                    xaxis=list(title=paste("<b>Cell Type</b>",sep="")),
                    yaxis=list(title=paste("<b>log2 Normalised Expression</b>",sep="")),showlegend=F)

            } else {

                ti<- paste("Split-violin plot of <b> Gene ID: <a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
                GeneNames(feature,what="Name"),"'>",
                GeneNames(feature,what="Name"),"</a> </b>",sep="")

                xgroups<-sort(unique(as.character(plotData[,my_colors$wgrouping])))
                vgroups<-sort(unique(as.character(plotData[,my_colors$wsplit])))

                yvals<-as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1],ncol(plotData)]))
                ll<-length(yvals[yvals>0])
                if(ll>0){
                    cols<-I(rep("#ee8572",length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1]]))))
                } else {
                    cols<-rep(NA,length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1]])))
                }



                    ppp2<-plot_ly(
                        x=plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1],my_colors$wgrouping],
                        y=as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1],ncol(plotData)])),
                        type="violin",
                        side="negative",
                        name=vgroups[1],
                        bandwidth=0.999,
                        spanmode="hard",
                        #jitter=0,
                        text=paste(rownames(plotData)[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1]],": ",
                            plotData$Condition[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1]]," ",
                            plotData$CellType[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[1]],sep=""),
                        box = list(visible = F),meanline = list(visible = T),
                        points = 'all',
                        marker=list(size=3),
                        color=cols) %>%
                                layout(title=list(text=ti,font=list(color="black",size=22),x=0.03,y=0.997),
                                yaxis=list(title=paste("<b>log2 Normalised Expression</b>",sep="")),
                                xaxis=list(title=paste("<b>Cell Type by Condition</b> <br> <br> (<i>Missing split-violins indicate that the gene is unexpressed in the cell type / condition</i>) ",sep="")),showlegend=T)

                yvals<-as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2],ncol(plotData)]))
                ll<-length(yvals[yvals>0])

                if(ll>0){
                    cols<-I(rep("#63b7af",length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2]]))))
                } else {
                    cols<-rep(NA,length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2]])))
                }

                    ppp2<-ppp2 %>%
                        add_trace(
                            x=plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2],my_colors$wgrouping],
                            y=as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2],ncol(plotData)])),
                            type="violin",
                            side="positive",
                            name=vgroups[2],
                             bandwidth=0.999,
                             spanmode="hard",
                            #jitter=0,
                            text=paste(rownames(plotData)[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2]],": ",
                                plotData$Condition[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2]]," ",
                                plotData$CellType[plotData[,my_colors$wgrouping]==xgroups[1] & plotData[,my_colors$wsplit]==vgroups[2]],sep=""),
                            box = list(visible = F),meanline = list(visible = T),
                            points = 'all',
                            marker=list(size=3),
                            color=cols)

                if(length(xgroups)>1){
                    for (i in 2:length(xgroups)){

                        yvals<-as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1],ncol(plotData)]))
                        ll<-length(yvals[yvals>0])
                        if(ll>0){
                            cols<-I(rep("#ee8572",length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1]]))))
                        } else {
                            cols<-rep(NA,length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1]])))
                        }

                            ppp2<-ppp2 %>%
                                add_trace(
                                    x=plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1],my_colors$wgrouping],
                                    y=as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1],ncol(plotData)])),
                                    type="violin",
                                    side="negative",
                                    name=vgroups[1],
                                     bandwidth=0.999,
                                     spanmode="hard",
                                    #jitter=0,
                                    text=paste(rownames(plotData)[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1]],": ",
                                        plotData$Condition[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1]]," ",
                                        plotData$CellType[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[1]],sep=""),
                                    box = list(visible = F),meanline = list(visible = T),
                                    points = 'all',
                                    marker=list(size=3),
                                    color=cols,showlegend=FALSE)


                        yvals<-as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2],ncol(plotData)]))
                        ll<-length(yvals[yvals>0])
                        if(ll>0){
                            cols<-I(rep("#63b7af",length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2]]))))
                        } else {
                            cols<-rep(NA,length(as.character(my_colors$colors[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2]])))
                        }

                            ppp2<-ppp2 %>%
                                add_trace(
                                    x=plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2],my_colors$wgrouping],
                                    y=as.numeric(as.character(plotData[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2],ncol(plotData)])),
                                    type="violin",
                                    side="positive",
                                    name=vgroups[2],
                                     bandwidth=0.99,
                                     spanmode="hard",
                                    #jitter=0,
                                    text=paste(rownames(plotData)[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2]],": ",
                                        plotData$Condition[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2]]," ",
                                        plotData$CellType[plotData[,my_colors$wgrouping]==xgroups[i] & plotData[,my_colors$wsplit]==vgroups[2]],sep=""),
                                    box = list(visible = F),meanline = list(visible = T),
                                    points = 'all',
                                    marker=list(size=3),
                                    color=cols,showlegend=FALSE)
                    }
                }
            }

    } else {


        # the non-interactive violin is done here
        if(is.null(split.by)){

            uc<-unique(cbind(as.character(plotData[,my_colors$wgrouping]),as.character(my_colors$colors)))
            uc<-uc[sort.list(uc[,1]),]

            ppp2<-ggplot(plotData,aes(
                x=as.character(plotData[,my_colors$wgrouping]),
                y=as.numeric(as.character(plotData[,ncol(plotData)])),
                fill=as.character(plotData[,my_colors$wgrouping]))) +
                geom_violin() +
                geom_jitter(shape=20, position=position_jitter(0.2),size=(dot.size/10)) +
                scale_fill_manual(breaks=as.character(uc[,1]),values=as.character(uc[,2]),name=grouping.by) +
                labs(x="",y="Expression",title=GeneNames(feature,what="Name")) +
                theme_light()

        } else {

            uc<-unique(cbind(as.character(plotData[,my_colors$wsplit]),as.character(my_colors$colors)))
            uc<-uc[sort.list(uc[,1]),]

            ppp2<-ggplot(plotData, aes(
                x=as.character(plotData[,my_colors$wgrouping]),
                y=as.numeric(as.character(plotData[,ncol(plotData)])),
                fill = as.character(plotData[,my_colors$wsplit]))) +
                geom_split_violin() +
                geom_jitter(shape=20, position=position_jitter(0.2),size=(dot.size/10)) +
                scale_fill_manual(breaks=as.character(uc[,1]),values=as.character(uc[,2]),name=split.by) +
                labs(x="",y="Expression",title=GeneNames(feature,what="Name")) +
                theme_light()

        }

    }

  return(list(Plot=ppp2,feature=feature))
}

#' Splits the violin into two components
#'
#' Helper function
#' @format NULL
#' @usage NULL
#' @return A data.frame
#' @import ggplot2 plyr scales
#' @importFrom ggplot2 ggproto
#'
GeomSplitViolin <- ggproto(`_class` = "GeomSplitViolin", `_inherit` = GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

#' Splits the violin into two components
#'
#' Helper function
#' @param ... other arguments
#'
#' @keywords geom_split_violin
#' @return A data.frame
#' @import ggplot2
#'
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



#' Generates the Volcano plot
#'
#' The function that generates the volcano plot
#' @param Data list. The outcome of the MarkerQuery().
#' @param de.comparison list. A filter specifying the comparison to be plotted.
#' @param logFC_cut numeric. A positive value of the logFC cut-off to be depicted.
#' @param FDR_cut numeric. A value in (0,1] of the FDR cut-off to be depicted.
#' @param show.more logical. If TRUE, an extensive set of gene names are shown.
#' @param dot.size numeric. The size of the dots to be plotted.
#' @param interactive logical. If TRUE, an interactive (plotly) plot is generated.
#' @keywords Volcano
#' @return Data for plotting
#' @import RColorBrewer plotly ggplot2
#'
Volcano<-function(Data,de.comparison,logFC_cut,FDR_cut,show.more,dot.size,interactive){

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    de.comparison<-as.expression(parse(text=de.comparison))
    plotData<-DEplot_entries(deData=Data$DEstats,filters=de.comparison)
    comparison<-rbind(plotData[[2]]$FactorFilters,plotData[[2]]$NumericFilters)
    comparison<-paste(apply(comparison,1,paste,collapse=":"),collapse=" & ")

    ppp3<-c()
    plotData<-plotData[[1]]
    plotData$PValue[plotData$PValue==0]<-min(as.numeric(as.character(plotData$PValue[plotData$PValue>0])))
    lp<- -log(as.numeric(as.character(plotData$PValue)),10)
    signf<-rep("Non-DE",nrow(plotData))
    signf[abs(as.numeric(as.character(plotData$logFC)))>=logFC_cut & as.numeric(as.character(plotData$FDR))<=FDR_cut]<-"DE"
    sel<-rep(" ",nrow(plotData))
    x<-as.character(plotData[,1])
    x<-cbind(x,GeneNames(x,what="Other"),GeneNames(x,what="Name"))
    vcat<-as.numeric(as.character(plotData$logPValue))
    vcat<-vcat[signf=="DE"]
    if(length(vcat)>0){
        vcat<-min(vcat)
    } else {
        vcat<-NA
    }

    if(show.more){
        gg<-match(Data$filteredData$Marker_List,as.character(plotData[,1]))
        gg<-plotData[gg,]
        gg<-gg[sort.list(gg$PValue),]
        gg<-unique(c(as.character(gg[1:min(70,nrow(gg)),1]),Data$selectedData$Marker_List))
        mm<-findGeneType(templ=x,genes=gg,multiple=FALSE)
        for(i in 1:length(mm)){
            signf[mm[i]]<-paste("Marker: ",GeneNames(gg[i],what="Name"),sep="")
        }
    } else {
        mm<-findGeneType(templ=x,genes=Data$selectedData$Marker_List,multiple=FALSE)
        signf[mm]<-paste("Marker: ",GeneNames(Data$selectedData$Marker_List,what="Name"),sep="")
    }
    sel[mm]<-"Selected"
    plotData<-cbind(plotData,logPValue=lp,Type=signf,Selected=sel)

    cols<-rep("#8ca0ca",nrow(plotData))
    cols[plotData$Type=="DE"]<-"#64ba9f"
    cols[plotData$Selected=="Selected"]<-col_vector[1:as.numeric(table(plotData$Selected)[2])]
    sh<-rep("circle",nrow(plotData))
    sh[plotData$Selected!=" "]<-"triangle"
    plotData<-cbind(plotData,Colors=cols,Shapes=sh)

    selData<-plotData[plotData$Selected=="Selected",]
    selData<-selData[sort.list(selData$Type),]
    uu<-cbind(as.character(selData$Type),sort(as.character(selData$Colors)))
    #uu<-uu[sort.list(uu[,1]),]

    complegend<-unlist(strsplit(comparison,":",fixed=T))
    complegend1<-paste(complegend[2:length(complegend)],collapse=": ")
    complegend2<-unlist(strsplit(complegend[length(complegend)],"[: -]"))
    complegend2<-complegend2[(length(complegend2)-1):length(complegend2)]
    if(interactive){
        ppp3<-plot_ly(type="scatter",mode="markers") %>%

        layout(title=list(text=paste("Volcano plot of <b>",complegend1,"</b> cells",sep=""),font=list(color="black",size=22),x=0.03,y=0.997),
        xaxis=list(title=paste("<b>log FC</b> <br> <br> (<i>log FC > 0 [log FC < 0] values indicate increased target expression in ",complegend2[1]," [",complegend2[2],"] cells</i>)",sep="")),
        yaxis=list(title="<b>-log PValue</b>"))

        ppp3<-ppp3 %>% add_trace(name="Non-DE",
            x=as.numeric(as.character(plotData$logFC))[plotData$Selected==" " & plotData$Type=="Non-DE"],
            y=as.numeric(as.character(plotData$logPValue))[plotData$Selected==" " & plotData$Type=="Non-DE"],
            marker=list(color="#82c4c3"),
            text = as.character(plotData[plotData$Selected==" " & plotData$Type=="Non-DE",1]))

        ppp3<-ppp3 %>% add_trace(name="DE",
            x=as.numeric(as.character(plotData$logFC))[plotData$Selected==" " & plotData$Type=="DE"],
            y=as.numeric(as.character(plotData$logPValue))[plotData$Selected==" " & plotData$Type=="DE"],
            marker=list(color="#fc8210"),
            text = as.character(plotData[plotData$Selected==" " & plotData$Type=="DE",1]))

        te<-"Non-significant"
        if(abs(as.numeric(as.character(selData$logFC)))>=logFC_cut & as.numeric(as.character(selData$FDR))<=FDR_cut){
            te<-"Significant"
        }
        ppp3<-ppp3 %>% add_trace(name=uu[,1],mode = 'text',textfont = list(color = 'transparent'),
            x=as.numeric(as.character(selData$logFC)),
            y=as.numeric(as.character(selData$logPValue)),
            marker=list(color="#6a097d",size=(dot.size*4)),
            symbol=selData$Selected,
            text = paste("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",GeneNames(uu[,1],what="Name"),"'>",uu[,1],"</a> (",te,")",sep="")) %>%
                add_lines(y = c(min(as.numeric(as.character(plotData$logPValue))), max(as.numeric(as.character(plotData$logPValue)))), x = c(logFC_cut, logFC_cut),color=I("gray"),size=I(1),showlegend=F) %>%
                add_lines(y = c(min(as.numeric(as.character(plotData$logPValue))), max(as.numeric(as.character(plotData$logPValue)))), x = c(-logFC_cut, -logFC_cut),color=I("gray"),size=I(1),showlegend=F) %>%
                add_lines(x = c(min(as.numeric(as.character(plotData$logFC))), max(as.numeric(as.character(plotData$logFC)))), y = c(vcat, vcat),color=I("gray"),size=I(1),showlegend=F)



    } else {

        if(length(which(names(table(signf))=="DE"))==0){
            uu<-rbind(uu,cbind(matrix("Non-DE",ncol=1),matrix("#8ca0ca",ncol=1)))
        }
        if(length(which(names(table(signf))=="Non-DE"))==0){
            uu<-rbind(uu,cbind(matrix("DE",ncol=1),matrix("#64ba9f",ncol=1)))
        }
        if(length(which(names(table(signf))=="DE"))>0 & length(which(names(table(signf))=="Non-DE"))>0){
            uu<-rbind(uu,cbind(matrix(c("DE","Non-DE"),ncol=1),matrix(c("#64ba9f","#8ca0ca"),ncol=1)))
        }
        uu<-uu[sort.list(uu[,1]),]

        ppp3<-ggplot(plotData,aes(
            x=as.numeric(as.character(logFC)),
            y=as.numeric(as.character(logPValue)),
            color=Type)) +
            geom_point(size=(dot.size/10)) +
            scale_color_manual(breaks=as.character(uu[,1]),values=as.character(uu[,2])) +
            xlim(range(as.numeric(as.character(plotData$logFC)))) +
            ylim(range(as.numeric(as.character(plotData$logPValue)))) +
            labs(x="logFC",y="log PValue",title=comparison) +
            theme_light() +
            geom_vline(xintercept=0,size=0.2)
        ppp3<-ppp3+geom_point(data=selData,aes(x=as.numeric(as.character(logFC)),y=as.numeric(as.character(logPValue))),size=((dot.size/10)*4),shape="triangle") +
                scale_fill_manual(breaks=uu[,1],values=as.character(uu[,2])) +
                geom_vline(xintercept=c(-logFC_cut,logFC_cut), color = "gray",size=0.4) +  geom_hline(yintercept= vcat, color = "gray",size=0.4)

    }

    return(list(Plot=ppp3,feature=paste(comparison,"_",as.character(Data$selectedData$Marker_List),sep="")))
}


#' Generates the MA plot
#'
#' The function that generates the MA plot
#' @param Data list. The outcome of the MarkerQuery().
#' @param de.comparison list. A filter specifying the comparison to be plotted.
#' @param logFC_cut numeric. A positive value of the logFC cut-off to be depicted.
#' @param FDR_cut numeric. A value in (0,1] of the FDR cut-off to be depicted.
#' @param show.more logical. If TRUE, an extensive set of gene names are shown.
#' @param dot.size numeric. The size of the dots to be plotted.
#' @param interactive logical. If TRUE, an interactive (plotly) plot is generated.
#' @keywords MA
#' @return Data for plotting
#' @import RColorBrewer plotly ggplot2
#'
MA<-function(Data,de.comparison,logFC_cut,FDR_cut,show.more,dot.size,interactive){

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    de.comparison<-as.expression(parse(text=de.comparison))
    plotData<-DEplot_entries(deData=Data$DEstats,filters=de.comparison)
    comparison<-rbind(plotData[[2]]$FactorFilters,plotData[[2]]$NumericFilters)
    comparison<-paste(apply(comparison,1,paste,collapse=":"),collapse=" & ")

    ppp4<-c()
    plotData<-plotData[[1]]
    signf<-rep("Non-DE",nrow(plotData))
    signf[abs(as.numeric(as.character(plotData$logFC)))>=logFC_cut & as.numeric(as.character(plotData$FDR))<=FDR_cut]<-"DE"
    sel<-rep(" ",nrow(plotData))
    x<-as.character(plotData[,1])
    x<-cbind(x,GeneNames(x,what="Other"),GeneNames(x,what="Name"))

    if(show.more){
        gg<-match(Data$filteredData$Marker_List,as.character(plotData[,1]))
        gg<-plotData[gg,]
        gg<-gg[sort.list(gg$PValue),]
        gg<-unique(c(as.character(gg[1:min(70,nrow(gg)),1]),Data$selectedData$Marker_List))
        mm<-findGeneType(templ=x,genes=gg,multiple=FALSE)
        for(i in 1:length(mm)){
            signf[mm[i]]<-paste("Marker: ",GeneNames(gg[i],what="Name"),sep="")
        }
    } else {
        mm<-findGeneType(templ=x,genes=Data$selectedData$Marker_List,multiple=FALSE)
        signf[mm]<-paste("Marker: ",GeneNames(Data$selectedData$Marker_List,what="Name"),sep="")
    }
    sel[mm]<-"Selected"
    plotData<-cbind(plotData,Type=signf,Selected=sel)

    cols<-rep("#8ca0ca",nrow(plotData))
    cols[plotData$Type=="DE"]<-"#64ba9f"
    cols[plotData$Selected=="Selected"]<-col_vector[1:as.numeric(table(plotData$Selected)[2])]
    sh<-rep("circle",nrow(plotData))
    sh[plotData$Selected!=" "]<-"triangle"
    plotData<-cbind(plotData,Colors=cols,Shapes=sh)

    selData<-plotData[plotData$Selected=="Selected",]
    selData<-selData[sort.list(selData$Type),]
    uu<-cbind(as.character(selData$Type),sort(as.character(selData$Colors)))

    complegend<-unlist(strsplit(comparison,":",fixed=T))
    complegend1<-paste(complegend[2:length(complegend)],collapse=": ")
    complegend2<-unlist(strsplit(complegend[length(complegend)],"[: -]"))
    complegend2<-complegend2[(length(complegend2)-1):length(complegend2)]

    if(interactive){
        ppp4<-plot_ly(type="scatter",mode="markers") %>%
        layout(title=list(text=paste("MA plot of <b>",complegend1,"</b> cells",sep=""),font=list(color="black",size=22),x=0.03,y=0.997),
            xaxis=list(title=paste("<b>log Average</b> <br> <br> (<i>log FC > 0 [log FC < 0] values indicate increased target expression in ",complegend2[1]," [",complegend2[2],"] cells</i>)",sep="")),
            yaxis=list(title="<b>log FC</b>"))

        ppp4<-ppp4 %>% add_trace(name="Non-DE",
                x=as.numeric(as.character(plotData$logCPM))[plotData$Selected==" " & plotData$Type=="Non-DE"],
                y=as.numeric(as.character(plotData$logFC))[plotData$Selected==" " & plotData$Type=="Non-DE"],
                marker=list(color="#82c4c3"),
                text = as.character(plotData[plotData$Selected==" " & plotData$Type=="Non-DE",1]))

        ppp4<-ppp4 %>% add_trace(name="DE",mode = 'text',textfont = list(color = 'transparent'),
                x=as.numeric(as.character(plotData$logCPM))[plotData$Selected==" " & plotData$Type=="DE"],
                y=as.numeric(as.character(plotData$logFC))[plotData$Selected==" " & plotData$Type=="DE"],
                marker=list(color="#fc8210"),
                text=as.character(plotData[plotData$Selected==" " & plotData$Type=="DE",1]))

        te<-"Non-significant"
        if(abs(as.numeric(as.character(selData$logFC)))>=logFC_cut & as.numeric(as.character(selData$FDR))<=FDR_cut){
            te<-"Significant"
        }
        ppp4<-ppp4 %>% add_trace(name=uu[,1],mode = 'text',textfont = list(color = 'transparent'),
                x=as.numeric(as.character(selData$logCPM)),
                y=as.numeric(as.character(selData$logFC)),
                marker=list(color="#6a097d",size=(dot.size*4)),
                symbol=selData$Selected,
                text = paste("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",GeneNames(uu[,1],what="Name"),"'>",uu[,1],"</a> (",te,")",sep="")) %>%
                    add_lines(x = c(min(as.numeric(as.character(plotData$logCPM))), max(as.numeric(as.character(plotData$logCPM)))), y = c(logFC_cut, logFC_cut),color=I("gray"),size=I(1),showlegend=F) %>%
                    add_lines(x = c(min(as.numeric(as.character(plotData$logCPM))), max(as.numeric(as.character(plotData$logCPM)))), y = c(-logFC_cut, -logFC_cut),color=I("gray"),size=I(1),showlegend=F)

         } else {

             if(length(which(names(table(signf))=="DE"))==0){
                 uu<-rbind(uu,cbind(matrix("Non-DE",ncol=1),matrix("#8ca0ca",ncol=1)))
             }
             if(length(which(names(table(signf))=="Non-DE"))==0){
                 uu<-rbind(uu,cbind(matrix("DE",ncol=1),matrix("#64ba9f",ncol=1)))
             }
             if(length(which(names(table(signf))=="DE"))>0 & length(which(names(table(signf))=="Non-DE"))>0){
                 uu<-rbind(uu,cbind(matrix(c("DE","Non-DE"),ncol=1),matrix(c("#64ba9f","#8ca0ca"),ncol=1)))
             }
             uu<-uu[sort.list(uu[,1]),]

             ppp4<-ggplot(plotData,aes(
                 x=as.numeric(as.character(logCPM)),
                 y=as.numeric(as.character(logFC)),
                 color=Type)) +
                 geom_point(size=(dot.size/10)) +
                 scale_color_manual(breaks=as.character(uu[,1]),values=as.character(uu[,2])) +
                 xlim(range(as.numeric(as.character(plotData$logCPM)))) +
                 ylim(range(as.numeric(as.character(plotData$logFC)))) +
                 labs(x="log Average",y="logFC",title=comparison) +
                 theme_light() +
                 geom_vline(xintercept=0,size=0.2)
             ppp4<-ppp4+geom_point(data=selData,aes(x=as.numeric(as.character(logCPM)),y=as.numeric(as.character(logFC))),size=((dot.size/10)*4),shape="triangle") +
                scale_fill_manual(breaks=uu[,1],values=as.character(uu[,2]))
                ppp4<-ppp4 +  geom_hline(yintercept=c(-logFC_cut,logFC_cut), color = "gray",size=0.4)

         }


         return(list(Plot=ppp4,feature=paste(comparison,"_",as.character(Data$selectedData$Marker_List),sep="")))
}






#' Plot generator function
#'
#' Helper function
#' @param myplot data.frame. The plot to be generated.
#' @param type character.  The type of plot.
#' @param output_folder character.  The output folder.
#' @param interactive logical.  If TRUE, the interative (plotly) plot is generated else a static (ggplot) is generated.
#' @param show.it logical. If TRUE, the plot is shown on screen.
#' @param save.it logical. If TRUE, the plot is stored in a file.
#' @param setDate character. The date stamp generated in CellsQuery().
#' @keywords plotSaver_helper
#' @return A data.frame.
#' @import RColorBrewer htmlwidgets grDevices
#'
plotSaver_helper<-function(myplot,type,output_folder,interactive,show.it,save.it,setDate){
    if(show.it){
       if(interactive){
           print(myplot$Plot)
       } else {
           #if(type!="Heatmap"){
               print(myplot$Plot)
            #} else {
            #   print(pheatmap(as.matrix(myplot$Plot$x),
            #   annotation_col = myplot$Plot$ann,
            #   show_colnames = TRUE,
            #   fontsize = myplot$font,
            #   color = colorRampPalette(c("#00aedb", "#ff4e50"))(n = 50),
            #   cluster_cols = FALSE,cluster_rows=FALSE,annotation_colors=list("",as.character(myplot$colors))))
            #}
       }
    }

    if(save.it){
       if(interactive){
        ff<-paste(unlist(strsplit(myplot$feature,":",fixed=T)),collapse="")
        dd<-unlist(strsplit(setDate,"[_.]"))[-1]
        dd<-paste(c(paste(dd[1:3],collapse=""),paste(dd[4:6],collapse="")),collapse="_")
        htmlwidgets::saveWidget(myplot$Plot,file=paste(output_folder,"InteractivePlots/",type,"_",ff,"_",dd,".html",sep=""))
       } else {
           #if(type!="Heatmap"){
               pdf(paste(output_folder,"NonInteractivePlots/",type,"_",myplot$feature,"_",setDate,".pdf",sep=""))
                plot(myplot$Plot)
               dev.off()
               #} else {
               #pdf(paste(output_folder,"NonInteractivePlots/Heatmap_",myplot$feature,"_",setDate,".pdf",sep=""))
               #   pheatmap(as.matrix(myplot$Plot$x),
               #       annotation_col = myplot$Plot$ann,
               #       show_colnames = TRUE,
               #       fontsize = myplot$font,
               #       color = colorRampPalette(c("#00aedb", "#ff4e50"))(n = 50),
               #       cluster_cols = FALSE,cluster_rows=FALSE,annotation_colors=list("",as.character(myplot$colors)))
               #dev.off()
               #}
        }
    }
}

#' Plot generator function
#'
#' Helper function
#' @param myplot data.frame. The plot to be generated.
#' @param type character. The type of plot.
#' @param output_folder character. The output folder.
#' @param interactive logical. If TRUE, the interative (plotly) plot is generated else a static (ggplot) is generated.
#' @param show.it logical. If TRUE, the plot is shown on screen.
#' @param save.it logical. If TRUE, the plot is stored in a file.
#' @param setDate character. The date stamp generated in CellsQuery().
#' @keywords plotSaver
#' @return A data.frame.
#'
plotSaver<-function(myplot,type,output_folder,interactive,show.it,save.it,setDate){

        if(interactive & save.it){
            if(!dir.exists(paste(output_folder,"InteractivePlots",sep=""))){
                dir.create(paste(output_folder,"InteractivePlots",sep=""))
            }
        }

        if(!interactive & save.it){
            if(!dir.exists(paste(output_folder,"NonInteractivePlots",sep=""))){
                dir.create(paste(output_folder,"NonInteractivePlots",sep=""))
            }
        }

        if(length(myplot)>0){
            plotSaver_helper(myplot=myplot,type=type,output_folder=output_folder,interactive=interactive,show.it=show.it,save.it=save.it,setDate=setDate)
        }
}


#' Table generator function
#'
#' Helper function
#' @param Data data.frame. The outcome og MarkerQuery().
#' @param save.table logical. If TRUE it stores the table in a file.
#' @keywords tableSaver
#' @return A list of tables
#'
tableSaver<-function(Data,save.table){

    # save the summary stats for the genes
    if(save.table){
        if(!dir.exists(paste(Data$Output_Folder,"Tables",sep=""))){
            dir.create(paste(Data$Output_Folder,"Tables",sep=""))
        }
    }

    if(save.table){
        summaries<-t(apply(Data$filteredData$Counts,1,quantile,seq(0,1,0.1)))
        summaries<-cbind(Marker=rownames(summaries),summaries)
        colnames(summaries)[c(2,12)]<-c("Min","Max")
        fi<-paste(Data$Output_Folder,"Tables/SummaryStats_Exact",Data$Exact_Marker_Match,"_",Data$Date_Stamp,".txt",sep="")
        if(!file.exists(fi)){
            write.table(summaries,fi,sep="\t",row.names=F)
        }
    } else {
        summaries<-c()
    }

    # save the values and the design
    if(save.table){
        data_to_store<-as.matrix(Data$filteredData$Counts)
        w1<-which(colnames(Data$filteredData$Design)=="Condition")
        w2<-which(colnames(Data$filteredData$Design)=="CellType")
        data_to_store<-rbind(as.matrix(t(Data$filteredData$Design[,c(w1,w2)])),data_to_store)
        data_to_store<-cbind(Marker=rownames(data_to_store),data_to_store)
        data_to_store[1:2,1]<-""
        fi<-paste(Data$Output_Folder,"Tables/MarkerData_Exact",Data$Exact_Marker_Match,"_",Data$Date_Stamp,".txt",sep="")
        if(!file.exists(fi)){
            write.table(data_to_store,fi,sep="\t",row.names=F)
        }
    } else {
        data_to_store<-c()
    }

    # save the DE stats
    if(save.table){
        if(!is.null(Data$filteredData$DEstats)){
            DE_genes<-data$filteredData$DEstats
            fi<-paste(Data$Output_Folder,"Tables/DEstats_Exact",Data$Exact_Marker_Match,"_",Data$Date_Stamp,".txt",sep="")
            if(!file.exists(fi)){
                write.table(DE_genes,fi,sep="\t",row.names=F)
            }
        } else {
            DE_genes<-c()
        }
    } else {
        DE_genes<-c()
    }

    res<-list(Summary_Stats=summaries,Marker_Data=data_to_store,DE_Stats=DE_genes)
 return(res)
}
