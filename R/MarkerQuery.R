#' Select the marker of interest for visualisation
#'
#' It selects the marker to be plotted and tabulated.
#' @param Data list. The outcome of readData().
#' @param marker character. The ID of the marker to be visualised.
#' @keywords MarkerQuery
#' @return data list. A list includying the outcome of readData() and an extra slot 'selectedData' with the data of the marker that has been selected for visualisation.
#' @export
#' @importFrom stringdist stringsim
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#' data<-MarkerQuery(Data=data,marker="Postn")
#'
MarkerQuery<-function(Data,marker){

    # remove previous selectedData folder
    Data<-DataTester(Data=Data,func="MarkerQuery_types")

    # define the selectedData from the filtered data without the genes
    cc<-setdiff(colnames(Data$filteredData$Design),colnames(Data$temp$Design))
    sel<-Data$filteredData
    sel$Design<-sel$Design[,-match(cc,colnames(Data$filteredData$Design))]

    # pick the gene of interest
    mm<-findGeneType(templ=sel$Annotation,genes=marker,correction=TRUE,multiple=FALSE)
    reportIt<-0
    if(mm[1]==0){
        if(length(grep("ENSMUSG",marker))==0){
            s <- strsplit(tolower(marker), " ")[[1]]
            marker<-paste(toupper(substring(s, 1,1)), substring(s, 2),sep="", collapse=" ")
        } else {
            marker<-unlist(strsplit(toupper(marker),":"))[1]
        }
        mm<-grep(marker,as.character(sel$Annotation[,1]))
        reportIt<-1
    }
    mm<-unique(mm)
    aa<-as.character(sel$Annotation[mm,1])
    sims<-rep(0,length(aa))

    for(i in 1:length(aa)){
        sims[i]<-stringsim(marker,aa[i])
    }
    aa<-aa[sort.list(sims,decreasing=TRUE)]
    mm<-findGeneType(templ=sel$Annotation,genes=aa[1],multiple=FALSE)

    # stop if the gene name does not exist
    if(mm==0){
        print("This gene ID does not match any of the existing IDs!")
      return(NULL)
    }

    # find the data of the particular gene
    sel$Counts<-sel$Counts[mm,]
    sel$Annotation<-sel$Annotation[mm,]
    sel$Marker_List<-rownames(sel$Counts)

    if(reportIt==0){
        print(paste("Gene ",sel$Marker_List," is selected for visualisation.",sep=""))
    } else {
        print(paste("Gene ",marker," does not exactly match any of the existing IDs! The top match ",sel$Marker_List," is selected for visualisation.",sep=""))
    }

    if(!is.null(sel$DEstats)){
        if(nrow(sel$DEstats)>0){
            x<-as.character(sel$DEstats[,1])
            x<-cbind(x,t(matrix(unlist(strsplit(x,":",fixed=T)),nrow=2)))
            mm<-findGeneType(templ=x,genes=rownames(sel$Counts),multiple=TRUE)
            sel$DEstats<-sel$DEstats[mm>0,]
        }
    }

    sel$Design<-cbind(sel$Design,t(sel$Counts))

    #if(moveit){
    #    Data$filteredData<-sel
    #}

    res<-c(Data,list(selectedData=sel))
  return(res)
}

