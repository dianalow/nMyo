#' Read the data
#'
#' It reads / loads the data of interest.
#' @param data_file character. The location of the nMyo data in the system (an RDS file). It is essentially a list of three components. The first component (CorrCounts) are the normalised marker expression profiles. It is a G x (N+1) data matrix of G markers (e.g. genes) and N samples (e.g. cells). The first column of the matrix is the Marker IDs and contains the marker names. The rest of the columns are named with a unique sampleID for each sample and should contain the normalised expression values. The second component (DesignTable) describes the experimental design. It is an N x M data matrix where M is a number of sample characteristics measured for each of the N samples (e.g. library size, cell type, condition etc). The table necessarily contains the following column names: (1) 'SampleID' the unique sampleID for each sample (should match exactly the column names of the Counts file), (2) 'Condition' is the experimental conditions (e.g. Normal, Control, Disease etc), (3) 'CellType' describes one or more cell types used in the experiment (if not relevant, a constant value such as 'Bulk_Cells' can be used)  and (4) at least 2 of "Dim1", "Dim2" and "Dim3" giving the values of the dimensionality reduction components (e.g. PCA). The third component contains the differential expression analysis results for each marker as produced by edgeR. The first column has the Marker IDs (in the same format as those in the Counts file). Each ID appears more than once since multiple comparisons have been done. Other necessary columns are: (1) 'logFC' with the log fold changes, (2) 'logCMP' with the log average values for each marker over all samples, (3) 'PValue' with the estimated P-values, (4) 'FDR' with the FDRs, (5) 'Comparison' with he relevant comparison of two conditions or other separated with '-' to have a meaningful logFC interpretation (e.g. Normal-Disease) and (6) 'CellType' with the cell type of the comparison (it should partly match the CellType levels of the design table).
#' @param Markers_file character. The location of the file with the marker names to be visualised. It should be a .txt file listing the marker IDs of interest (preferably in the same format as those of the Counts file). The contents of this file will be contrasted to the marker IDs of the Counts file to obtain the marker data for visualisation. If NULL (default) the user is prompted to key the marker ID of interest in the MarkerQuery().
#' @param is.Exact logical. If TRUE, the exact marker IDs of the Markers_file are retrieved from the Counts_file. Otherwise a pattern algorithm is used to return all pattern matching names. Default is TRUE.
#' @param logFC_cutoff numeric. A positive value specifying the logFC criterion of the differentially expressed genes of the DE_file. If the Markers_file is missing and the DE_file is present, the algorithm will retrieve the markers that satisfy the logFC and FDR criteria. Default is 0.
#' @param FDR_cutoff numeric. A value in (0,1] specifying the FDR criterion of the differentially expressed genes of the DE_file. If the Markers_file is missing and the DE_file is present, the algorithm will retrieve the markers that satisfy the logFC and FDR criteria. Default is 1.
#' @keywords readData
#' @return data list. A list with the 'Counts' data, the 'Design' data, the 'Annotation' data (marker IDs), the 'DEstats' data (differential expression results), the 'Output_Folder' where the final results (plots and tables) will be stored, the 'Dimensions' to be plotted, the 'Exact_Marker_Match' specifying the marker finding pattern, the 'logFC' and 'FDR' cutoffs and a sublist 'filteredData' containing only the data of the markers of interest for visualisation. Slot filteredData contains an extra component with the marker IDs to be visualised (filteredData$Marker_List). If both the DE_file and the Markers_files are missing then Marker_List contains all the marker names of the Counts_file.
#' @export
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#'
readData<-function(data_file,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1){

    # stopping rule if the names are not correct or the required files do not exist
    #fs<-firstStop(Counts_file=Counts_file,Design_file=Design_file)

    # read the required data: counts, design; Fix the annotation; Read the optional data: DE, Markers
    data<-RequiredDataReader(data_file)
    data<-testDataValidity(data=data)
    genes<-ReadMarkers(Markers_file=Markers_file)
    if(!is.null(genes)){
        data<-c(data,list(Status="Marker file",Exact_Marker_Match=is.Exact,logFC=logFC_cutoff,FDR=FDR_cutoff))
    } else {
        data<-c(data,list(Status="DE file",Exact_Marker_Match=is.Exact,logFC=logFC_cutoff,FDR=FDR_cutoff))
    }

    # find the final gene list through a list of rules and existing files
    sigDE_genes_status<-finaliseMarkerList(data=data,genes=genes)

    # find existing PCA/tSNE/UMAP dimensions
    w<-which(colnames(data$Design)=="Dim1" | colnames(data$Design)=="Dim2" | colnames(data$Design)=="Dim3")
    if(length(w)>0){
        dims<-cbind(colnames(data$Design)[w],w)
    }

    # define the output folder
    folder<-system.file('', package = 'nMyo')
    data_folder<-paste(paste(folder[1:(length(folder)-1)],collapse="/"),"/",sep="")
    output_folder<-paste(data_folder,"output/",sep="")
    if (!file.exists(output_folder)){
         dir.create(output_folder)
    }

    #summarise the data
    data<-c(data,list(Output_Folder=output_folder,Dimensions=dims))

    # summarise the filtered data
    redData<-fixFilteredData(data=data,sigDE.genes.status=sigDE_genes_status)

    res<-c(data,list(temp=redData,filteredData=redData))

    setDate<-unlist(strsplit(date()," ",fixed=T))
    setDate<-setDate[c(1,2,3,5,4)]
    setDate[5]<-chartr(old = ":", new = ".", setDate[5])
    setDate<-paste(setDate,collapse="_")
    res$filteredData$Design<-cbind(res$filteredData$Design,t(res$filteredData$Counts))
    res$Date_Stamp<-setDate

    return(res)
}
