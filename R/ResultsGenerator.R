#' Generates the dimensionality reduction scatterplot
#'
#' It creates a dimensionality reduction plot (e.g. PCA or t-SNE) highlighting the samples of interest
#' @param Data list. The outcome of MarkerQuery().
#' @param grouping.by character. A grouping variable to be depicted in the plot. It groups (color code) the samples with common values for this variable. It can take any of the factor variables depicted at the column names of Data$Design. Default is NULL.
#' @param highlight.by character. A second grouping variable to highlight certain point on the INTERACTIVE plot upon double click. It highlights the samples with common values for this variable and shades the rest. It can take any of the factor variables depicted at the column names of Data$Design. Default is CellType.
#' @param showDims character vector. A character vector of size 2 or 3 specifing the dimensions to be plotted (Dim1, Dim2 and Dim3). If all three are specified then a 3D plot is generated otherwise a 2D plot of the specified dimensions is returned. Default is c("Dim1","Dim2"), i.e. the first 2 dimensionality reduction components.
#' @param size numeric. It specifies the plotted dot size. Default is 5.
#' @param show.plot logical. If TRUE, it shows the plot on the screen. Default is TRUE.
#' @param save.plot logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/<InteractivePlots or NonInteractivePlots>. Default is FALSE.
#' @param save.table logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/Tables. Default if FALSE.
#' @param interactive logical. If TRUE, it generates an interactive plot with plotly otherwise a static plot with ggplot. Default is TRUE.
#' @keywords doScatterDR
#' @return data list. The plots and the tables of interest stored in automatically generated <data folder>/Output/ folders. The tables include a summary table of quantiles of the marker expression, the differential expression analysis results of the markers of interest and the normalised data for the markers of interest.
#' @export
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#'
#' data<-MarkerQuery(Data=data,marker="Postn")
#'
#' data<-doScatterDR(Data=data,
#'          highlight.by="CellType",
#'          show.plot=TRUE,
#'          save.plot=FALSE)
#'
doScatterDR<-function(Data,
                        grouping.by=NULL,
                        highlight.by="CellType",
                        showDims=c("Dim1","Dim2"),
                        size=5,
                        show.plot=TRUE,
                        save.plot=FALSE,
                        save.table=FALSE,
                        interactive=TRUE){

    # do the plot
    my_plot<-Scatter(Data=Data,
                    feature=Data$selectedData$Marker_List,
                    grouping.by=grouping.by,
                    highlight.by=highlight.by,
                    showDims=showDims,
                    dot.size=size,
                    interactive=interactive)

    gg<-"CT"
    if(highlight.by!="CellType"){
        gg<-"CO"
    }
    # save the plot
    plotSaver(myplot=my_plot,
                type=paste("DRScatterplot_by_",gg,sep=""),
                output_folder=Data$Output_Folder,
                interactive=interactive,
                show.it=show.plot,
                save.it=save.plot,
                setDate=Data$Date_Stamp)

    # save the tabulated information: summary stats, expression and design, DE stats
    tab<-tableSaver(Data=Data,save.table=save.table)

    # update existing entries
    Data<-DataTester(Data=Data,func="DRScatterplot")

    res<-c(Data,tab,list(DRScatterplot=my_plot))

 return(res)
}




#' Generates the Violin plot
#'
#' It creates the Violin plot of a marker for the cell types of interest
#' @param Data list. The outcome of MarkerQuery().
#' @param grouping.by character. A grouping variable to be depicted in the plot. It groups (color code) the samples with common values for this variable. It can take any of the factor variables depicted at the column names of Data$Design. Default is CellType.
#' @param grouping2.by character. A grouping variable to split each violin plot. Default is Condition.
#' @param size numeric. It specifies the plotted dot size. Default is 5.
#' @param show.plot logical. If TRUE, it shows the plot on the screen. Default is TRUE.
#' @param save.plot logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/<InteractivePlots or NonInteractivePlots>. Default is FALSE.
#' @param save.table logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/Tables. Default if FALSE.
#' @param interactive logical. If TRUE, it generates an interactive plot with plotly otherwise a static plot with ggplot. Default is TRUE.
#' @keywords doViolin
#' @return data list. The plots and the tables of interest stored in automatically generated <data folder>/Output/ folders. The tables include a summary table of quantiles of the marker expression, the differential expression analysis results of the markers of interest and the normalised data for the markers of interest.
#' @export
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#'
#' data<-MarkerQuery(Data=data,marker="Postn")
#'
#' data<-doViolin(Data=data,
#'          show.plot=TRUE,
#'          save.plot=FALSE)
#'
doViolin<-function(Data,
                    grouping.by="CellType",
                    grouping2.by="Condition",
                    size=5,
                    show.plot=TRUE,
                    save.plot=FALSE,
                    save.table=FALSE,
                    interactive=TRUE){

    # do the plot
    my_plot<-Violin(Data=Data,
                    feature=Data$selectedData$Marker_List,
                    grouping.by=grouping.by,
                    split.by=grouping2.by,
                    dot.size=size,
                    interactive=interactive)

    gg<-"CT"
    if(grouping2.by=="Condition"){
        gg<-"CO"
    }
    # save the plot
    plotSaver(myplot=my_plot,
                type=paste("Violin_by_",gg,sep=""),
                output_folder=Data$Output_Folder,
                interactive=interactive,
                show.it=show.plot,
                save.it=save.plot,
                setDate=Data$Date_Stamp)

    # save the tabulated information: summary stats, expression and design, DE stats
    tab<-tableSaver(Data=Data,save.table=save.table)

    # update existing entries
    # update existing entries
    Data<-DataTester(Data=Data,func="Violin")

    res<-c(Data,tab,list(Violin=my_plot))

 return(res)
}



#' Generates the MA plot
#'
#' It creates the MA plot for a particular cell type and condition
#' @param Data list. The outcome of MarkerQuery().
#' @param filters character. A set of conditions that specifies which differential expression results should be plotted. It should ALWAYS be a level of the Comparison variable of the DEstats component.
#' @param logFC numeric. A positive value specifying the symmetric logFC cut-off to be depicted. Default is 1.
#' @param FDR numeric. A value between (0,1] specifying the FDR cutoff to be depicted. Default is 0.01.
#' @param show.more logical. If TRUE, all markers of Data$filteredData$Marker_List are highlighted, otherwise only the marker of Data$selectedData$Marker_List is highlighted. Default is FALSE.
#' @param size numeric. It specifies the plotted dot size. Default is 5.
#' @param show.plot logical. If TRUE, it shows the plot on the screen. Default is TRUE.
#' @param save.plot logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/<InteractivePlots or NonInteractivePlots>. Default is FALSE.
#' @param save.table logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/Tables. Default if FALSE.
#' @param interactive logical. If TRUE, it generates an interactive plot with plotly otherwise a static plot with ggplot. Default is TRUE.
#' @keywords doMA
#' @return data list. The plots and the tables of interest stored in automatically generated <data folder>/Output/ folders. The tables include a summary table of quantiles of the marker expression, the differential expression analysis results of the markers of interest and the normalised data for the markers of interest.
#' @export
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#'
#' data<-MarkerQuery(Data=data,marker="Postn")
#'
#' data<-doMA(Data=data,
#'          filters="Comparison=='Fibroblast: Sham-MI'",
#'          logFC=1,
#'          FDR=0.01,
#'          show.plot=TRUE,
#'          save.plot=FALSE)
#'
doMA<-function(Data,
                filters,
                logFC=1,
                FDR=0.01,
                show.more=FALSE,
                size=5,
                show.plot=TRUE,
                save.plot=FALSE,
                save.table=FALSE,
                interactive=TRUE){

    # do the plot
    if(!is.null(Data$DEstats)){
        my_plot<-MA(Data=Data,
                    de.comparison=filters,
                    logFC_cut=logFC,
                    FDR_cut=FDR,
                    show.more=show.more,
                    dot.size=size,
                    interactive=interactive)
    } else {
        stop("
        ***** The MA plot cannot be generated because the DE file does not exist! *****
        ")
    }

    #gg<-unlist(strsplit(filters,"[= : ']"))[4]
    # save the plot
    plotSaver(myplot=my_plot,
                #type=paste("MA_",gg,sep=""),
                type="MA",
                output_folder=Data$Output_Folder,
                interactive=interactive,
                show.it=show.plot,
                save.it=save.plot,
                setDate=Data$Date_Stamp)

    # save the tabulated information: summary stats, expression and design, DE stats
    tab<-tableSaver(Data=Data,save.table=save.table)

    # update existing entries
    # update existing entries
    Data<-DataTester(Data=Data,func="MA")

    res<-c(Data,tab,list(MA=my_plot))

 return(res)
}





#' Generates the Volcano plot
#'
#' It creates the Volcano plot for a particular cell type and condition
#' @param Data list. The outcome of MarkerQuery().
#' @param filters character. A set of conditions that specifies which differential expression results should be plotted. It should ALWAYS be a level of the Comparison variable of the DEstats component.
#' @param logFC numeric. A positive value specifying the symmetric logFC cut-off to be depicted. Default is 1.
#' @param FDR numeric. A value between (0,1] specifying the FDR cutoff to be depicted. Default is 0.01.
#' @param show.more logical. If TRUE, all markers of Data$filteredData$Marker_List are highlighted, otherwise only the marker of Data$selectedData$Marker_List is highlighted. Default is FALSE.
#' @param size numeric. It specifies the plotted dot size. Default is 5.
#' @param show.plot logical. If TRUE, it shows the plot on the screen. Default is TRUE.
#' @param save.plot logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/<InteractivePlots or NonInteractivePlots>. Default is FALSE.
#' @param save.table logical. If TRUE, it saves the plot in an automatically generated folder <data folder>/Output/Tables. Default if FALSE.
#' @param interactive logical. If TRUE, it generates an interactive plot with plotly otherwise a static plot with ggplot. Default is TRUE.
#' @keywords doVolcano
#' @return data list. The plots and the tables of interest stored in automatically generated <data folder>/Output/ folders. The tables include a summary table of quantiles of the marker expression, the differential expression analysis results of the markers of interest and the normalised data for the markers of interest.
#' @export
#' @examples
#' data(nMyo_Data)
#' data<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)
#'
#' data<-MarkerQuery(Data=data,marker="Postn")
#'
#' data<-doVolcano(Data=data,
#'          filters="Comparison=='Fibroblast: Sham-MI'",
#'          logFC=1,
#'          FDR=0.01,
#'          show.plot=TRUE,
#'          save.plot=FALSE)
#'
doVolcano<-function(Data,
                filters,
                logFC=1,
                FDR=0.01,
                show.more=FALSE,
                size=5,
                show.plot=TRUE,
                save.plot=FALSE,
                save.table=FALSE,
                interactive=TRUE){

    # do the plot
    if(!is.null(Data$DEstats)){
        my_plot<-Volcano(Data=Data,
                         de.comparison=filters,
                         logFC_cut=logFC,
                         FDR_cut=FDR,
                         show.more=show.more,
                         dot.size=size,
                         interactive=interactive)
    } else {
        stop("
        ***** The Volcano plot cannot be generated because the DE file does not exist! *****
        ")
    }

    #gg<-unlist(strsplit(filters,"[= : ']"))[4]
    # save the plot
    plotSaver(myplot=my_plot,
                #type=paste("Volcano_",gg,sep=""),
                type="Volcano",
                output_folder=Data$Output_Folder,
                interactive=interactive,
                show.it=show.plot,
                save.it=save.plot,
                setDate=Data$Date_Stamp)

    # save the tabulated information: summary stats, expression and design, DE stats
    tab<-tableSaver(Data=Data,save.table=save.table)

    # update existing entries
    Data<-DataTester(Data=Data,func="Volcano")

    res<-c(Data,tab,list(Volcano=my_plot))

 return(res)
}


