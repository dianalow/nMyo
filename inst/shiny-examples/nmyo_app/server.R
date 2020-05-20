options(shiny.maxRequestSize=30*1024^2)

shinyServer(function(input, output,session) {

  #### Auto stop app when browser is closed ####
  session$onSessionEnded(stopApp)

  # #### Libraries & functions ####
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(plotly)
  library(RColorBrewer)
  library(stringdist)
  library(plyr)
  library(viridis)
  library(stringr)

  source_folder<-getwd()
  genenameselected<-NULL

  values <- reactiveValues(dataloaded = 0, datasaved = 0)

  #### Load data file : Run analysis of a single dataset ####
  data(nMyo_Data)
  thedata<-readData(nMyo_Data,Markers_file=NULL,is.Exact=TRUE,logFC_cutoff=0,FDR_cutoff=1)

  hide(id = "loading-content", anim = TRUE, animType = "fade")
  shinyjs::show("app-content")

  #### Define variables ####

  #### Button toggles ####

  # Toggle load gene button
  observe({
    shinyjs::toggleState("loadGeneButton", !is.null(thedata) && length(input$genename)!=0)
  })

  observe({
    shinyjs::toggleState("showdataButton", values$dataloaded!=0)
  })

  observe({
    shinyjs::toggleState("savedataButton", values$dataloaded!=0)
  })

  observe({
    shinyjs::toggleState("downloadsavedData", values$datasaved!=0)
  })

  #### Button clicks : run/load/download ####

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("nmyo_data", "zip", sep=".")
    },
    content = function(fname) {
      withBusyIndicatorServer("downloadData",{
        shiny::withProgress(
          message = "Preparing zip file....",
          value = 0,
          {
            shiny::incProgress(1/5)
            thedir=dirname(system.file("data", package = "nMyo"))
            fs = list.files(thedir, full.names = TRUE)
            shiny::incProgress(2/5)
            zip(zipfile=fname, files=fs)
            shiny::incProgress(1/5)
            if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
          }
        )
      })
    },
    contentType = "application/zip"
  )

  output$downloadsavedData <- downloadHandler(
    filename = function() {
      paste("generated_data", "zip", sep=".")
    },
    content = function(fname) {
      withBusyIndicatorServer("downloadsavedData",{
        shiny::withProgress(
          message = "Preparing zip file....",
          value = 0,
          {
            currdir=getwd()
            setwd(thedata$Output_Folder)
            shiny::incProgress(1/5)
            fs = list.files(".", full.names = TRUE)
            shiny::incProgress(2/5)
            zip(zipfile=fname, files=fs)
            shiny::incProgress(1/5)
            if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
            setwd(currdir)
          }
        )
      })
    },
    contentType = "application/zip"
  )


  ##### > load gene #####
  observeEvent(input$loadGeneButton, {
    geneinput<-input$genename
    withBusyIndicatorServer("loadGeneButton", {
        output$geneloadstatus <- renderPrint({genenameselected<<-MarkerQuery(Data=thedata,marker=geneinput)})
        values$dataloaded=1
    })
  })

  observeEvent(input$showdataButton, {
    #filter1=strsplit(input$celltypeSelect," & ")
    filter1=input$celltypeSelect
    #filter1=paste("CellType=='",filter1[[1]][1],"' & ","Comparison=='",filter1[[1]][2],"'",sep="")
    filter1=paste("Comparison=='",filter1,"'",sep="")
    switch(input$plottypeSelect,
           "t-SNE"={

             if(input$scatterCheckbox==TRUE){
               scattercheck="Condition"
             } else {
               scattercheck="CellType"
             }
             output$plotgraph <- renderPlotly({
               ddata<-doScatterDR(Data=genenameselected,
                               highlight.by=scattercheck,
                               show.plot=TRUE,
                               save.plot=FALSE)
               ddata$DRScatterplot$Plot
               })

           },
           "MA"={
             output$plotgraph <- renderPlotly({
               ddata<-doMA(Data=genenameselected,
                        filters=filter1,
                        logFC=input$logfc,
                        FDR=input$fdr,
                        show.plot=TRUE,
                        save.plot=FALSE)
               ddata$MA$Plot
               })
           },
           "Volcano"={
             output$plotgraph <- renderPlotly({
               ddata<-doVolcano(Data=genenameselected,
                             filters=filter1,
                             logFC=input$logfc,
                             FDR=input$fdr,
                             show.plot=TRUE,
                             save.plot=FALSE)
               ddata$Volcano$Plot
             })
           },
           "Violin"={
             if(input$scatterCheckbox==TRUE){
               viocheck="Condition"
             } else {
               viocheck=NULL
             }
             output$plotgraph <- renderPlotly({
               ddata<-doViolin(Data=genenameselected,
                               grouping2.by=viocheck,
                            show.plot=TRUE,
                            save.plot=FALSE)
               ddata$Violin$Plot
             })
           }
           )
  })

  observeEvent(input$savedataButton, {
    values$datasaved=1
    filter1=paste("Comparison=='",input$celltypeSelect,"'",sep="")
    switch(input$plottypeSelect,
           "t-SNE"={

             if(input$scatterCheckbox==TRUE){
               scattercheck="CellType"
             } else {
               scattercheck=NULL
             }
               ddata<-doScatterDR(Data=genenameselected,
                                  highlight.by=scattercheck,
                                  show.plot=FALSE,
                                  save.plot=TRUE)


           },
           "MA"={
               ddata<-doMA(Data=genenameselected,
                           filters=filter1,
                           logFC=input$logfc,
                           FDR=input$fdr,
                           show.plot=FALSE,
                           save.plot=TRUE)
           },
           "Volcano"={
               ddata<-doVolcano(Data=genenameselected,
                                filters=filter1,
                                logFC=input$logfc,
                                FDR=input$fdr,
                                show.plot=FALSE,
                                save.plot=TRUE)
           },
           "Violin"={
             if(input$scatterCheckbox==TRUE){
               viocheck="Condition"
             } else {
               viocheck=NULL
             }
               ddata<-doViolin(Data=genenameselected,
                               grouping2.by=viocheck,
                               show.plot=FALSE,
                               save.plot=TRUE)
           }
    )
  })
  #### Button clicks : changing tabs ####

  observeEvent(input$gotostart, {
    updateTabItems(session, "tabs", "quickTab")
  })


  observeEvent(input$gotofaq, {
    updateTabItems(session, "tabs", "faq")
  })


})
