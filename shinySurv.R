##### Shiny App for Survival Analyses #####
#K. Nixon - March 24, 2022

#Begin by loading required data and functions
message("Reading in gene-level TPM data...")
exp<-readRDS("glTPM.RDS")
#exp<-readRDS("batchPANCAN.RDS")
message("Reading in survival data...")
survival<-readRDS("survival.RDS")

survival<-survival[which(rownames(survival) %in% colnames(exp)),]
survival<-survival[match(colnames(exp), rownames(survival)),]
all.equal(colnames(exp), rownames(survival))

source("R/TCGA_Survival Functions.R")

ui<-shiny::fluidPage(
  #App title
  shiny::titlePanel("TCGA PANCAN Survival Analysis"),
  #Sidebar layout
  shiny::sidebarLayout(
    #Pandel for inputs
    shiny::sidebarPanel(
      #Dropdown menu for selecting cancer type
      shiny::selectInput(inputId = "type",
                         label="Select Cancer Type:",
                         choices=c("NULL",unique(survival$cancer.type.abbreviation)),
                         selected="NULL"),
      #Dropdown menu for selecting gender
      shiny::selectInput(inputId = "gender",
                         label="Select Gender:",
                         choices=c("NULL",unique(survival$gender)),
                         selected="NULL"),
      #Dropdown menu for selecting race
      shiny::selectInput(inputId = "race",
                         label="Select Race:",
                         choices=c("NULL", unique(survival$race)),
                         selected="NULL"),
      #Option to select sample types
      shiny::radioButtons(inputId = "norm",
                          label = "Remove Normal Samples:",
                          choices=c("No","Yes"),
                          selected="No"),
      shiny::radioButtons(inputId = "primary",
                          label = "Use only primary tumor samples:",
                          choices=c("No","Yes"),
                          selected = "No"),
      #Option to show confidence interval
      shiny::radioButtons(inputId = "conf",
                          label = "Show Confidence Intervals:",
                          choices=c("No","Yes"),
                          selected="No"),
      #Dropdown menu for selecting the primary gene
      shiny::selectInput(inputId = "x",
                         label="Select Primary Gene:",
                         choices=rownames(exp),
                         selected=rownames(exp)[1]),
      #Dropdown menu for selecting the secondary gene
      shiny::selectInput(inputId = "y",
                         label="Secondary Gene (optional):",
                         choices=c("NULL", rownames(exp)),
                         selected="NULL"),
      #Radio buttons for selecting method
      shiny::radioButtons(inputId = "whichFun",
                          label="Sample Stratification Method:",
                          choices=c("Optimized","Median","Quartile"),
                          selected="Optimized"),
      #Radio buttons for using gene ratio
      shiny::radioButtons(inputId="ratio",
                          label="Use gene expression ratio (primary/secondary):",
                          choices=c("Yes","No"),
                          selected="No"),
      #Dropdown for selecting time
      shiny::selectInput(inputId = "time",
                         label="Survival type:",
                         choices=c("Overall Survival","Disease-Specific Survival","Disease-Free Interval","Progression-Free Interval"),
                         selected="Overall Survival"),
      #Selection for which strata to plot:
      shiny::uiOutput('stratOpt'),
      #Slider for selecting maximum time
      shiny::sliderInput(inputId = "timeLim",
                         label="Select Time Cutoff:",
                         min=0,
                         max=max(survival$OS.time, na.rm=T),
                         value=max(survival$OS.time, na.rm=T)),
      # #Input for a filename
      # shiny::textInput(inputId="outfile",
      #                  label="Filename to save plot:",
      #                  value="SurvivalPlot.png"),
      #Set up a download button for the plot
      shiny::downloadButton(outputId = "downloadPlot",
                            label="Download Plot (survival curve only)")
    ),
    #Main Panel display
    shiny::mainPanel(
      shiny::fluidRow(
        shiny::column(12, shiny::plotOutput(outputId="plot")),
        shiny::column(12, shiny::plotOutput(outputId="violin"))
      )
    )
  )
)

#Set up the server

server<-function(input, output){
  source("R/TCGA_Survival Functions.R")
  makePlot<-shiny::reactive({
    type<-input$type
    types<-"all"
    if(type!="NULL"){
      types<-type
    }
    gender<-input$gender
    race<-input$race
    y<-input$y
    norm<-ifelse(input$norm=="Yes",TRUE, FALSE)
    prim<-ifelse(input$primary=="Yes",TRUE,FALSE)
    conf<-ifelse(input$conf=="Yes",TRUE,FALSE)
    ratio<-FALSE
    strat<-input$strat
    if(type=="NULL"){
      type<-NULL
    }
    if(gender=="NULL"){
      gender<-NULL
    }
    if(race=="NULL"){
      race<-NULL
    }
    if(y=="NULL"){
      y<-NULL
    }
    if(input$ratio=="Yes"){
      ratio<-TRUE
    }
    if(strat[1]=="NULL"){
      strat<-NULL
    }
    #print(strat)
    #print(class(strat))
    #retrieve the subsetted data
    newDat<-subDat(exp, survival, type, gender, race, norm, prim)
    #set up the time/groupings according to the subsetted data and inputs
    time<-newDat$survival$OS.time
    group=newDat$survival$OS
    if(input$time!="Overall Survival"){
      if(input$time=="Disease-Specific Survival"){
        time<-newDat$survival$DSS.time
        group<-newDat$survival$DSS
      }
      if(input$time=="Disease-Free Interval"){
        time<-newDat$survival$DFI.time
        group<-newDat$survival$DFI
      }
      if(input$time=="Progression-Free Interval"){
        time<-newDat$survival$PFI.time
        group<-newDat$survival$PFI
      }
    }
    main<-paste(input$time,"for",types,"samples -",input$x)
    if(!is.null(y)){
      if(isTRUE(ratio)){
        main<-paste0(input$time," for ",types,"samples - log2(",input$x,"/",y,")")
      } else{
        main<-paste(input$time,"for",types,"samples -",input$x,"&",y)
      }
    }
    #Now make the figure
    if(input$whichFun=="Optimized"){
      main<-paste(main,"- Optimized Threshold")
      return(optStrat(x=input$x, y, data=newDat$exp, time, event=group, title=main, useRatio=ratio, plotStrat=strat, retPlot=T, timeLim=input$timeLim, show.conf=conf))
    }
    if(input$whichFun=="Median"){
      ydat<-NULL
      main<-paste(main,"- Median as Threshold")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(plotSurv(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                        method="median", useRatio=ratio), time, status=group, title=main, plotStrat=strat, retPlot=T, timeLim=input$timeLim, show.conf=conf))
    }
    if(input$whichFun=="Quartile"){
      ydat<-NULL
      main<-paste(main,"- Quartiles")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(plotSurv(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                               method="median", useQuartile=T, useRatio=ratio), time, status=group, title=main, plotStrat=strat, retPlot=T, timeLim=input$timeLim, show.conf=conf))
    }
  })
  getDat<-shiny::reactive({
    type<-input$type
    types<-"all"
    if(type!="NULL"){
      types<-type
    }
    gender<-input$gender
    race<-input$race
    y<-input$y
    norm<-ifelse(input$norm=="Yes",TRUE, FALSE)
    prim<-ifelse(input$primary=="Yes",TRUE,FALSE)
    ratio<-FALSE
    if(type=="NULL"){
      type<-NULL
    }
    if(gender=="NULL"){
      gender<-NULL
    }
    if(race=="NULL"){
      race<-NULL
    }
    if(y=="NULL"){
      y<-NULL
    }
    if(input$ratio=="Yes"){
      ratio<-TRUE
    }

    #retrieve the subsetted data
    newDat<-subDat(exp, survival, type, gender, race, norm, prim)
    #set up the time/groupings according to the subsetted data and inputs
    time<-newDat$survival$OS.time
    group=newDat$survival$OS
    if(input$time!="Overall Survival"){
      if(input$time=="Disease-Specific Survival"){
        time<-newDat$survival$DSS.time
        group<-newDat$survival$DSS
      }
      if(input$time=="Disease-Free Interval"){
        time<-newDat$survival$DFI.time
        group<-newDat$survival$DFI
      }
      if(input$time=="Progression-Free Interval"){
        time<-newDat$survival$PFI.time
        group<-newDat$survival$PFI
      }
    }
    main<-paste(input$time,"for",types,"samples -",input$x)
    if(!is.null(y)){
      if(isTRUE(ratio)){
        main<-paste0(input$time," for ",types,"samples - log2(",input$x,"/",y,")")
      } else{
        main<-paste(input$time,"for",types,"samples -",input$x,"&",y)
      }
    }
    #Now make the figure
    if(input$whichFun=="Optimized"){
      main<-paste(main,"- Optimized Threshold")
      return(optStrat(x=input$x, y, data=newDat$exp, time, event=group, title=main, useRatio=ratio, retDat=T, timeLim=input$timeLim)$data)
    }
    if(input$whichFun=="Median"){
      ydat<-NULL
      main<-paste(main,"- Median as Threshold")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                               method="median", useRatio=ratio, retDat=T))
    }
    if(input$whichFun=="Quartile"){
      ydat<-NULL
      main<-paste(main,"- Quartiles")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                      method="median", useQuartile=T, useRatio=ratio, retDat=T))
    }
  })
  getStrat<-shiny::reactive({
    type<-input$type
    types<-"all"
    if(type!="NULL"){
      types<-type
    }
    gender<-input$gender
    race<-input$race
    y<-input$y
    norm<-ifelse(input$norm=="Yes",TRUE, FALSE)
    prim<-ifelse(input$primary=="Yes",TRUE,FALSE)
    ratio<-FALSE
    if(type=="NULL"){
      type<-NULL
    }
    if(gender=="NULL"){
      gender<-NULL
    }
    if(race=="NULL"){
      race<-NULL
    }
    if(y=="NULL"){
      y<-NULL
    }
    if(input$ratio=="Yes"){
      ratio<-TRUE
    }

    #retrieve the subsetted data
    newDat<-subDat(exp, survival, type, gender, race, norm, prim)
    #set up the time/groupings according to the subsetted data and inputs
    time<-newDat$survival$OS.time
    group=newDat$survival$OS
    if(input$time!="Overall Survival"){
      if(input$time=="Disease-Specific Survival"){
        time<-newDat$survival$DSS.time
        group<-newDat$survival$DSS
      }
      if(input$time=="Disease-Free Interval"){
        time<-newDat$survival$DFI.time
        group<-newDat$survival$DFI
      }
      if(input$time=="Progression-Free Interval"){
        time<-newDat$survival$PFI.time
        group<-newDat$survival$PFI
      }
    }
    main<-paste(input$time,"for",types,"samples -",input$x)
    if(!is.null(y)){
      if(isTRUE(ratio)){
        main<-paste0(input$time," for ",types,"samples - log2(",input$x,"/",y,")")
      } else{
        main<-paste(input$time,"for",types,"samples -",input$x,"&",y)
      }
    }
    #Now make the figure
    if(input$whichFun=="Optimized"){
      main<-paste(main,"- Optimized Threshold")
      stratDat<-optStrat(x=input$x, y, data=newDat$exp, time, event=group, title=main, useRatio=ratio, retDat=T, timeLim=input$timeLim)$data
      if(length(grep("category", colnames(stratDat)))==2){
        return(unique(paste(stratDat[,grep("category",colnames(stratDat))[1]],stratDat[,grep("category",colnames(stratDat))[2]], sep=".")))
      } else {
        return(unique(stratDat[,grep("category", colnames(stratDat))]))
      }
    }
    if(input$whichFun=="Median"){
      ydat<-NULL
      main<-paste(main,"- Median as Threshold")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(unique(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                             method="median", useRatio=ratio, retStrat=T)))
    }
    if(input$whichFun=="Quartile"){
      ydat<-NULL
      main<-paste(main,"- Quartiles")
      if(!is.null(y)){
        ydat<-newDat$exp[which(rownames(newDat$exp) %in% y),]
      }
      return(unique(ExpStrat(x=newDat$exp[which(rownames(newDat$exp) %in% input$x),], ydat,
                             method="median", useQuartile=T, useRatio=ratio, retStrat=T)))
    }
  })
  output$stratOpt<-shiny::renderUI({
    opt<-getStrat()
    shiny::selectInput(inputId = "strat",
                       label = "Strata to plot:",
                       choices = c("NULL",opt),
                       multiple = TRUE,
                       selected="NULL")
  })
  output$plot<-shiny::renderPlot({
    makePlot()
  })
  output$violin<-shiny::renderPlot({
    plotStrat(getDat())
  })
  output$downloadPlot<-shiny::downloadHandler(
    filename=function(){
      "plot.png"
    },
    content=function(file){
      ggplot2::ggsave(file, makePlot()$plot, device=strsplit(file,".",T)[[1]][2])
      #ggplot2::ggsave(paste0(file,"_riskTable.png"), makePlot()$table, device=strsplit(file,".",T)[[1]][2])
    },
    contentType="image/png"
  )
}
shiny::shinyApp(ui=ui, server=server)
