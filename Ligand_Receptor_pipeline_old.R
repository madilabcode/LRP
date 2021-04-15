library(dplyr)
library(Seurat)
library(knitr)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyverse)  
library(reshape2)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(grid)
library(igraph)
library(MAST)
library(circlize)
library(viridis)
library(ComplexHeatmap)
library(gridBase)
library(formattable)
library(kableExtra)
library(reticulate)
library(hash)
source_python('Codes/Tf_Graph.py') 
source_python("Codes/utils.py")

#title: "Ligand Recptor pipline algorithm"
#author: "Ron Sheinin"
#date: "14 4 2020"
#output: pdf_document 

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

th = theme_bw() +
  theme(
    plot.title = element_text(size = 24),
    axis.text.x = element_text(size = 24, color='black'),
    axis.text.y = element_text(size = 24, color='black'),
    strip.text.x = element_text(size = 24, color='black',angle=90),
    strip.text.y = element_text(size = 24, color='black',angle=90),
    axis.title.x = element_text(size = 24, color='black'),
    axis.title.y = element_text(size = 24,angle=90),
    legend.text=element_text(size=12))

createLRtable = function(seuratObj,LigandReceptorTable,fromName = fromName,toName=toName,assyaName,thrshold = 0.1,num2compare = "de"){
  
  Expression = GetAssayData(object = seuratObj, slot = assyaName) %>% as.data.frame()
  toN = names(seuratObj@active.ident[seuratObj@active.ident %in% toName])
  fromN = names(seuratObj@active.ident[seuratObj@active.ident %in% fromName])
  toExpression = Expression[,names(Expression) %in% toN]
  fromExpression = Expression[,names(Expression) %in% fromN]
  remove(Expression)
  
  
  if (assyaName == "counts"){
    normel = function(x){
      x = x / sum(x)
      return(x)
    }
    
    for (i in 1:length(fromExpression)){
      fromExpression[i] = normel(fromExpression[i])
    }
    
    for(i in 1:length(toExpression)){
      toExpression[i] = normel(toExpression[i])
    }
    
  }
  
  #subset ligands in "from" and receptors in "to"
  fromLigands = fromExpression[row.names(fromExpression) %in%  LigandReceptorTable$from,] 
  
  toReceptors = toExpression[row.names(toExpression) %in% LigandReceptorTable$to,]
  
  #get expression mean from cells
  fromMean = as.data.frame(apply(fromLigands,1, mean))  
  colnames(fromMean) = "meanExp"
  
  toMean = as.data.frame(apply(toReceptors,1, mean))
  colnames(toMean) = "meanExp"
  
  #take only count expressed Genes:
  fromMean = subset.data.frame(fromMean, meanExp>0)
  toMean = subset.data.frame(toMean, meanExp>0)
  
  #few cutoff options for comparing ligands and receptors
  if(num2compare %>% is.numeric()){
    fromtop = as.data.frame(fromMean[order(-fromMean$meanExp),] %>% head(num2compare))
    totop = as.data.frame(toMean[order(-toMean$meanExp),] %>% head(num2compare))
    
  }else if(num2compare == "Median"){ #above median
    fromMedian = fromMean$meanExp %>% median()
    fromtop = fromMean %>% subset(., meanExp > fromMedian)
    toMedian = toMean$meanExp %>% median()
    totop = toMean %>% subset(., meanExp > toMedian)
    
  }else if(num2compare == "topQ"){ #top Quantile
    fromTQ = fromMean$meanExp %>% quantile() %>% .["75%"]
    fromtop = subset(fromMean, meanExp > fromTQ)
    toTQ = toMean$meanExp %>% quantile() %>% .["75%"]
    totop = subset(toMean, meanExp > toTQ)
    
  }else if(num2compare == "topDecile"){ #top Decile
    fromTD = fromMean$meanExp %>% quantile(prob = seq(0, 1, length = 11)) %>% .["90%"]
    fromtop = subset(fromMean, meanExp > fromTD)
    toTD = toMean$meanExp %>% quantile(prob = seq(0, 1, length = 11)) %>% .["90%"]
    totop = subset(toMean, meanExp > toTD)
    
  }else if (num2compare ==  "de"){
    de_genes_redction = function(geneT,id){
      markerall = FindMarkers(seuratObj, ident.1 = id, only.pos = TRUE,
                              min.pct = 0.2, logfc.threshold = 0.2,test.use = "MAST")
      
      geneT$geneName = row.names(geneT)
      geneT = geneT[which(row.names(geneT) %in% row.names(markerall)),]
      return (list(markerall,geneT))
    }
    ft = de_genes_redction(fromMean,fromName)
    tt = de_genes_redction(toMean,toName)
    
    fromtop = ft[[2]]
    markerallL = ft[[1]]
    totop = tt[[2]]
    markerallR = tt[[1]]
  }
  else{
    stop("invaled num2compare")
  }
  
  LigandReceptorTableSub = LigandReceptorTable[,c(1,2)] %>% unique()
  LigandReceptorTableSub = LigandReceptorTableSub[which(LigandReceptorTableSub$from %in% rownames(fromtop) & LigandReceptorTableSub$to %in% rownames(totop)),]
  #LigandReceptorTableSub = LigandReceptorTableSub[which(LigandReceptorTableSub$to %in% rownames(totop)),]
  if (dim(LigandReceptorTableSub)[1] == 0){
    return("ex")
  }
  
  for(i in unique(LigandReceptorTableSub$from)){
    expVal = fromtop[which(rownames(fromtop) == i),"meanExp"]
    LigandReceptorTableSub[which(LigandReceptorTableSub$from == i),"fromEXP"] = expVal}
  
  for(i in unique(LigandReceptorTableSub$to)){
    expVal = totop[which(rownames(totop) == i),"meanExp"]
    LigandReceptorTableSub[which(LigandReceptorTableSub$to == i),"toEXP"] = expVal}
  
  
  colnames(LigandReceptorTableSub) = c("Ligand","Receptor","LigandExp","ReceptorExp")
  #remove overlaps between ligand and receptor
  LR = LigandReceptorTableSub[!(LigandReceptorTableSub$Ligand %in% intersect(LigandReceptorTableSub$Ligand,LigandReceptorTableSub$Receptor)),] 
  LR = LR[!(LR$Receptor %in% intersect(LR$Ligand,LR$Receptor)),] #reciprocally
  
  if (assyaName == "counts"){
    LR[3] = log10(LR[3])
    LR[4] = log10(LR[4])
  }
  if (num2compare ==  "de"){
    return(list(LR,markerallL,markerallR))
  }else{
    return(LR)
  }
}

python_phyper = function(q,m,n,k){
  return = phyper(q,m,n,k,lower.tail = FALSE)
}

DElegenedNcolor = function(seuratObj, fromName, toName, LR,markerallL = NULL,markerallR = NULL, use_location_capacity = FALSE){ 
  #***create the color function and legend for differential expression data***#
  
  if (is.null(markerallL) | is.null(markerallR)){
    #*#* DE Ligands calculation #*#*
    markerallL = FindMarkers(seuratObj, ident.1 = fromName, only.pos = TRUE,
                             min.pct = 0.1, logfc.threshold = 0.1,test.use = "MAST")
    #*#* DE Receptors calculation #*#*
    markerallR = FindMarkers(seuratObj, ident.1 = toName, only.pos = TRUE,
                             min.pct = 0.1, logfc.threshold = 0.1,test.use = "MAST")
  }
  
  #create the color scheme and legend DE pval
  DE_DataL = markerallL[which(markerallL %>% row.names() %in% LR$Ligand),] 
  DE_DataL$gene = DE_DataL %>% row.names()
  DE_DataL_Sub = DE_DataL[,c("p_val_adj","gene")]
  DE_DataL_Sub$p_val_adj = -log10(DE_DataL_Sub$p_val_adj) 
  
  
  #create the color scheme and legend DE pval
  DE_DataR = markerallR[which(markerallR %>% row.names() %in% LR$Receptor),] 
  DE_DataR$gene = DE_DataR %>% row.names()
  DE_DataR_Sub = DE_DataR[,c("p_val_adj","gene")]
  DE_DataR_Sub$p_val_adj = -log10(DE_DataR_Sub$p_val_adj) 
  
  DE_Data_Sub = rbind(DE_DataL_Sub,DE_DataR_Sub) %>% unique()
  DE_Data_Sub_No_Inf = DE_Data_Sub[DE_Data_Sub$p_val_adj < Inf,]
  
  abs_maxDE = quantile(abs(c(DE_Data_Sub_No_Inf$p_val_adj) - 0.5), 0.95, na.rm = TRUE)
  DEcol_fun = colorRamp2(c(0, 0.5 + abs_maxDE), c("white", "#800080"))
  lgd_DE = Legend(at = c(0, (DE_Data_Sub_No_Inf$p_val_adj %>% max/2) %>% round,DE_Data_Sub_No_Inf$p_val_adj %>% max %>% round),
                  col_fun = DEcol_fun, title_position = "leftcenter-rot", title = "-log(DE Pval)",legend_height = unit(2, "cm"))
  return(list(DEcol_fun,lgd_DE,DE_Data_Sub))
}

TF_DSA_FUNC_dircted  = function(ProtNamesInfo,ProtActions,toExpression,legRet ,TfDict, path = "C:/Users/Ron/Desktop/10X_Eilam" ){ 
  setwd(path)
  DSA = DSA_anaylsis(toExpression,legRet$Receptor,ProtActions,ProtNamesInfo,tfs = TfDict)
  DSA_Table = DSA[[1]]
  DSA_grpahs = DSA[[2]]
  return(list(DSA_Table,DSA_grpahs))
}
 
DSAlegenedNcolor = function(DSAFunck, ProtNamesInfo,ProtActions,toExpression,legRet,TfDict,path = "C:/Users/Ron/Desktop/10X_Eilam" ){
  DSA = DSAFunck(ProtNamesInfo,ProtActions,toExpression,legRet,TfDict,path)
  allReceptorDSA = DSA[[1]]
  DSA_grpahs = DSA[[2]]
  #allReceptorDSASacle = allReceptorDSA
  allReceptorDSA$DSA = allReceptorDSA$DSA
  allReceptorDSA = allReceptorDSA[allReceptorDSA$DSA != Inf &  allReceptorDSA$DSA > 0,]
  
  if (length(allReceptorDSA$DSA) == 0){
    return("ex")
  }
  
  abs_maxmean_weight = quantile(abs(c(allReceptorDSA$DSA) - 0.5),0.85, na.rm = TRUE)
  DSAcol_fun = colorRamp2(c(0, (0.5 + abs_maxmean_weight)), c("white", "#28a10a"))
  lgd_DSA = Legend(at = c(0, (allReceptorDSA$DSA %>% max()/2) %>% round(digits = 2),
                          allReceptorDSA$DSA %>% max %>% round(digits = 2)),
                   col_fun = DSAcol_fun, title_position = "leftcenter-rot", title = "DSA",legend_height = unit(2, "cm"))
  return(list(DSAcol_fun,lgd_DSA,allReceptorDSA,DSA_grpahs))
}

createCircosPlots = function(toExpression,LR,DE ,DSA,fromName  ,toName , num2compare, seuratObj,de_recptors = NULL,path = "/plotOut"){
  #***integrate everything and create the circos plot***#
  TFR = hash(DSA[[4]]$Tfs_Per_Recp_in_max_flow())
  maxTf = 0
  for (key in keys(TFR)){
    if(TFR[[key]] > maxTf){
      maxTf =  as.numeric(TFR[[key]])
    }
  } 
  # """LR = createLRtable(seuratObj,fromName ,toName , num2compare )
  #LR = LR[order(LR$Ligand),] 
  
  #DE = DElegenedNcolor(seuratObj)
  DEcol_fun = DE[[1]]
  lgd_DE = DE[[2]]
  DE_Data_Sub = DE[[3]]
  
  #  DSA = DSAlegenedNcolor(LR,threshold,toName,seuratObj)
  DSAcol_fun = DSA[[1]]
  lgd_DSA = DSA[[2]]
  DSA_Table = DSA[[3]]
  
  #*#* create the expression color scheme and legend #*#*
  GeneExpression = data.frame(gene = c(LR[,1] %>% as.character,LR[,2]%>% as.character), mean = c(LR[,3]%>% as.numeric(),LR[,4]%>% as.numeric)) %>% unique()
  GeneExpression$gene = GeneExpression$gene %>% as.character()
  rownames(GeneExpression) = GeneExpression$gene %>% as.character()
  
  EXPRcol_fun = colorRamp2(c(min(GeneExpression$mean), (max(GeneExpression$mean)+min(GeneExpression$mean))/2,
                             max(GeneExpression$mean)), c("#FFFF00", "#FF9966", "red"))
  lgd_EXPR = Legend(at = c(min(GeneExpression$mean) %>% round(digits = 2),
                           ((max(GeneExpression$mean)+min(GeneExpression$mean))/2) %>% round(digits = 2),
                           max(GeneExpression$mean)%>% round(digits = 2)),
                    col_fun = EXPRcol_fun, title_position = "leftcenter-rot", title = " log(Relative Expression)",legend_height = unit(4, "cm"))
  
  
  # EXPRcol_fun = colorRamp2(c(min(GeneExpression$mean), mean(GeneExpression$mean), mean(GeneExpression$mean)) + sd(GeneExpression$mean)
  #            , c("blue", "white", "red"))
  #lgd_EXPR = Legend(at = c(min(GeneExpression$mean) %>% round(digits = 20),
  #                        (mean(GeneExpression$mean)) %>% round(digits = 20),
  #                      (mean(GeneExpression$mean)) + sd(GeneExpression$mean)%>% round(digits = 2)),
  #               col_fun = EXPRcol_fun, title_position = "leftcenter-rot", title = "Expression",legend_height = unit(2, "cm"))
  
  
  
  #combine legends
  lgd_list_vertical = packLegend(lgd_EXPR,lgd_DE ,lgd_DSA, gap = unit(0.7, "cm"))
  
  #receprors rectangles are grey
  
  if (!is.null(de_recptors)){
    grid.colPRE = cbind(as.character(LR$Receptor), ifelse(LR$Receptor %in% de_recptors, "#33CCCC","grey")) %>% as.data.frame()
  }else{
    grid.colPRE = cbind(as.character(LR$Receptor), c("grey")) %>% as.data.frame()
    
  }
  #naming the sectors (slices) of the circle
  grid.col = setNames(as.character(grid.colPRE$V2), grid.colPRE$V1)
  #vector of features that are in the plot
  factors =  unique(c(as.character(LR$Ligand),as.character(LR$Receptor)))
  #subset LR to just names, otherwize it affect the width of the arrows
  LR = LR[,c("Ligand","Receptor")]

  #the plot itself
  pdf(paste0(paste(sep = "" ,getwd(),path),"/",fromName,"_",toName,"Circos.pdf"),width=10,height=7)
  circos.clear()
  # circos.clear()
  #making the middle of the circle vertical as opposed to horizontale
  circosPlot = circos.par(start.degree = 90, clock.wise = F, points.overflow.warning=FALSE)
  
  chordDiagram(LR,scale = F,
               big.gap = 20,
               grid.col = grid.col,
               directional = 1,
               annotationTrack = "grid",
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               link.sort = TRUE,
               link.visible = ifelse(is.null(de_recptors) | LR$Receptor %in% de_recptors,TRUE,FALSE),
               preAllocateTracks = list(track.height = LR %>% dimnames() %>% unlist() %>% strwidth() %>% max())) 
  
  #Expression loop  
  for(i in intersect(GeneExpression$gene, factors)) {
      circos.rect(xleft = 0,
       xright = c(grep(paste0("^",i,"$"),LR$Ligand),
                  grep(paste0("^",i,"$"),LR$Receptor)) %>% length(),
       ybottom = 0.7, ytop = 0.9, 
       col = ifelse(i %in% LR$Ligand |is.null(de_recptors) | i %in% de_recptors ,EXPRcol_fun(GeneExpression[i,"mean"]),"gray"), 
       border = ifelse(i %in% LR$Ligand |is.null(de_recptors) | i %in% de_recptors ,EXPRcol_fun(GeneExpression[i,"mean"]),"gray"),
       sector.index = i, track.index = 1)
  }
  #DE loop  
  for(a in intersect(DE_Data_Sub$gene, factors)) {
      circos.rect(xleft = 0,
      xright = c(grep(paste0("^",a,"$"),LR$Ligand),
                 grep(paste0("^",a,"$"),LR$Receptor)) %>% length(),
      ybottom = 0.95, ytop = 1.15, 
      col = ifelse(a %in% LR$Ligand | is.null(de_recptors) | a %in% de_recptors , DEcol_fun(DE_Data_Sub[a,"p_val_adj"]),"gray"), 
      border = ifelse(a %in% LR$Ligand |is.null(de_recptors) | a %in% de_recptors , DEcol_fun(DE_Data_Sub[a,"p_val_adj"]),"gray"),
      sector.index = a, track.index = 1)
  }
  #DSA loop
  for(receptor in intersect (DSA_Table$Recp, factors)){
      circos.rect(xleft = 0,
      xright = grep(paste0("^",receptor,"$"),LR$Receptor) %>% length,
      ybottom = 1.2, ytop = 1.2 + 2*max(0.25, (as.numeric(TFR[[receptor]])/maxTf)), 
      col = ifelse(is.null(de_recptors) | receptor %in% de_recptors,DSAcol_fun(DSA_Table[receptor,"DSA"]),"gray"), 
      border = ifelse(is.null(de_recptors) | receptor %in% de_recptors,DSAcol_fun(DSA_Table[receptor,"DSA"]),"gray"),
      sector.index = receptor, track.index = 1)
  }
  #gene names track
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(cex = 0.7,CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))},  bg.border = NA)
  
  circle_size = unit(1, "snpc")
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  upViewport()
  draw(lgd_list_vertical,just = c("left", "bottom"),x = unit(1, "cm"), y = unit(1, "cm"))
  dev.off()
  print(paste0("Done making ",fromName," to ",toName))

}

dsa_score_per_cell_all_cluster = function(obj,toExpression, DSA_Table, sacle_factor = 1){
  toExpression = toExpression[row.names(toExpression) %in% unlist(DSA_Table$Recp),]
  for(i in 1:length(row.names(toExpression))){
    toExpression[i,] = (log1p(toExpression[i,]) *  DSA_Table[unlist(DSA_Table$Recp) == row.names(toExpression)[i],]$DSA)
  }
  
  dsa_per_cell = sapply(toExpression, sum)
  dsa_per_cell = dsa_per_cell * sacle_factor
  names(dsa_per_cell) = names(toExpression)
  obj[["DSA_SCORE"]] = dsa_per_cell
  
  return(obj)
}

dsa_score_per_cell = function(obj,toExpression, DSA_Table, sacle_factor = 1){
  toExpression = toExpression[row.names(toExpression) %in% unlist(DSA_Table$Recp),]
  for(i in 1:length(row.names(toExpression))){
    toExpression[i,] = (log1p(toExpression[i,]) *  DSA_Table[unlist(DSA_Table$Recp) == row.names(toExpression)[i],]$DSA)
  }
  
  if (!"DSA_SCORE" %in% names(obj@meta.data)){
    dsa_per_cell = sapply(toExpression, sum)
    dsa_per_cell = dsa_per_cell * sacle_factor
    names(dsa_per_cell) = names(toExpression)
    obj[["DSA_SCORE"]] = dsa_per_cell
  }else{
    dsa_per_cell = sapply(toExpression, sum)
    dsa_per_cell = dsa_per_cell * sacle_factor
    names(dsa_per_cell) = names(toExpression)
    dsa_per_cell = as.data.frame(dsa_per_cell)
    dsa_per_cell$cell = row.names(dsa_per_cell)
    DSA_temp = obj[["DSA_SCORE"]]
    DSA_temp$cell = row.names(DSA_temp)
    DSA_temp = merge(DSA_temp,dsa_per_cell,by.x = "cell", by.y = "cell",all.x = TRUE)
    DSA_temp$DSA_SCORE = ifelse(is.na(DSA_temp$DSA_SCORE),DSA_temp$dsa_per_cell,DSA_temp$DSA_SCORE)
    dsa_factor = DSA_temp$DSA_SCORE
    names(dsa_factor) = DSA_temp$cell
    obj[["DSA_SCORE"]] = dsa_factor
  }
  return(obj)
}

dsa_score_mean_per_cluster = function(obj,to,from,legRet,DSA_Table, sacle_factor = 1){
  toGenes  = subset(obj, idents = to)
  toExpression  = GetAssayData(object = toGenes , slot = "counts") %>% as.data.frame()
  
  fromGenes = subset(obj, idents = from )
  fromGExpression  = GetAssayData(object = fromGenes , slot = "counts") %>% as.data.frame()
  
  return(mean(dsa_with_lig(toExpression,fromGExpression,DSA_Table,legRet)))
}

build_or_use_classfier = function(obj1,obj2,from1,from2,to1,to2,lr1,lr2=NULL,modle = NULL){
  toGenes_1  = subset(obj1, idents = to1)
  toExpression_1  = GetAssayData(object = toGenes_1 , slot = "counts") %>% as.data.frame()
  
  fromGenes_1 = subset(obj1, idents = from1 )
  fromGExpression_1  = GetAssayData(object = fromGenes_1 , slot = "counts") %>% as.data.frame()
  
  
  toGenes_2 = subset(obj2, idents = to2 )
  toExpression_2 = GetAssayData(object = toGenes_2, slot = "counts") %>% as.data.frame()
  
  fromGenes_2 = subset(obj2, idents = from2 )
  fromGExpression_2 = GetAssayData(object = fromGenes_2, slot = "counts") %>% as.data.frame()
  
  if(!is.null(lr2)){
    legRet = rbind(lr1,lr2)
  }else{
    legRet = lr1
  }
  
  ProtNamesInfo = read.csv("files/humanProtinInfo.csv")
  ProtActions = read.csv("files/humanProtinAction.csv")
  
  #Expression_scale = GetAssayData(object = obj1, slot = "scale.data") %>% as.data.frame()
  #toExpression = GetAssayData(object = pbmc_unite, slot = "counts") %>% as.data.frame()
  
  
  source_python('Codes/pyDictToR.py')
  lst = main_py_to_R(ProtNamesInfo, ProtActions)
  source_python('Codes/Tf_Graph.py') 
  
  TfDict = lst[[1]]
  pro_action = lst[[3]]
  pro_name = lst[[2]]
  
  #distDic = hash(py$dist_dickt)
  
  
  DSA_Table_1 = TF_DSA_FUNC_dircted(ProtNamesInfo,ProtActions,toExpression_1,legRet,TfDict,getwd())
  DSA_Table_2 = TF_DSA_FUNC_dircted(ProtNamesInfo,ProtActions,toExpression_2,legRet,TfDict,getwd())
  args1 = list(toExpression_1,fromGExpression_1,DSA_Table_1[[1]],legRet)
  args2 = list(toExpression_2,fromGExpression_2,DSA_Table_2[[1]],legRet)
  if (!is.null(modle)){
    lst = DSA_Classfier(args1,args2,modle = modle,plot_name = paste("roc",from1,from2,sep = "_"))
  }else{
    lst = DSA_Classfier(args1,args2,return_modle = TRUE,plot_name = paste("roc",from1,from2,sep = "_"))
  }
  return (lst)
}

draw_flow_graph = function(DSA_graphs,recp){
  df = DSA_graphs$make_df_flow(recp)
  g = graph_from_data_frame(df,directed = TRUE)
  V(g)$color = "#00CCCC"
  tkplot(g,layout = layout_as_tree(g,recp),edge.label=round(as.numeric(E(g)$flow), 4))
  
}

draw_flow_graph_from_df = function(df,root = "t"){
  g = graph_from_data_frame(df,directed = TRUE)
  V(g)$color = "#00CCCC"
  tkplot(g,layout = layout_as_tree(g,root),edge.label=round(as.numeric(E(g)$flow), 4))
}

deCircosPlot = function(obj,fromName, toName,legRet,DSA_Table,DSA_graph,de_recptors,pathPlot){
  markerallL = FindMarkers(obj, ident.1 = fromName,ident.2 = toName,features = unlist(legRet["Ligand"]), only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.1,test.use = "MAST")
  markerallR = FindMarkers(obj, ident.1 = toName,ident.2 = fromName,features = unlist(legRet["Receptor"]), only.pos = TRUE,
                           min.pct = 0.1, logfc.threshold = 0.1,test.use = "MAST")
  
  DE = DElegenedNcolor(obj,fromName,toName,legRet,markerallL,markerallR)
  
  DEcol_fun = DE[[1]]
  lgd_DE = DE[[2]]
  DE_Data_Sub = DE[[3]]
  
  DSA_Table = DSA_Table[DSA_Table$DSA != Inf,]
  abs_maxmean_weight = quantile(abs(c(DSA_Table$DSA) - 0.5),0.85, na.rm = TRUE)
  DSAcol_fun = colorRamp2(c(0, (0.5 + abs_maxmean_weight)), c("white", "#28a10a"))
  lgd_DSA = Legend(at = c(0, (DSA_Table$DSA %>% max()/2) %>% round(digits = 2),
                          DSA_Table$DSA %>% max %>% round(digits = 2)),
                   col_fun = DSAcol_fun, title_position = "leftcenter-rot", title = "DSA",legend_height = unit(2, "cm"))
  
  createCircosPlots(toExpression,legRet,DE,list(DSAcol_fun,lgd_DSA,DSA_Table,DSA_graph),fromName,toName,"de",obj,de_recptors,path = pathPlot) 
}

LigandReceptorPipline = function(obj,fromName, toName,ProtActions,ProtNamesInfo,TfDict,stw,thrshold = 0.1,plot_path = NULL,per_claster = FALSE){
  setwd(stw)
  #source_python("Tf_Graph.py") 
  LigandReceptorTable = read.delim("files/LigandReceptorTableMouse.tsv",sep = "\t", header = T, quote = "" )  
  #remove overlaps between ligand and receptor
  LigandReceptorTable = LigandReceptorTable[!LigandReceptorTable$from %in% LigandReceptorTable$to,]
  LigandReceptorTable = LigandReceptorTable[!LigandReceptorTable$to %in% LigandReceptorTable$from,]
  
  LR  = createLRtable(obj,LigandReceptorTable,fromName,toName, "counts",thrshold = thrshold )
  if (length(LR) == 1){
    return("ex")
  }
  legRet = LR[[1]]
  markerallL = LR[[2]]
  markerallR = LR[[3]]
  #  return(legRet)
  #saveRDS(legRet,file="legRet.Robj")
  DE = DElegenedNcolor(obj,fromName,toName,legRet,markerallL,markerallR)
  
  DEcol_fun = DE[[1]]
  lgd_DE = DE[[2]]
  DE_Data_Sub = DE[[3]]
  
  toGenes = subset(obj, idents = toName )
  toExpression.scale = GetAssayData(object = toGenes, slot = "scale.data") %>% as.data.frame()
  toExpression = GetAssayData(object = toGenes, slot = "counts") %>% as.data.frame()
  
  DSA = DSAlegenedNcolor(TF_DSA_FUNC_dircted, ProtNamesInfo,ProtActions,toExpression,legRet,TfDict,stw)
  
  if (length(DSA) == 1){
    return("ex")
  }
  
  DSAcol_fun = DSA[[1]]
  lgd_DSA = DSA[[2]]
  DSA_Table = DSA[[3]]
  DSA_grpahs  = DSA[[4]]
  DSA_Table = DSA_Table[DSA_Table$DSA != Inf,]
  
  print(plot_path)
  if(is.null(plot_path)){
    createCircosPlots(toExpression,legRet,DE,DSA,fromName,toName,"de",obj)
  }else{
    createCircosPlots(toExpression,legRet,DE,DSA,fromName,toName,"de",obj,path = plot_path )
    
  }

  ##################3
  
  if(per_claster == TRUE ){
    obj = dsa_score_per_cell_all_cluster(obj,toExpression,DSA_Table)
    
    DSAPlot = FeaturePlot(obj, features = c("DSA_SCORE"),pt.size=1.5)  +
      scale_colour_gradientn(colours = jet.colors(10)) + th
    
    
  }else{
    obj = dsa_score_per_cell(obj,toExpression,DSA_Table)
    
    DSAPlot = FeaturePlot(obj, features = c("DSA_SCORE"),pt.size=1.5)  +
      scale_colour_gradientn(colours = jet.colors(10)) + th
    
  }
  
  if(is.null(plot_path)){
    ggsave( plot = DSAPlot, filename =  paste0(fromName,"_",toName,"_","DSAPlot.pdf")                ,path = paste(sep = "" ,getwd(),"/plotOut"), device = "pdf")
  }else{
    print(paste(sep = "" ,getwd(),plot_path))
    ggsave( plot = DSAPlot, filename =  paste0(fromName,"_",toName,"_","DSAPlot.pdf")               ,path = paste(sep = "" ,getwd(),plot_path), device = "pdf")
  }
  return(list(legRet,DSA_Table,DSA_grpahs,obj))
}
