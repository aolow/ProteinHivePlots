### Hive Plot Function with two colors
## A. Olow, last update 7/29/2015

#load packages (install before use if necessary)
library(HiveR); library(grid); library(RColorBrewer)

# set working directory to where you stor the files on nodes/connections

# data annotating protein ids and hive node parameters (id, axis, radius, size, color)
load("nodes.RData")
# data with connections - full names (for searching) and numbered ids id.x (Kinase), id.y(Peptide), id (Substrate)
load("tst.RData")


# Hive Plot for protein acting both as kinase and substrate
myHive <- function(protein, color, color2){
    #find id of your protein
    kin <- unique(tst3[tst3$Kinase==paste0("kin_", protein),4])
    sub <- unique(tst3[tst3$Substrate==paste0("sub_", protein),6])
    query <- unique(c(kin,sub))
    # make edges for proteins in query only
    load("tst3.RData")
    tst <- tst3[(tst3$id %in% query | tst3$id.x %in% query),]
    total_kin <- length(unique(tst$Kinase))
    total_pep <- length(unique(tst$Peptide))-1
    total_sub <- length(unique(tst$Substrate))
    tmp1 <- tst[, 4:5]
    if(length(row.names(tmp1))>0){
        tmp1$color <- color
        names(tmp1) <- c("id1", "id2", "color")}
    tmp2 <- tst[, 5:6]
    if(length(row.names(tmp2))>0){
        tmp2$color <- color
        names(tmp2) <- c("id1", "id2", "color")}
    tmp3 <- tst[is.na(tst$Peptide),]
    tmp3 <- tmp3[, c(4,6)]
    if(length(row.names(tmp3))>0){
        tmp3$color <- color2
        names(tmp3) <- c("id1", "id2", "color")
    }
    edges <- rbind(tmp1, tmp2, tmp3)
    dead_edges <-  which(apply(edges, 1, function(x) sum(is.na(x))) > 0)
    if(length(dead_edges) > 0){
        edges <- edges[-dead_edges,]
    }
    edges <- data.frame(edges, weight = rep(1.5, nrow(edges)), stringsAsFactors = F)
    total_edges <- length(edges$id1)
    data <- list()
    data$nodes <- nodes
    data$edges <- edges
    data$type <- "2D"
    data$desc <- "proteins"
    data$axis.cols <- rep("#00000000", 3)
    attr(data, "class") <- "HivePlotData"
    # check if dataset is correct
    chkHPD(data, confirm=T)
    # plot
    
    plotHive(data, bkgnd = NULL, ch=100, axLabs=c(paste0("Kinase\n", total_kin), paste0("Peptide\n", total_pep), paste0("Substrate\n", total_sub) ),
    axLab.gpar = gpar(col = "black", fontsize = 14),
    axLab.pos = c(450, 450, 450),
    )
    grid.text(paste0("Total edges: ", total_edges), x = 1200, y = -2700, default.units = "native", gp = gpar(fontsize = 12))
    grid.text(protein, x = 0, y = 2000, default.units = "native", gp = gpar(fontsize = 20))
}




### Hive plot for protein acting as kinase only

myHive_kin <- function(protein, color, color2){
    #find id of your protein
    kin <- unique(tst3[tst3$Kinase==paste0("kin_", protein),4])
    query <- kin
    # make edges
    load("tst3.RData")
    tst <- tst3[(tst3$id %in% query | tst3$id.x %in% query),]
    total_kin <- length(unique(tst$Kinase))
    total_pep <- length(unique(tst$Peptide))-1
    total_sub <- length(unique(tst$Substrate))
    tmp1 <- tst[, 4:5]
    if(length(row.names(tmp1))>0){
        tmp1$color <- color
        names(tmp1) <- c("id1", "id2", "color")}
    tmp2 <- tst[, 5:6]
    if(length(row.names(tmp2))>0){
        tmp2$color <- color
        names(tmp2) <- c("id1", "id2", "color")}
    tmp3 <- tst[is.na(tst$Peptide),]
    tmp3 <- tmp3[, c(4,6)]
    if(length(row.names(tmp3))>0){
        tmp3$color <- color2
        names(tmp3) <- c("id1", "id2", "color")
    }
    edges <- rbind(tmp1, tmp2, tmp3)
    dead_edges <-  which(apply(edges, 1, function(x) sum(is.na(x))) > 0)
    if(length(dead_edges) > 0){
        edges <- edges[-dead_edges,]
    }
    edges <- data.frame(edges, weight = rep(1.5, nrow(edges)), stringsAsFactors = F)
    total_edges <- length(edges$id1)
    data <- list()
    data$nodes <- nodes
    data$edges <- edges
    data$type <- "2D"
    data$desc <- "proteins"
    data$axis.cols <- rep("#00000000", 3)
    attr(data, "class") <- "HivePlotData"
    
    # check if dataset is correct
    chkHPD(data, confirm=T)
    # plot
    plotHive(data, bkgnd = NULL, ch=100, axLabs=c(paste0("Kinase\n", total_kin), paste0("Peptide\n", total_pep), paste0("Substrate\n", total_sub) ),
    axLab.gpar = gpar(col = "black", fontsize = 14),
    axLab.pos = c(450, 450, 450),
    )
    grid.text(paste0("Total edges: ", total_edges), x = 1200, y = -2700, default.units = "native", gp = gpar(fontsize = 12))
    grid.text(paste0(protein, " - as kinase"), x = 0, y = 2000, default.units = "native", gp = gpar(fontsize = 20))
}

# Hive plot for protein acting as substrate only

myHive_sub <- function(protein, color, color2){
    #find id of your protein
    sub <- unique(tst3[tst3$Substrate==paste0("sub_", protein),6])
    query <- sub
    load("tst3.RData")
    tst <- tst3[(tst3$id %in% query | tst3$id.x %in% query),]
    total_kin <- length(unique(tst$Kinase))
    total_pep <- length(unique(tst$Peptide))-1
    total_sub <- length(unique(tst$Substrate))
    tmp1 <- tst[, 4:5]
    if(length(row.names(tmp1))>0){
        tmp1$color <- color
        names(tmp1) <- c("id1", "id2", "color")}
    tmp2 <- tst[, 5:6]
    if(length(row.names(tmp2))>0){
        tmp2$color <- color
        names(tmp2) <- c("id1", "id2", "color")}
    tmp3 <- tst[is.na(tst$Peptide),]
    tmp3 <- tmp3[, c(4,6)]
    if(length(row.names(tmp3))>0){
        tmp3$color <- color2
        names(tmp3) <- c("id1", "id2", "color")
    }
    edges <- rbind(tmp1, tmp2, tmp3)
    dead_edges <-  which(apply(edges, 1, function(x) sum(is.na(x))) > 0)
    if(length(dead_edges) > 0){
        edges <- edges[-dead_edges,]
    }
    edges <- data.frame(edges, weight = rep(1.5, nrow(edges)), stringsAsFactors = F)
    total_edges <- length(edges$id1)
    data <- list()
    data$nodes <- nodes
    data$edges <- edges
    data$type <- "2D"
    data$desc <- "proteins"
    data$axis.cols <- rep("#00000000", 3)
    attr(data, "class") <- "HivePlotData"
    
    # check if dataset is correct
    chkHPD(data, confirm=T)
    # plot
    plotHive(data, bkgnd = NULL, ch=100, axLabs=c(paste0("Kinase\n", total_kin), paste0("Peptide\n", total_pep), paste0("Substrate\n", total_sub) ),
    axLab.gpar = gpar(col = "black", fontsize = 14),
    axLab.pos = c(450, 450, 450),
    )
    grid.text(paste0("Total edges: ", total_edges), x = 1200, y = -2700, default.units = "native", gp = gpar(fontsize = 12))
    grid.text(paste0(protein, " - as substrate"), x = 0, y = 2000, default.units = "native", gp = gpar(fontsize = 20))
}


### EXAMPLE USE

myHive("AKT1", "grey", "yellow")
myHive_kin("AKT1", "grey", "yellow")
myHive_sub("AKT1", "grey", "yellow")


