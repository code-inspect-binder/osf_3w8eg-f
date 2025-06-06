
# Packages ----------------------------------------------------------------

library(lavaan)
library(semPlot)
library(qgraph)
library(huge)


# Extra function ----------------------------------------------------------

## MM to inchconversion for plot ##

MMtoInch <- function(MM) {
  z <- MM / 25.4
  return(z)
}


# Import datasets ---------------------------------------------------------

# SF-36
datafile_NKI_C <- read.delim("SF_36_NKI_CANCER_1.txt", na.strings = -9, header = TRUE)
datafile_NKI_H <- read.delim("SF_36_NKI_HEALTHY.txt", na.strings = -9, header = TRUE)
DS <- c(rep(0, times = nrow(datafile_NKI_C)), rep(1, times = nrow(datafile_NKI_H)))

data <- rbind(datafile_NKI_C, datafile_NKI_H)
data <- cbind(data, DS)

datafile_NKI_C <- data[1:485, 1:36] # delete the dummy variable
datafile_NKI_H <- data[486:2227, 1:36] # delete the dummy variable

#Subscale analysis

Sub_C <- read.delim("NKI_C_SubScale.txt", na.strings = -9, header = TRUE)
Sub_H <- read.delim("NKI_H_SubScale.txt", na.strings = -9, header = TRUE)
DS <- c(rep(0, times = nrow(Sub_C)), rep(1, times = nrow(Sub_H))) # group status

data <- rbind(Sub_C, Sub_H)
data <- cbind(data, DS)

Sub_C <- data[1:485, 1:8] # delete the dummy variable
Sub_H <- data[486:2227, 1:8] # delete the dummy variable

# Figure 1 ----------------------------------------------------------------


#reflective measurement model
model2 <- 'MH =~ NP + DC + CP + DB + HP'
data2 <- simulateData(model2)
fit2 <- lavaan::cfa(model2, data = data2, std.lv = TRUE)

#formative measurement model
model3 <- 'MH ~ NP + DC + CP + DB + HP'
data3 <- simulateData(model3)
fit3 <- lavaan::cfa(model3, data = data3)

#network model
model4 <- matrix(c(0, 0, 1, 0, 0,
                   0, 0, 0, 1, 1,
                   1, 0, 0, 1, 1,
                   0, 1, 1, 0, 0,
                   0, 1, 1, 0, 0), 5, 5, byrow = TRUE)

layout = matrix(c(1,2,
                  4,1,
                  3,3,
                  2,1,
                  5,2), 5, 2, byrow = TRUE)

h <- MMtoInch(234)
w <- MMtoInch(84)

pdf("Fig1.pdf", width = w, height = h)
layout(1:3)
#par(omd = c(.05, .95, .05, .95))
par("mar"=c(10, 5, 10, 0))

semPaths (fit3, title = FALSE, intercepts = FALSE, style = "lisrel", layout = "tree", residuals = TRUE, thresholds = FALSE, sizeMan = 25, sizeLat = 25, esize = 5, asize = 10, arrowAngle = pi/4, nCharNodes = 5, curve = 2.5)
mtext("(a)", cex = 1.25, side = 1, line = 10.5, las = 1)

semPaths (fit2, what = "path", title = FALSE, intercepts = FALSE, style = "lisrel", layout = "tree", residuals = TRUE, thresholds = FALSE, exoCov = TRUE, exoVar = FALSE, structural = FALSE, ask = FALSE, sizeMan = 25, sizeLat = 25, esize = 5, asize = 10, arrowAngle = pi/4, nCharNodes = 5)
mtext("(b)", cex = 1.25, side = 1, line = 10, las = 1)

qgraph(model4, labels = c("NP", "DC", "CP", "DB", "HP"), vsize = 25, esize = 5, asize = 7, mar = c(5.5, 5.5, 5.5, 5.5), arrowAngle = pi/4, layout = layout)
mtext("(c)", cex = 1.25, side = 1, line = 8, las = 1)
dev.off()


# Figure 2 ----------------------------------------------------------------

labels1 <- c("GH 01", "02", "PF 03", "PF 04", "PF 05", "PF 06", "PF 07", "PF 08", "PF 09", "PF 10", "PF 11", "PF 12", "RP 13", "RP 14", "RP 15", "RP 16", "RE 17", "RE 18", "RE 19", "SF 20", "BP 21", "BP 22", "VT 23", "MH 24", "MH 25", "MH 26", "VT 27", "MH 28", "VT 29", "MH 30", "VT 31", "SF 32", "GH 33", "GH 34", "GH 35", "GH 36", "DS")

labels2 <- labels1[-37]

names(data) <- labels1
names(datafile_NKI_C) <- labels2
names(datafile_NKI_H) <- labels2

groups1 <- list("General health" = c(1, 33, 34, 35, 36),
                "Physical functioning" = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                "Mental health" = c(24, 25, 26, 28, 30),
                "Role limitations physical" = c(13, 14, 15, 16),
                "Role limitations emotional" = c(17, 18, 19),
                "Bodily pain" = c(21, 22),
                "Social functioning" = c(20, 32),
                "Vitality" = c(23, 27, 29, 31),
                "No Dimension" = c(2, 37))

groups2 <- groups1
groups2$'No Dimension' <- groups2$'No Dimension'[-2]

## Order data as ordinal ##

for(i in 1:ncol(datafile_NKI_C))
{
  datafile_NKI_C[,i] <- as.ordered(datafile_NKI_C[,i])
}

for(i in 1:ncol(datafile_NKI_H))
{
  datafile_NKI_H[,i] <- as.ordered(datafile_NKI_H[,i])
}

for(i in 1:ncol(data))
{
  data[,i] <- as.ordered(data[,i])
}

NPN_graph_NKI_C <- cor_auto(datafile_NKI_C)
NPN_graph_NKI_H <- cor_auto(datafile_NKI_H)
NPN_graph_NKI <- cor_auto(data)

row.names(NPN_graph_NKI_C) <- labels2
colnames(NPN_graph_NKI_C) <- labels2

row.names(NPN_graph_NKI_H) <- labels2
colnames(NPN_graph_NKI_H) <- labels2

Layout_Graph <- qgraph(NPN_graph_NKI, graph = "glasso", layout = "spring", sampleSize = nrow(data))
L1 <- averageLayout(Layout_Graph)
L2 <- L1[1:36,]

h <- MMtoInch(234)
w <- MMtoInch(84)

pdf("Fig2.pdf", width = w, height = h)
layout(1:3)
par("mar"=c(10, 2, 10, 2))


NPN_Graph_NKI_C <- qgraph(NPN_graph_NKI_C, graph = "glasso", layout = L2, sampleSize = nrow(datafile_NKI_C), esize = 20, cut = 0.1, minimum = 0, maximum = 1, groups = groups2, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels2, legend = FALSE)
mtext("(a)", cex = 1.25, side = 1, line = 8.7, las = 1)

NPN_Graph_NKI_H <- qgraph(NPN_graph_NKI_H, graph = "glasso", layout = L2, sampleSize = nrow(datafile_NKI_H), esize = 20, cut = 0.1, minimum = 0, maximum = 1, groups = groups2, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels2, legend = FALSE)
mtext("(b)", cex = 1.25, side = 1, line = 8.7, las = 1)

NPN_Graph_NKI <- qgraph(NPN_graph_NKI, graph = "glasso", layout = L1, sampleSize = nrow(data), esize = 20, cut = 0.1, minimum= 0, maximum = 1, groups = groups1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels1, legend = FALSE)
mtext("(c)", cex = 1.25, side = 1, line = 8.7, las = 1)

dev.off()


# Figure 3 ----------------------------------------------------------------

labels <- c(names(Sub_C), "DS")

graph_C <- cor_auto(Sub_C)
graph_H <- cor_auto(Sub_H)
graph <- cor_auto(data)

Layout_Graph <- qgraph(graph, graph = "glasso", layout = "spring", sampleSize = nrow(data))

L1 <- averageLayout(Layout_Graph)
L2 <- L1[1:8,]

h <- MMtoInch(234)
w <- MMtoInch(84)

colours <- c("yellow", "cornflowerblue", "purple", "grey", "orange", "green", "maroon1", "red", "chocolate3")

pdf("Fig3.pdf", width = w, height = h)
layout(1:3)
par("mar"=c(10, 2, 10, 2))

NPN_Graph_NKI_C <- qgraph(graph_C, graph = "glasso", layout = L2, sampleSize = nrow(Sub_C), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = colours[-9], borders = FALSE, vsize= 10, legend = FALSE)
mtext("(a)", cex = 1.25, side = 1, line = 8.7, las = 1)

NPN_Graph_NKI_H <-qgraph(graph_H, graph = "glasso", layout = L2, sampleSize = nrow(Sub_H), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = colours[-9], borders = FALSE, vsize= 10, legend = FALSE)
mtext("(b)", cex = 1.25, side = 1, line = 8.7, las = 1)

NPN_Graph_NKI <- qgraph(graph, graph = "glasso", layout = L1, sampleSize = nrow(data), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = colours, borders = FALSE, vsize= 10, legend = FALSE)
mtext("(c)", cex = 1.25, side = 1, line = 8.7, las = 1)

dev.off()


# Figure 4 ----------------------------------------------------------------

labels1 <- c("Item 01", "Item 02", "Item 03a", "Item 03b", "Item 03c", "Item 03d", "Item 03e", "Item 03f", "Item 03g", "Item 03h", "Item 03i", "Item 03j", "Item 04a", "Item 04b", "Item 04c", "Item 04d", "Item 05a", "Item 05b", "Item 05c", "Item 06", "Item 07", "Item 08", "Item 09a", "Item 09b", "Item 09c", "Item 09d", "Item 09e", "Item 09f", "Item 09g", "Item 09h", "Item 09i", "Item 10", "Item 11a", "Item 11b", "Item 11c", "Item 11d", "HS")

labels2 <- labels1[-37]

names(data) <- labels1
names(datafile_NKI_C) <- labels2
names(datafile_NKI_H) <- labels2

groups1 <- list("General health" = c(1, 33, 34, 35, 36),
                "Physical functioning" = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                "Mental health" = c(24, 25, 26, 28, 30),
                "Role limitations physical" = c(13, 14, 15, 16),
                "Role limitations emotional" = c(17, 18, 19),
                "Bodily pain" = c(21, 22),
                "Social functioning" = c(20, 32),
                "Vitality" = c(23, 27, 29, 31),
                "No Dimension" = c(2, 37))

groups2 <- groups1
groups2[[9]] <- groups1[[9]][-2]

## Order data as ordinal ##

for(i in 1:ncol(datafile_NKI_C))
{
  datafile_NKI_C[,i] <- as.ordered(datafile_NKI_C[,i])
}

for(i in 1:ncol(datafile_NKI_H))
{
  datafile_NKI_H[,i] <- as.ordered(datafile_NKI_H[,i])
}

for(i in 1:ncol(data))
{
  data[,i] <- as.ordered(data[,i])
}

NPN_graph_NKI_C <- cor_auto(datafile_NKI_C)
NPN_graph_NKI_H <- cor_auto(datafile_NKI_H)
NPN_graph_NKI <- cor_auto(data)


NPN_graph_NKI_C <- rbind(NPN_graph_NKI_C,rep(0, times = nrow(NPN_graph_NKI_C)))
NPN_graph_NKI_C <- cbind(NPN_graph_NKI_C,rep(0, times = nrow(NPN_graph_NKI_C)))
NPN_graph_NKI_C[37,37] <- 1
row.names(NPN_graph_NKI_C) <- labels1
colnames(NPN_graph_NKI_C) <- labels1

NPN_graph_NKI_H <- rbind(NPN_graph_NKI_H,rep(0, times = nrow(NPN_graph_NKI_H)))
NPN_graph_NKI_H <- cbind(NPN_graph_NKI_H,rep(0, times = nrow(NPN_graph_NKI_H)))
NPN_graph_NKI_H[37,37] <- 1
row.names(NPN_graph_NKI_H) <- labels1
colnames(NPN_graph_NKI_H) <- labels1

Layout_Graph <- qgraph(NPN_graph_NKI, graph = "glasso", layout = "spring", sampleSize = nrow(data), DoNotPlot = TRUE)

L <- averageLayout(Layout_Graph)

NPN_Graph_NKI_C <- qgraph(NPN_graph_NKI_C, graph = "glasso", layout = L, sampleSize = nrow(datafile_NKI_C), esize = 10, cut = 0.1, minimum= 0, maximum = 1, groups = groups1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels1, legend = FALSE)

NPN_Graph_NKI_H <- qgraph(NPN_graph_NKI_H, graph = "glasso", layout = L, sampleSize = nrow(datafile_NKI_H), esize = 10, cut = 0.1, minimum= 0, maximum = 1, groups = groups1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels1, legend = FALSE)

NPN_Graph_NKI <- qgraph(NPN_graph_NKI, graph = "glasso", layout = L, sampleSize = nrow(data), esize = 10, cut = 0.1, minimum= 0, maximum = 1, groups = groups1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, labels = labels1, legend = FALSE)

graphs <- list(NPN_Graph_NKI_C, NPN_Graph_NKI_H, NPN_Graph_NKI)
names(graphs) <- c("Cancer patient sample", "National sample", "Combined sample")
Long <- centralityTable(graphs, standardized=TRUE, labels=labels1, relative=FALSE)

Long <- subset(Long, measure %in% "Closeness")

# Ordering by node name to make nice paths:
Long <- Long[gtools::mixedorder(Long$node),] 
Long$node <- factor(as.character(Long$node), levels = unique(gtools::mixedsort(as.character(Long$node))))

table <- data.frame("node" = Long$node, "type" = Long$type, "value" = round(Long$value, 3))

table <- reshape(table, direction = "wide", idvar = "node", timevar = "type")

# find highest and lowest centrality measures
table[order(table[,2]),c(1,2)]
table[order(table[,3]),c(1,3)]
table[order(table[,4]),c(1,4)]

h <- MMtoInch(84)
w <- MMtoInch(174)

pdf("Fig4.pdf", height = h, width = w)
g <- ggplot(Long, 
            aes(x = value, 
                y = node, 
                group = type, 
                colour = type))  +  
  geom_path() +  
  xlab("") + 
  ylab("") + 
  geom_point() + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1, family = "ArialMT"), 
    axis.text.y = element_text(size = 7, family = "ArialMT"), 
    legend.text =  element_text(size = 7, family = "ArialMT"),
    legend.position=c(0.66,0.8),
    legend.background = element_rect(fill=NA),
    legend.title=element_blank()
  ) +  
  scale_x_continuous(breaks = (-10):10, minor = seq(-10,10,by=0.5)) + 
  coord_flip()

print(g)

dev.off()


# Figure 5 ----------------------------------------------------------------

labels <- c(names(Sub_C), "DS")

graph_C <- cor_auto(Sub_C)
graph_H <- cor_auto(Sub_H)
graph <- cor_auto(data)

graph_C <- rbind(graph_C,rep(0, times = nrow(graph_C)))
graph_C <- cbind(graph_C,rep(0, times = nrow(graph_C)))
graph_C[9,9] <- 1
row.names(graph_C)[9] <- "DS"
colnames(graph_C)[9] <- "DS"

graph_H <- rbind(graph_H,rep(0, times = nrow(graph_H)))
graph_H <- cbind(graph_H,rep(0, times = nrow(graph_H)))
graph_H[9,9] <- 1
row.names(graph_H)[9] <- "DS"
colnames(graph_H)[9] <- "DS"

Layout_Graph <- qgraph(graph, graph = "glasso", layout = "spring", sampleSize = nrow(data))

NPN_Graph_NKI_C <- qgraph(graph_C, graph = "glasso", layout = Layout_Graph, sampleSize = nrow(Sub_C), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, legend = FALSE)

NPN_Graph_NKI_H <-qgraph(graph_H, graph = "glasso", layout = Layout_Graph, sampleSize = nrow(Sub_H), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, legend = FALSE)

NPN_Graph_NKI <- qgraph(graph, graph = "glasso", layout = Layout_Graph, sampleSize = nrow(data), esize = 20, cut = 0.1, minimum = 0, maximum = 1, color = c("red", "yellow", "orange", "cornflowerblue", "green", "purple", "grey", "maroon1", "chocolate3"), borders = FALSE, vsize= 10, legend = FALSE)

graphs <- list(NPN_Graph_NKI_C, NPN_Graph_NKI_H, NPN_Graph_NKI)
names(graphs) <- c("Cancer patient sample", "National sample", "Combined sample")
Long <- centralityTable(graphs, standardized=TRUE, labels=labels, relative=FALSE)

Long <- subset(Long, measure %in% "Closeness")

# Ordering by node name to make nice paths:
Long <- Long[gtools::mixedorder(Long$node),] 
Long$node <- factor(as.character(Long$node), levels = unique(gtools::mixedsort(as.character(Long$node))))
Long <- rbind(Long[4:6,], Long[1:3,], Long[7:27,])
table <- data.frame("node" = Long$node, "type" = Long$type, "value" = round(Long$value, 3))
table <- reshape(table, direction = "wide", idvar = "node", timevar = "type")

# find highest and lowest centrality measures
table[order(table[,2]),c(1,2)]
table[order(table[,3]),c(1,3)]
table[order(table[,4]),c(1,4)]

h <- MMtoInch(84)
w <- MMtoInch(174)

pdf("Fig5.pdf", height = h, width = w)

g <- ggplot(Long, aes(x = value, 
                      y = node, 
                      group = type, 
                      colour = type))  +  
  geom_path() +  
  xlab("") + 
  ylab("") + 
  geom_point() + 
  theme_bw() 


g$data$node <- c(g$data$node[4:6], g$data$node[1:3], g$data$node[7:27])

g <- g +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1, family = "ArialMT"), 
    axis.text.y = element_text(size = 7, family = "ArialMT"), 
    legend.text =  element_text(size = 7, family = "ArialMT"),
    legend.position=c(0.66,0.2),
    legend.background = element_rect(fill=NA),
    legend.title=element_blank()
  ) + scale_x_continuous(breaks = (-10):10, minor = seq(-10,10,by=0.5)) +
  scale_y_discrete(limits = c("DS", "BP", "GH", "MH", "PF", "RE", "RP", "SF", "VT")) + coord_flip()

print(g)

dev.off()


