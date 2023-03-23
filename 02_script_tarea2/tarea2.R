library(BoolNet)


#A)


redjoo <- loadNetwork("01_RawData/redjoo.txt")
redjoo
atrajoo <- getAttractors(redjoo)

plotAttractors(atrajoo)

plotStateGraph(atrajoo, drawLabels = T)


#######Discusion

########Grafica 2


immu <- loadNetwork("01_RawData/immuno.txt") #Agregar un renglon mas para el txt
immu
atraimm <- getAttractors(immu)

plotAttractors(atraimm)




###########STARWARS#############

library(igraph)

starwars <- read.csv("01_RawData/star-wars-network-edges.csv")
View(starwars)




starwarsmat <- as.matrix(starwars)

redstarwars <- graph_from_adj_list(starwarsmat, mode = "all")


plot(redstarwars)


