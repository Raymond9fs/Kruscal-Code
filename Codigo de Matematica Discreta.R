# Cargar librerías necesarias
library(igraph)
library(leaflet)
library(geosphere)
library(optrees)


# Limpiar el ambiente
rm(list=ls())

# Definir funciones útiles
# Función para verificar si vectores están dentro de una secuencia de vectores
vectorsInVectorSeq <- function(testVectorMatrix, seqVectorMatrix) {
  r <- rep(FALSE,nrow(seqVectorMatrix))
  for (k in 1:nrow(testVectorMatrix)) {
    testElement = testVectorMatrix[k,]
    r <- r | apply(seqVectorMatrix, 1, function(x, test) isTRUE(all.equal(x, test)), testElement)
  }
  return(r)
}

# Función para crear una matriz de distancias
mat2list <- function(D) {
  n = dim(D)[1]
  k <- 1
  e <- matrix(ncol = 3,nrow = n*(n-1)/2)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      e[k,] = c(i,j,D[i,j])
      k<-k+1
    }
  }
  return(e)
}

# Coordenadas de algunas ciudades en Perú
lonlat <- read.table(textConnection(
  "lon, lat 
  -77.028240, -12.046374
  -71.537451, -16.409047
  -79.028240, -8.109051
  -79.930922, -6.771372
  -80.632691, -3.749229
  -73.247785, -3.743673
  -71.978694, -13.535761
  -75.208367, -11.123606
  -70.238679, -15.839819
  -74.535416, -13.163141
  "), header=TRUE, strip.white=TRUE, sep=",")

# Nombres de las ciudades
nname <- c("Lima", "Arequipa", "Trujillo", "Chiclayo", "Piura",
           "Iquitos", "Cusco", "Huancayo", "Tacna", "Pucallpa")

# Crear un data frame con IDs, nombres, y coordenadas
v <- data.frame(ids = 1:10, 
                name = nname, 
                x = lonlat$lon, 
                y = lonlat$lat)

# Mostrar las ubicaciones en un mapa interactivo usando leaflet
leaflet(data=v) %>% 
  addTiles() %>%
  addMarkers(lng=v$x, lat=v$y, popup=v$name) %>%
  addPolylines(lng=v$x, lat=v$y, popup=v$name)


#---------------------------------------------------------
#---------------------------------------------------------

leaflet(data=v[1:10,]) %>% addTiles() %>%
  addMarkers(lng=v$x, lat=v$y, popup=v$name) %>%
  addPolylines(lng=v$x, lat=v$y, popup=v$name)

#---------------------------------------------------------
#---------------------------------------------------------

# Calcular la matriz de distancias y crear el árbol de expansión mínima (MST)
D <- distm(lonlat, lonlat, fun=distVincentyEllipsoid) # en metros
eD <- mat2list(D/1000) # en kilómetros
eD <- eD[order(eD[,3]),] # ordenar aristas por distancia

# Suponiendo que eD está estructurada como V1, V2, V3
colnames(eD) <- c("ept1", "ept2", "weight")

# Redondear los valores de peso en la lista de aristas
eD[,3] <- round(eD[,3])

#---------------------------------------------------------
#---------------------------------------------------------

# Ordenar la tabla eD por la tercera columna (distancia)
eD <- eD[order(eD[,3]),]

# Crear un grafo con los datos de eD
net <- graph.data.frame(eD[,1:2], directed = FALSE)

# Asignar nombres a los vértices
V(net)$name <- v$name

# Asignar pesos a las aristas
E(net)$weight <- eD[,3]

# Obtener el MST usando el algoritmo de Prim
mst <- minimum.spanning.tree(net)

# Ordenar los vértices en un círculo para una mejor visualización
coords <- layout.circle(net)

# Graficar el grafo con el MST destacado en rojo y etiquetas de peso
plot(net, layout=coords, vertex.label=V(net)$name, vertex.size=5, edge.label=E(net)$weight,
     edge.width=1, edge.color=ifelse(E(net) %in% E(mst), "red", "grey"),
     main="Árbol de Expansión Mínima (MST)")

# Agregar el MST en rojo
plot(mst, layout=coords, add=TRUE, edge.color="red", edge.width=2)


#---------------------------------------------------------
#---------------------------------------------------------

# 4.) create a graph of MST
net <- graph.data.frame(eD[,1:2],directed = FALSE, vertices = v)
E(net)$weight <- eD[,3]
mstgraph <- minimum.spanning.tree(net)
par(mfrow=c(1,2), mar=c(0,1,0.75,0)) 
plot(net,vertex.label=nname)
plot(mstgraph,vertex.shape="none",edge.label=round(E(mstgraph)$weight))

#---------------------------------------------------------
#---------------------------------------------------------

# Graficar el MST resaltando las aristas del cluster
k <- 3
e.mst  = get.edges(mstgraph, 1:ecount(mstgraph)) # matriz de aristas del MST
e.clust = get.edges(mstgraph, 1:(ecount(mstgraph)-k+1)) # matriz de aristas del cluster k
clust_idx = vectorsInVectorSeq(e.clust, e.mst) # índices de filas del cluster en el MST

ecol <- rep("grey80", ecount(mstgraph)) # color predeterminado de las aristas
ecol[clust_idx] <- "green" # color para el cluster

ew <- rep(2, ecount(mstgraph)) # ancho predeterminado de las aristas
ew[clust_idx] <- 5 # ancho para las aristas del cluster

# Graficar el MST con colores y etiquetas
plot(mstgraph, vertex.shape="none", vertex.label.cex=0.85,
     edge.color=ecol, edge.width=ew, edge.label.cex=0.85,
     edge.label=round(E(mstgraph)$weight))

print("Fin del análisis")