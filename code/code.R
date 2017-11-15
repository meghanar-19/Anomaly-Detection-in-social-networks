library(igraph)

getwd()

normalize_vector <- function(x) {
  minx = min(x)
  print(minx)
  maxx = max(x)
  print(maxx)
  x <- (x - minx)/(maxx - minx)
}
 
get_graph <- function() {
  graph <- read.graph("./dataset/facebook_combined.txt", format="edgelist", directed=FALSE)
     
  dataset <- read.csv("./dataset/features.txt", header = FALSE, sep = ' ')
  num_features = ncol(dataset)
  V <- vcount(graph)
  E <- ecount(graph)
   
  features <- list()
  for(i in 1:V) {
    flag = 1
    for(j in 1:num_features) {
      if(dataset[i,j] == 1) {
           if(flag == 1) {
               features[[i]] <- c(j)
               flag = 0
              }
           else {
                features[[i]] <- c(features[[i]], j)
               }
           }
       }
   }
 
  ends_graph <- ends(graph, E(graph))
  for(i in 1:ecount(graph)) {
     E(graph)[i]$weight <- length(intersect(features[[ends_graph[i,1]]],features[[ends_graph[i,2]]]))
   }
  graph
}

get_outlierness_EDPL <- function(g) {
  V <- vcount(g)
  E <- ecount(g)
       
  # ego_vertexsize is an integer vector of size V, and stores no.of vertices in egonet of each node
  ego_vertexsize <- ego_size(g, 1)
  
  # ego_graph is a list of ego graphs of each node
  ego_graph <- make_ego_graph(g, 1)
  
  # ego_edgesize is an integer vector of size V, and stores no.of edges in egonet of each node
  ego_edgesize = rep(0, V)
  for (i in 1:V) {
    ego_edgesize[i] <- ecount(ego_graph[[i]])
  }
  
  #print(ego_vertexsize)
  #print(ego_edgesize)
  
  #plot(ego_vertexsize, ego_edgesize, type = "o", xlab = "Vertex Size", ylab = "Edge Size")
  C = 1
  alpha = 0.5
  #oulierness is calculated according to the Edge Density Power Law of the Oddball Algorithm
  outlierness = rep(0, V)
  for(i in 1:V) {
    outlierness[i] <- (max(ego_edgesize[i], C*(ego_vertexsize[i]**alpha)) / min(ego_edgesize[i], C*(ego_vertexsize[i]**alpha))) * (log10(abs(ego_edgesize[i] - C*(ego_vertexsize[i]**alpha)) + 1))
  }
  outlierness <- normalize_vector(outlierness)
  outlierness
}

get_outlierness_EWPL <- function(g) {
  V <- vcount(g)
  E <- ecount(g)
  
  # Ego_graph is a list of ego graphs of each node
  ego_graph <- make_ego_graph(g, 1)
  
  # ego_edgesize is an integer vector of size V, and stores no.of edges in egonet of each node
  ego_edgesize = rep(0, V)
  for (i in 1:V) {
    ego_edgesize[i] <- ecount(ego_graph[[i]])
  }
  
  # ego_weight is an integer vector of size V, and stores the sum of the weight of egonet of each node
  ego_weight = rep(0, V)
  for(i in 1:V) {
    ego_weight[i] <- sum(E(ego_graph[[i]])$weight) + 1
  } 
  
  C = 1
  beta = 1.15
  #oulierness is calculated according to the Edge Density Power Law of the Oddball Algorithm
  outlierness = rep(0, V)
  for(i in 1:V) {
    outlierness[i] <- (max(ego_weight[i], C*(ego_edgesize[i]**beta)) / min(ego_weight[i], C*(ego_edgesize[i]**beta))) * (log10(abs(ego_weight[i] - C*(ego_edgesize[i]**beta)) + 1))
  }
  outlierness <- normalize_vector(outlierness)
  outlierness
}

get_outlierness_ELWPL <- function(g) {
  V <- vcount(g)
  E <- ecount(g)
  
  # Ego_graph is a list of ego graphs of each node
  ego_graph <- make_ego_graph(g, 1)
  
  
  # ego_eigen_value is the eigen value of the weighted adjacency matrix of the egonet of each node
  ego_eigen_value = rep(0, V)
  for(i in 1:V) {
  A <- as_adjacency_matrix(ego_graph[[i]], attr="weight")
    ego_eigen_value[i] <- eigen(A)$values[1] + 1
  }
  
  # ego_weight is an integer vector of size V, and stores the sum of the weight of egonet of each node
  ego_weight = rep(0, V)
  for(i in 1:V) {
    ego_weight[i] <- sum(E(ego_graph[[i]])$weight) + 1
  }
  
  C = 1
  gamma = 0.75
  #oulierness is calculated according to the Edge Density Power Law of the Oddball Algorithm
  outlierness = rep(0, V)
  for(i in 1:V) {
  outlierness[i] <- (max(ego_eigen_value[i], C*(ego_weight[i]**gamma)) / min(ego_eigen_value[i], C*(ego_weight[i]**gamma))) * (log10(abs(ego_eigen_value[i] - C*(ego_weight[i]**gamma)) + 1))
  }
  outlierness <- normalize_vector(outlierness)
  outlierness
}

graph = get_graph()
#E(graph)$weight
outlierness_EDPL <- get_outlierness_EDPL(graph)
outlierness_EWPL <- get_outlierness_EWPL(graph)
outlierness_ELWPL <- get_outlierness_ELWPL(graph)
outlierness <- (outlierness_EDPL + outlierness_EWPL + outlierness_ELWPL)/3
par(mar = rep(2, 4))
plot(outlierness, type = "h", main = "Total Outlierness Score")

print(outlierness)

cnt <- 0
for(i in 1:vcount(graph)) {
  if (outlierness[i] >= 0.5) {
    print (c(outlierness[i], as.integer(i)))
    cnt <- cnt + 1;
  }
}

print(min(outlierness))
print(max(outlierness))
print(cnt)
