# method="GSC"; clustering.method="average"; plot=TRUE
tree.weight=function (cor.mat, method="GSC", clustering.method="average", plot=TRUE, orientation=c("vertical","horizontal"), label.cex=.3, ...) {
  
  orientation=match.arg(orientation)
  dist.mat=1-abs(cor.mat) # negative values need to be dealt with
  
  #    paid = calcPairwiseIdentity (getDat(dat.str), dissimilarity=TRUE)
  hclust1= hclust (as.dist(dist.mat, diag=FALSE, upper=FALSE), method=clustering.method)#average works, mcquitty works better, single works too
  merg=hclust1$merge
  height=hclust1$height
  n=length(height) # number of internal nodes, 1 less than sample size
  
  par(cex=label.cex)
  if (plot) plot(as.dendrogram(hclust1), horiz=orientation=="horizontal", xaxt="n", yaxt="n", ...)
  
  #    # get total weight
  #    total.weight=sum(height>38)
  
  # each element in nodes is the collection of sequence under that internal node
  nodes=list()
  for (i in 1:n) {
    members=c()
    for (j in 1:2) {
      item=merg[i,j]
      if (item<0) members=append (members, -item)
      else members=append (members, nodes[[item]])
    }
    nodes[[i]]=members
  }
  
  # get relative weights
  relative.weights=array (0, dim=n+1)
  for (i in 1:n) 
    if (height[i]!=0) 
      for (j in 1:2) {
        node=merg[i,j]
        if (node<0) 
          relative.weights[-node]=height[i]
        else {
          branch.weight = height[i]-height[node]
          proportion = relative.weights[ nodes[[node]] ]
          total=sum(proportion)
          if (total==0) proportion=rep(1/length(proportion),length(proportion))
          else proportion = proportion/total
          relative.weights[ nodes[[node]] ] = relative.weights[ nodes[[node]] ] + branch.weight * proportion
        }
      }
  relative.weights = relative.weights/sum(relative.weights) 
  #    weight = total.weight * relative.weights
  names(relative.weights)=hclust1$labels
  relative.weights
}
