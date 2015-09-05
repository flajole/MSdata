PCA.Anal <- function(){
    pca <- prcomp(dataSet$norm, center=F, scale=F);
    
    # obtain variance explained
    sum.pca <- summary(pca);
    imp.pca <- sum.pca$importance;
    std.pca <- imp.pca[1,]; # standard devietation
    var.pca <- imp.pca[2,]; # variance explained by each PC
    cum.pca <- imp.pca[3,]; # cummulated variance explained
    
    # store the item to the pca object
    analSet$pca <<- append(pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
    
    write.csv(signif(analSet$pca$x, 5), file="pca_score.csv");
    write.csv(signif(analSet$pca$rotation, 5), file="pca_loadings.csv");
}

