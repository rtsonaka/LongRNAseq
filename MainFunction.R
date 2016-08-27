
longRNAseq <- function(formula.lmer, data, Time, ID, gene.nams, dmat.voom){
  
  n.genes <- length(gene.nams)
  n.times <- length(unique(data[, Time]))
  data.time <- split(data, data[, Time])
  RNAseq.time <- split(data[, gene.nams], unique(data[, Time]))
  X.time <- split(as.data.frame(dmat.voom), data[, Time])
  voom.time <- lapply(seq_len(n.times), function(x){
    voom(t(RNAseq.time[[x]]), design = X.time[[x]])
  } )
  weights.time <- lapply(seq_len(n.times), function(x){
    data.frame(t(voom.time[[x]]$weights), id = data.time[[x]][, ID], time = data.time[[x]][, Time])
  })
  w <- do.call(rbind, weights.time)
  w <- w[order(w[, ID], w[, Time]), -(n.genes + 1:2)]
  
  RNAseq <- data[, gene.nams]
  test <- try(lmer(formula = as.formula(paste("log(RNAseq[, 1] + 0.5, base = 2)", formula.lmer, sep = "~")), 
                              data = data, REML = FALSE))
  n.betas <- length(fixef(test))
    
  coefs.lmer <- matrix(NA, nrow = n.genes, ncol = n.betas)
  var.lmer <- matrix(NA, nrow = n.genes, ncol = 1)
  pvals.lmer <- matrix(NA, nrow = n.genes, ncol = n.betas)
  rownames(pvals.lmer) <- rownames(coefs.lmer) <- rownames(var.lmer) <- gene.nams
  colnames(pvals.lmer) <- colnames(coefs.lmer) <- names(fixef(test))
  
  for(k in 1:n.genes){
    form <- as.formula(paste("log(RNAseq[, k] + 0.5, base = 2)", formula.lmer, sep = "~"))
    
    model.lmer.1 <- try(lmer(form, 
                             data = data, weights = w[, k], REML = FALSE))
    if(!inherits(model.lmer.1, "try-error")){
      s <- summary( model.lmer.1, ddf = "Satterthwaite" )
      coefs.lmer[k, ] <- fixef(model.lmer.1)
      var.lmer[k,] <- unlist(VarCorr(model.lmer.1))
      pvals.lmer[k, ] <- s$coefficients[ , "Pr(>|t|)" ]
    }
  }
  list(Gene = gene.nams, coefs = coefs.lmer, pvals = pvals.lmer, disp = var.lmer)
}




