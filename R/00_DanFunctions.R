
# 00_DanFunctions.R ####

## function for I2 for rma.mv
i2 = function(model){
  
  ## metafor site code for I2
  W = diag(1/model$vi)
  X = model.matrix(model)
  P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 = 100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2 = round(I2,2)
  
  ## summarize by each variance component
  allI2 = 100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2 = round(allI2,3)
  return(list(I2 = I2,allI2 = allI2))
}

## add R2

r2f = function(x){
  mod = x
  base = update(mod,mods = ~1)
  r2 = (sum(base$sigma2) - sum(mod$sigma2)) / sum(base$sigma2)
  r2 = ifelse(r2<0,0,r2)
  return(r2)
}
