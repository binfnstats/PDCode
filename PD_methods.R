# Required Packages

library('sigFeature')
library('caret')
library('tictoc')
library('Boruta')

##########

check_data = function(x, 
                      y, 
                      B, 
                      k, 
                      repB, 
                      setSeed = 1234567, 
                      nCores = parallel::detectCores() - 1, 
                      methods, topN = 50, 
                      compareMethod = "BIC"){
  
  p = ncol(x)
  df = base::data.frame(x, y)
  nMethods = length(methods)
  kFolds = list()
  
  output = list()
  
  # K-folds
  set.seed(1234)
  
  for(i in 1:repB){
    kFold = caret::createFolds(y, 
                        k = k, 
                        list = TRUE, 
                        returnTrain = FALSE)
    kFolds[[i]] = kFold
  }
  
  saveRDS(kFolds,paste0("folds_seed_",setSeed,".rds",collapse=""))
  
  
  for(l in 1:nMethods){
    accVals = matrix(data = 0,
                     nrow = k,
                     ncol = repB)
    aucVals = matrix(data = 0,
                     nrow = k,
                     ncol = repB)
    timeVals = matrix(data = 0,
                      nrow = k,
                      ncol = repB)
    features = list()
    
    for(i in 1:repB){
      kFold = kFolds[[i]]
      
      features2 = list()
      
      for(j in 1:k){
        trainTest = unlist(kFolds[[i]][-j])
        xTrain = x[trainTest,]
        yTrain = y[trainTest]
        xTest = x[-trainTest,]
        yTest = y[-trainTest]
        preProcValues = caret::preProcess(xTrain, method = c("center", "scale"))
        xTrain = stats::predict(preProcValues, xTrain)
        xTest = stats::predict(preProcValues, xTest)
        
        if(methods[l] == "vivid"){
          tic()
          vivid = VIVID::vivid(x = xTrain,
                               y = yTrain,
                               bootstraps = B,
                               cores = nCores,
                               seed = setSeed,
                               compareMethod = compareMethod)
          time = toc()
          saveRDS(vivid,paste0(methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,unlist(vivid$optModel) == 1]
          vividGLM = glmnet::cv.glmnet(newX, 
                                       yTrain, 
                                       alpha = 0, 
                                       family = "binomial")
          newXTest = xTest[,unlist(vivid$optModel) == 1]
          vividPred = predict(vividGLM,
                              s = "lambda.1se",
                              newx = newXTest)
          fitted = c(exp(vividPred)/(1+exp(vividPred)))
          features2[[j]] = vivid$optFeatures
        }
        
        if(methods[l] == "rf"){
          tic()
          rf = randomForest::randomForest(x = xTrain,
                                          y = yTrain,
                                          ntree = 500)
          time = toc()
          saveRDS(rf,paste0(methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          fitted = predict(rf, 
                           newdata = xTest, 
                           type = "prob")[,2]
        }
        
        if(methods[l] == "boruta"){
          tic()
          boruta = Boruta::Boruta(x = xTrain,
                                  y = yTrain,
                                  doTrace = 0,
                                  maxRuns = 500)
          time = toc()
          
          saveRDS(boruta,paste0(methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,boruta$finalDecision == "Confirmed"]
          borutaGLM = glmnet::cv.glmnet(newX, 
                                        yTrain, 
                                        alpha = 0, 
                                        family = "binomial")
          newXTest = xTest[,boruta$finalDecision == "Confirmed"]
          borutaPred = predict(borutaGLM,
                               s = "lambda.1se",
                               newx = newXTest)
          fitted = c(exp(borutaPred)/(1+exp(borutaPred)))
          features2[[j]] = colnames(x)[which(boruta$finalDecision == "Confirmed")]
        }
        
        if(methods[l] == "RFE"){
          tic()
          RFE = sigFeature::sigFeature(X = xTrain,
                                       Y = yTrain)
          time = toc()
          saveRDS(RFE,paste0(methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,RFE > (p - topN)]
          RFEGLM = glmnet::cv.glmnet(newX, 
                                     yTrain, 
                                     alpha = 0, 
                                     family = "binomial")
          newXTest = xTest[,RFE > (p - topN)]
          RFEPred = predict(RFEGLM,
                            s = "lambda.1se",
                            newx = newXTest)
          fitted = c(exp(RFEPred)/(1+exp(RFEPred)))
          features2[[j]] = base::colnames(x)[RFE > (p - topN)]
          
        }
        
        modelPred = 1*(fitted > 0.5)
        binary = 1*(yTest == levels(y)[2])
        accVals[j, i] = mean(1*(modelPred == binary))
        roc = pROC::roc(binary, fitted)
        aucVals[j, i] = roc$auc
        timeVals[j, i] = time$toc - time$tic
        
      }
      features[[i]] = features2
    }
    
    output[[l]] = list(
      acc = accVals,
      auc = aucVals,
      time = timeVals,
      features = features
    )
  }
  names(output) = methods
  
  return(output)
  
}