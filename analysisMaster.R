#### File Information #####################################################################################################################
#File Name: analysisMaster.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file can act as a master file to analyze a given dataset with all modeling techniques and evaluation metrics.

### Functions #############################################################################################################################

## Function 1: analysisMaster = function(survivalDataset, numberOfFolds =5,
#                                        CoxKP = T,CoxKPEN = T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T, GBMModel = T, BayesianNetModel = T,#Models
#                                        DCal = T, OneCal = T, Concor = T, L1Measure = T, BrierInt = T, BrierSingle = T, #Evaluations
#                                        DCalBins = 10, OneCalTime = NULL,  concordanceTies = "Risk", #Evaluation args
#                                        SingleBrierTime = NULL, IntegratedBrierTimes = NULL, numBrierPoints = 1000, Ltype = "Margin", 
#                                        Llog = F, typeOneCal = "DN", oneCalBuckets = 10, survivalPredictionMethod = "Median", 
#                                        AFTDistribution = "weibull", #Model args,
#                                        FS = T, imputeZero=T, verbose = T, # Misc args
#                                        foldIndex = NULL, useAllData=F)

#Inputs: 
# survivalDataset - This is the dataset one wishes to analyze. This must include 'time', 'delta', and at least 1 more feature. No default.
# numberOfFolds - The number of desired cross-validation folds. Default is 5.
# CoxKP, CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel: Booleans specifying whether or not to run that model. Default is TRUE.
# DCal, OneCal, Concor, L1Measure, BrierSingle, BrierInt: Booleans specifying whether or not to run that evaluation metric. Default is TRUE.
# DCalBins: Number of bins for D-Calibration. Default is 10.
# OneCalTime: An int specifying the time to evaluate 1-Calibration. If left as NULL but OneCal = TRUE, then the 10th, 25th, 50th, 75th,
#             and 90th percentiless of all event times are used. Default is NULL.
# concordanceTies: A string ("None", "Time", "Risk","All") indicating how to handle ties in concordance. Default is "Risk".
# SingleBrierTime: The time to evaluate the Brier Score. If left as null, the 50th percentile of all event times is used. Default is NULL.
# IntegratedBrierTimes: A 2 length vector (e.g. c(0,100)) specifying the lower and upper bounds on the integrated Brier score. If NULL then
#                 the default is 0 as a lower bound and the max event time of the entire dataset is used as an upper bound. Default is NULL.
# numBrierPoints: The number of points to evaluate the integrated Brier score. A simple trapezoidal numerical approximation is used. Default
#                 is 1000 points.
# Ltype: The type of L1-loss. Must be one of "Uncensored","Hinge", or "Margin". Default is "Margin".
# Llog: A boolean specifying whether or not to use log-L1 metric. Default is FALSE.
# typeOneCal: A string indicating the type of 1-Calibrtion to use. Must be one of "DN" or "Uncensored". Default is "DN".
# oneCalBuckets: An int specifying number of bins for 1-Calibration. Default is 10.
# survivalPredictionMethod: The way in which to estimate average survival times. Must be one of "Mean" or "Median". Default is "Median".
# AFTDistribution: The distribution to use for AFT, default is "weibull". Must be one of "weibull","exponential","lognormal","gaussian",
#                   "loglogistic","logistic".
# FS: A boolean specifying whether or not to use feature selection. Default is TRUE.
# imputeZero: A boolean specifying whether 0 valued times should be imputed (AFT breaks for 0 valued times). If TRUE then 0 valued times are
# imputed to half the minimum non-zero time. Default is TRUE. 
# verbose: A boolean specifying whether or not to return results and progress information.
# foldIndex: Define each cross-validation fold by predifined fold index.
# useAllData: Use all the data for training. For drawing survival curves for exaplain the survival model. 


#Output: A list of (3) items:
#(1) datasetUsed: This is the dataset that is actually used post feature selection but pre normalization and imputation. datasetUsed
#will have all the patients who had acceptable time and delta values and the features that were selected.
#(2) survivalCurves: This is a list containing the survival curves for all patients for each model that was tested. 
#(3) results: This is a dataframe containing all the evaluation results with specified model and fold number. Additionally the sample size
#feature size, and censoring percnetage are returned. Notice that the feature sizes before and after one hot encoding are returned. 
#If none of the features were factors then NumFeatures should equal NumFeaturesOneHot.

#Note that survivalCurves can be plotted by plotSurvivalCurves().

## Function 2: getSurvivalCurves()

# coxTimes, coxENTimes, kmTimes, aftTimes, rsfTimes, mtlrTimes, gbmTimes, bayesianNetTimes - The times used for prediction of each model.
# CoxKP = T,CoxKPEN=T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T: The models used in analysisMaster.
# combinedTestResults: A List containing all model survival curves. 
# numberOfFolds: Number of folds for cross validation.
# originalIndexing: The original indexing prior to CV folds.

#Output: The survival curves of all survival models for all test patients.

#Usage: This is a helper function for analysisMaster(). This is used to get the survival curves for each model for each patient.

### Code ##################################################################################################################################
#Data processing files:
source("ValidateCleanCV/validateAndClean.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")

#Modeling files:
source("Models/CoxPH_KP.R")
source("Models/KaplanMeier.R")
source("Models/RandomSurvivalForests.R")
source("Models/AcceleratedFailureTime.R")
source("Models/MTLR.R")
source("Models/BayesianNetUpper.R")
source("Models/GBMCox.R")

#Evaluation files:
source("Evaluations/DCalibration.R")
source("Evaluations/OneCalibration.R")
source("Evaluations/Concordance.R")
source("Evaluations/L1Measures.R")
source("Evaluations/BrierScore.R")


#Misc files:
source("FeatureSelection/FeatureSelection.R")
source("Plotting/plotSurvivalCurves.R")

analysisMaster = function(survivalDataset, numberOfFolds =5, BayesianC1 = NULL,
                          CoxKP = T,CoxKPEN = T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T, GBMModel = T, BayesianNetModel = T,#Models
                          DCal = T, OneCal = T, Concor = T,ConcorCurve = T, L1Measure = T, BrierInt = T, BrierSingle = T, #Evaluations
                          DCalBins = 10, OneCalTime = NULL,  concordanceTies = "All", #Evaluation args
                          SingleBrierTime = NULL, IntegratedBrierTimes = NULL, numBrierPoints = 1000, Ltype = "Margin", #Evaluation args
                          Llog = F, typeOneCal = "DN", oneCalBuckets = 10, survivalPredictionMethod = "Median", #Evaluation args
                          AFTDistribution = "weibull", #Model args,
                          FS = T, imputeZero=T, verbose = T, # Misc args
                          foldIndex = NULL, useAllData=F
                          ){
  validatedData = validateAndClean(survivalDataset, imputeZero)
  if(FS) {validatedData = FeatureSelection(validatedData, type = "UniCox")}
  if(is.null(foldIndex)) {foldsAndNormalizedData = createFoldsAndNormalize(validatedData, numberOfFolds)}
  else if(!is.null(foldIndex)) {foldsAndNormalizedData = createFoldsAndNormalize(validatedData, numberOfFolds, T, foldIndex)}
  originalIndexing = foldsAndNormalizedData[[1]]
  normalizedData = foldsAndNormalizedData[[2]]
  evaluationResults = data.frame()
  combinedTestResults = list(Cox = list(),CoxEN = list(), KM = list(), AFT = list(), RSF = list(), MTLR = list(), GBM = list(), BayesianNet = list())
  coxTimes = NULL;coxENTimes = NULL; kmTimes = NULL; rsfTimes = NULL; aftTimes = NULL; mtlrTimes = NULL; gbmTimes = NULL; bayesianNetTimes = NULL;
  ConcordanceCurve = NULL
  for(i in 1:numberOfFolds){
    if(verbose){
      print(Sys.time())
      print(paste("Starting fold",i,"of", numberOfFolds, "total folds."))
    }
    #Models - We evaluate values to NULL so we can pass them to evaluations, regardless if the models were ran or not.
    coxMod = NULL;coxENMod =NULL; kmMod = NULL; rsfMod = NULL; aftMod = NULL; mtlrMod = NULL; gbmMod = NULL; bayesianNetMod = NULL;
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]
    if(useAllData) {training = rbind(training,testing)}
    if(verbose){
      print(paste("Beginning model training."))
    }
    if(CoxKP){
      if(verbose){
        print("Starting Cox Proportional Hazards.")
      }
      coxMod = CoxPH_KP(training, testing)
      if(length(coxMod) ==1){
        combinedTestResults$Cox = list()
        coxTimes = NULL
        CoxKP = F
        if(i > 1)
          evaluationResults = with(evaluationResults,evaluationResults[-which(Model == "CoxKP"),])
      }
      else{
        combinedTestResults$Cox[[i]] = coxMod
        coxTimes = c(coxTimes,coxMod[[1]]$time)
      }
    }
    if(CoxKPEN){
      if(verbose){
        print("Starting Cox Proportional Hazards - Elastic Net.")
      }
      coxENMod = CoxPH_KP(training, testing,ElasticNet = T)
      combinedTestResults$CoxEN[[i]] = coxENMod
      coxENTimes = c(coxENTimes,coxENMod[[1]]$time)
    }
    if(KaplanMeier){
      if(verbose){
        print("Starting Kaplan Meier.")
      }
      kmMod = KM(training, testing)
      combinedTestResults$KM[[i]] = kmMod
      kmTimes = c(kmTimes,kmMod[[1]]$time)
    }
    if(RSFModel){
      if(verbose){
        print("Starting Random Survival Forests.")
      }
      rsfMod = RSF(training, testing)
      combinedTestResults$RSF[[i]] = rsfMod
      rsfTimes = c(rsfTimes,rsfMod[[1]]$time)
    }
    if(AFTModel){
      if(verbose){
        print("Starting Accelerated Failure Time.")
      }
      aftMod = AFT(training, testing, AFTDistribution)
      if(length(aftMod)==1){
        combinedTestResults$AFT = list()
        aftTimes = NULL
        AFTModel = F
        if(i >1)
          evaluationResults = with(evaluationResults,evaluationResults[-which(Model == "AFT"),])
      }
      else{
        combinedTestResults$AFT[[i]] = aftMod
        aftTimes = c(aftTimes,aftMod[[1]]$time)
      }
    }
    if(MTLRModel){
      if(verbose){
        print("Starting Multi-task Logistic Regression (PSSP).")
      }
      mtlrMod = MTLR(training, testing)
      combinedTestResults$MTLR[[i]] = mtlrMod
      mtlrTimes = c(mtlrTimes,mtlrMod[[1]]$time)
    }
    if(GBMModel){
      if(verbose){
        print("Starting GBM Cox.")
      }
      gbmMod = GBMCox_KP(training, testing)
      combinedTestResults$GBM[[i]] = gbmMod
      gbmTimes = c(gbmTimes,gbmMod[[1]]$time)
    }
    if(BayesianNetModel){
      if(verbose){
        print("Starting Bayesian Net.")
      }
      bayesianNetMod = BayesianNetUpper(training, testing)
      combinedTestResults$BayesianNet[[i]] = bayesianNetMod
      bayesianNetTimes = c(bayesianNetTimes,bayesianNetMod[[1]]$time)
    }
    #Evaluations - Note that if evaluations are passed a NULL value they return a NULL.
    DCalResults = NULL;OneCalResults = NULL;ConcordanceResults = NULL;
    BrierResultsInt = NULL;BrierResultsSingle = NULL;L1Results = NULL; L2Results = NULL; 
    if(Concor){
      if(verbose){
        print("Staring Evaluation: Concordance")
      }
      coxConc = Concordance(coxMod, concordanceTies,survivalPredictionMethod)
      coxENConc = Concordance(coxENMod, concordanceTies,survivalPredictionMethod)
      kmConc = Concordance(kmMod, concordanceTies,survivalPredictionMethod)
      rsfConc = Concordance(rsfMod, concordanceTies,survivalPredictionMethod)
      aftConc = Concordance(aftMod, concordanceTies,survivalPredictionMethod)
      mtlrConc = Concordance(mtlrMod, concordanceTies,survivalPredictionMethod)
      gbmConc = Concordance(gbmMod, concordanceTies,survivalPredictionMethod)
      bayesConc = Concordance(bayesianNetMod, concordanceTies,survivalPredictionMethod)
      ConcordanceResults = rbind(coxConc,coxENConc, kmConc, rsfConc, aftConc, mtlrConc, gbmConc, bayesConc)
    }
    if(F){
      if(verbose){
        print("Staring Evaluation: Concordance Curve")
      }
      coxConCurve = Concordance(coxMod, concordanceTies,survivalPredictionMethod)
      coxENConCurve = Concordance(coxENMod, concordanceTies,survivalPredictionMethod)
      kmConCurve = Concordance(kmMod, concordanceTies,survivalPredictionMethod)
      rsfConCurve = Concordance(rsfMod, concordanceTies,survivalPredictionMethod)
      aftConCurve = Concordance(aftMod, concordanceTies,survivalPredictionMethod)
      mtlrConCurve = Concordance(mtlrMod, concordanceTies,survivalPredictionMethod)
      gbmConCurve = Concordance(gbmMod, concordanceTies,survivalPredictionMethod)
      bayesConCurve = Concordance(bayesianNetMod, concordanceTies,survivalPredictionMethod)
      
      ConcordanceCurveResults = rbind(coxConCurve, coxENConCurve, kmConCurve, rsfConCurve, aftConCurve, mtlrConCurve, gbmConCurve, bayesConCurve)
      if(is.null(ConcordanceCurve)) {ConcordanceCurve = rbind(timeOfInterest,ConcordanceCurveResults)}
      else{ConcordanceCurve = rbind(ConcordanceCurve,ConcordanceCurveResults)}
      print(ConcordanceCurve)
      
    }
    if(BrierInt){
      if(verbose){
        print("Staring Evaluation: Brier Score- Integrated")
      }
      coxBrierInt = BrierScore(coxMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      coxENBrierInt = BrierScore(coxENMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      kmBrierInt = BrierScore(kmMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      rsfBrierInt = BrierScore(rsfMod, type = "Integrated",numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      aftBrierInt = BrierScore(aftMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      mtlrBrierInt = BrierScore(mtlrMod, type = "Integrated", numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      gbmBrierInt = BrierScore(gbmMod, type = "Integrated", numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      bayesBrierInt = BrierScore(bayesianNetMod, type = "Integrated", numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      
      BrierResultsInt = rbind(coxBrierInt,coxENBrierInt, kmBrierInt, rsfBrierInt, aftBrierInt, mtlrBrierInt, gbmBrierInt, bayesBrierInt)
      
    }
    if(BrierSingle){
      if(verbose){
        print("Staring Evaluation: Brier Score - Single")
      }
      coxBrierSingle = BrierScore(coxMod, type = "Single", singleBrierTime =SingleBrierTime )
      coxENBrierSingle = BrierScore(coxENMod, type = "Single", singleBrierTime =SingleBrierTime )
      kmBrierSingle = BrierScore(kmMod, type = "Single", singleBrierTime =SingleBrierTime )
      rsfBrierSingle = BrierScore(rsfMod, type = "Single", singleBrierTime =SingleBrierTime )
      aftBrierSingle = BrierScore(aftMod, type = "Single", singleBrierTime =SingleBrierTime )
      mtlrBrierSingle = BrierScore(mtlrMod, type = "Single", singleBrierTime =SingleBrierTime )
      gbmBrierSingle = BrierScore(gbmMod, type = "Single", singleBrierTime =SingleBrierTime )
      bayesBrierSingle = BrierScore(bayesianNetMod, type = "Single", singleBrierTime =SingleBrierTime )
      
      BrierResultsSingle = rbind(coxBrierSingle,coxENBrierSingle, kmBrierSingle, rsfBrierSingle, aftBrierSingle, mtlrBrierSingle, gbmBrierSingle, bayesBrierSingle)
      
    }
    if(L1Measure){
      if(verbose){
        print("Staring Evaluation: L1 Loss")
      }
      coxL1 = L1(coxMod, Ltype, Llog,survivalPredictionMethod)
      coxENL1 = L1(coxENMod, Ltype, Llog,survivalPredictionMethod)
      kmL1 = L1(kmMod, Ltype, Llog,survivalPredictionMethod)
      rsfL1 = L1(rsfMod, Ltype, Llog,survivalPredictionMethod)
      aftL1 = L1(aftMod, Ltype, Llog,survivalPredictionMethod)
      mtlrL1 = L1(mtlrMod, Ltype, Llog,survivalPredictionMethod)
      gbmL1 = L1(gbmMod, Ltype, Llog,survivalPredictionMethod)
      bayesL1 = L1(bayesianNetMod, Ltype, Llog,survivalPredictionMethod)
      
      L1Results = rbind(coxL1,coxENL1,kmL1,rsfL1,aftL1,mtlrL1,gbmL1,bayesL1)
    }
    
    toAdd = as.data.frame(cbind(ConcordanceResults,
                                BrierResultsInt, BrierResultsSingle,L1Results))
    metricsRan = c(Concor,BrierInt,BrierSingle, L1Measure)
    names(toAdd) = c("Concordance",
                     "BrierInt","BrierSingle", "L1Loss")[metricsRan]
    modelsRan = c(CoxKP,CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel, GBMModel, BayesianNetModel)
    models = c("CoxKP","CoxKPEN","Kaplan-Meier","RSF","AFT", "MTLR", "GBM", "Bayes")[modelsRan]
    if(any(metricsRan)){
      toAdd = cbind.data.frame(Model = models,FoldNumer = i, toAdd)
    }else{
      toAdd = cbind.data.frame(Model = models,FoldNumer = i)
    }
    evaluationResults = rbind.data.frame(evaluationResults, toAdd)
    if(verbose){
      print(evaluationResults)
    }
  }
  if(DCal){
    if(verbose){
      print("Staring Evaluation: Cumulative D-Calibration")
    }
    coxDcal = DCalibrationCumulative(combinedTestResults$Cox,DCalBins)
    coxENDcal = DCalibrationCumulative(combinedTestResults$CoxEN,DCalBins)
    kmDcal = DCalibrationCumulative(combinedTestResults$KM,DCalBins)
    rsfDcal = DCalibrationCumulative(combinedTestResults$RSF,DCalBins)
    aftDcal = DCalibrationCumulative(combinedTestResults$AFT,DCalBins)
    mtlrDcal = DCalibrationCumulative(combinedTestResults$MTLR,DCalBins)
    gbmDcal = DCalibrationCumulative(combinedTestResults$GBM,DCalBins)
    bayesDcal = DCalibrationCumulative(combinedTestResults$BayesianNet,DCalBins)
    
    DCalResults = c(coxDcal,coxENDcal, kmDcal, rsfDcal, aftDcal, mtlrDcal, gbmDcal, bayesDcal)
    evaluationResults$DCalibration = rep(DCalResults, numberOfFolds)
  }
  if(OneCal){
    if(verbose){
      print("Staring Evaluation: Cumulative One-Calibration")
    }
    cox1cal = OneCalibrationCumulative(combinedTestResults$Cox, OneCalTime, typeOneCal, oneCalBuckets)
    coxEN1cal = OneCalibrationCumulative(combinedTestResults$CoxEN, OneCalTime, typeOneCal, oneCalBuckets)
    km1cal = OneCalibrationCumulative(combinedTestResults$KM, OneCalTime, typeOneCal, oneCalBuckets)
    rsf1cal = OneCalibrationCumulative(combinedTestResults$RSF, OneCalTime, typeOneCal, oneCalBuckets)
    aft1cal = OneCalibrationCumulative(combinedTestResults$AFT, OneCalTime, typeOneCal, oneCalBuckets)
    mtlr1cal = OneCalibrationCumulative(combinedTestResults$MTLR, OneCalTime, typeOneCal, oneCalBuckets)
    gbm1cal = OneCalibrationCumulative(combinedTestResults$GBM, OneCalTime, typeOneCal, oneCalBuckets)
    bayes1cal = OneCalibrationCumulative(combinedTestResults$BayesianNet, OneCalTime, typeOneCal, oneCalBuckets)
    
    numTimes = max(sapply(list(cox1cal,coxEN1cal, km1cal, rsf1cal, aft1cal, mtlr1cal, gbm1cal, bayes1cal),length))
    
    for(times in 1:numTimes){
      varName = paste("OneCalibration_",times, sep="")
      assign(varName,c(cox1cal[times],coxEN1cal[times], km1cal[times], rsf1cal[times],aft1cal[times], mtlr1cal[times], gbm1cal[times], bayes1cal[times]))
      evaluationResults[varName] = rep(eval(parse(text=varName)), numberOfFolds)
    }
    if(verbose){
      print(evaluationResults)
    }
  }
  #We will add some basic information about the dataset.
  evaluationResults$N = nrow(validatedData)
  #Note we subtract 2 to not count `time` and `delta`.
  evaluationResults$NumFeatures = ncol(training) - 2
  evaluationResults$PercentCensored = sum(!validatedData$delta)/nrow(validatedData)
  survivalCurves = getSurvivalCurves(coxTimes,coxENTimes, kmTimes, aftTimes, rsfTimes, mtlrTimes, gbmTimes, bayesianNetTimes,
                                     CoxKP,CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel, GBMModel, BayesianNetModel,
                                     combinedTestResults, numberOfFolds,originalIndexing)
  names(survivalCurves) = c("Cox","CoxEN","KM","AFT","RSF","MTLR","GBM","Bayes")[c(CoxKP,CoxKPEN, KaplanMeier, AFTModel,RSFModel, MTLRModel, GBMModel, BayesianNetModel)]
  rownames(evaluationResults) = NULL
  
  combinedBins = list(Cox=NULL, CoxEN=NULL, KM=NULL, RSF=NULL, AFT=NULL, MTLR=NULL, GBM=NULL, BayesianNet=NULL)
  combinedBins$Cox =colSums(ldply(lapply(seq_along(combinedTestResults$Cox), function(x) getBinned(combinedTestResults$Cox[[x]], DCalBins)), rbind))
  combinedBins$CoxEN =colSums(ldply(lapply(seq_along(combinedTestResults$CoxEN), function(x) getBinned(combinedTestResults$CoxEN[[x]], DCalBins)), rbind))
  combinedBins$RSF =colSums(ldply(lapply(seq_along(combinedTestResults$RSF), function(x) getBinned(combinedTestResults$RSF[[x]], DCalBins)), rbind))
  combinedBins$AFT =colSums(ldply(lapply(seq_along(combinedTestResults$AFT), function(x) getBinned(combinedTestResults$AFT[[x]], DCalBins)), rbind))
  combinedBins$MTLR =colSums(ldply(lapply(seq_along(combinedTestResults$MTLR), function(x) getBinned(combinedTestResults$MTLR[[x]], DCalBins)), rbind))
  combinedBins$GBM =colSums(ldply(lapply(seq_along(combinedTestResults$GBM), function(x) getBinned(combinedTestResults$GBM[[x]], DCalBins)), rbind))
  combinedBins$BayesianNet =colSums(ldply(lapply(seq_along(combinedTestResults$BayesianNet), function(x) getBinned(combinedTestResults$BayesianNet[[x]], DCalBins)), rbind))
  
  return(list(datasetUsed = validatedData, survivalCurves = survivalCurves, results = evaluationResults, DcalHistogram = combinedBins, ConCurve = ConcordanceCurve))
}



#This function combines survival curves across the folds into one dataframe (we must get predictions for all
#the times across all folds otherwise we cannot combine patients from different folds into a dataframe.)
getSurvivalCurves = function(coxTimes,coxENTimes, kmTimes, aftTimes, rsfTimes, mtlrTimes, gbmTimes, bayesianNetTimes,
                             CoxKP = T,CoxKPEN=T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel = T, GBMModel = T, BayesianNetModel = T,
                             combinedTestResults, numberOfFolds, originalIndexing){
  originalIndexOrder = order(unname(unlist(originalIndexing)))
  if(!is.null(coxTimes))
    coxTimes = sort(unique(coxTimes))
  if(!is.null(coxENTimes))
    coxENTimes = sort(unique(coxENTimes))
  if(!is.null(kmTimes))
    kmTimes = sort(unique(kmTimes))
  if(!is.null(rsfTimes))
    rsfTimes = sort(unique(rsfTimes))
  if(!is.null(aftTimes))
    aftTimes = sort(unique(aftTimes))
  if(!is.null(mtlrTimes))
    mtlrTimes = sort(unique(mtlrTimes))
  if(!is.null(gbmTimes))
    gbmTimes = sort(unique(gbmTimes))
  if(!is.null(bayesianNetTimes))
    bayesianNetTimes = sort(unique(bayesianNetTimes))
  models = c(CoxKP,CoxKPEN, KaplanMeier, AFTModel,RSFModel,MTLRModel, GBMModel, BayesianNetModel)
  allTimes = list(coxTimes,coxENTimes,kmTimes,aftTimes,rsfTimes,mtlrTimes,gbmTimes,bayesianNetTimes)
  survivalCurves = list()
  count = 0
  for(j in which(models)){
    count =count+1
    fullCurves = data.frame(row.names = 1:length(allTimes[[j]]))
    for(i in 1:numberOfFolds){
      #Index method -> fold -> survival curves
      times = combinedTestResults[[j]][[i]][[1]]$time
      maxTime = max(times)
      curves  = combinedTestResults[[j]][[i]][[1]][,-1]
      timesToEvaluate = setdiff(allTimes[[j]],times)
      #Here we are going to combine the times from all folds and fit a spline so all patients have predictions for all times
      #across all folds.
      fullCurves = cbind.data.frame(fullCurves,sapply(curves,
                                                      function(x){
                                                        curveSpline = splinefun(times,x,method='hyman')
                                                        maxSpline = curveSpline(maxTime)
                                                        curveSplineConstant = function(time){
                                                          timeToEval = ifelse(time > maxTime, maxTime,time)
                                                          toReturn = rep(NA,length(time))
                                                          toReturn[timeToEval== maxTime] = max(maxSpline,0)
                                                          toReturn[timeToEval !=maxTime] = curveSpline(timeToEval[timeToEval!=maxTime])
                                                          return(toReturn)
                                                        }
                                                        extraPoints =curveSplineConstant(timesToEvaluate)
                                                        toReturn = rep(NA, length(allTimes[[j]]))
                                                        originalIndex = which(!allTimes[[j]] %in% timesToEvaluate)
                                                        newIndex = which(allTimes[[j]] %in% timesToEvaluate)
                                                        toReturn[originalIndex] = x
                                                        toReturn[newIndex] = extraPoints
                                                        return(toReturn)
                                                      }
      ))
    }
    fullCurves =  fullCurves[originalIndexOrder]
    fullCurves = cbind.data.frame(allTimes[j], fullCurves)
    colnames(fullCurves) = c("time",1:(ncol(fullCurves)-1))
    survivalCurves[[count]] = fullCurves
  }
  return(survivalCurves)
}
