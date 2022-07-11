# Hospitalization Rates from National Covariates
myWB = read.csv("data/HospitalizationRates_World_Bank_Indicators.csv")
myWB = myWB[1:12]

# Calculate share of urban population
myWB$UrbanPop = myWB$Urban.population  / myWB$PopCNTY #problem with values higher than one
myWB = myWB[-which(names(myWB) %in% c("Urban.population","PopCNTY"))]
myWB$Physicians..per.1000.people. = as.numeric(myWB$Physicians..per.1000.people.)

myCovs = myWB[3:ncol(myWB)]
ExclusionDF = vifstep(myCovs, th = 5)
myCovsI = myCovs[-(which(names(myCovs) %in% ExclusionDF@excluded))]

# do imputation of missing covariates 
myCovsImp = complete(mice(myCovsI, action = "all"))

myWB = cbind(myWB[1:2],myCovsImp)

# use logistic regression to calculate missing hospitalization rates from covariates 
myForm = as.formula(paste0("HRNTYR~",paste0(names(myCovsI), collapse = "+")))
myGLM = glm(myForm, data = myWB, family = "quasibinomial")
myWB$HRNTYR_p = predict(myGLM, newdata =  myCovsImp, type = "response")

