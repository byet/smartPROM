source("inference.R")
source("information.R")

infer_test = function(){
  bn <- loadHuginNet("models/bety_test.net")
  bety_factor <- read.csv("data/bety_factor.csv", sep = ";")
  bety_factor <- na.omit(bety_factor)
  bety_factor <- apply(bety_factor, 2, as.factor)
  bety_factor <- bety_factor[, -1]
  bety_factor <- as.data.frame(bety_factor)
  response <- c("Factor_1","Factor_2","Factor_3","Factor_4")
  input <- colnames(bety_factor)
  input <- input[which(!input %in% response)]
  n = nrow(bety_factor)
  
  res <- matrix(NA,ncol=12,nrow=n)
  
  for(i in 1:2){
    print(i)
    test_data <- bety_factor[i,]
    test_data <- apply(test_data,2,as.character)
    this_res <- informative_calc(data = test_data,target = response,query = input,bn_model = bn,num_iter = 3)
    res[i,] <- unlist(this_res)
  }
  print(res)
}


inf_test <- function(){
  
  print(entropy(c(0.1451612903225806,
                  0.1754032258064515,
                  0.6794354838709679)))
  
  file_path <- "data/bety_factor.csv"
  bety_factor <- read.csv(file_path, sep = ";")
  bety_factor <- na.omit(bety_factor)
  bety_factor <- apply(bety_factor, 2, as.factor)
  bety_factor <- bety_factor[, -1]
  bety_factor <- as.data.frame(bety_factor)
  response <- c("Factor_1","Factor_2","Factor_3","Factor_4")
  input <- colnames(bety_factor)
  input <- input[which(!input %in% response)]
  
  #bn <- tabu(bety_factor)
  #bn <- bn.fit(bn, bety_factor,method = "bayes")
  #write.net(bn,"bety_c1.net")
  #bn <- as.grain(bn)
  #saveHuginNet(bn,"bety_c.net")
  bn <- loadHuginNet("models/bety_test.net")
  cond_ent(target = response,evidence_e = "B22_0",bn_model = bn)
  mut_inf(target = response, evidence = input, bn_model = bn)
}
inf_test()
infer_test()
