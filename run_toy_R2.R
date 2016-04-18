# 2016 04 18 I.Zliobaite 
# for plotting R2 dependence of n

rep <- 1

compute_R2 <- function(true_fun,pred_fun)
{
  R2 <- 1 - sum((true_fun - pred_fun)^2)/sum((true_fun - mean(true_fun))^2)
  return(R2)
}

R2 <- c()
for (n in 2:50){
  data_now <- runif(n, 0, 19)
  predict_now <- c()
  for (sk in 1:n){
    test_now <- data_now[sk]
    predict_now <- c(predict_now,mean(data_now[-sk]))
  }
  R2 <- rbind(R2,c(n,compute_R2(data_now,predict_now),(1 - (n/(n-1))^2)))
}

