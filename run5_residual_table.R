# 2015 11 09 I.Zliobaite
# print residual table

input_file_predictions_t1_cv <- 'working_data/out_predictions_t1_cv.csv'
input_file_predictions_t1_fit <- 'working_data/out_predictions_t1_fit.csv'
input_file_predictions_t3_cv <- 'working_data/out_predictions_t3_cv.csv'
input_file_predictions_t3_fit <- 'working_data/out_predictions_t3_fit.csv'
input_file_predictions_t5_cv <- 'working_data/out_predictions_t5_cv.csv'
input_file_predictions_t5_fit <- 'working_data/out_predictions_t5_fit.csv'

out_file_t1 <- 'results/table3.txt'
out_file_t3 <- 'results/table4.txt'
out_file_t5 <- 'results/table5.txt'

######################################

predictions_t1_cv <-read.csv(input_file_predictions_t1_cv, header = TRUE, sep = '\t')
predictions_t1_fit <-read.csv(input_file_predictions_t1_fit, header = TRUE, sep = '\t')
predictions_t3_cv <-read.csv(input_file_predictions_t3_cv, header = TRUE, sep = '\t')
predictions_t3_fit <-read.csv(input_file_predictions_t3_fit, header = TRUE, sep = '\t')
predictions_t5_cv <-read.csv(input_file_predictions_t5_cv, header = TRUE, sep = '\t')
predictions_t5_fit <-read.csv(input_file_predictions_t5_fit, header = TRUE, sep = '\t')

######################################

plot_name_res <- 'results/figure5.pdf'

no_model_t1 <- 2 # 1 - true, 2 ... model predictions
no_model_t3 <- 12 # 1 - true, 2 - naive, 3 ... model predictions
no_model_t5 <- 5 # 1 - true, 2 - naive, 3 ... model predictions

print_res_table <- function(predictions,no_model){
  true <- predictions[,1]  
  pred <- predictions[,no_model]
  nv <- mean(true)
  pred_nv <- pred*0 + nv
  res_main <- true - pred
  res_nv <- true - pred_nv
  res_per <- res_main/true
  res_rel <- res_main/res_nv
  #pred_nv <- round(pred_nv)
  #res_nv <- round(res_nv)
  #res_per <- round(res_per,digits = 3)
  #res_rel <- round(res_rel,digits = 3)
  res <- predictions[,'SITE']
  res <- as.data.frame(res)
  res <- cbind(res,true,pred,pred_nv,res_main,res_nv,res_per,res_rel)
  colnames(res) <- c('Site','observed','model','naive','res. model','res. naive','res. perc.','res. rel.')
  return(res)
}

#####
fit_per <- lm(predictions_t1_fit[,'DPER'] ~ predictions_t1_fit[,'true'], data=predictions_t1_fit)
out_per <- predict(fit_per)
res_per <- predictions_t1_fit[,'DPER'] - out_per

fit_ses <- lm(predictions_t1_fit[,'DTEM'] ~ predictions_t1_fit[,'true'], data=predictions_t1_fit)
out_ses <- predict(fit_ses)
res_ses <- predictions_t1_fit[,'DTEM'] - out_ses

fit_riv <- lm(predictions_t1_fit[,'DRIV'] ~ predictions_t1_fit[,'true'], data=predictions_t1_fit)
out_riv <- predict(fit_riv)
res_riv <- predictions_t1_fit[,'DRIV'] - out_riv

compute_R2 <- function(true_fun,pred_fun)
{
  R2 <- 1 - sum((true_fun - pred_fun)^2)/sum((true_fun - mean(true_fun))^2)
  return(R2)
}

res_table_PREC <- print_res_table(predictions_t1_cv,no_model_t1)
res_table_PREC[,c(2:6)] <- round(res_table_PREC[,c(2:6)])
res_table_PREC[,c(7:8)] <- round(res_table_PREC[,c(7:8)],digits = 2)
write.table(res_table_PREC,file = out_file_t1,quote = FALSE,row.names = FALSE,sep='\t')

res_table_NDVI <- print_res_table(predictions_t3_cv,no_model_t3)
res_table_NDVI[,c(2:8)] <- round(res_table_NDVI[,c(2:8)],digits = 2)
write.table(res_table_NDVI,file = out_file_t3,quote = FALSE,row.names = FALSE,sep='\t')

res_table_NPP <- print_res_table(predictions_t5_cv,no_model_t5)
res_table_NPP[,c(2:6)] <- round(res_table_NPP[,c(2:6)])
res_table_NPP[,c(2:8)] <- round(res_table_NPP[,c(2:8)],digits = 2)
write.table(res_table_NPP,file = out_file_t5,quote = FALSE,row.names = FALSE,sep='\t')


sites <- as.vector(predictions_t1_fit[,'SITE'])

res <- res_table_PREC[,5] #residuals
data_res <- as.data.frame(cbind(res,res_per,res_ses))

pdf(plot_name_res,width = 7,height = 4)
par(mfrow=c(1,2))
plot(data_res[,'res'],data_res[,'res_per'], main="Permanent rivers",xlab = 'precipitation residual, mm',ylab = 'river index residual',xlim = c(-400,400),ylim=c(-85,85))
text(data_res[,'res'],data_res[,'res_per'],sites,cex=0.7,pos=3,col='#d95319')
reg1 <- lm(res_per~res,data=data_res)
abline(reg1,col='orange')
cc <- cor(data_res[,'res'],data_res[,'res_per'])
legend("bottomright", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
plot(data_res[,'res'],data_res[,'res_ses'],main="Seasonal rivers",xlab = 'precipitation residual, mm',ylab = 'river index residual',xlim = c(-400,400),ylim = c(-85,85))
text(data_res[,'res'],data_res[,'res_ses'],sites,cex=0.7, pos = 3, col='#d95319')
reg1 <- lm(res_ses~res,data=data_res)
abline(reg1,col='orange')
cc <- cor(data_res[,'res'],data_res[,'res_ses'])
legend("bottomright", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
dev.off()

#######
if (sum(as.vector(res_table_PREC[,'Site'])!= sites)>0){
  print('!!!site mismatch')
}