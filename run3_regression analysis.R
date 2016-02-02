# 31.1.2015 I.Zliobaite (Ilaret)
# model elev, precip and temp targets

library(lars)
library(pls)

#Net primary productivity (NPP) is computed from temperature and precipitation using the classic formula from Leith (1972), as cited in Liu et al (2012). 

#read
input_file_name <- 'working_data/data_all.csv'

output_file_R2_pls_fit <- 'working_data/out_R2_pls_fit.csv'
output_file_R2_pls_cv <- 'working_data/out_R2_pls_cv.csv'
output_file_R2_lars_fit <- 'working_data/out_R2_lars_fit.csv'
output_file_R2_lars_cv <- 'working_data/out_R2_lars_cv.csv'
output_file_predictions_t1_cv <- 'working_data/out_predictions_t1_cv.csv'
output_file_predictions_t1_fit <- 'working_data/out_predictions_t1_fit.csv'
output_file_predictions_t3_cv <- 'working_data/out_predictions_t3_cv.csv'
output_file_predictions_t3_fit <- 'working_data/out_predictions_t3_fit.csv'
output_file_predictions_t5_cv <- 'working_data/out_predictions_t5_cv.csv'
output_file_predictions_t5_fit <- 'working_data/out_predictions_t5_fit.csv'
output_file_models_lars <- 'working_data/out_models_lars.csv'
#output_file_features <- 'working_data/out_feature_selection.csv'
output_file_features_t1 <- 'working_data/out_feature_selection_t1.csv'
output_file_features_t3 <- 'working_data/out_feature_selection_t3.csv'
output_file_features_t5 <- 'working_data/out_feature_selection_t5.csv'
#output_file_fet_coefficients <- 'working_data/out_feature_coefficients.csv'
output_file_fet_coefficients_t1 <- 'working_data/out_feature_coefficients_t1.csv'
output_file_fet_coefficients_t3 <- 'working_data/out_feature_coefficients_t3.csv'
output_file_fet_coefficients_t5 <- 'working_data/out_feature_coefficients_t5.csv'
output_model1 <- 'working_data/out_model1_fit.txt'
output_model3 <- 'working_data/out_model3_fit.txt'
output_model5 <- 'working_data/out_model5_fit.txt'

input_file_features_t1 <- output_file_features_t1
input_file_features_t3 <- output_file_features_t3
input_file_features_t5 <- output_file_features_t5
input_file_fet_coefficients_t1 <- output_file_fet_coefficients_t1
input_file_fet_coefficients_t3 <- output_file_fet_coefficients_t3
input_file_fet_coefficients_t5 <- output_file_fet_coefficients_t5

output_file_table2 <- 'results/table6.txt'
output_file_table3 <- 'results/table7.txt'
output_file_table4 <- 'results/table8.txt'

plot_name_fig2 <- 'results/figure2.pdf'
plot_name_NPP <- 'results/figureNPP.pdf'

R2_files <- c(output_file_R2_lars_fit,output_file_R2_pls_fit,output_file_R2_lars_cv,output_file_R2_pls_cv)

fet_targets <- c('ELEV','TEMP','TEMP_MIN','PREC','PREC_MIN','NPP','NPP_MIN','NDVI','NDVI_MIN')

data_sites_all <- read.csv(input_file_name, header = TRUE)
#print('removing elevation')
#data_sites_all[,fet_targets] <- data_sites_all[,fet_targets] - data_sites_all[,'ELEV']
p <- dim(data_sites_all)[2]

#0s removed due to linear dependence
fet_inputs <- c('HYP','WCT_CS','WCT_AL','WCT_OL','WCT_CP','WCT_SF','WCT_CM',"HYP_2","HYP_3","WCT_CS_4","WCT_CS_5","WCT_CS_6",'WCT_CS_7',"WCT_AL_1","WCT_AL_2","WCT_OL_2","WCT_OL_3","WCT_OL_5","WCT_CP_1","WCT_CP_2","WCT_SF_1","WCT_SF_2",'MASS_log_mean')

# without proportions
#fet_inputs <- c('HYP','WCT_CS','WCT_AL','WCT_OL','WCT_CP','WCT_SF','WCT_CM','MASS_log_mean')

do_figure_a1 <- FALSE #select TRUE for doing Figure A1 in Appendix A.2

param_pls_method <- 'svdpc' #pca, works almost the same

param_steps <- 10 #buvo13
param_pls_comps <- 10
param_round_digits <- 3

#for model reporting
param_model_select_lars <- 10+1 #the first model is 0
param_model_select_pls <- 1

#for recording models for analysis of variable importance
param_target1 <- 'PREC_MIN'
param_select_t1 <- 1 #which model
dg1 <- 0
param_target3 <- 'NDVI_MIN'
param_select_t3 <- 10
dg3 = 3 #digits for rounding
param_target5 <- 'NPP_MIN'
param_select_t5 <- 3
dg5 = 0 #digits for rounding


if (do_figure_a1){
  fet_inputs <- c('HYP','WCT_CS','WCT_AL','WCT_OL','WCT_CP','WCT_SF','WCT_CM','MASS_log_mean')
  plot_name_fig2 <- 'results/figureA1.pdf'
  param_steps <- 8 
  param_pls_comps <- 8
  param_select_t3 <- 8 
  param_model_select_lars <- 1
}

##############################################################
#functions
compute_R2_matrix <- function(predictions_fun)
{  
  R2_all <- c()
  for (sk in 2:dim(predictions_fun)[2])
  {
    R2_all <- c(R2_all,compute_R2(predictions_fun[,1],predictions_fun[,sk]))  
  }
  return(R2_all)
}

compute_R2 <- function(true_fun,pred_fun)
{
  R2 <- 1 - sum((true_fun - pred_fun)^2)/sum((true_fun - mean(true_fun))^2)
  return(R2)
}

convert_predictions <- function(true_fun,pred_fun)
{
  comp_fun <- dim(pred_fun)[3]
  predictions_fun <- true_fun
  for (sk in 1:comp_fun)
  {
    predictions_fun <- cbind(predictions_fun,pred_fun[,,sk])
  }
  return(predictions_fun)
}

make_write_down_model <- function(model_now,mnames)
{
  k <- length(mnames)
  mnames <- process_feature_names(mnames)
  fml <- as.character(model_now[k+1])
  for(sk in 1:k)
  {
    if (model_now[sk]!=0)
    {
      if (model_now[sk]>0)
      {
        sgn = '+'
      }else
      {
        sgn = '-'
      }
      fml = paste(fml,sgn,as.character(abs(model_now[sk])),mnames[sk],sep = ' ') 
    }
  }
  return(fml)
}

convert_models <- function(coef_fun)
{
  comp_fun <- dim(coef_fun)[3]
  coefficients_fun <- c()
  for (sk in 1:comp_fun)
  {
    coefficients_fun <- rbind(coefficients_fun,coef_fun[,,sk])
  }
  return(coefficients_fun)
}

recover_coefficients <- function(fit_pls,comp,data_sites)
{
  model_now <- coef(fit_pls,ncomp = comp)
  n <- length(model_now)
  coefs <- model_now
  data <-data_sites*0
  data <- data[1,]
  for (sk in 1:n)
  {
    fet_now <- rownames(model_now)[sk]
    data1 <- data
    data1[fet_now] <- 1
    pred <- predict(fit_pls, newdata = data1)[,,comp] - predict(fit_pls, newdata = data)[,,comp]
    model_now[sk] <- pred
  }
  model_now <- c(model_now,predict(fit_pls, newdata = data)[,,comp])
  return(model_now)
}

process_feature_names <- function(Feature)
{
  Feature <- gsub('MASS_log_mean','log(MASS)',Feature)
  Feature <- gsub('WCT_', '', Feature)
  Feature <- gsub('new', '', Feature)
  Feature <- gsub('_', '=', Feature)
  for (sk in 1:length(Feature))
  {
    if (grepl('=',Feature[sk]))
    {
      Feature[sk] <- paste('prop(',Feature[sk],')',sep='')
    }
    else
    {
      Feature[sk] <- paste('mean(',Feature[sk],')',sep='')
    }
    if (grepl('med',Feature[sk]))
    {
      Feature[sk] <- gsub('med', '', Feature[sk])
      Feature[sk] <- gsub('mean', 'median', Feature[sk])
    }
  }
  return(Feature)
}

##############################################################

# manual cross validation
N <- dim(data_sites_all)[1]
all_R2_pls_cv <- c(1:param_pls_comps)
all_R2_lars_cv <- c(0:param_steps)
all_R2_pls_fit <- c(1:param_pls_comps)
all_R2_lars_fit <- c(0:param_steps)

print('target now:')
for (cik in 1:length(fet_targets))
{
  target_now <- fet_targets[cik]
  print(target_now)
  
  data_sites <- data_sites_all[,c(target_now,fet_inputs)] #for pls filter data
  fml <- as.formula(paste(target_now,'  ~.')) #regression variables, formula for pls
  
  predictions_lars_cv <- c()
  predictions_lars_fit <- c()
  predictions_pls_cv <- c()
  predictions_pls_fit <- c()
  
  model_collection_lars <- c()
  model_collection_pls <- c()
  
  fet_selection_lars <- data_sites*0
  fet_coefs_lars <- data_sites*0
  
  #cross validation
  for (ind_test in 1:N)
  {
    #training data
    ind_train <- setdiff(c(1:N),ind_test)
    
    #fit models
    fit_lars <- lars(as.matrix(data_sites[ind_train,fet_inputs]),as.matrix(data_sites[ind_train,target_now]), normalize = TRUE,max.steps = param_steps)
    fit_pls <- mvr(fml, param_pls_comps, method = param_pls_method, data = data_sites[ind_train,],scale = TRUE)
    
    #extract models (coefficients) plain regression
    model_lars <- coef(fit_lars, mode = 'step')
    model_pls <- convert_models(fit_pls$coefficients)
    intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
    intercept_pls <- convert_models(predict(fit_pls, newdata = data_sites[1,]*0))
    model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
    model_collection_pls <- rbind(model_collection_pls,c(model_pls[param_model_select_pls,],intercept_pls[param_model_select_pls]))
    
    #selection order
    for (sk6 in 1:param_steps)
    {
      if (length(fit_lars$actions[[sk6]])>1)
      {
        fta <- fit_lars$actions[[sk6]][1] 
        print(fit_lars$actions[[sk6]])
      }else
      {
        fta <- fit_lars$actions[[sk6]]
      }
      ind_fet_now <- which(colnames(fet_selection_lars)==names(fta))
      ind_coef_now <- which(as.data.frame(colnames(model_lars))==names(fta))
      if (fit_lars$actions[[sk6]][1]>0)
      {
        fet_selection_lars[ind_test,ind_fet_now] <- param_steps + 1 - sk6
        fet_coefs_lars[ind_test,ind_fet_now] <- sign(model_lars[sk6+1,ind_coef_now])
      }
    }
    
    #predictions together with ground truth in the first column
    pp_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs], type="fit")
    predictions_lars_cv <- rbind(predictions_lars_cv,cbind(data_sites[ind_test,target_now],t(pp_lars$fit)))
    pp_pls <- predict(fit_pls, newdata = data_sites[ind_test,])
    predictions_pls_cv <- rbind(predictions_pls_cv,convert_predictions(data_sites[ind_test,target_now],pp_pls))
    
    #manual prediction (model lars)
    #sum(model[13,]*data_sites[ind_test,fet_inputs]) + intercept$fit[13])
    #manual predictions pls
    #sum(fit_pls$coefficients[,,5]*data_sites[4,2:32]) + intercepts_pls[5]
  }
  
  R2_lars_cv <- compute_R2_matrix(predictions_lars_cv)
  R2_pls_cv <- compute_R2_matrix(predictions_pls_cv)
  
  all_R2_pls_cv <- cbind(all_R2_pls_cv,round(R2_pls_cv,digits = param_round_digits))
  all_R2_lars_cv <- cbind(all_R2_lars_cv,round(R2_lars_cv,digits = param_round_digits))
  
  #full fit and add to models as the last row
  fit_lars <- lars(as.matrix(data_sites[,fet_inputs]),as.matrix(data_sites[,target_now]), normalize = TRUE,max.steps = param_steps)
  fit_pls <- mvr(fml, param_pls_comps, method = param_pls_method, data = data_sites, scale = TRUE)
  model_lars <- coef(fit_lars, mode = 'step')
  model_pls <- convert_models(fit_pls$coefficients)
  intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
  intercept_pls <- convert_models(predict(fit_pls, newdata = data_sites[1,]*0))
  model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
  model_collection_pls <- rbind(model_collection_pls,c(model_pls[param_model_select_pls,],intercept_pls[param_model_select_pls]))
  
  
  #predictions all
  pp_lars <- predict.lars(fit_lars, data_sites[,fet_inputs], type="fit")
  predictions_lars_fit <- cbind(data_sites[,target_now],pp_lars$fit)
  pp_pls <- predict(fit_pls, newdata = data_sites)
  predictions_pls_fit <- convert_predictions(data_sites[,target_now],pp_pls)
  
  R2_lars_fit <- compute_R2_matrix(predictions_lars_fit)
  R2_pls_fit <- compute_R2_matrix(predictions_pls_fit)
  
  all_R2_lars_fit <- cbind(all_R2_lars_fit,round(R2_lars_fit,digits = param_round_digits))
  all_R2_pls_fit <- cbind(all_R2_pls_fit,round(R2_pls_fit,digits = param_round_digits))
  
  if (target_now == param_target1)
  {  
    model_now <- recover_coefficients(fit_pls,param_select_t1,data_sites)
    mnames <- rownames(coef(fit_pls,ncomp=param_select_t1))
    md <- make_write_down_model(round(model_now,digits = dg1),mnames)
    write.table(md,file = output_model1,quote = FALSE,row.names = FALSE,sep=',')
    pp_pls <- predict(fit_pls, newdata = data_sites)
    predictions_target1 <- convert_predictions(data_sites[,target_now],pp_pls)
    #sum(model_now[1:18]*data_sites_all[1,names(model_now)[1:18]]) + model_now[19]
    predictions_t1_cv <- cbind(as.data.frame(round(predictions_pls_cv,digits = dg1)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_pls_cv)[param_select_t1],digits = 2))
    R2vec <- rep(round(compute_R2_matrix(predictions_pls_fit)[param_select_t1],digits = 2),dim(data_sites_all)[1])
    predictions_t1_fit <- cbind(as.data.frame(round(predictions_pls_fit,digits = dg1)),data_sites_all[,c('SITE','LAT','LON')],R2vec,data_sites_all[,c('den_perm_riv','den_temp_riv','den_riv','ELEV')])
    colnames(predictions_t1_cv) <- c('true',rep('pr',(dim( predictions_t1_cv)[2]-5)),'SITE','LAT','LON','R2')
    colnames(predictions_t1_fit) <- c('true',rep('pr',(dim( predictions_t1_fit)[2]-9)),'SITE','LAT','LON','R2','DPER','DTEM','DRIV','ELEV')
    fet_selection_t1 <- fet_selection_lars
    fet_coefs_t1 <- fet_coefs_lars
  }
  
  if (target_now == param_target3){
    mod <- coef(fit_lars, mode = 'step')
    int <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
    model_now <- c(mod[param_select_t3+1,],int$fit[param_select_t3+1]) #+1 because 1st model 0
    mnames <- names(model_now)[1:(length(model_now)-1)]
    model_now <- as.vector(model_now)
    md <- make_write_down_model(round(model_now,digits = dg3),mnames)
    write.table(md,file = output_model3,quote = FALSE,row.names = FALSE,sep=',')
    # predictions_lars_cv <- rbind(predictions_lars_cv,cbind(data_sites[ind_test,target_now],t(pp_lars$fit)))
    predictions_t3_cv <- cbind(as.data.frame(round(predictions_lars_cv,digits = dg3)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_cv)[param_select_t3+1],digits = 3))
    predictions_t3_fit <- cbind(as.data.frame(round(predictions_lars_fit,digits = dg3)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_fit)[param_select_t3+1],digits = 3))
    colnames(predictions_t3_cv) <- c('true',rep('pr',(dim( predictions_t3_cv)[2]-5)),'SITE','LAT','LON','R2')
    colnames(predictions_t3_fit) <- c('true',rep('pr',(dim( predictions_t3_fit)[2]-5)),'SITE','LAT','LON','R2')
    fet_selection_t3 <- fet_selection_lars
    fet_coefs_t3 <- fet_coefs_lars
  }
  
  if (target_now == param_target5){
    mod <- coef(fit_lars, mode = 'step')
    int <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
    model_now <- c(mod[param_select_t5+1,],int$fit[param_select_t5+1]) #+1 because 1st model 0
    mnames <- names(model_now)[1:(length(model_now)-1)]
    model_now <- as.vector(model_now)
    md <- make_write_down_model(round(model_now,digits = dg5),mnames)
    write.table(md,file = output_model5,quote = FALSE,row.names = FALSE,sep=',')
    # predictions_lars_cv <- rbind(predictions_lars_cv,cbind(data_sites[ind_test,target_now],t(pp_lars$fit)))
    predictions_t5_cv <- cbind(as.data.frame(round(predictions_lars_cv,digits = dg5)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_cv)[param_select_t5+1],digits = 3))
    predictions_t5_fit <- cbind(as.data.frame(round(predictions_lars_fit,digits = dg5)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_fit)[param_select_t5+1],digits = 3))
    colnames(predictions_t5_cv) <- c('true',rep('pr',(dim( predictions_t5_cv)[2]-5)),'SITE','LAT','LON','R2')
    colnames(predictions_t5_fit) <- c('true',rep('pr',(dim( predictions_t5_fit)[2]-5)),'SITE','LAT','LON','R2')
    fet_selection_t5 <- fet_selection_lars
    fet_coefs_t5 <- fet_coefs_lars
  }
  
  #jacknife
  #jack = sqrt(apply((model_collection[1:13,] - matrix(1,13,1)%*%model_collection[14,])^2,2,sum)*(N-1)/N)
}

colnames(all_R2_pls_cv) <- c('step',fet_targets)
colnames(all_R2_pls_fit) <- c('step',fet_targets)
colnames(all_R2_lars_cv) <- c('step',fet_targets)
colnames(all_R2_lars_fit) <- c('step',fet_targets)

write.table(all_R2_pls_cv,file = output_file_R2_pls_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_pls_fit,file = output_file_R2_pls_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_lars_cv,file = output_file_R2_lars_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_lars_fit,file = output_file_R2_lars_fit,quote = FALSE,row.names = FALSE,sep='\t')

write.table(predictions_t1_fit,file = output_file_predictions_t1_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t1_cv,file = output_file_predictions_t1_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t3_fit,file = output_file_predictions_t3_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t3_cv,file = output_file_predictions_t3_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t5_fit,file = output_file_predictions_t5_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t5_cv,file = output_file_predictions_t5_cv,quote = FALSE,row.names = FALSE,sep='\t')


#write LARS models
#write.table(model_collection_lars,file = output_file_models_lars,quote = FALSE,row.names = FALSE,sep='\t')
#write feature selections
write.table(fet_selection_t1,file = output_file_features_t1,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_coefs_t1,file = output_file_fet_coefficients_t1,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_selection_t3,file = output_file_features_t3,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_coefs_t3,file = output_file_fet_coefficients_t3,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_selection_t5,file = output_file_features_t5,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_coefs_t5,file = output_file_fet_coefficients_t5,quote = FALSE,row.names = FALSE,sep='\t')


############################################################################################
#plotting 

xlabels <- c(NA,NA,'selection step','components')
ylabels <- c('R2',NA,'R2',NA)
titles <- c('LARS fit (modeling)','PLS fit (modeling)','LARS CV (testing)','PLS CV (testing)')
#parula
mycolors <-c('#0072bd','#d95319','#edb120','#7e2f8e','#77ac30','#4dbeee','#a2142f','#808080','#000000')
mylwd <- 2.2

# figure 2
pdf(plot_name_fig2,width = 6, height = 6)
#png(plot_name)
par(mfrow = c(2,2), mar= c(4, 4, 1, 1) + 0.2) #matrix of subplots, mar is for margins
for (plt in 1:4)
{
  #read file
  R2 <- read.csv(R2_files[plt], header = TRUE, sep = '\t')
  n <- dim(R2)[1]
  p <- dim(R2)[2]
  #setting up plotting area
  plot(NA,xlab = xlabels[plt],ylab = ylabels[plt],main = titles[plt],xlim=c(R2[1,'step'],R2[n,'step']),ylim=c(-1,1))
  #,ylim=c(min(R2[,c(2:p)]),max(R2[,c(2:p)])))
  #xaxs="i" reduces white space within axes
  for (sk in 1:length(fet_targets))
  {
    points(R2[,'step'],R2[,fet_targets[sk]],type="o",col = mycolors[sk],lwd=mylwd,pch=sk,cex=0.6)
  }  
  lines(c(R2[1,'step'],R2[n,'step']),c(0,0),col = 1,lwd=1, lty=2)
  if (plt == 1 )
  {
    legend("bottom",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.6)
  }
}
#legend("right",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.7)
dev.off()

#################################################
#feature selection tables

make_feature_table <- function(input_file_features,input_file_fet_coefficients,out_file,out_file_latex)
{
  feature_selection <-read.csv(input_file_features, header = TRUE, sep = '\t')
  feature_selection <- feature_selection[,c(2:dim(feature_selection)[2])]
  mean_rating <- apply(feature_selection,2,mean)
  ind <- order(-mean_rating)
  Feature <- names(mean_rating[ind])
  
  Feature <- process_feature_names(Feature)
  
  Importance <-  round(mean_rating[ind],digits = 1)
  Frequency <- round(apply(feature_selection[,ind]>0,2,sum)/dim(feature_selection)[1],digits=2)
  feature_table <- cbind(Feature,Importance,Frequency)
  
  feature_coefficients <-read.csv(input_file_fet_coefficients, header = TRUE, sep = '\t')
  #assume the same order as features
  feature_coefficients <- feature_coefficients[,c(2:dim(feature_coefficients)[2])]
  feature_coefficients <- feature_coefficients[,ind]
  Sign <- c()
  Consistency <- c()
  for (sk in 1:dim(feature_coefficients)[2])
  {
    ind_pl <- which(feature_coefficients[,sk]==1)
    ind_mn <- which(feature_coefficients[,sk]==-1)
    lpl <- length(ind_pl)
    lmn <- length(ind_mn)
    if (lpl>lmn)
    {
      Sign <- rbind(Sign,'+')
      Consistency <- rbind(Consistency,lpl/(lpl+lmn))
    }
    else
    {
      if (lpl<lmn)
      {
        Sign <- rbind(Sign,'-')
        Consistency <- rbind(Consistency,lmn/(lpl+lmn))
      }
      else
      {
        Sign <- rbind(Sign,'.')
        Consistency <- rbind(Consistency,NA)
      }
    }  
  }
  Consistency <- round(Consistency,digits = 2)
  colnames(Sign) <- 'Sign'
  colnames(Consistency) <- 'Consistency'
  
  feature_table <- cbind(feature_table,Sign,Consistency)
  
  write.table(feature_table,file = out_file,quote = FALSE,row.names = FALSE,sep='\t')
}

make_feature_table(input_file_features_t1,input_file_fet_coefficients_t1,output_file_table2)
make_feature_table(input_file_features_t3,input_file_fet_coefficients_t3,output_file_table3)
make_feature_table(input_file_features_t5,input_file_fet_coefficients_t5,output_file_table4)

#plot NPP

pdf(plot_name_NPP,width = 6,height = 6)
par(mfrow = c(2,2), mar= c(4, 4, 1, 1) + 0.2)
plot(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP_MIN'], main="(a) min. NPP",xlab = 'prop(SF=2)',ylab = 'NPP_MIN',xlim = c(-0.01,0.12),ylim=c(400,1800))
text(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP_MIN'],data_sites_all[,'SITE'],cex=0.7,pos=3,col='#d95319')
reg1 <- lm(NPP_MIN~WCT_SF_2,data=data_sites_all)
#abline(reg1,col='orange')
cc <- cor(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP_MIN'])
legend("bottomleft", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
##
plot(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP'], main="(b) NPP",xlab = 'prop(SF=2)',ylab = 'NPP',xlim = c(-0.01,0.12),ylim=c(600,2200))
text(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP'],data_sites_all[,'SITE'],cex=0.7,pos=3,col='#d95319')
reg1 <- lm(NPP~WCT_SF_2,data=data_sites_all)
cc <- cor(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NPP'])
legend("bottomleft", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
##
plot(data_sites_all[,'WCT_SF_2'],data_sites_all[,'PREC_MIN'], main="(c) min. PREC",xlab = 'prop(SF=2)',ylab = 'PREC_MIN',xlim = c(-0.01,0.12),ylim=c(300,1700))
text(data_sites_all[,'WCT_SF_2'],data_sites_all[,'PREC_MIN'],data_sites_all[,'SITE'],cex=0.7,pos=3,col='#d95319')
reg1 <- lm(PREC_MIN~WCT_SF_2,data=data_sites_all)
#abline(reg1,col='orange')
cc <- cor(data_sites_all[,'WCT_SF_2'],data_sites_all[,'PREC_MIN'])
legend("bottomleft", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
##
plot(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NDVI_MIN'], main="(d) min. NDVI",xlab = 'prop(SF=2)',ylab = 'NDVI_MIN',xlim = c(-0.01,0.12),ylim=c(0.1,0.7))
text(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NDVI_MIN'],data_sites_all[,'SITE'],cex=0.7,pos=3,col='#d95319')
reg1 <- lm(NDVI_MIN~WCT_SF_2,data=data_sites_all)
#abline(reg1,col='orange')
cc <- cor(data_sites_all[,'WCT_SF_2'],data_sites_all[,'NDVI_MIN'])
legend("bottomleft", bty="n", legend=paste("cor=",round(cc,digits=2),'\n',"R2=",round(summary(reg1)$adj.r.squared, digits=2)),text.col='orange')
dev.off()
