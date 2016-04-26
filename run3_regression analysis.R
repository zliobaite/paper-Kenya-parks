# 31.1.2015 I.Zliobaite (Ilaret)
# model elev, precip and temp targets

library(lars)

#Net primary productivity (NPP) is computed from temperature and precipitation using the classic formula from Leith (1972), as cited in Liu et al (2012). 

#read
input_file_name <- 'working_data/data_all.csv'

output_file_R2_lars_fit <- 'working_data/out_R2_lars_fit.csv'
output_file_R2_lars_cv <- 'working_data/out_R2_lars_cv.csv'
output_file_R2_one_cv <- 'working_data/out_R2_one_cv.csv'
output_file_models_lars <- 'working_data/out_models_lars.csv'

output_file_table <- 'results/table11.txt'
output_file_results_all <- 'results/table13.txt'

do_plot_fig2 <- FALSE

R2_files <- c(output_file_R2_lars_fit,output_file_R2_lars_cv,output_file_R2_one_cv)

fet_targets <- c('TEMP','TEMP_MIN','TEMPmin_MIN','TEMPmin','TEMP_MAX','TEMPmax_MAX','TEMPmax','PREC','PREC_MIN','PRECsp_MIN','PREC_MAX','PRECsp_MAX','NPP','NPP_MIN','NDVI','NDVI_MIN')

data_sites_all <- read.csv(input_file_name, header = TRUE)
p <- dim(data_sites_all)[2]

fet_inputs <- c('HYP','HOR','AL','OL','SF','OT','CM')
fet_inputs_all <- c('HYP','HOR','AL','OL','SF','OT','CM','MASS_log_mean')
fet_inputs_one <- c('MASS_log_mean','no_species_fact')

param_steps <- 6 #buvo13
param_round_digits <- 3

#for model reporting
param_model_select_lars <- 3+1 #the first model is 0

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
  #R2 <-sqrt(sum((true_fun - pred_fun)^2)/length(true_fun))/15
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

process_feature_names <- function(Feature)
{
  Feature <- gsub('MASS_log_mean','log(MASS)',Feature)
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
all_R2_lars_cv <- c(0:param_steps)
all_R2_lars_fit <- c(0:param_steps)
all_R2_ols_all_cv <- c()
all_R2_ols_mass_cv <- c()
all_R2_ols_nspec_cv <- c()

print('target now:')
for (cik in 1:length(fet_targets))
{
  target_now <- fet_targets[cik]
  print(target_now)
  
  data_sites <- data_sites_all[,c(target_now,union(union(fet_inputs,fet_inputs_one),fet_inputs_all))] #for pls filter data
  fml <- as.formula(paste(target_now,'  ~.')) #regression variables, formula for pls
  
  predictions_lars_cv <- c()
  predictions_lars_fit <- c()
  predictions_ols_all_cv <- c()
  predictions_ols_mass_cv <- c()
  predictions_ols_nspec_cv <- c()
  
  model_collection_lars <- c()
  
  fet_selection_lars <- data_sites*0
  fet_coefs_lars <- data_sites*0
  
  #cross validation
  for (ind_test in 1:N)
  {
    #training data
    ind_train <- setdiff(c(1:N),ind_test)
    
    #fit models
    fit_lars <- lars(as.matrix(data_sites[ind_train,fet_inputs]),as.matrix(data_sites[ind_train,target_now]), normalize = TRUE,max.steps = param_steps)
    
    fit_ols_all <- lm(as.formula(paste(target_now,'~.')), data=data_sites[ind_train,c(fet_inputs_all,target_now)])
    fit_ols_mass <- lm(as.formula(paste(target_now,'~ MASS_log_mean')), data=data_sites[ind_train,c(fet_inputs_one,target_now)])
    fit_ols_nspec <- lm(as.formula(paste(target_now,'~ no_species_fact')), data=data_sites[ind_train,c(fet_inputs_one,target_now)])
    
    #extract models (coefficients) plain regression
    model_lars <- coef(fit_lars, mode = 'step')

    intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
  
    model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
    
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
    predictions_ols_all_cv <- rbind(predictions_ols_all_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_all,data_sites[ind_test,c(fet_inputs_all,target_now)])))
    predictions_ols_mass_cv <- rbind(predictions_ols_mass_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_mass,data_sites[ind_test,c(fet_inputs_one,target_now)])))
    predictions_ols_nspec_cv <- rbind(predictions_ols_nspec_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_nspec,data_sites[ind_test,c(fet_inputs_one,target_now)])))
    
    #manual prediction (model lars)
    #sum(model[13,]*data_sites[ind_test,fet_inputs]) + intercept$fit[13])
    #manual predictions pls
    #sum(fit_pls$coefficients[,,5]*data_sites[4,2:32]) + intercepts_pls[5]
  }
  
  R2_lars_cv <- compute_R2_matrix(predictions_lars_cv)
  R2_ols_all_cv <- compute_R2_matrix(predictions_ols_all_cv)
  R2_ols_mass_cv <- compute_R2_matrix(predictions_ols_mass_cv)
  R2_ols_nspec_cv <- compute_R2_matrix(predictions_ols_nspec_cv)
  
  all_R2_lars_cv <- cbind(all_R2_lars_cv,round(R2_lars_cv,digits = param_round_digits))
  all_R2_ols_all_cv <- cbind(all_R2_ols_all_cv,round(R2_ols_all_cv,digits = param_round_digits))
  all_R2_ols_mass_cv <- cbind(all_R2_ols_mass_cv,round(R2_ols_mass_cv,digits = param_round_digits))
  all_R2_ols_nspec_cv <- cbind(all_R2_ols_nspec_cv,round(R2_ols_nspec_cv,digits = param_round_digits))
  
  #full fit and add to models as the last row
  fit_lars <- lars(as.matrix(data_sites[,fet_inputs]),as.matrix(data_sites[,target_now]), normalize = TRUE,max.steps = param_steps)
 
  model_lars <- coef(fit_lars, mode = 'step')
 
  intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
  
  model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
  
  #predictions all
  pp_lars <- predict.lars(fit_lars, data_sites[,fet_inputs], type="fit")
  predictions_lars_fit <- cbind(data_sites[,target_now],pp_lars$fit)
  
  R2_lars_fit <- compute_R2_matrix(predictions_lars_fit)
  
  all_R2_lars_fit <- cbind(all_R2_lars_fit,round(R2_lars_fit,digits = param_round_digits))
  
  #jacknife
  #jack = sqrt(apply((model_collection[1:13,] - matrix(1,13,1)%*%model_collection[14,])^2,2,sum)*(N-1)/N)
}

colnames(all_R2_lars_cv) <- c('step',fet_targets)
colnames(all_R2_lars_fit) <- c('step',fet_targets)
colnames(all_R2_ols_all_cv) <- fet_targets
colnames(all_R2_ols_mass_cv) <- fet_targets
colnames(all_R2_ols_nspec_cv) <- fet_targets
step <- c('all','mass','nspec')
all_R2_one_cv <- cbind(step,rbind(all_R2_ols_all_cv,all_R2_ols_mass_cv,all_R2_ols_nspec_cv))


modelno <- c()
for (sk in 2:dim(all_R2_lars_cv)[2]){
  ind <- which(all_R2_lars_cv[,sk] == max(all_R2_lars_cv[,sk]))
  modelno <- c(modelno,all_R2_lars_cv[ind,'step'])
}
ind <- which(all_R2_lars_cv[,'step']==3)
results_all <- cbind(t(all_R2_ols_mass_cv),t(all_R2_ols_nspec_cv),all_R2_lars_cv[ind,2:dim(all_R2_lars_cv)[2]],apply(all_R2_lars_cv[,2:dim(all_R2_lars_cv)[2]],2,max),t(all_R2_ols_all_cv))
colnames(results_all) <- c('OLS mass','OLS no. spec.','LARS(3)','LARS best','OLS all')

results_all <- (results_all - (-0.174))/(1 - (-0.174))
results_all <- round(results_all,digits = 2)
rownames(all_R2_ols_all_cv) <- names(all_R2_ols_all_cv)
results_all <- cbind(results_all,modelno)


write.table(results_all,file = output_file_results_all,quote = FALSE,row.names = TRUE,sep='\t')

write.table(all_R2_lars_cv,file = output_file_R2_lars_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_lars_fit,file = output_file_R2_lars_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_one_cv,file = output_file_R2_one_cv,quote = FALSE,row.names = FALSE,sep='\t')

if (do_plot_fig2){
  ############################################################################################
  #plotting 
  
  xlabels <- c('selection step','selection step')
  ylabels <- c('R2',NA,NA)
  naive <- c(0,-0.174,-0.174)
  titles <- c('LARS fit (modeling)','LARS CV (testing)', 'OLS CV')
  #parula
  mycolors <-c('#0072bd','#d95319','#edb120','#7e2f8e','#77ac30','#4dbeee','#a2142f','#808080','#000000')
  mylwd <- 1.5
  
  # figure 2
  pdf(plot_name_fig2,width = 6, height = 2.7)
  #png(plot_name)
  #par(mfrow = c(1,3), mar= c(4, 4, 1, 1) + 0.2) #matrix of subplots, mar is for margins
  layout(t(c(1,2,3)),widths = c(2.5,2.5,1.5), heights = c(2.5,2.5,2.5))
  for (plt in 1:3)
  {
    #read file
    R2 <- read.csv(R2_files[plt], header = TRUE, sep = '\t')
    n <- dim(R2)[1]
    p <- dim(R2)[2]
    #setting up plotting area
    if (plt==3){
      plot(NA,xlab = xlabels[plt],ylab = ylabels[plt],main = titles[plt],xlim=c(0,4),ylim=c(-1,1))
    }else{
      plot(NA,xlab = xlabels[plt],ylab = ylabels[plt],main = titles[plt],xlim=c(R2[1,'step'],R2[n,'step']),ylim=c(-1,1))
    }
    #,ylim=c(min(R2[,c(2:p)]),max(R2[,c(2:p)])))
    #xaxs="i" reduces white space within axes
    for (sk in 1:length(fet_targets))
    {
      if (plt==3){
        points(R2[,'step'],R2[,fet_targets[sk]],type="b",col = mycolors[sk],lwd=mylwd,pch=sk,cex=0.6)
      }else{
        points(R2[,'step'],R2[,fet_targets[sk]],type="o",col = mycolors[sk],lwd=mylwd,pch=sk,cex=0.6)
      }
    }  
    lines(c(R2[1,'step'],R2[n,'step']),c(naive[plt],naive[plt]),col = 1,lwd=1, lty=2)
    if (plt == 1 )
    {
      legend("bottom",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.6)
    }
  }
  #legend("right",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.7)
  dev.off()
}

