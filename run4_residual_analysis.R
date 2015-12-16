# 6.11.2015 I.Zliobaite
# residual analysis

input_file_predictions_t1_cv <- 'working_data/out_predictions_t1_cv.csv'
input_file_predictions_t1_fit <- 'working_data/out_predictions_t1_fit.csv'
input_file_predictions_t3_cv <- 'working_data/out_predictions_t3_cv.csv'
input_file_predictions_t3_fit <- 'working_data/out_predictions_t3_fit.csv'

plot_name_preds <- 'results/figure3.pdf'

######################################
predictions_t1_cv <-read.csv(input_file_predictions_t1_cv, header = TRUE, sep = '\t')
predictions_t1_fit <-read.csv(input_file_predictions_t1_fit, header = TRUE, sep = '\t')
predictions_t3_cv <-read.csv(input_file_predictions_t3_cv, header = TRUE, sep = '\t')
predictions_t3_fit <-read.csv(input_file_predictions_t3_fit, header = TRUE, sep = '\t')
wdt <- 5.5
hgt <- 8.3
txtcex = 0.7
txtoff <- 0.3
xylm <- c(2,22)
selected_model <- 8 #actually 7, but 1 for true, and 1 for 0
pdf(plot_name_preds,width = wdt, height = hgt)
#png(plot_name_preds)
# one plot (3.7 x 4.2)
par(mfrow = c(3,2), mar= c(4, 4, 1, 1) + 0.2) #matrix of subplots, mar is for margins
####
xylm <- c(200,1600)
selected_model <- 2 #actually 1, but 1 for true
plot(NA,xlab = 'modeled PREC_MIN',ylab = 'true PREC_MIN', main = 'PLS fit (model)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t1_fit[,selected_model],predictions_t1_fit[,'true'])
text(predictions_t1_fit[,selected_model],predictions_t1_fit[,'true'], predictions_t1_fit[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(predictions_t1_fit[1,'R2']),sep=' = ')
text(1300,310,tt, col = '#0072bd')
text(1300,210,paste('model:',as.character(selected_model-1),'comp.',sep=' '),col = '#0072bd')
plot(NA,xlab = 'predicted PREC_MIN',ylab = NA, main = 'PLS CV (testing)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t1_cv[,selected_model],predictions_t1_cv[,'true'])
text(predictions_t1_cv[,selected_model],predictions_t1_cv[,'true'], predictions_t1_cv[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(predictions_t1_cv[1,'R2']),sep=' = ')
text(1300,310,tt, col = '#0072bd')
text(1300,210,paste('model:',as.character(selected_model-1),'comp.',sep=' '),col = '#0072bd')
###
xylm <- c(0.1,0.7)
selected_model <- 12# actually 10 + 1 true + zero stage
plot(NA,xlab = 'modeled NDVI_MIN',ylab = 'true NDVI_MIN', main = 'LARS fit (model)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t3_fit[,selected_model],predictions_t3_fit[,'true'])
text(predictions_t3_fit[,selected_model],predictions_t3_fit[,'true'], predictions_t3_fit[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(round(predictions_t3_fit[1,'R2'],digits = 2)),sep=' = ')
text(0.6,0.15,tt, col = '#0072bd')
text(0.55,0.1,paste('model:',as.character(selected_model-2),'steps',sep=' '),col = '#0072bd')
plot(NA,xlab = 'predicted NDVI_MIN',ylab = NA, main = 'LARS CV (testing)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t3_cv[,selected_model],predictions_t3_cv[,'true'])
text(predictions_t3_cv[,selected_model],predictions_t3_cv[,'true'], predictions_t3_cv[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(round(predictions_t3_cv[1,'R2'],digits = 2)),sep=' = ')
text(0.6,0.15,tt, col = '#0072bd')
text(0.55,0.1,paste('model:',as.character(selected_model-2),'steps',sep=' '),col = '#0072bd')
###
xylm <- c(400,1800)
selected_model <- 5# actually 3 + 1 true
plot(NA,xlab = 'modeled NPP_MIN',ylab = 'true NPP_MIN', main = 'LARS fit (model)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t5_fit[,selected_model],predictions_t5_fit[,'true'])
text(predictions_t5_fit[,selected_model],predictions_t5_fit[,'true'], predictions_t5_fit[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(round(predictions_t5_fit[1,'R2'],digits = 2)),sep=' = ')
text(1550,500,tt, col = '#0072bd')
text(1550,400,paste('model:',as.character(selected_model-2),'steps',sep=' '),col = '#0072bd')
plot(NA,xlab = 'predicted NPP_MIN',ylab = NA, main = 'LARS CV (testing)',xlim=xylm,ylim=xylm)
abline(a=0,b=1)
points(predictions_t5_cv[,selected_model],predictions_t5_cv[,'true'])
text(predictions_t5_cv[,selected_model],predictions_t5_cv[,'true'], predictions_t5_cv[,'SITE'], cex=txtcex, pos = 3, offset = txtoff, col='#d95319')
tt = paste('R2',as.character(round(predictions_t5_cv[1,'R2'],digits = 2)),sep=' = ')
text(1550,500,tt, col = '#0072bd')
text(1550,400,paste('model:',as.character(selected_model-2),'steps',sep=' '),col = '#0072bd')
dev.off()