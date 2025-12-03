library(MTS)
library(vars)
library(tseries)
library(xts)
library(zoo)
library(dplyr)
library(VARtests)
library(readxl)
library(rmgarch)
library(rugarch)
library(psych)

data <- read_excel("C:/Users/jessi/OneDrive/Documents/Skripsi/Data/Data Emas, IHSG, Nilai Tukar (2014 - Juni 2025).xlsx")

# Mengubah kolom 'Tanggal' menjadi format tanggal
dataTanggal <- as.Date(data$Tanggal)

# Mengatur data sebagai objek time series dengan tanggal (xts)
data_xts <- xts(data[, -1], order.by = data$Tanggal)

# Menggunakan log return: rt = log(Pt) - log(Pt-1)
log_return_xts <- diff(log(data_xts)) %>%
  na.omit()

# Nama kolom untuk log return
colnames(log_return_xts) <- c("rEmas", "rIHSG", "rUSD")

train_data_xts <- log_return_xts[1:2772]

# ADF Test
print(adf.test(train_data_xts$rEmas))
print(adf.test(train_data_xts$rIHSG))
print(adf.test(train_data_xts$rUSD))

# Konversi ke matriks/data.frame untuk model VAR/VARMA
train_data_mat <- as.matrix(train_data_xts)

# Menggunakan informasi kriteria untuk memilih orde VAR (p)
VARselect_results <- VARselect(train_data_mat, lag.max = 10, type = "const")
print(VARselect_results$selection)
p_optimal_var <- VARselect_results$selection[1]

# Membangun model VAR
p <- p_optimal_var
var_model <- VAR(train_data_mat, p = p, type = "const")
print(summary(var_model))
Acoef(var_model)

# Intercept µ (vector 3x1)
mu <- matrix(c(
  0.0003725,   # μ1 untuk rEmas
  0.0001340,   # μ2 untuk rIHSG
  1.520e-04    # μ3 untuk rUSD
), ncol = 1)


# Koefisien Φ1 ... Φ9 (list of matrices)
Phi_list <- Acoef(var_model)

# FITTED VALUE
returns <- train_data_mat
N <- nrow(returns)
fitted_manual <- matrix(NA, nrow=N, ncol=3)
colnames(fitted_manual) <- colnames(returns)

for (t in (p+1):N) {
  pred <- mu
  
  for (lag in 1:p) {
    pred <- pred + Phi_list[[lag]] %*% as.numeric(returns[t-lag,])
  }
  
  fitted_manual[t,] <- pred
}

idx <- as.Date(index(train_data_xts)) 
fitted_return_xts <- xts(fitted_manual, order.by = idx)

##### BUAT RMSE DLL
pred_returns <- fitted_return_xts
actual_returns <- log_return_xts

common_index <- complete.cases(pred_returns)

pred_returns   <- pred_returns[common_index, ]
actual_returns <- actual_returns[common_index, ]

print(dim(actual_returns))
print(dim(pred_returns))

############## RMSE
err_norm <- sqrt(
  (actual_returns[,1] - pred_returns[,1])^2 +
    (actual_returns[,2] - pred_returns[,2])^2 +
    (actual_returns[,3] - pred_returns[,3])^2
)


rmse <- function(actual_returns, pred_returns){
  sqrt(mean((actual_returns - pred_returns)^2, na.rm = TRUE))
}

# RMSE per variabel
rmse_emas <- rmse(actual_returns[,1], pred_returns[,1])
rmse_ihsg <- rmse(actual_returns[,2], pred_returns[,2])
rmse_usd  <- rmse(actual_returns[,3], pred_returns[,3])

# RMSE multivariat (norm-based)
err_norm <- sqrt(
  (actual_returns[,1] - pred_returns[,1])^2 +
    (actual_returns[,2] - pred_returns[,2])^2 +
    (actual_returns[,3] - pred_returns[,3])^2
)

rmse_norm <- sqrt(mean(err_norm^2, na.rm = TRUE))

rmse_table <- data.frame(
  Aset = c("Emas", "IHSG", "USD/IDR", "Multivariat (Norm)"),
  RMSE = c(rmse_emas, rmse_ihsg, rmse_usd, rmse_norm)
)

print(rmse_table)

###### SMAPE
smape <- function(actual_returns, pred_returns){
  mean(abs(pred_returns - actual_returns) / ((abs(actual_returns) + abs(pred_returns)) / 2),
       na.rm = TRUE) * 100
}

smape_emas <- smape(actual_returns[,1], pred_returns[,1])
smape_ihsg <- smape(actual_returns[,2], pred_returns[,2])
smape_usd  <- smape(actual_returns[,3], pred_returns[,3])

act_norm<-sqrt(train_data_xts$rIHSG^2+train_data_xts$rUSD^2+train_data_xts$rEmas^2)

smape_norm <- mean(
  err_norm / act_norm,
  na.rm = TRUE
) * 100

smape_table <- data.frame(
  Aset = c("Emas", "IHSG", "USD/IDR", "Multivariat (Norm)"),
  SMAPE = c(smape_emas, smape_ihsg, smape_usd, smape_norm)
)

print(smape_table)

#MAE
mae <- function(actual_returns, pred_returns){
  mean(abs(pred_returns - actual_returns),
       na.rm = TRUE)
}

mae_emas <- mae(actual_returns[,1], pred_returns[,1])
mae_ihsg <- mae(actual_returns[,2], pred_returns[,2])
mae_usd  <- mae(actual_returns[,3], pred_returns[,3])
mae_norm <- mean(abs(err_norm), na.rm = TRUE)

mae_table <- data.frame(
  Aset = c("Emas", "IHSG", "USD/IDR", "Multivariat (Norm)"),
  MAE = c(mae_emas, mae_ihsg, mae_usd, mae_norm)
)

print(mae_table)

annualized_mae_emas <- mae(actual_returns[,1], pred_returns[,1])*252
annualized_mae_ihsg <- mae(actual_returns[,2], pred_returns[,2])*252
annualized_mae_usd  <- mae(actual_returns[,3], pred_returns[,3])*252
annualized_mae_norm <- mae_norm * 252

annualized_mae_table <- data.frame(
  Aset = c("Emas", "IHSG", "USD/IDR", "Multivariat (Norm)"),
  MAE = c(annualized_mae_emas, annualized_mae_ihsg, annualized_mae_usd, annualized_mae_norm)
)

print(annualized_mae_table)

### Uji Asumsi VAR
# Uji Normalitas (Multivariate Skewness and Kurtosis)
cat("\n--- Uji Normalitas Sisaan VAR ---\n")
print(normality.test(var_model))

# Uji Autokorelasi (Portmanteau Test)
cat("\n--- Uji Autokorelasi Sisaan VAR ---\n")
# Nilai lags.pt > p_var (misal lags.pt=p_var+1, atau 10)
p_var=9
print(serial.test(var_model, lags.pt = p_var + 5, type = "PT.asymptotic"))

# Uji Heteroskedastisitas (ARCH-LM Test)
cat("\n--- Uji Heteroskedastisitas Sisaan VAR (ARCH-LM) ---\n")
print(arch.test(var_model, lags.multi = 10, multivariate.only = TRUE))

# Spesifikasi GARCH(1,1)
uspec <- ugarchspec(
  mean.model = list(armaOrder = c(0,0)), 
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  distribution.model = "norm"
)

# Spesifikasi DCC + VAR(9)
spec_dcc_var9 <- dccspec(
  uspec = multispec(replicate(3, uspec)),  
  dccOrder = c(1,1),
  VAR = TRUE,    
  lag = 9,      
  distribution = "mvt" 
)

########################################

Pmax <- 10
Qmax <- 2
Eccm(train_data_mat, maxp=10, maxq=5)
best <- data.frame(p = integer(), q = integer(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)

for(p in 0:Pmax){
  for(q in 0:Qmax){
    try({
      fit <- VARMA(train_data_mat, p = p, q = q, prelim = TRUE)   # prelim=TRUE untuk estimasi awal
      # fit$aic dan fit$bic tersedia menurut dokumentasi MTS
      best <- rbind(best, data.frame(p = p, q = q, AIC = fit$aic, BIC = fit$bic))
    }, silent = TRUE)
  }
}

best <- best[order(best$AIC), ]
print(head(best, 10))   # 10 kombinasi teratas menurut AIC
# pilih p,q dari baris pertama
p_opt <- 9
q_opt <- 0
cat("Pilihan (p,q) terbaik menurut AIC:", p_opt, q_opt, "\n")

varma_model <- VARMA(train_data_mat, p = p_opt, q = q_opt, include.mean = TRUE, prelim = FALSE, details = TRUE)

# Residual
resid_varma <- residuals(varma_model)

cat("\n--- Uji Diagnostik VARMA ---\n")

# 1. Uji White Noise / Autokorelasi Multivariat
cat("\n--- Portmanteau Test (mq) ---\n")
mq(resid_varma, lag = 20)

# 2. Uji Normalitas Multivariat
cat("\n--- Uji Normalitas (Mardia JB) ---\n")
res_varma <- residuals(varma_model)
mardia(res_varma)

# 3. Uji Heteroskedastisitas (ARCH effect)
cat("\n--- Uji Heteroskedastisitas (mq pada residual^2) ---\n")
mq(resid_varma^2, lag = 10)


# Granger Causality Test
causality(var_model, cause = "rEmas")
causality(var_model, cause = "rIHSG")
causality(var_model, cause = "rUSD")


# Fit model langsung
dcc_var9_model <- dccfit(spec_dcc_var9, data = train_data_xts)

irf_var <- irf(
  var_model, 
  impulse = c("rEmas", "rIHSG", "rUSD"),  # variabel yang diberi shock
  response = c("rEmas", "rIHSG", "rUSD"), # variabel yang merespons
  n.ahead = 30,   # jumlah periode ke depan yang diamati
  boot = TRUE,    # bootstrap untuk interval kepercayaan
  ci = 0.95,       # level confidence interval
  runs = 1000
)

# Menampilkan hasil numeriknya
print(irf_var)

# Visualisasi IRF
plot(irf_var)

######### PLOT
par(mfrow=c(3,1), mar=c(4,4,2,1))

# Plot rEmas
plot(actual_returns$rEmas, type="l", col="black", lwd=1,
     main="Actual vs Fitted - rEmas", ylab="Return")
lines(pred_returns$rEmas, col="red", lwd=1.5)
legend("topleft", legend=c("Actual","Fitted"), col=c("black","red"), lwd=2)

# Plot rIHSG
plot(actual_returns$rIHSG, type="l", col="black", lwd=1,
     main="Actual vs Fitted - rIHSG", ylab="Return")
lines(pred_returns$rIHSG, col="blue", lwd=1.5)
legend("topleft", legend=c("Actual","Fitted"), col=c("black","blue"), lwd=2)

# Plot rUSD
plot(actual_returns$rUSD, type="l", col="black", lwd=1,
     main="Actual vs Fitted - rUSD", ylab="Return")
lines(pred_returns$rUSD, col="green", lwd=1.5)
legend("topleft", legend=c("Actual","Fitted"), col=c("black","green"), lwd=2)
