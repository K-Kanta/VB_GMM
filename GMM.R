# *必要なライブラリのインストール ####
library(mvnfast)
library(gganimate)
library(dplyr)



# *データの読み込みと可視化 ####
data <- read.csv("GMM/data/label.csv", header = FALSE) %>%
  rename(x1 = "V1", x2 = "V2")


x_nd <- as.matrix(read.csv("GMM/data/train.csv", header = FALSE) %>% 
                    rename(x1 = "V1", x2 = "V2")) #データを行列に型変換





# *事前分布の設定 ####
# 次元数
D <- 2

# クラスター数
K <- 4

# データ数
N <- nrow(x_nd) 

# muの事前分布のパラメータを指定
beta <- 1
m_d <- rep(0, D)

# lambdaの事前分布のパラメータを指定
w_dd <- diag(D)
nu <- 100

# piの事前分布のパラメータを指定
alpha_k <- rep(1, K)



# *初期値の設定 ####
# 潜在変数の近似事後分布の期待値を初期化
E_s_nk <- runif(n = N * K, min = 0, max = 1) %>% 
  matrix(nrow = N, ncol = K)
E_s_nk <- E_s_nk / rowSums(E_s_nk) # 正規化

# 初期値によるmuの近似事後分布のパラメータを計算:式(4.114)
beta_hat_k <- colSums(E_s_nk) + beta
m_hat_kd <- t(t(x_nd) %*% E_s_nk + beta * m_d) / beta_hat_k

# 初期値によるlambdaの近似事後分布のパラメータを計算:式(4.118)
w_hat_ddk <- array(0, dim = c(D, D, K))
term_m_dd <- beta * matrix(m_d) %*% t(m_d)
for(k in 1:K) {
  term_x_dd <- t(E_s_nk[, k] * x_nd) %*% x_nd
  term_m_hat_dd <- beta_hat_k[k] * matrix(m_hat_kd[k, ]) %*% t(m_hat_kd[k, ])
  w_hat_ddk[, , k] <- solve(term_x_dd + term_m_dd - term_m_hat_dd + solve(w_dd))
}
nu_hat_k <- colSums(E_s_nk) + nu

# 初期値によるpiの近似事後分布のパラメータを計算:式(4.58)
alpha_hat_k <- colSums(E_s_nk) + alpha_k



# *初期値による混合分布####
# 軸の作成
x1_vec <- seq(-2, 2, length.out = 40)
x2_vec <- seq(-2, 2, length.out = 40)
x_point_mat <- cbind(
  rep(x1_vec, times = length(x2_vec)),
  rep(x2_vec, each = length(x1_vec))
) 

# 確率密度の計算
init_density <- 0
for(k in 1:K) {
  # クラスタkの分布の確率密度を計算
  tmp_density <- mvnfast::dmvn(
    X = x_point_mat, mu = m_hat_kd[k, ], sigma = solve(nu_hat_k[k] * w_hat_ddk[, , k])
  )
  
  # K個の分布の加重平均を計算
  init_density <- init_density + alpha_hat_k[k] / sum(alpha_hat_k) * tmp_density
}

# 初期値による分布をデータフレームに格納
init_df <- tibble(
  x_1 = x_point_mat[, 1], 
  x_2 = x_point_mat[, 2], 
  density = init_density
)

# 初期値による分布を作図
ggplot() + 
  geom_contour(data = init_df, aes(x = x_1, y = x_2, z = density, color = ..level..)) + # 期待値による分布
  labs(title = "Gaussian Mixture Model", 
       subtitle = paste0("iter:", 0, ", K=", K), 
       x = expression(x[1]), y = expression(x[2]))

# 観測データと確率が最大のクラスタ番号をデータフレームに格納
s_df <- tibble(
  x_n1 = x_nd[, 1], 
  x_n2 = x_nd[, 2], 
  cluster = as.factor(max.col(E_s_nk)), # 確率が最大のクラスタ番号
  prob = E_s_nk[cbind(1:N, max.col(E_s_nk))] # 確率の最大値
)

# クラスタの初期値を作図
ggplot() + 
  geom_contour_filled(data = init_df, aes(x = x_1, y = x_2, z = density, fill = ..level..), 
                      alpha = 0.6) + # 期待値による分布
  geom_point(data = s_df, aes(x = x_n1, y = x_n2, color = cluster), 
             alpha = s_df[["prob"]]) + # 確率値によるクラスタ
  labs(title = "Gaussian Mixture Model", 
       subtitle = paste0("iter:", 0, ", N=", N, ", K=", K), 
       x = expression(x[1]), y = expression(x[2]))



# *推論 ####
# 試行回数を指定
MaxIter <- 100

# 途中計算に用いる変数を作成
ln_eta_nk <- matrix(0, nrow = N, ncol = K)

# 推移の確認用の受け皿を作成
trace_E_s_ink  <- array(0, dim = c(MaxIter + 1, N, K))
trace_beta_ik  <- matrix(0, nrow = MaxIter + 1, ncol = K)
trace_m_ikd    <- array(0, dim = c(MaxIter + 1, K, D))
trace_w_iddk   <- array(0, dim = c(MaxIter + 1, D, D, K))
trace_nu_ik    <- matrix(0, nrow = MaxIter + 1, ncol = K)
trace_alpha_ik <- matrix(0, nrow = MaxIter + 1, ncol = K)

# 初期値を記録
trace_E_s_ink[1, , ]  <- E_s_nk
trace_beta_ik[1, ]    <- beta_hat_k
trace_m_ikd[1, , ]    <- m_hat_kd
trace_w_iddk[1, , , ] <- w_hat_ddk
trace_nu_ik[1, ]      <- nu_hat_k
trace_alpha_ik[1, ]   <- alpha_hat_k

# 変分推論
for(i in 1:MaxIter) {
  # 潜在変数の近似事後分布のパラメータを計算:式(4.109)
  for(k in 1:K) {
    # クラスタkの中間変数を計算:式(4.119-4.122,4.62)
    E_lmd_dd <- nu_hat_k[k] * w_hat_ddk[, , k]
    E_ln_det_lmd <- sum(digamma(0.5 * (nu_hat_k[k] + 1 - 1:D))) + D * log(2) + log(det(w_hat_ddk[, , k]))
    E_lmd_mu_d1 <- E_lmd_dd %*% matrix(m_hat_kd[k, ])
    E_mu_lmd_mu <- (t(m_hat_kd[k, ]) %*% E_lmd_mu_d1 + D / beta_hat_k[k]) %>% 
      as.vector()
    E_ln_pi <- digamma(alpha_hat_k[k]) - digamma(sum(alpha_hat_k))
    term_x1_n <- - 0.5 * x_nd %*% E_lmd_dd %*% t(x_nd) %>% 
      diag()
    term_x2_n <- x_nd %*% E_lmd_mu_d1 %>% 
      as.vector()
    ln_eta_nk[, k] <- term_x1_n + term_x2_n - 0.5 * E_mu_lmd_mu + 0.5 * E_ln_det_lmd + E_ln_pi
  }
  tmp_eta_nk <- exp(ln_eta_nk)
  eta_nk <- (tmp_eta_nk + 1e-7) / rowSums(tmp_eta_nk + 1e-7) # 正規化
  
  # 潜在変数の近似事後分布の期待値を計算:式(4.59)
  E_s_nk <- eta_nk
  
  # 初期値によるmuの近似事後分布のパラメータを計算:式(4.114)
  beta_hat_k <- colSums(E_s_nk) + beta
  m_hat_kd <- t(t(x_nd) %*% E_s_nk + beta * m_d) / beta_hat_k
  
  # 初期値によるlambdaの近似事後分布のパラメータを計算:式(4.118)
  w_hat_ddk <- array(0, dim = c(D, D, K))
  term_m_dd <- beta * matrix(m_d) %*% t(m_d)
  for(k in 1:K) {
    term_x_dd <- t(E_s_nk[, k] * x_nd) %*% x_nd
    term_m_hat_dd <- beta_hat_k[k] * matrix(m_hat_kd[k, ]) %*% t(m_hat_kd[k, ])
    w_hat_ddk[, , k] <- solve(term_x_dd + term_m_dd - term_m_hat_dd + solve(w_dd))
  }
  nu_hat_k <- colSums(E_s_nk) + nu
  
  # piの近似事後分布のパラメータを計算:式(4.58)
  alpha_hat_k <- colSums(E_s_nk) + alpha_k
  
  # i回目のパラメータを記録
  trace_E_s_ink[i + 1, , ]  <- E_s_nk
  trace_beta_ik[i + 1, ]    <- beta_hat_k
  trace_m_ikd[i + 1, , ]    <- m_hat_kd
  trace_w_iddk[i + 1, , , ] <- w_hat_ddk
  trace_nu_ik[i + 1, ]      <- nu_hat_k
  trace_alpha_ik[i + 1, ]   <- alpha_hat_k
  
  # 動作確認
  print(paste0(i, ' (', round(i / MaxIter * 100, 1), '%)'))
}

# *結果 ####
# 最後に更新したパラメータの期待値による混合分布を計算
res_density <- 0
for(k in 1:K) {
  # クラスタkの分布の確率密度を計算
  tmp_density <- mvnfast::dmvn(
    X = x_point_mat, mu = m_hat_kd[k, ], sigma = solve(nu_hat_k[k] * w_hat_ddk[, , k])
  )
  
  # K個の分布の加重平均を計算
  res_density <- res_density + alpha_hat_k[k] / sum(alpha_hat_k) * tmp_density
}

# 最終的な分布をデータフレームに格納
res_df <- tibble(
  x_1 = x_point_mat[, 1], 
  x_2 = x_point_mat[, 2], 
  density = res_density
)

s_df <- tibble(
  x_n1 = x_nd[, 1], 
  x_n2 = x_nd[, 2], 
  cluster = as.factor(max.col(E_s_nk)), # 確率が最大のクラスタ番号
  prob = E_s_nk[cbind(1:N, max.col(E_s_nk))] # 確率の最大値
)

# 最終的なクラスタを作図
ggplot() + 
  geom_contour_filled(data = res_df, aes(x = x_1, y = x_2, z = density, fill = ..level..), 
                      alpha = 0.6) + # 期待値による分布
  geom_point(data = s_df, aes(x = x_n1, y = x_n2, color = cluster), 
             alpha = s_df[["prob"]]) + # 確率値によるクラスタ
  labs(title = "Gaussian Mixture Model:Variational Inference", 
       subtitle = paste0("iter:", MaxIter, ", N=", N, 
                         ", E[pi]=(", paste0(round(alpha_hat_k / sum(alpha_hat_k), 3), collapse = ", "), ")"), 
       x = expression(x[1]), y = expression(x[2]))

# *クラスターの推移 ####
# 作図用のデータフレームを作成
trace_model_df <- tibble()
trace_cluster_df <- tibble()
for(i in 1:(MaxIter + 1)) {
  # i回目の混合分布を計算
  res_density <- 0
  for(k in 1:K) {
    # クラスタkの分布の確率密度を計算
    tmp_density <- mvnfast::dmvn(
      X = x_point_mat, 
      mu = trace_m_ikd[i, k, ], 
      sigma = solve(trace_nu_ik[i, k] * trace_w_iddk[i, , , k])
    )
    
    # K個の分布の加重平均を計算
    res_density <- res_density + trace_alpha_ik[i, k] / sum(trace_alpha_ik[i, ]) * tmp_density
  }
  
  # i回目の分布をデータフレームに格納
  res_df <- tibble(
    x_1 = x_point_mat[, 1], 
    x_2 = x_point_mat[, 2], 
    density = res_density, 
    iteration = as.factor(i - 1), 
    label = paste0(
      "iter:", i - 1, ", N=", N, 
      ", pi=(", paste0(round(trace_alpha_ik[i, ] / sum(trace_alpha_ik[i, ]), 3), collapse = ", "), ")"
    ) %>% 
      as.factor()
  )
  
  # 結果を結合
  trace_model_df <- rbind(trace_model_df, res_df)
  
  # 観測データとi回目のクラスタ番号をデータフレームに格納
  tmp_E_s_nk <- trace_E_s_ink[i, , ] # i回目の結果
  s_df <- tibble(
    x_n1 = x_nd[, 1], 
    x_n2 = x_nd[, 2], 
    cluster = as.factor(max.col(tmp_E_s_nk)), # 確率が最大のクラスタ番号
    prob = tmp_E_s_nk[cbind(1:N, max.col(tmp_E_s_nk))], # 確率の最大値
    iteration = as.factor(i - 1), 
    label = paste0(
      "iter:", i - 1, ", N=", N, 
      ", pi=(", paste0(round(trace_alpha_ik[i, ] / sum(trace_alpha_ik[i, ]), 3), collapse = ", "), ")"
    ) %>% 
      as.factor()
  )
  
  # 結果を結合
  trace_cluster_df <- rbind(trace_cluster_df, s_df)
  
  # 動作確認
  print(paste0(i - 1, ' (', round((i - 1) / MaxIter * 100, 1), '%)'))
}

# クラスタの推移を作図
trace_cluster_graph <- ggplot() + 
  geom_contour_filled(data = trace_model_df, aes(x = x_1, y = x_2, z = density), 
                      alpha = 0.6) + # 期待値による分布
  geom_point(data = trace_cluster_df, aes(x = x_n1, y = x_n2, color = cluster), 
             alpha = trace_cluster_df[["prob"]]) + # 確率が最大のクラスタ
  gganimate::transition_manual(label) + # フレーム
  labs(title = "Gaussian Mixture Model:Variational Inference", 
       subtitle = "{current_frame}", 
       x = expression(x[1]), y = expression(x[2]))

# gif画像を作成


model_animation <- gganimate::animate(trace_cluster_graph, nframes = MaxIter + 1, fps = 10)
anim_save(model_animation, filename = "GMM/result/VB_GMM.mp4")



