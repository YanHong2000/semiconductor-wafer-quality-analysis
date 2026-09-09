# ============================================================
# 半導體晶圓製程變異與品質分析
# Semiconductor Wafer Process Variation and Quality Analysis
#
# 整理自原始統計品管課程專案程式。
# 內容保留：EDA、製程能力、管制圖、孿生變數、OGA、LASSO、
#          Two-stage OGA + HDIC + Trim。
# ============================================================


# ============================================================
# 0. 套件與基本設定
# ============================================================

library(dplyr)
library(ggplot2)
library(qcc)
library(igraph)
library(glmnet)

set.seed(42)


# ============================================================
# 1. 資料讀取與整理
# ============================================================

file_path <- file.path(
  "data",
  "new_wafer_data.csv"
)

df <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 移除 CSV 匯出時產生的流水號欄位
if (
  names(df)[1] %in% c("", "X", "X.1") &&
  is.numeric(df[[1]]) &&
  all(df[[1]] == seq_len(nrow(df)))
) {
  df <- df[, -1, drop = FALSE]
}

required_cols <- c("LOT", "WAFER", "OBSERVATION")

if (!all(required_cols %in% names(df))) {
  stop("資料缺少 LOT、WAFER 或 OBSERVATION 欄位。")
}

# LOT 依數字順序排列
lot_levels <- unique(df$LOT)
lot_number <- suppressWarnings(as.integer(sub("^LOT", "", lot_levels)))

if (all(!is.na(lot_number))) {
  lot_levels <- lot_levels[order(lot_number)]
}

df$LOT <- factor(df$LOT, levels = lot_levels)

# 找出 Txxxx 製程變數
process_cols <- grep(
  "^T[0-9]{4}$",
  names(df),
  value = TRUE
)

# 基本資料檢查
cat("資料維度：", nrow(df), "x", ncol(df), "\n")
cat("製程變數數量：", length(process_cols), "\n")
cat("缺失值數量：", sum(is.na(df)), "\n")
cat("重複資料列數：", sum(duplicated(df)), "\n")
print(table(df$LOT))

# 檢查 T 變數是否皆為 0/1
is_binary <- sapply(
  df[, process_cols, drop = FALSE],
  function(x) all(x %in% c(0, 1, NA))
)

if (!all(is_binary)) {
  warning(
    "以下製程變數不是純 0/1：",
    paste(process_cols[!is_binary], collapse = ", ")
  )
}

# 找出常數變數
constant_cols <- process_cols[
  sapply(
    df[, process_cols, drop = FALSE],
    function(x) length(unique(x[!is.na(x)])) <= 1
  )
]

cat("常數製程變數數量：", length(constant_cols), "\n")
print(constant_cols)


# ============================================================
# 2. 探索性資料分析
# ============================================================

# OBSERVATION 描述統計
print(summary(df$OBSERVATION))
cat("標準差：", sd(df$OBSERVATION), "\n")
cat("IQR：", IQR(df$OBSERVATION), "\n")

# 直方圖
p_hist <- ggplot(df, aes(x = OBSERVATION)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = mean(df$OBSERVATION),
    linetype = "dashed"
  ) +
  labs(
    title = "OBSERVATION 整體分布",
    x = "OBSERVATION",
    y = "樣本數"
  ) +
  theme_minimal()

print(p_hist)

# Q-Q plot
p_qq <- ggplot(df, aes(sample = OBSERVATION)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "OBSERVATION Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()

print(p_qq)

# Shapiro-Wilk 常態性檢定
shapiro_result <- shapiro.test(df$OBSERVATION)
print(shapiro_result)

# 各 LOT 描述統計
lot_summary <- df %>%
  group_by(LOT) %>%
  summarise(
    n = n(),
    mean = mean(OBSERVATION),
    median = median(OBSERVATION),
    sd = sd(OBSERVATION),
    min = min(OBSERVATION),
    max = max(OBSERVATION),
    .groups = "drop"
  )

print(lot_summary)

# 各 LOT 分布
p_lot_box <- ggplot(df, aes(x = LOT, y = OBSERVATION)) +
  geom_boxplot(fill = "lightblue") +
  labs(
    title = "各 LOT 的 OBSERVATION 分布",
    x = "LOT",
    y = "OBSERVATION"
  ) +
  theme_minimal()

print(p_lot_box)

# LOT 內 wafer sequence
p_wafer <- ggplot(
  df,
  aes(x = WAFER, y = OBSERVATION, group = 1)
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ LOT, scales = "free_x") +
  labs(
    title = "各 LOT 內 Wafer 的 OBSERVATION 變化",
    x = "WAFER",
    y = "OBSERVATION"
  ) +
  theme_minimal()

print(p_wafer)

# 檢查不同 LOT 的分布位置是否一致
kruskal_result <- kruskal.test(
  OBSERVATION ~ LOT,
  data = df
)

print(kruskal_result)


# ============================================================
# 3. 製程能力與規格
# ============================================================

obs <- df$OBSERVATION
mean_obs <- mean(obs)
sd_obs <- sd(obs)

cat("OBSERVATION 平均值：", round(mean_obs, 4), "\n")
cat("OBSERVATION 標準差：", round(sd_obs, 4), "\n")

# ------------------------------------------------------------
# 3.1 課堂假設規格：LSL = 12、USL = 20
# ------------------------------------------------------------

LSL <- 12
USL <- 20

Cp <- (USL - LSL) / (6 * sd_obs)

Cpk <- min(
  (USL - mean_obs) / (3 * sd_obs),
  (mean_obs - LSL) / (3 * sd_obs)
)

out_of_spec_rate <- mean(obs < LSL | obs > USL)

capability_assumed <- data.frame(
  LSL = LSL,
  USL = USL,
  Mean = mean_obs,
  SD = sd_obs,
  Cp = Cp,
  Cpk = Cpk,
  Out_of_spec_rate = out_of_spec_rate
)

print(capability_assumed)

# ------------------------------------------------------------
# 3.2 以 1% / 99% 分位數作為規格代理
# ------------------------------------------------------------

LSL_est <- as.numeric(quantile(obs, 0.01))
USL_est <- as.numeric(quantile(obs, 0.99))

Cp_quantile <- (USL_est - LSL_est) / (6 * sd_obs)

cat("1% / 99% quantile Cp：", round(Cp_quantile, 3), "\n")

# ------------------------------------------------------------
# 3.3 Mean ± 3 SD 自然界限
# ------------------------------------------------------------

LSL_natural <- mean_obs - 3 * sd_obs
USL_natural <- mean_obs + 3 * sd_obs

natural_out_rate <- mean(
  obs < LSL_natural |
    obs > USL_natural
)

cat(
  "Mean ± 3 SD 範圍外比例：",
  round(natural_out_rate * 100, 4),
  "%\n"
)


# ============================================================
# 4. 管制圖
# ============================================================

# ------------------------------------------------------------
# 4.1 Global Individuals Chart
# ------------------------------------------------------------

i_chart <- qcc(
  df$OBSERVATION,
  type = "xbar.one",
  title = "Individuals Chart：OBSERVATION"
)

cat("I Chart Center：", i_chart$center, "\n")
cat("I Chart Std.Dev：", i_chart$std.dev, "\n")
print(i_chart$limits)

# ------------------------------------------------------------
# 4.2 各 LOT 每 5 片 wafer 的平均變化
# ------------------------------------------------------------

df_avg <- df %>%
  arrange(LOT, WAFER) %>%
  group_by(LOT) %>%
  mutate(group_id = ceiling(row_number() / 5)) %>%
  group_by(LOT, group_id) %>%
  summarise(
    avg_obs = mean(OBSERVATION),
    n = n(),
    .groups = "drop"
  )

p_group_avg <- ggplot(
  df_avg,
  aes(
    x = factor(group_id),
    y = avg_obs,
    group = LOT,
    color = LOT
  )
) +
  geom_line() +
  geom_point() +
  labs(
    title = "各 LOT 內每 5 片 Wafer 的平均 OBSERVATION",
    x = "每 5 片分組",
    y = "平均 OBSERVATION"
  ) +
  theme_minimal()

print(p_group_avg)

# ------------------------------------------------------------
# 4.3 LOT 層級 X-bar、R、S Chart
# ------------------------------------------------------------

lot_groups <- split(
  df$OBSERVATION,
  df$LOT
)

lot_groups <- lot_groups[
  lengths(lot_groups) > 1
]

subgroup_sizes <- lengths(lot_groups)
cat("各 LOT subgroup size：", subgroup_sizes, "\n")

# 不同 LOT 樣本數不足處補 NA
max_size <- max(subgroup_sizes)

group_matrix <- t(
  sapply(
    lot_groups,
    function(x) {
      length(x) <- max_size
      x
    }
  )
)

rownames(group_matrix) <- names(lot_groups)

xbar_chart <- qcc(
  group_matrix,
  type = "xbar",
  title = "X-bar Chart for LOT Means"
)

r_chart <- qcc(
  group_matrix,
  type = "R",
  title = "R Chart for LOT Ranges"
)

s_chart <- qcc(
  group_matrix,
  type = "S",
  title = "S Chart for LOT Standard Deviations"
)


# ============================================================
# 5. 孿生製程變數
# ============================================================

# 移除常數欄位後計算 Pearson correlation
param_data <- df[
  ,
  setdiff(process_cols, constant_cols),
  drop = FALSE
]

cor_matrix <- cor(
  param_data,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# 課堂專案使用 0.8 < |r| < 1 定義高度相關變數
twin_threshold <- 0.8

high_cor_pairs <- which(
  abs(cor_matrix) > twin_threshold &
    abs(cor_matrix) < 1,
  arr.ind = TRUE
)

high_cor_df <- data.frame(
  Param1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
  Param2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
  Correlation = cor_matrix[high_cor_pairs],
  stringsAsFactors = FALSE
)

cat("高度相關 pair 數量：", nrow(high_cor_df), "\n")
print(head(high_cor_df, 10))

if (nrow(high_cor_df) == 0) {
  stop("未找到高度相關製程變數，無法建立孿生變數群集。")
}

# 建立孿生變數群集
g <- graph_from_data_frame(
  high_cor_df[, c("Param1", "Param2")],
  directed = FALSE
)

clusters <- components(g)

cat("孿生變數群集數量：", clusters$no, "\n")
print(summary(clusters$csize))

cluster_membership <- clusters$membership

# 建立完整 X 與 y
X <- as.matrix(
  df[, process_cols, drop = FALSE]
)

y <- df$OBSERVATION

# 計算各製程變數與 OBSERVATION 的相關性
cor_with_y <- apply(
  X,
  2,
  function(x) {
    suppressWarnings(
      cor(
        x,
        y,
        use = "pairwise.complete.obs",
        method = "pearson"
      )
    )
  }
)

names(cor_with_y) <- colnames(X)

# 每個孿生群集選擇與 OBSERVATION 絕對相關最高的代表變數
representative_vars <- character(clusters$no)

for (i in seq_len(clusters$no)) {
  cluster_params <- names(
    cluster_membership[
      cluster_membership == i
    ]
  )

  cluster_cor <- abs(
    cor_with_y[cluster_params]
  )

  representative_vars[i] <- cluster_params[
    which.max(cluster_cor)
  ]
}

# 加入不屬於孿生群集的變數
independent_vars <- setdiff(
  process_cols,
  names(cluster_membership)
)

selected_vars <- unique(
  c(
    representative_vars,
    independent_vars
  )
)

X_selected <- as.matrix(
  df[, selected_vars, drop = FALSE]
)

cat(
  "孿生群集代表變數：",
  length(representative_vars),
  "\n"
)

cat(
  "獨立變數：",
  length(independent_vars),
  "\n"
)

cat(
  "群集整理後變數總數：",
  length(selected_vars),
  "\n"
)


# ============================================================
# 6. OGA
# ============================================================

# OGA 依序挑選與目前 residual 關聯最大的變數
OGA <- function(X, y, Kn = NULL, c1 = 5) {

  if (!is.vector(y)) stop("y should be a vector")
  if (!is.matrix(X)) stop("X should be a matrix")

  n <- nrow(X)
  p <- ncol(X)

  if (n != length(y)) {
    stop("X 與 y 的樣本數不一致。")
  }

  if (is.null(Kn)) {
    K <- max(
      1,
      min(
        floor(c1 * sqrt(n / log(p))),
        p
      )
    )
  } else {
    K <- Kn
  }

  dy <- y - mean(y)

  dX <- apply(
    X,
    2,
    function(x) x - mean(x)
  )

  Jhat <- rep(0, K)
  sigma2hat <- rep(0, K)
  XJhat <- matrix(0, n, K)

  u <- as.matrix(dy)
  xnorms <- sqrt(colSums(dX^2))

  aSSE <- abs(t(u) %*% dX) / xnorms
  aSSE[!is.finite(aSSE)] <- -Inf

  Jhat[1] <- which.max(aSSE)

  XJhat[, 1] <- dX[, Jhat[1]] /
    sqrt(sum(dX[, Jhat[1]]^2))

  u <- u -
    XJhat[, 1] %*%
    t(XJhat[, 1]) %*%
    u

  sigma2hat[1] <- mean(u^2)

  if (K > 1) {
    for (k in 2:K) {

      aSSE <- abs(t(u) %*% dX) / xnorms
      aSSE[!is.finite(aSSE)] <- -Inf
      aSSE[Jhat[1:(k - 1)]] <- -Inf

      Jhat[k] <- which.max(aSSE)

      rq <- dX[, Jhat[k]] -
        XJhat[, 1:(k - 1), drop = FALSE] %*%
        t(XJhat[, 1:(k - 1), drop = FALSE]) %*%
        dX[, Jhat[k]]

      rq_norm <- sqrt(sum(rq^2))

      if (!is.finite(rq_norm) || rq_norm == 0) {
        break
      }

      XJhat[, k] <- rq / rq_norm

      u <- u -
        XJhat[, k] %*%
        t(XJhat[, k]) %*%
        u

      sigma2hat[k] <- mean(u^2)
    }
  }

  valid <- Jhat > 0

  list(
    n = n,
    p = p,
    Kn = sum(valid),
    J_OGA = Jhat[valid],
    sigma2hat = sigma2hat[valid]
  )
}

# ------------------------------------------------------------
# 6.1 OGA：不考慮孿生變數
# ------------------------------------------------------------

oga_raw <- OGA(
  X,
  y
)

oga_raw_vars <- colnames(X)[
  oga_raw$J_OGA
]

cat(
  "Raw OGA 選出變數：\n",
  paste(oga_raw_vars, collapse = ", "),
  "\n"
)

# ------------------------------------------------------------
# 6.2 OGA：考慮孿生變數
# ------------------------------------------------------------

oga_twin <- OGA(
  X_selected,
  y
)

oga_twin_vars <- colnames(X_selected)[
  oga_twin$J_OGA
]

cat(
  "Twin-aware OGA 選出變數：\n",
  paste(oga_twin_vars, collapse = ", "),
  "\n"
)


# ============================================================
# 7. LASSO
# ============================================================

# 使用 10-fold cross-validation 選擇 lambda
set.seed(42)

cvfit <- cv.glmnet(
  x = X,
  y = y,
  alpha = 1,
  nfolds = 10
)

# 繪製交叉驗證曲線
# 增加上方邊界，避免標題與非零係數數量重疊
old_par <- par(no.readonly = TRUE)

par(
  mar = c(5.1, 5.1, 7.0, 2.1)
)

plot(
  cvfit,
  main = ""
)

title(
  main = "LASSO 10-fold Cross-Validation",
  line = 4.5,
  cex.main = 1.2
)

par(old_par)


# 取得 lambda.min 與 lambda.1se
lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se


# 取得兩個 lambda 下的模型係數
coef_min <- coef(
  cvfit,
  s = "lambda.min"
)

coef_1se <- coef(
  cvfit,
  s = "lambda.1se"
)


# 找出非零係數變數
vars_min <- rownames(coef_min)[-1][
  as.vector(
    coef_min[-1, 1] != 0
  )
]

vars_1se <- rownames(coef_1se)[-1][
  as.vector(
    coef_1se[-1, 1] != 0
  )
]


# 輸出 lambda
cat(
  "lambda.min：",
  lambda_min,
  "\n"
)

cat(
  "lambda.1se：",
  lambda_1se,
  "\n"
)


# 輸出選取變數
cat(
  "lambda.min 選到",
  length(vars_min),
  "個變數：\n",
  paste(
    vars_min,
    collapse = ", "
  ),
  "\n"
)

cat(
  "lambda.1se 選到",
  length(vars_1se),
  "個變數：\n",
  paste(
    vars_1se,
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 8. Two-stage OGA + HDIC + Trim
# ============================================================

# ------------------------------------------------------------
# 8.1 OGA + HDIC + Trim
# ------------------------------------------------------------

Ohit2011 <- function(
  X,
  y,
  Kn = NULL,
  c1 = 5,
  HDIC_Type = "HDHQ",
  const2 = seq(0.5, 2, 0.5),
  const3 = seq(0.2, 1, 0.2)
) {

  n <- nrow(X)
  p <- ncol(X)

  K <- ifelse(
    is.null(Kn),
    max(1, min(floor(c1 * sqrt(n / log(p))), p)),
    Kn
  )

  dy <- y - mean(y)

  dX <- apply(
    X,
    2,
    function(x) x - mean(x)
  )

  Jhat <- rep(0, K)
  sigma2hat <- rep(0, K)
  XJhat <- matrix(0, n, K)

  u <- as.matrix(dy)
  xnorms <- sqrt(colSums(dX^2))

  aSSE <- abs(t(u) %*% dX) / xnorms
  aSSE[!is.finite(aSSE)] <- -Inf

  Jhat[1] <- which.max(aSSE)

  XJhat[, 1] <- dX[, Jhat[1]] /
    sqrt(sum(dX[, Jhat[1]]^2))

  u <- u -
    XJhat[, 1] %*%
    t(XJhat[, 1]) %*%
    u

  sigma2hat[1] <- mean(u^2)

  if (K > 1) {
    for (k in 2:K) {

      aSSE <- abs(t(u) %*% dX) / xnorms
      aSSE[!is.finite(aSSE)] <- -Inf
      aSSE[Jhat[1:(k - 1)]] <- -Inf

      Jhat[k] <- which.max(aSSE)

      rq <- dX[, Jhat[k]] -
        XJhat[, 1:(k - 1), drop = FALSE] %*%
        t(XJhat[, 1:(k - 1), drop = FALSE]) %*%
        dX[, Jhat[k]]

      rq_norm <- sqrt(sum(rq^2))

      if (!is.finite(rq_norm) || rq_norm == 0) {
        break
      }

      XJhat[, k] <- rq / rq_norm

      u <- u -
        XJhat[, k] %*%
        t(XJhat[, k]) %*%
        u

      sigma2hat[k] <- mean(u^2)
    }
  }

  valid <- Jhat > 0

  Jhat <- Jhat[valid]
  sigma2hat <- sigma2hat[valid]
  K <- length(Jhat)

  if (HDIC_Type == "HDAIC") penalty <- 1
  if (HDIC_Type == "HDBIC") penalty <- log(n)
  if (HDIC_Type == "HDHQ") penalty <- log(log(n))

  kn_hat_now <- 0
  result_all <- NULL

  for (c2 in const2) {

    omega_n <- c2 * penalty

    hdic <- n * log(sigma2hat) +
      (1:K) * omega_n * log(p)

    kn_hat <- which.min(hdic)

    if (kn_hat_now != kn_hat) {

      for (c3 in const3) {

        benchmark <- n * log(sigma2hat[kn_hat]) +
          kn_hat * c3 * penalty * log(p)

        J_Trim <- Jhat[1:kn_hat]
        trim_pos <- rep(0, kn_hat)

        if (kn_hat > 1) {

          for (l in seq_len(kn_hat)) {

            JDrop1 <- J_Trim[-l]

            fit <- lm(
              dy ~ . - 1,
              data = data.frame(
                dX[, JDrop1, drop = FALSE]
              )
            )

            uDrop1 <- fit$residuals

            HDICDrop1 <- n * log(mean(uDrop1^2)) +
              (kn_hat - 1) * c3 * penalty * log(p)

            if (HDICDrop1 > benchmark) {
              trim_pos[l] <- 1
            }
          }

          J_Trim <- J_Trim[trim_pos == 1]
        }

        result <- c(
          sort(J_Trim),
          rep(0, K - length(J_Trim))
        )

        result_all <- rbind(
          result_all,
          result
        )
      }

      kn_hat_now <- kn_hat
    }
  }

  list(
    OGA = Jhat,
    Trim = result_all
  )
}


# ------------------------------------------------------------
# 8.2 Two-stage mean / dispersion model
# ------------------------------------------------------------

Twohit <- function(
  X,
  y,
  Kn = NULL,
  c1 = 5,
  HDIC_Type = "HDHQ",
  const2 = seq(0.5, 2, 0.5),
  const3 = seq(0.2, 1, 0.2)
) {

  n <- nrow(X)
  p <- ncol(X)

  benchmark <- Inf
  Trim_reg <- NULL
  Trim_dis <- NULL
  Beta <- NULL
  Alpha <- NULL

  vs_reg <- Ohit2011(
    X,
    y,
    Kn,
    c1,
    HDIC_Type,
    const2,
    const3
  )

  variable_all_reg <- vs_reg$Trim

  for (i in seq_len(nrow(variable_all_reg))) {

    variable_reg <- variable_all_reg[
      i,
      variable_all_reg[i, ] != 0
    ]

    nowX_reg <- cbind(
      1,
      X[, variable_reg, drop = FALSE]
    )

    u <- lm(
      y ~ nowX_reg - 1
    )$residuals

    pos <- which(u^2 < 1e-8)

    if (length(pos) > 0) {
      u[pos] <- sign(u[pos]) * 1e-4
    }

    newU <- log(u^2)

    vs_dis <- Ohit2011(
      X,
      newU,
      Kn,
      c1,
      HDIC_Type,
      const2,
      const3
    )

    variable_all_dis <- vs_dis$Trim

    for (j in seq_len(nrow(variable_all_dis))) {

      variable_dis <- variable_all_dis[
        j,
        variable_all_dis[j, ] != 0
      ]

      nowX_dis <- cbind(
        1,
        X[, variable_dis, drop = FALSE]
      )

      if (
        length(variable_reg) == 0 &&
        length(variable_dis) == 0
      ) {

        L <- function(theta) {
          A <- rep(1, n)
          B <- rep(1, n)

          sum(B * theta[2]) +
            sum(
              (y - A * theta[1])^2 /
                exp(B * theta[2])
            )
        }

        initial_estimate <- c(
          mean(y),
          mean(newU)
        )
      }

      if (
        length(variable_reg) > 0 &&
        length(variable_dis) == 0
      ) {

        L <- function(theta) {
          A <- nowX_reg
          B <- rep(1, n)

          beta_index <- 1:(length(variable_reg) + 1)
          alpha_index <- length(variable_reg) + 2

          sum(B * theta[alpha_index]) +
            sum(
              (y - A %*% theta[beta_index])^2 /
                exp(B * theta[alpha_index])
            )
        }

        initial_estimate <- c(
          as.vector(
            lm(y ~ nowX_reg - 1)$coefficients
          ),
          mean(newU)
        )
      }

      if (
        length(variable_reg) == 0 &&
        length(variable_dis) > 0
      ) {

        L <- function(theta) {
          A <- rep(1, n)
          B <- nowX_dis

          sum(B %*% theta[2:length(theta)]) +
            sum(
              (y - A * theta[1])^2 /
                exp(B %*% theta[2:length(theta)])
            )
        }

        initial_estimate <- c(
          mean(y),
          as.vector(
            lm(newU ~ nowX_dis - 1)$coefficients
          )
        )
      }

      if (
        length(variable_reg) > 0 &&
        length(variable_dis) > 0
      ) {

        L <- function(theta) {
          A <- nowX_reg
          B <- nowX_dis

          beta_index <- 1:(length(variable_reg) + 1)
          alpha_index <- (length(variable_reg) + 2):length(theta)

          sum(B %*% theta[alpha_index]) +
            sum(
              (y - A %*% theta[beta_index])^2 /
                exp(B %*% theta[alpha_index])
            )
        }

        initial_estimate <- c(
          as.vector(
            lm(y ~ nowX_reg - 1)$coefficients
          ),
          as.vector(
            lm(newU ~ nowX_dis - 1)$coefficients
          )
        )
      }

      paraEst <- try(
        optim(
          initial_estimate,
          L,
          method = "L-BFGS-B"
        ),
        silent = TRUE
      )

      if ("try-error" %in% class(paraEst)) {
        next
      }

      hdicNow <- L(paraEst$par) +
        (
          length(variable_reg) +
            length(variable_dis)
        ) *
        log(2 * p) *
        log(log(n)) *
        2

      if (hdicNow < benchmark) {

        Trim_reg <- variable_reg
        Trim_dis <- variable_dis

        Beta <- paraEst$par[
          1:(length(variable_reg) + 1)
        ]

        Alpha <- paraEst$par[
          (length(variable_reg) + 2):
            (
              length(variable_reg) +
                length(variable_dis) +
                2
            )
        ]

        benchmark <- hdicNow
      }
    }
  }

  list(
    Trim_reg = Trim_reg,
    Trim_dis = Trim_dis,
    Beta = Beta,
    Alpha = Alpha
  )
}


# ------------------------------------------------------------
# 8.3 執行 Two-stage OGA + HDIC + Trim
# ------------------------------------------------------------

twohit_result <- Twohit(
  X,
  y
)

mean_vars <- colnames(X)[
  twohit_result$Trim_reg
]

dispersion_vars <- colnames(X)[
  twohit_result$Trim_dis
]

cat(
  "Mean model 重要變數：",
  paste(mean_vars, collapse = ", "),
  "\n"
)

cat(
  "Dispersion model 重要變數：",
  paste(dispersion_vars, collapse = ", "),
  "\n"
)

cat("Mean model coefficients (Beta)：\n")
print(twohit_result$Beta)

cat("Dispersion model coefficients (Alpha)：\n")
print(twohit_result$Alpha)


# ------------------------------------------------------------
# 8.4 Two-stage model 圖
# ------------------------------------------------------------

if (
  length(mean_vars) == 2 &&
  length(dispersion_vars) == 1
) {

  df_plot <- df

  df_plot$index <- seq_len(
    nrow(df_plot)
  )

  df_plot$y_hat <-
    twohit_result$Beta[1] +
    twohit_result$Beta[2] *
    df_plot[[mean_vars[1]]] +
    twohit_result$Beta[3] *
    df_plot[[mean_vars[2]]]

  df_plot$sigma_hat <- exp(
    (
      twohit_result$Alpha[1] +
      twohit_result$Alpha[2] *
      df_plot[[dispersion_vars[1]]]
    ) / 2
  )

  df_plot$upper <-
    df_plot$y_hat +
    2 * df_plot$sigma_hat

  df_plot$lower <-
    df_plot$y_hat -
    2 * df_plot$sigma_hat

  p_twohit <- ggplot(
    df_plot,
    aes(
      x = index,
      y = OBSERVATION,
      color = LOT
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.8
    ) +
    geom_step(
      aes(y = y_hat),
      color = "black",
      linewidth = 0.9
    ) +
    geom_step(
      aes(y = upper),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_step(
      aes(y = lower),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    labs(
      title = "Two-stage OGA + HDIC + Trim",
      x = "Wafer Sequence",
      y = "OBSERVATION"
    ) +
    theme_minimal()

  print(p_twohit)
}


# ============================================================
# 9. 主要結果摘要
# ============================================================

cat("\n================ 主要結果摘要 ================\n")
cat("Raw OGA 變數數量：", length(oga_raw_vars), "\n")
cat("Twin-aware OGA 變數數量：", length(oga_twin_vars), "\n")
cat("LASSO lambda.min 變數數量：", length(vars_min), "\n")
cat("LASSO lambda.1se 變數數量：", length(vars_1se), "\n")
cat("Two-stage Mean variables：", paste(mean_vars, collapse = ", "), "\n")
cat(
  "Two-stage Dispersion variables：",
  paste(dispersion_vars, collapse = ", "),
  "\n"
)
cat("================================================\n")
