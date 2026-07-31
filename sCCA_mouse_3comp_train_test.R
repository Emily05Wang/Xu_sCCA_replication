#!/usr/bin/env Rscript

# ============================================================
# Nutrimumouse sCCA 教学复现
# 版本B：划分 train/test
# 软件环境：mixOmics 6.32.0 + nscancor
#
# 分析目标：
# 1. 使用 gene 作为 X 数据块、lipid 作为 Y 数据块；
# 2. 将小鼠划分为 train 和 test；
# 3. 只在 train 中学习 canonical weights；
# 4. test 使用 train weights 计算 canonical scores；
# 5. 依次拟合 3 个 sCCA component；
# 6. 输出 train/test canonical r、weights、scores、
#    loadings 和 cross-loadings。
#
# 重要说明：
# - 标准化参数只从 train 估计，再应用到 test；
# - 为尽量接近论文公开代码，后续 component 的 test deflation
#   仍在 test 矩阵上调用 acor()，但使用 train weights；
# - 因此它是“paper-like 教学版”，不是完全冻结的
#   train-only projection mapping；
# - Nutrimumouse 只有 40 只小鼠，test canonical r 波动会较大。
# ============================================================


# ------------------------------------------------------------
# 0. 加载分析包
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(mixOmics)  # 提供 wrapper.sgcca()
  library(nscancor)  # 提供 acor()，用于 projection deflation
})


# ------------------------------------------------------------
# 1. 设置随机种子和主要参数
# ------------------------------------------------------------
# 固定随机种子，使 train/test 划分可以重复。
set.seed(20260729)

# 最终依次拟合 3 个 component。
ncomp <- 3

# 75% 样本进入 train，25% 进入 test。
train_fraction <- 0.75

# 每个 component 在 X 端保留 10 个 gene variables。
keepX.X <- 10

# 每个 component 在 Y 端保留 5 个 lipid variables。
keepX.Y <- 5

# 设置输出文件夹。
output_dir <- "sCCA_nutrimouse_train_test_output"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 2. 读取 Nutrimumouse 数据并构建 X/Y
# ------------------------------------------------------------
data(nutrimouse)

# X 数据块：40 只小鼠 × 120 个 gene variables。
X_raw <- as.matrix(nutrimouse$gene)

# Y 数据块：同样 40 只小鼠 × 21 个 lipid variables。
Y_raw <- as.matrix(nutrimouse$lipid)

storage.mode(X_raw) <- "double"
storage.mode(Y_raw) <- "double"

# 检查 X/Y 是否样本一一对应，并确认没有 NA。
stopifnot(
  nrow(X_raw) == nrow(Y_raw),
  !anyNA(X_raw),
  !anyNA(Y_raw)
)


# ------------------------------------------------------------
# 3. 随机划分 train/test
# ------------------------------------------------------------
# 总样本数。
n <- nrow(X_raw)

# 随机抽取约 75% 样本作为 train。
train_id <- sort(
  sample(
    seq_len(n),
    size = floor(train_fraction * n),
    replace = FALSE
  )
)

# 剩余样本作为 test。
test_id <- setdiff(
  seq_len(n),
  train_id
)

# 根据相同的行号同步拆分 X 和 Y，
# 保证每一行仍对应同一只小鼠。
X_train_raw <- X_raw[train_id, , drop = FALSE]
Y_train_raw <- Y_raw[train_id, , drop = FALSE]

X_test_raw <- X_raw[test_id, , drop = FALSE]
Y_test_raw <- Y_raw[test_id, , drop = FALSE]


# ------------------------------------------------------------
# 4. 只用 train 估计标准化参数
# ------------------------------------------------------------
# X block 的 train 均值和标准差。
X_center <- colMeans(X_train_raw)
X_scale <- apply(X_train_raw, 2, sd)

# Y block 的 train 均值和标准差。
Y_center <- colMeans(Y_train_raw)
Y_scale <- apply(Y_train_raw, 2, sd)

# 检查 train 中是否存在标准差为 0 或无效的变量。
if (any(!is.finite(X_scale) | X_scale == 0)) {
  stop("Train 中至少有一个 X variable 的标准差为 0 或无效。")
}

if (any(!is.finite(Y_scale) | Y_scale == 0)) {
  stop("Train 中至少有一个 Y variable 的标准差为 0 或无效。")
}


# ------------------------------------------------------------
# 5. 用同一套 train 参数标准化 train 和 test
# ------------------------------------------------------------
# train 根据自身均值和标准差进行标准化。
train <- list(
  X = scale(
    X_train_raw,
    center = X_center,
    scale = X_scale
  ),
  Y = scale(
    Y_train_raw,
    center = Y_center,
    scale = Y_scale
  )
)

# test 必须使用 train 的均值和标准差，
# 不能使用 test 自己的均值和标准差。
test <- list(
  X = scale(
    X_test_raw,
    center = X_center,
    scale = X_scale
  ),
  Y = scale(
    Y_test_raw,
    center = Y_center,
    scale = Y_scale
  )
)


# ------------------------------------------------------------
# 6. 显式定义 X/Y 两个数据块的连接结构
# ------------------------------------------------------------
design <- matrix(
  c(
    0, 1,
    1, 0
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("X", "Y"),
    c("X", "Y")
  )
)


# ------------------------------------------------------------
# 7. 初始化结果保存对象
# ------------------------------------------------------------
# xt/yt 保存每个 component 使用的 train/test X/Y 数据矩阵。
xt <- list(
  train = list(),
  test = list()
)

yt <- list(
  train = list(),
  test = list()
)

# 保存每个 component 的 weights。
# test weights 不重新学习，而是复制 train weights。
weight <- list(
  trainX = list(),
  trainY = list(),
  testX = list(),
  testY = list()
)

# 保存每只小鼠的 canonical scores。
score <- list(
  trainX = list(),
  trainY = list(),
  testX = list(),
  testY = list()
)

# 保存 train/test canonical correlations。
canonical_r <- list(
  train = numeric(ncomp),
  test = numeric(ncomp)
)

# 保存 train/test loadings。
loading <- list(
  trainX = list(),
  trainY = list(),
  testX = list(),
  testY = list()
)

# 保存 train/test cross-loadings。
crossloading <- list(
  trainX_to_Yscore = list(),
  trainY_to_Xscore = list(),
  testX_to_Yscore = list(),
  testY_to_Xscore = list()
)

# Component 1 使用标准化后的原始 train/test 矩阵。
xt$train[[1]] <- train$X
yt$train[[1]] <- train$Y

xt$test[[1]] <- test$X
yt$test[[1]] <- test$Y


# ------------------------------------------------------------
# 8. 依次拟合 Component 1、2、3
# ------------------------------------------------------------
for (i in seq_len(ncomp)) {

  # ----------------------------------------------------------
  # 8.1 Component 2 和 3：先做 projection deflation
  # ----------------------------------------------------------
  if (i > 1) {

    # 仅使用 train 当前矩阵和上一 component 的 train weights
    # 对 train 做 deflation。
    ns_train <- nscancor::acor(
      xt$train[[i - 1]],
      matrix(weight$trainX[[i - 1]], ncol = 1),
      yt$train[[i - 1]],
      matrix(weight$trainY[[i - 1]], ncol = 1),
      xscale = TRUE,
      yscale = TRUE
    )

    # 保存 train deflation 后的矩阵。
    xt$train[[i]] <- ns_train$xp
    yt$train[[i]] <- ns_train$yp

    # 恢复行列名。
    dimnames(xt$train[[i]]) <- dimnames(train$X)
    dimnames(yt$train[[i]]) <- dimnames(train$Y)


    # --------------------------------------------------------
    # Paper-like test deflation
    # --------------------------------------------------------
    # 为贴近论文公开代码，在 test 当前矩阵上调用 acor()，
    # 但传入的 weights 来自 train。
    #
    # 注意：
    # 这不是完全冻结的 train-only projection mapping；
    # 正式方法学研究应进一步实现 train-fitted projection
    # 原样应用到 test。
    ns_test <- nscancor::acor(
      xt$test[[i - 1]],
      matrix(weight$testX[[i - 1]], ncol = 1),
      yt$test[[i - 1]],
      matrix(weight$testY[[i - 1]], ncol = 1),
      xscale = TRUE,
      yscale = TRUE
    )

    # 保存 test deflation 后的矩阵。
    xt$test[[i]] <- ns_test$xp
    yt$test[[i]] <- ns_test$yp

    # 恢复行列名。
    dimnames(xt$test[[i]]) <- dimnames(test$X)
    dimnames(yt$test[[i]]) <- dimnames(test$Y)
  }


  # ----------------------------------------------------------
  # 8.2 仅在 train 中拟合当前 sCCA component
  # ----------------------------------------------------------
  result_train <- wrapper.sgcca(

    # 当前 component 的 train X/Y 数据块。
    X = list(
      X = xt$train[[i]],
      Y = yt$train[[i]]
    ),

    # 显式指定 X/Y 相互连接。
    design = design,

    # 当前 6.32.0 教学实现使用 keepX 控制稀疏度。
    penalty = NULL,

    # 每轮只拟合 1 个 component。
    ncomp = 1,

    # 分别规定 X/Y 两端保留的变量数。
    keepX = list(
      X = keepX.X,
      Y = keepX.Y
    ),

    # canonical 模式：
    # 使 X/Y 两端 canonical scores 的相关尽可能强。
    mode = "canonical",

    # 已在函数外用 train 参数完成标准化。
    scale = FALSE,

    tol = .Machine$double.eps,
    max.iter = 1000,
    near.zero.var = FALSE,
    all.outputs = TRUE
  )


  # ----------------------------------------------------------
  # 8.3 提取 train 学到的 X/Y weights
  # ----------------------------------------------------------
  weight$trainX[[i]] <- drop(
    result_train$loadings$X[, 1]
  )

  weight$trainY[[i]] <- drop(
    result_train$loadings$Y[, 1]
  )

  # 补回变量名。
  names(weight$trainX[[i]]) <- colnames(xt$train[[i]])
  names(weight$trainY[[i]]) <- colnames(yt$train[[i]])

  # test 不重新拟合 weights，
  # 而是直接使用 train 中学到的同一组 weights。
  weight$testX[[i]] <- weight$trainX[[i]]
  weight$testY[[i]] <- weight$trainY[[i]]


  # ----------------------------------------------------------
  # 8.4 提取 train 中每只小鼠的 canonical scores
  # ----------------------------------------------------------
  score$trainX[[i]] <- drop(
    result_train$variates$X[, 1]
  )

  score$trainY[[i]] <- drop(
    result_train$variates$Y[, 1]
  )

  names(score$trainX[[i]]) <- rownames(xt$train[[i]])
  names(score$trainY[[i]]) <- rownames(yt$train[[i]])


  # ----------------------------------------------------------
  # 8.5 用 train weights 计算 test canonical scores
  # ----------------------------------------------------------
  # X-side test score = test X matrix × train X weight。
  score$testX[[i]] <- drop(
    xt$test[[i]] %*% weight$trainX[[i]]
  )

  # Y-side test score = test Y matrix × train Y weight。
  score$testY[[i]] <- drop(
    yt$test[[i]] %*% weight$trainY[[i]]
  )

  names(score$testX[[i]]) <- rownames(xt$test[[i]])
  names(score$testY[[i]]) <- rownames(yt$test[[i]])


  # ----------------------------------------------------------
  # 8.6 分别计算 train/test canonical r
  # ----------------------------------------------------------
  canonical_r$train[i] <- cor(
    score$trainX[[i]],
    score$trainY[[i]]
  )

  canonical_r$test[i] <- cor(
    score$testX[[i]],
    score$testY[[i]]
  )


  # ----------------------------------------------------------
  # 8.7 计算 train/test loadings
  # ----------------------------------------------------------
  # loading = 原始变量与本端 canonical score 的相关。
  loading$trainX[[i]] <- drop(
    cor(xt$train[[i]], score$trainX[[i]])
  )

  loading$trainY[[i]] <- drop(
    cor(yt$train[[i]], score$trainY[[i]])
  )

  loading$testX[[i]] <- drop(
    cor(xt$test[[i]], score$testX[[i]])
  )

  loading$testY[[i]] <- drop(
    cor(yt$test[[i]], score$testY[[i]])
  )


  # ----------------------------------------------------------
  # 8.8 计算 train/test cross-loadings
  # ----------------------------------------------------------
  # X variable 与 Y-side score 的相关。
  crossloading$trainX_to_Yscore[[i]] <- drop(
    cor(xt$train[[i]], score$trainY[[i]])
  )

  # Y variable 与 X-side score 的相关。
  crossloading$trainY_to_Xscore[[i]] <- drop(
    cor(yt$train[[i]], score$trainX[[i]])
  )

  # test X variable 与 test Y-side score 的相关。
  crossloading$testX_to_Yscore[[i]] <- drop(
    cor(xt$test[[i]], score$testY[[i]])
  )

  # test Y variable 与 test X-side score 的相关。
  crossloading$testY_to_Xscore[[i]] <- drop(
    cor(yt$test[[i]], score$testX[[i]])
  )
}


# ------------------------------------------------------------
# 9. 将 list 形式结果整理为矩阵
# ------------------------------------------------------------
component_names <- paste0(
  "Component",
  seq_len(ncomp)
)

# 辅助函数：
# 将三个 component 的结果按列合并。
bind_components <- function(object_list) {
  output_matrix <- do.call(
    cbind,
    object_list
  )
  colnames(output_matrix) <- component_names
  output_matrix
}

# 每行一个变量、每列一个 component。
weights_X <- bind_components(weight$trainX)
weights_Y <- bind_components(weight$trainY)

# 每行一只 train/test 小鼠、每列一个 component。
scores_train_X <- bind_components(score$trainX)
scores_train_Y <- bind_components(score$trainY)

scores_test_X <- bind_components(score$testX)
scores_test_Y <- bind_components(score$testY)

# Train/test loadings。
loadings_train_X <- bind_components(loading$trainX)
loadings_train_Y <- bind_components(loading$trainY)

loadings_test_X <- bind_components(loading$testX)
loadings_test_Y <- bind_components(loading$testY)

# Train/test cross-loadings。
crossloadings_train_X <- bind_components(
  crossloading$trainX_to_Yscore
)
crossloadings_train_Y <- bind_components(
  crossloading$trainY_to_Xscore
)

crossloadings_test_X <- bind_components(
  crossloading$testX_to_Yscore
)
crossloadings_test_Y <- bind_components(
  crossloading$testY_to_Xscore
)


# ------------------------------------------------------------
# 10. 整理 train/test canonical r
# ------------------------------------------------------------
canonical_r_table <- data.frame(
  Component = component_names,
  Train_canonical_r = canonical_r$train,
  Test_canonical_r = canonical_r$test
)


# ------------------------------------------------------------
# 11. 提取每个 component 中非零权重的变量
# ------------------------------------------------------------
selected_X <- do.call(
  rbind,
  lapply(seq_len(ncomp), function(i) {

    current_weights <- weights_X[, i]

    data.frame(
      Component = component_names[i],
      Variable = names(current_weights)[current_weights != 0],
      Weight = unname(current_weights[current_weights != 0]),
      row.names = NULL
    )
  })
)

selected_Y <- do.call(
  rbind,
  lapply(seq_len(ncomp), function(i) {

    current_weights <- weights_Y[, i]

    data.frame(
      Component = component_names[i],
      Variable = names(current_weights)[current_weights != 0],
      Weight = unname(current_weights[current_weights != 0]),
      row.names = NULL
    )
  })
)


# ------------------------------------------------------------
# 12. 保存每只小鼠属于 train 还是 test
# ------------------------------------------------------------
split_table <- data.frame(
  Sample = rownames(X_raw),

  # 根据样本行号标记 Train/Test。
  Set = ifelse(
    seq_len(n) %in% train_id,
    "Train",
    "Test"
  ),

  # diet 和 genotype 仅作为元数据保存，
  # 不进入本次 sCCA 的 X/Y。
  Diet = nutrimouse$diet,
  Genotype = nutrimouse$genotype,

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 13. 保存全部输出结果
# ------------------------------------------------------------

# 13.1 Train 学到的 X/Y weights。
# Test 使用的是同一组 weights，因此无需重复保存。
write.csv(
  weights_X,
  file.path(output_dir, "weights_X_gene.csv")
)

write.csv(
  weights_Y,
  file.path(output_dir, "weights_Y_lipid.csv")
)

# 13.2 每个 component 中非零权重变量。
write.csv(
  selected_X,
  file.path(output_dir, "selected_variables_X_gene.csv"),
  row.names = FALSE
)

write.csv(
  selected_Y,
  file.path(output_dir, "selected_variables_Y_lipid.csv"),
  row.names = FALSE
)

# 13.3 Train/test 每只小鼠的 canonical scores。
write.csv(
  scores_train_X,
  file.path(output_dir, "canonical_scores_train_X_gene.csv")
)

write.csv(
  scores_train_Y,
  file.path(output_dir, "canonical_scores_train_Y_lipid.csv")
)

write.csv(
  scores_test_X,
  file.path(output_dir, "canonical_scores_test_X_gene.csv")
)

write.csv(
  scores_test_Y,
  file.path(output_dir, "canonical_scores_test_Y_lipid.csv")
)

# 13.4 Train/test canonical correlations。
write.csv(
  canonical_r_table,
  file.path(output_dir, "canonical_correlations_train_test.csv"),
  row.names = FALSE
)

# 13.5 Train/test loadings。
write.csv(
  loadings_train_X,
  file.path(output_dir, "loadings_train_X_gene.csv")
)

write.csv(
  loadings_train_Y,
  file.path(output_dir, "loadings_train_Y_lipid.csv")
)

write.csv(
  loadings_test_X,
  file.path(output_dir, "loadings_test_X_gene.csv")
)

write.csv(
  loadings_test_Y,
  file.path(output_dir, "loadings_test_Y_lipid.csv")
)

# 13.6 Train/test cross-loadings。
write.csv(
  crossloadings_train_X,
  file.path(output_dir, "crossloadings_train_X_to_Yscore.csv")
)

write.csv(
  crossloadings_train_Y,
  file.path(output_dir, "crossloadings_train_Y_to_Xscore.csv")
)

write.csv(
  crossloadings_test_X,
  file.path(output_dir, "crossloadings_test_X_to_Yscore.csv")
)

write.csv(
  crossloadings_test_Y,
  file.path(output_dir, "crossloadings_test_Y_to_Xscore.csv")
)

# 13.7 样本拆分和元数据。
write.csv(
  split_table,
  file.path(output_dir, "sample_split_and_metadata.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. 保存 train 标准化参数
# ------------------------------------------------------------
# 保存为 RDS，便于将来对新数据使用完全相同的预处理。
saveRDS(
  list(
    X_center = X_center,
    X_scale = X_scale,
    Y_center = Y_center,
    Y_scale = Y_scale
  ),
  file.path(
    output_dir,
    "training_scaling_parameters.rds"
  )
)


# ------------------------------------------------------------
# 15. 保存软件版本和会话信息
# ------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)


# ------------------------------------------------------------
# 16. 在 Console 中显示运行完成信息
# ------------------------------------------------------------
cat(
  "\n分析完成。\n输出文件夹：",
  normalizePath(output_dir),
  "\n\n"
)

# 打印三个 component 的 train/test canonical r。
print(canonical_r_table)
