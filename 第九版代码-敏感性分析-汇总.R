## ############################################################################
## 第九版代码-敏感性分析-汇总.R
## 
## 核心任务：
##   汇总 4 次分析（基准 + 3 次敏感性分析）的结果
##   生成跨真实模型的对比表
## 
## 设计说明：
##   - 此脚本由主控脚本调用，在所有敏感性分析完成后运行
##   - 也可以单独运行（会读取已保存的结果）
## 
## 输出：
##   - outputs9_SA/sensitivity_summary.rds
##   - outputs9_SA/sensitivity_summary.csv
##   - outputs9_SA/sensitivity_report.txt
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析汇总")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 0. 路径设置
## ############################################################################

PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")
OUT_ROOT_SA <- file.path(PROJ_ROOT, "outputs9_SA")

## ############################################################################
## 1. 读取基准分析结果（zang2021）
## ############################################################################

message("\n--- 1. 读取基准分析结果 ---")

baseline_results <- readRDS(file.path(OUT_DIR_9, "evaluation_results.rds"))

baseline_summary <- tibble(
  true_model = "zang2021",
  analysis_type = "基准",
  
  ## 完全一致率
  agreement_ind = baseline_results$cross_model_consistency$agreement_ind,
  agreement_pop = baseline_results$cross_model_consistency$agreement_pop,
  
  ## Fleiss' Kappa
  fleiss_kappa_ind = baseline_results$cross_model_consistency$fleiss_kappa_ind,
  fleiss_kappa_pop = baseline_results$cross_model_consistency$fleiss_kappa_pop,
  
  ## 全对率
  all_correct_ind = baseline_results$multilevel_accuracy$ind_pct[1],
  all_correct_pop = baseline_results$multilevel_accuracy$pop_pct[1],
  
  ## 平均准确率
  mean_accuracy_ind = mean(baseline_results$accuracy_by_model$accuracy_ind),
  mean_accuracy_pop = mean(baseline_results$accuracy_by_model$accuracy_pop)
)

message("✅ 基准分析结果已读取")

## ############################################################################
## 2. 读取敏感性分析结果
## ############################################################################

message("\n--- 2. 读取敏感性分析结果 ---")

SENSITIVITY_MODELS <- c("li2018", "yin2016", "sun2021")

sa_summaries <- list()

for (model_name in SENSITIVITY_MODELS) {
  
  result_path <- file.path(OUT_ROOT_SA, model_name, "evaluation_results.rds")
  
  if (file.exists(result_path)) {
    
    sa_result <- readRDS(result_path)
    
    sa_summaries[[model_name]] <- tibble(
      true_model = model_name,
      analysis_type = "敏感性分析",
      
      agreement_ind = sa_result$cross_model_consistency$agreement_ind,
      agreement_pop = sa_result$cross_model_consistency$agreement_pop,
      
      fleiss_kappa_ind = sa_result$cross_model_consistency$fleiss_kappa_ind,
      fleiss_kappa_pop = sa_result$cross_model_consistency$fleiss_kappa_pop,
      
      all_correct_ind = sa_result$multilevel_accuracy$ind_pct[1],
      all_correct_pop = sa_result$multilevel_accuracy$pop_pct[1],
      
      mean_accuracy_ind = mean(sa_result$accuracy_by_model$accuracy_ind),
      mean_accuracy_pop = mean(sa_result$accuracy_by_model$accuracy_pop)
    )
    
    message("  ✓ ", model_name, " 结果已读取")
    
  } else {
    message("  ⚠️ ", model_name, " 结果未找到: ", result_path)
  }
}

## ############################################################################
## 3. 合并汇总表
## ############################################################################

message("\n--- 3. 合并汇总表 ---")

sensitivity_summary <- bind_rows(
  baseline_summary,
  bind_rows(sa_summaries)
) %>%
  mutate(
    ## 计算 Δκ（IND - POP）
    delta_kappa = fleiss_kappa_ind - fleiss_kappa_pop,
    
    ## 计算 Δ 完全一致率
    delta_agreement = agreement_ind - agreement_pop,
    
    ## 计算 Δ 全对率
    delta_all_correct = all_correct_ind - all_correct_pop,
    
    ## 计算 Δ 平均准确率
    delta_mean_accuracy = mean_accuracy_ind - mean_accuracy_pop
  )

## ############################################################################
## 4. 打印汇总表
## ############################################################################

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("                     敏感性分析汇总表\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n【核心指标对比】\n\n")

print_table <- sensitivity_summary %>%
  select(
    `真实模型` = true_model,
    `Fleiss κ (IND)` = fleiss_kappa_ind,
    `Fleiss κ (POP)` = fleiss_kappa_pop,
    `Δκ` = delta_kappa,
    `完全一致率 IND` = agreement_ind,
    `完全一致率 POP` = agreement_pop
  ) %>%
  mutate(
    `Fleiss κ (IND)` = round(`Fleiss κ (IND)`, 3),
    `Fleiss κ (POP)` = round(`Fleiss κ (POP)`, 3),
    `Δκ` = sprintf("%+.3f", `Δκ`),
    `完全一致率 IND` = sprintf("%.1f%%", `完全一致率 IND`),
    `完全一致率 POP` = sprintf("%.1f%%", `完全一致率 POP`)
  )

print(print_table, n = Inf)

## ############################################################################
## 5. 结论判定
## ############################################################################

cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("结论判定\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

## 检查所有 Δκ 是否 > 0
all_delta_kappa_positive <- all(sensitivity_summary$delta_kappa > 0)

## 检查所有 Δ 完全一致率是否 > 0
all_delta_agreement_positive <- all(sensitivity_summary$delta_agreement > 0)

## 计算 κ 的变异系数
kappa_cv_ind <- sd(sensitivity_summary$fleiss_kappa_ind) / mean(sensitivity_summary$fleiss_kappa_ind) * 100
kappa_cv_pop <- sd(sensitivity_summary$fleiss_kappa_pop) / mean(sensitivity_summary$fleiss_kappa_pop) * 100

cat("1. IND 模式一致性优于 POP 模式：\n")
if (all_delta_kappa_positive) {
  cat("   ✅ 在所有真实模型下，Δκ > 0\n")
} else {
  cat("   ❌ 存在 Δκ ≤ 0 的情况\n")
}

cat("\n2. 结论稳健性：\n")
cat(sprintf("   κ (IND) 变异系数 = %.1f%%\n", kappa_cv_ind))
cat(sprintf("   κ (POP) 变异系数 = %.1f%%\n", kappa_cv_pop))

if (kappa_cv_ind < 30) {
  cat("   ✅ IND 模式结果在不同真实模型间稳定\n")
} else {
  cat("   ⚠️ IND 模式结果存在一定变异\n")
}

cat("\n3. 总体结论：\n")
if (all_delta_kappa_positive && all_delta_agreement_positive) {
  cat("   ✅ 敏感性分析支持原结论：\n")
  cat("      MAPB 个体化参数使不同模型对同一患者的判别结果趋于一致，\n")
  cat("      此结论不依赖于特定真实模型的选择。\n")
} else {
  cat("   ⚠️ 敏感性分析结果需要进一步分析\n")
}

## ############################################################################
## 6. 保存输出
## ############################################################################

message("\n--- 6. 保存输出 ---")

## 保存汇总表
saveRDS(sensitivity_summary, file.path(OUT_ROOT_SA, "sensitivity_summary.rds"))
write_csv(sensitivity_summary, file.path(OUT_ROOT_SA, "sensitivity_summary.csv"))
message("  ✓ sensitivity_summary.rds / .csv 已保存")

## 生成文本报告
report_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "第九版敏感性分析报告",
  paste(rep("=", 80), collapse = ""),
  "",
  paste("生成时间:", Sys.time()),
  "",
  "【分析概述】",
  "  - 基准分析：zang2021 作为真实模型（一室模型）",
  "  - 敏感性分析 1：li2018 作为真实模型（两室模型）",
  "  - 敏感性分析 2：yin2016 作为真实模型（两室模型）",
  "  - 敏感性分析 3：sun2021 作为真实模型（两室模型）",
  "",
  "【核心指标汇总】",
  ""
)

for (i in 1:nrow(sensitivity_summary)) {
  row <- sensitivity_summary[i, ]
  report_lines <- c(report_lines,
    sprintf("  %s:", row$true_model),
    sprintf("    Fleiss' κ: IND=%.3f, POP=%.3f, Δ=%+.3f",
            row$fleiss_kappa_ind, row$fleiss_kappa_pop, row$delta_kappa),
    sprintf("    完全一致率: IND=%.1f%%, POP=%.1f%%",
            row$agreement_ind, row$agreement_pop),
    sprintf("    全对率: IND=%.1f%%, POP=%.1f%%",
            row$all_correct_ind, row$all_correct_pop),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  "【结论】",
  sprintf("  所有 Δκ > 0: %s", ifelse(all_delta_kappa_positive, "是", "否")),
  sprintf("  κ (IND) 变异系数: %.1f%%", kappa_cv_ind),
  sprintf("  结论稳健性: %s", 
          ifelse(all_delta_kappa_positive && kappa_cv_ind < 30, "通过", "需进一步分析")),
  "",
  paste(rep("=", 80), collapse = "")
)

writeLines(report_lines, file.path(OUT_ROOT_SA, "sensitivity_report.txt"))
message("  ✓ sensitivity_report.txt 已保存")

## ############################################################################
## 7. 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 敏感性分析汇总完成")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_ROOT_SA)
message("\n📄 输出文件：")
message("  - sensitivity_summary.rds / .csv  （汇总表）")
message("  - sensitivity_report.txt          （文本报告）")

message("\n📊 核心结论：")
message(sprintf("   所有 Δκ > 0: %s", ifelse(all_delta_kappa_positive, "✅ 是", "❌ 否")))
message(sprintf("   结论稳健性: %s", 
                ifelse(all_delta_kappa_positive && kappa_cv_ind < 30, "✅ 通过", "⚠️ 需进一步分析")))