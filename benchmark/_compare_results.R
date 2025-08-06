pr_files <- list.files("pr-benchmark-results", pattern = "*.bench.Rds",
                       full.names = TRUE)

output <- paste0(
  "report<<EOF\n### Performance benchmark results\n\n",
  "| Call     | Status | Change | Time (ms) |\n",
  "|----------|--------|--------|-----------|\n"
  )

regressions <- FALSE

for (pr_file in pr_files) {
  file_name <- basename(pr_file)
  replicate_file <- file.path("pr2-benchmark-results", file_name)
  main_file <- file.path("main-benchmark-results", file_name)
  if (!file.exists(main_file)) next;
  
  # Load the results
  rep_exists <- file.exists(replicate_file)
  pr1 <- readRDS(pr_file)
  pr2 <- if (rep_exists) readRDS(replicate_file) else pr1
  main <- readRDS(main_file)
  
  # Prepare a report
  report <- list()
  
  # Iterate over each function benchmarked
  for (fn_name in unique(as.character(unlist(pr1[["expression"]])))) {
    pr1_times <-  as.numeric(pr1[["time"]][[1]] * 1e9)
    pr2_times <-  as.numeric(pr2[["time"]][[1]] * 1e9)
    pr_times <- if (rep_exists) c(pr1_times, pr2_times) else pr1_times
    main_times <- as.numeric(main[["time"]][[1]] * 1e9)
    matched <- if (length(main_times)) {
      TRUE
    } else {
      main_times <- main[, "time"]
      FALSE
    }
    
    median_pr <- median(pr_times)
    median_main <- median(main_times)
    percentage_change <- ((median_main - median_pr) / median_main) * 100
    
    q <- 0.25
    main_iqr <- quantile(main_times, c(q, 1 - q))
    
    is_faster <- median_pr < main_iqr[[1]] && matched
    is_slower <- median_pr > main_iqr[[2]] && matched
    
    report[[fn_name]] <- list(
      matched = matched,
      slower = is_slower,
      faster = is_faster,
      p_value = worse_result$p.value,
      median_pr = median_pr1,
      median_cf = median_pr2,
      median_main = median_main,
      change = percentage_change
    )
  }
  
  # Create a markdown-formatted message
  has_significant_regression <- FALSE
  
  for (fn_name in names(report)) {
    res <- report[[fn_name]]
    status <- ifelse(res$matched,
                     ifelse(res$slower, "\U1F7E0 Slower \U1F641",
                            ifelse(res$faster, "\U1F7E2 Faster!",
                                   ifelse(res$p_value < 0.05,
                                          "\U1F7E1 ?Slower",
                                          "\U26AA NSD")
                            )
                     ),
                     "\U1F7E4 ?Mismatch")
    
    if (res$slower) {
      has_significant_regression <- TRUE
    }
    
    bold <- ifelse(res$faster | res$slower, "**", "")
    
    message <- paste0(
      "| `", fn_name, "` | ", status, " | ", 
      bold, round(res$change, 2), "%", bold, " | ", 
      signif(res$median_main * 1e-6, 3), " \u2192<br />",
      signif(res$median_pr   * 1e-6, 3), ",  ",
      signif(res$median_cf   * 1e-6, 3), " |\n"
    )
  }
  
  if (has_significant_regression) {
    regressions <- TRUE
  }
  
  cat(message)
  output <- paste0(output, message)
}

cat(paste0(output, "\nEOF"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
}

