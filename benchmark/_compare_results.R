pr_files <- list.files("pr-benchmark-results", pattern = "*.bench.Rds",
                       full.names = TRUE)

output <- "report<<EOF\n ### Performance benchmark results\n"
regressions <- FALSE

for (pr_file in pr_files) {
  file_name <- basename(pr_file)
  main_file <- file.path("main-benchmark-results", file_name)
  if (!file.exists(main_file)) next;
  
  # Load the results
  pr_results <- readRDS(pr_file)
  main_results <- readRDS(main_file)
  
  
  # A simple comparison function using t-test
  # You can make this much more sophisticated.
  compare_timings <- function(pr_results, main_results) {
    # Get the data frames of timings
    pr_df <- as.data.frame(pr_results)
    main_df <- as.data.frame(main_results)
    
    # Prepare a report
    report <- list()
    
    # Iterate over each function benchmarked
    for (fn_name in unique(pr_df$expr)) {
      pr_times <- pr_df[pr_df$expr == fn_name, "time"]
      main_times <- main_df[main_df$expr == fn_name, "time"]
      
      better_result <- t.test(pr_times, main_times, alternative = "less")
      worse_result <- t.test(pr_times, main_times, alternative = "greater")
      
      # The p-value tells us if the PR's performance is significantly slower
      # A small p-value (e.g., < 0.05) suggests it is.
      is_faster <- better_result$p.value < 0.01
      is_slower <- worse_result$p.value < 0.01
      median_pr <- median(pr_times)
      median_main <- median(main_times)
      percentage_change <- ((median_pr - median_main) / median_main) * 100
      
      report[[fn_name]] <- list(
        slower = is_slower,
        faster = is_faster,
        p_value = worse_result$p.value,
        median_pr = median_pr,
        median_main = median_main,
        change = percentage_change
      )
    }
    
    return(report)
  }
  
  # Generate the report
  report <- compare_timings(pr_results, main_results)
  
  # Create a markdown-formatted message
  has_significant_regression <- FALSE
  
  for (fn_name in names(report)) {
    res <- report[[fn_name]]
    status <- ifelse(res$slower, "\U1F7E0 Slower", 
                     ifelse(res$faster, "\U1F7E2 Faster!",
                            ifelse(res$p_value < 0.05,
                                   "\U1F7E1 A little slower (0.01 < p < 0.05)",
                                   "\U26AA No significant change")
                     )
    )
    if (res$slower) {
      has_significant_regression <- TRUE
    }
    
    message <- paste0(
      "#### `", fn_name, "`: ", status, "\n",
      "Change: **", round(res$change, 2), "%** (p = ", 
      format.pval(res$p_value), "): ",
      round(res$median_main / 1e6, 2), " \U2192 ",
      round(res$median_pr / 1e6, 2), " ms\n\n"
    )
  }
  
  if (has_significant_regression) {
    regressions <- TRUE
  }
  
  cat(message)
  output <- paste0(output, message)
}

cat(paste0(output, "\nEOF"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)

message(readLines(Sys.getenv("GITHUB_OUTPUT")))

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
}
