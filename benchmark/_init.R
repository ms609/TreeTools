library("TreeTools")

Benchmark <- function(..., times = 1L, min_iterations = NULL, min_time = NULL) {
  # Capture expressions unevaluated, preserving non-standard evaluation.
  env <- parent.frame()
  exprs <- as.list(substitute(list(...)))[-1]
  labels <- vapply(exprs, deparse1, "")
  if (times > 1L) {
    # Each timed call should take several ms, so that CI's 0.2 ms floor
    # (benchmark/_compare_results.R) still exposes a change of a few percent.
    # Repeat fast calls `times` times in a function that is byte-compiled
    # once: a bare `for` loop is JIT-compiled at every evaluation, adding ~2 ms.
    labels <- paste(times, "*", labels)
    exprs <- lapply(exprs, function(e) {
      f <- compiler::cmpfun(eval(call("function", NULL, bquote(
        for (.i in seq_len(.(as.integer(times)))) .(e))), env))
      as.call(list(f))
    })
  }
  names(exprs) <- labels
  result <- if (is.null(min_time)) {
    bench::mark(exprs = exprs, env = env,
                min_iterations = min_iterations %||% 3, time_unit = "us")
  } else {
    bench::mark(exprs = exprs, env = env, min_iterations = min_iterations %||% 3,
                min_time = min_time, time_unit = "us")
  }
  if (interactive()) {
    print(result)
  } else {
    fileroot <- gsub("[\"']", "",
                     gsub("[\\(\\):, /*]", "_", labels))
    .FileName <- function(fileRoot, i) {
      paste0(c(fileroot, i, "bench.Rds"), collapse = ".")
    }
    i <- double(0)
    while(file.exists(.FileName(fileroot, i))) {
      if (length(i) == 0) i <- 0
      i <- 1 + i
    }
    saveRDS(result, .FileName(fileroot, i))
  }
}
