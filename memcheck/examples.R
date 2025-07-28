# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
topics <- tools::Rd_db("TreeTools")

cat("Running examples for", length(topics), "topics\n")

for (topic in names(topics)) {
  cat("\n>>> Example:", topic, "\n")
  tryCatch({
    tools::runExample("TreeTools", topic = topic)
  }, error = function(e) {
    cat("Error in topic:", topic, "\n", conditionMessage(e), "\n")
  })
}
