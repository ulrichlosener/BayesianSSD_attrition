# Function to translate informative hypotheses in string format into a list with two contrast matrices

hyp2mat <- function(hypothesis) {
  
  # Split the hypothesis into parts
  parts <- unlist(strsplit(hypothesis, "(?=[<>=])|(?<=[<>=])", perl = TRUE))
  
  # Extract parameters
  parameters <- unique(parts[!parts %in% c("<", ">", "=")])
  num_params <- length(parameters)
  
  # storage for restriction matrices
  ineq_constraints <- list()
  eq_constraints <- list()
  
  # Parse the hypothesis
  for (i in seq(1, length(parts) - 1, by = 2)) {
    param1 <- parts[i]
    operator <- parts[i + 1]
    param2 <- parts[i + 2]
    
    index1 <- match(param1, parameters)
    index2 <- match(param2, parameters)
    
    if (operator == "<" || operator == ">") {
      constraint <- rep(0, num_params)
      if (operator == "<") {
        constraint[c(index1, index2)] <- c(-1, 1)
      } else {
        constraint[c(index1, index2)] <- c(1, -1)
      }
      ineq_constraints[[length(ineq_constraints) + 1]] <- constraint
    } else if (operator == "=") {
      constraint <- rep(0, num_params)
      constraint[c(index1, index2)] <- c(1, -1)
      eq_constraints[[length(eq_constraints) + 1]] <- constraint
    }
  }
  
  # Convert lists to matrices
  if (length(ineq_constraints) > 0) {
    Ineq.mat <- do.call(rbind, ineq_constraints)
  } else {
    Ineq.mat <- matrix(0, 0, num_params)
  }
  
  if (length(eq_constraints) > 0) {
    Eq.mat <- do.call(rbind, eq_constraints)
  } else {
    Eq.mat <- matrix(0, 0, num_params)
  }
  
  list(Ineq.mat = Ineq.mat, Eq.mat = Eq.mat)
}

# Example usage
hypothesis1 <- "a>b=c<d"
hypothesis2 <- "a=b=c"
hypothesis3 <- "a1>a2<a3"
hypothesis4 <- "a>b=c"

result1 <- hyp2mat(hypothesis1)
result2 <- hyp2mat(hypothesis2)
result3 <- hyp2mat(hypothesis3)
result4 <- hyp2mat(hypothesis4)

