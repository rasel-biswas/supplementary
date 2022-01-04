library(mitml)
library(mice)
# Care-Survival Data (taken from Little and Rubin (1987, Table 9.8))

counts_full_part <- data.frame(
  survived = gl(n = 2, k = 4, labels = 0:1, length = 8),
  cared = gl(n=2, k=2, labels = 0:1, length = 8),
  clinic_B = gl(n=2, k=1, labels = 0:1, length = 8),
  Freq = c(3, 17, 4, 2, 176, 197, 293, 23)
)
xtabs(Freq ~ cared+survived+clinic_B, data=counts_full_part) # contingency table

counts_partial_part <- data.frame(
  survived = gl(n=2, k=2, labels = 0:1, length = 4),
  cared = gl(n=2, k=1, labels = 0:1, length = 4),
  clinic_B = NA,
  Freq = c(10, 5, 150, 90)
)
xtabs(Freq ~ cared+survived, data=counts_partial_part)


countsToCases <- function(x, count_column = "Freq") {
  idx <- rep.int(x = seq_len(nrow(x)),
                 times = x[[count_column]])# Get the row_indices to pull from x
  x[[count_column]] <- NULL                # Drop count column
  x[idx,]                                  # Get the rows from x
}


data1 <- countsToCases(x = counts_full_part, count_column ="Freq")
data2 <- countsToCases(x = counts_partial_part, count_column ="Freq")
data3 <- rbind(data1, data2)

#m=50
#d <- complete(imputed_datasets, action = 1)
#as.data.frame(table(data.frame(d$survived, d$cared, d$clinic_B)))


sim_func <- function(m=m){

  imputed_datasets <<- mice::mice(data = data3, m = m, maxit = 100, print = FALSE, seed = 2021)

  fit_independence <- with(
    imputed_datasets,
    glm(
      Freq ~ survived + cared + clinic_B,
      family = poisson(link = "log"),
      data = as.data.frame(table(data.frame(survived, cared, clinic_B)))
    )
  )

  fit_conditional_ind <- with(
    imputed_datasets,
    glm(
      Freq ~ survived + cared + clinic_B +
        survived * clinic_B + clinic_B * cared,
      family = poisson(link = "log"),
      data = as.data.frame(table(data.frame(survived, cared, clinic_B)))
    )
  )

  fit_saturated <- with(
    imputed_datasets,
    glm(
      Freq ~ survived * cared * clinic_B,
      family = poisson(link = "log"),
      data = as.data.frame(table(data.frame(survived, cared, clinic_B)))
    )
  )

  out = c(
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses)$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses)$test[, c(1,4,5)],
    
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses, df.com = 970-4-1)$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses, df.com = 970-2-1)$test[, c(1,4,5)],
    
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses, 'D2')$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses, 'D2')$test[, c(1,4,5)],
    
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses, method = 'D2', use = 'likelihood')$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses, method = 'D2', use = 'likelihood')$test[, c(1,4,5)],
    
    D3(fit_saturated, fit_conditional_ind)$result[c(1, 4, 5)],
    D3(fit_saturated, fit_independence)$result[c(1, 4, 5)],
    
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses, 'D4')$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses, 'D4')$test[, c(1,4,5)],
    
    testModels(fit_saturated$analyses, fit_conditional_ind$analyses, 'D4', ariv = "robust")$test[, c(1,4,5)],
    testModels(fit_saturated$analyses, fit_independence$analyses, 'D4', ariv = "robust")$test[, c(1,4,5)]
  )
  out
}

result_matrix <- matrix(
  data = c(sim_func(m = 200)),
  byrow = TRUE,
  ncol = 6
)

rownames(result_matrix) <- c("D1", "D1-Reiter", "D2-Wald", "D2-Lik", "D3", "D4", "D4-robust")
colnames(result_matrix) <- c("D", "p-value", "r", "D", "p-value", "r")

knitr::kable(result_matrix, digits = 3)
