# uses A standardized diagnostic approach and ongoing feedback improves outcome in psychogenic nonepileptic seizures by Drane et al., 2016 Table 2

library(pwr)

# 8 week event frequency using standard practice
standard_practice_mean <- 2.5 
standard_practice_sd <- 1.0
standard_practice_n <- 12

# 8 week event frequency using structured feedback
structured_feedback_mean <- 1.7
structured_feedback_sd <- 0.5
structured_feedback_n <- 10

# calculate pooled standard deviation
pooled_sd <- sqrt(((standard_practice_n-1)*standard_practice_sd^2+(structured_feedback_n-1)*structured_feedback_sd^2)/(standard_practice_n+structured_feedback_n-2))

# Calculate effect size
cohen_d <- (standard_practice_mean-structured_feedback_mean)/pooled_sd

# do power calculation
pwr.t.test(n = NULL, d = cohen_d, sig.level = 0.05, power = 0.8, type = "two.sample")

# 17 October 2019
#24 September 2019
#19 August 2019