#### Import
library(readr)
library(dplyr)
library(tidyr)

#### Import
data <- read_delim("inst/extdata/data.csv", escape_double = FALSE, trim_ws = TRUE)

###### 2.1 Select variables for methods and materials // DATA0
data0 <- data %>%
  select(ends_with("A"), gender, age) %>%
  filter(gender > 0) %>%  #67 with 0,  1=Male, 2=Female, 3=Other
  filter(age < 95) %>%
  drop_na()

data0$gender <- as.factor(data0$gender)


###### 2.2 Select & prepare dataset for further analysis // DATA1
data_0_rename <- data0 %>%
  select(Q1A:Q42A)

colnames(data_0_rename) <- c(paste0("Q", 1:(ncol(data_0_rename))), "label")

### Select items from Anxiety subscale
data_test <- data_0_rename %>%  # Anxiety subscale
  select(Q2,Q7,Q9,Q15,Q19, Q20, Q23, Q25, Q28, Q30, Q36, Q40, Q41)

usethis::use_data(data_test, overwrite = TRUE)
