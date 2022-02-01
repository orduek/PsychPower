#### Import
library(readr)
library(dplyr)
library(tidyr)

#### Import
data <- read_delim("data-raw/data.csv", escape_double = FALSE, trim_ws = TRUE)

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

#' DASS - Anxiety
#'
#' A dataset containing A39,775 participant ratings of the items of the
#' anxiety subscale of the Depression Anxiety Stress Scales (DASS).
#'
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{Q2}{I was aware of dryness of my mouth.}
#'   \item{Q7}{I had a feeling of shakiness (eg, legs going to give way).}
#'   \item{Q9}{I found myself in situations that made me so anxious I was most relieved when they ended.}
#'   \item{Q15}{I had a feeling of faintness.}
#'   \item{Q19}{I perspired noticeably (eg, hands sweaty) in the absence of high temperatures or physical exertion.}
#'   \item{Q20}{I felt scared without any good reason.}
#'   \item{Q23}{I had difficulty in swallowing.}
#'   \item{Q25}{I was aware of the action of my heart in the absence of physical exertion (eg, sense of heart rate increase, heart missing a beat).}
#'   \item{Q28}{I felt I was close to panic.}
#'   \item{Q30}{I feared that I would be &quot;thrown&quot; by some trivial but unfamiliar task.}
#'   \item{Q36}{I felt terrified.}
#'   \item{Q40}{I was worried about situations in which I might panic and make a fool of myself.}
#'   \item{Q41}{I experienced trembling (eg, in the hands).}
#'
#' }
#' @source \url{https://openpsychometrics.org/_rawdata/}
"data_test"

usethis::use_data(data_test, overwrite = TRUE)
