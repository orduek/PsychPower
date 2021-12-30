library(foreign)
dfPCL5 <- read.spss('/home/or/Documents/pcl5_vaData/phq_pcl_19/PCLpct19.sav', to.data.frame = TRUE)
#### Prepare dataset
data_subset <- dfPCL5[,5:24]
# omitting data just based on PCL NAs
dfPCL5 <- dfPCL5[complete.cases(data_subset), ]

dB <- binarize_data(dfPCL5[,5:24], 2)
db_counted <- get_freq(dB)

# building the package
library(devtools)
library(roxygen2)
load_all()

roxygenise()

# installing from github
library(devtools) # Make sure that the devtools library is loaded
install_github("orduek/PsychPower")
library(PsychPower)


