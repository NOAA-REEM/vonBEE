# Load raw data from .csv file
# exampleData <- read.csv("data-raw/simulated-data.csv")
load("data-raw/LWAdat.rda")
LWA <- LWAdat
# Apply preprocessing...
# Save the cleaned data in the required R package location
#usethis::use_data(exampleData)
usethis::use_data(LWA)

load("data-raw/diet_p4CEATTLE.Rdata")

usethis::use_data(diet_p)
