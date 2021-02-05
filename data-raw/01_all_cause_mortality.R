### code to prepare `all_cause_mortality` dataset goes here

test_all <- readHMDweb(CNTRY = "USA", 
                       item = "fltper_1x1", 
                       username = "feralaes@gmail.com", 
                       password = "sip667", 
                       fixup = TRUE)
test_2015 <- test_all %>% 
  filter(Year == 2015)
test_2018 <- test_all %>% 
  filter(Year == 2018)


file.init <- "data-raw/01_all_cause_mortality.csv"

all_cause_mortality  <- read.csv(file = file.init, stringsAsFactors = F)

# Create .rda object for initial set of parameters and store it in 'data' folder
usethis::use_data(all_cause_mortality, overwrite = TRUE)