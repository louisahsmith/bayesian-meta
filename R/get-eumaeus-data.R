library(DatabaseConnector)
library(tidyverse)

method_name <- 'SCCS' # CaseControl, CohortMethod, HistoricalComparator, SCCS
analysis_ids <- 1 # 1:30 (depending on method)
exposure_ids <- 21215 # c(21184, 21185, 21214, 21215, 211981, 211982, 211983, 211831, 211832, 211833)
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = keyring::key_get("eumaeusDriverPath"))

eumaeusConnectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"
  ),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword")
)

connection <- connect(eumaeusConnectionDetails)

# all_tbls_sql <- "
# SELECT *
# FROM pg_catalog.pg_tables
# WHERE schemaname != 'pg_catalog' AND 
#     schemaname != 'information_schema';
# "
# 
# eumaeus_tbls <- querySql(connection, all_tbls_sql) %>% 
#   janitor::clean_names() %>% 
#   filter(schemaname == "eumaeus")

get_tbl <- function(tbl_name) {
  sql <- str_glue("
  SELECT *
  FROM eumaeus.{tbl_name};
  ")
  
  tbl <- querySql(connection, sql, snakeCaseToCamelCase = TRUE)
  
  tbl
}

database <- get_tbl("database")
analysis <- get_tbl("analysis")
exposure <- get_tbl("exposure")
negative_control_outcome <- get_tbl("negative_control_outcome")
positive_control_outcome <- get_tbl("positive_control_outcome") %>% 
  select(-exposureId)
time_period <- get_tbl("time_period")

for (analysis_id in analysis_ids) {
  for (exposure_id in exposure_ids) {
    select_est <- str_glue("
                           SELECT *
                           FROM eumaeus.estimate
                           NATURAL JOIN eumaeus.likelihood_profile
                           WHERE estimate.method = '{method_name}'
                             AND estimate.analysis_id = {analysis_id}
                             AND estimate.exposure_id = {exposure_id};
                           ")
    
    est <- querySql(connection, select_est, snakeCaseToCamelCase = TRUE)
    
    dat <- est %>% 
      left_join(database, by = "databaseId") %>% 
      left_join(analysis, by = c("analysisId", "method"), suffix = c("Database", "Analysis")) %>% 
      left_join(exposure, by = c("exposureId")) %>% 
      left_join(time_period, by = c("exposureId", "periodId"), suffix = c("Exposure", "Period")) %>% 
      left_join(negative_control_outcome, by = "outcomeId", suffix = c("", "Negative")) %>% 
      left_join(positive_control_outcome, by = "outcomeId", suffix = c("", "Positive")) %>% 
      mutate(outcomeName = ifelse(is.na(outcomeName), outcomeNamePositive, outcomeName)) %>% 
      select(-outcomeNamePositive)
    
    get_ll <- function(point, value, ...) {
      tibble(point = as.numeric(strsplit(point, ";")[[1]]),
             value = as.numeric(strsplit(value, ";")[[1]]))
    }
    
    with_ll <- dat %>% 
      mutate(ll = pmap(list(point, value), get_ll)) %>% 
      select(-point, -value)
    
    write_rds(with_ll, here::here("data", str_glue("eumaeus_{method_name}_{analysis_id}_{exposure_id}.rds")))
    
  }
}

disconnect(connection)
