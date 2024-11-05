#' Zotero example
#'
#' A dataset of scientific journal abstracts exported directly from Zotero. The
#' necessary recodes are included here so that it is easier to replicate on your
#' own Zotero library.
#'
#' @name abstracts
#'
#' @docType data
#'
#' @usage data(abstracts)
#'
#' @return A dataframe named 'abstracts' with the minimum required variables for making the 'topicspace' package work correctly, namely: ID, documents, title, start.year. It also returns 'participants' dataframe with additional variables categorizing authors, which can be projected into the space of topics.
#'
#' @examples
#'
#' library(dplyr)
#' library(purrr)
#'
#' data(abstracts)
#' abstracts$nchar <- nchar(abstracts$`Abstract Note`)
#' abstracts       <- abstracts %>% filter(!is.na(nchar)) %>% filter(nchar > 500)
#' nmiss <- abstracts %>% map_dbl(~sum(is.na(.x))/length(.x))
#' abstracts     <- abstracts[, nmiss < 0.90]
#'
#' # Required columns -----
#' abstracts$ID           <- abstracts$Key
#' abstracts$documents     <- abstracts$`Abstract Note`
#' abstracts$title        <- abstracts$Title
#' abstracts$start.year   <- abstracts$`Date Added` %>% lubridate::year()
#'
#' # Distinct ----
#' abstracts    <- abstracts %>% distinct(documents, .keep_all = TRUE)
#' abstracts    <- abstracts %>% distinct(DOI, .keep_all = TRUE)
#' abstracts    <- abstracts %>% distinct(title, .keep_all = TRUE)
#'
#' # Authors ----
#' authors <- abstracts %>% select(ID, Author = Author)
#' authors <- authors %>%
#' tidyr::separate_rows(Author, sep = ";") %>%
#' mutate(Author = trimws(Author))
#'
"abstracts"
