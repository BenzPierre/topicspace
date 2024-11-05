# Section 1: Text Processing Functions
# -------------------------------------

#' Make tidyterms from raw text
#'
#' This function prepares raw text data for topic modeling by tokenizing, cleaning, and transforming it into a tidy format. It also calculates term frequency (tf), inverse document frequency (idf), and term frequency-inverse document frequency (tf-idf).
#'
#' @param abstracts A dataframe including the following minimum required columns: ID, documents, title, start.year.
#' @param long.label Logical indicating whether to include long labels for word stems. Default is FALSE.
#' @param add.titles Logical indicating whether to add titles to the abstract text. Default is TRUE.
#' @return A tibble with the tidy terms ready for topic modeling.
#' @export
#'
#' @examples
#' # Load example abstracts data
#' data(abstracts)
#' # Set minimum requirements for abstract data
#' example(abstracts)
#' # Prepare tidy terms for topic modeling
#' topic.model.ready.tidy.terms(abstracts)
#' @seealso
#' \code{\link{documents.to.topics}}
topic.model.ready.tidy.terms <- function(abstracts, long.label = FALSE, add.titles = TRUE){
  obj    <- abstracts$documents
  if(identical(add.titles, TRUE)) obj <- paste(obj, abstracts$title)

  x       <- tibble(text = obj, ID = abstracts$ID)
  x.tidy  <- x %>% tidytext::unnest_tokens(word, text)
  x.tidy  <- x.tidy %>% anti_join(tidytext::stop_words)
  x.tidy  <- x.tidy[is.na(as.numeric(gsub("[^\\d]+", "", x.tidy$word, perl=TRUE))), ] # Here we remove all numbers!
  x.tidy  <- x.tidy[nchar(x.tidy$word, type = "bytes")>2, ]

  x.tidy$stem    <- textstem::lemmatize_words(x.tidy$word)
  x.tidy$wordset <- x.tidy$stem
  x.tidy$wordset[is.na(x.tidy$wordset)] <- x.tidy$word[is.na(x.tidy$wordset)]

  # Alias in i stem
  if (identical(long.label, TRUE)){
    lab.tab           <- x.tidy %>% group_by(wordset) %>% distinct() %>% summarise(strings = word %>% unique %>% sort %>% paste(collapse = ", "))
    x.tidy            <- left_join(x.tidy, lab.tab)
    x.tidy$stem       <- paste(x.tidy$wordset, x.tidy$strings, sep = ": ")
  }
  x.freq            <- x.tidy %>% count(ID, stem, sort = TRUE) %>% tidytext::bind_tf_idf(ID, stem, n)
  x.freq
}


#' Prune incidence matrix
#'
#' This function prunes an incidence matrix by removing rows and columns with too few nonzero elements, based on specified thresholds. The pruning occurs in multiple iterations.
#'
#' @param incidence An incidence matrix to be pruned.
#' @param min.col Minimum number of nonzero elements (columns) required for retention. Default is 5 percent of the number of rows in the original matrix.
#' @param min.row Minimum number of nonzero elements (rows) required for retention. Default is 5.
#' @param silent Logical indicating whether to suppress output messages. Default is FALSE.
#'
#' @return A pruned incidence matrix.
#' @export
#'
#' @examples
#' # Example usage:
#' # incidence_matrix <- ...
#' # pruned_matrix <- prune.incidence(incidence_matrix)
#'
#' @seealso
#' This function is used within other functions. See \code{\link{documents.to.topics}} and \code{\link{mca.from.topics}}.
prune.incidence <- function(incidence, min.col = nrow(incidence) * 0.05, min.row = 5, silent = FALSE){
  inc <- incidence


  level.up <- function(inc, k, j = 3){
    while(
      all(!is.null(dim(inc)),
          any(
            c(
              any(rowSums(inc) < j), # Is there any individuals with less than j positions
              any(colSums(inc) < k)  # Is there any affiliations with less than k members
            )
          )
      )
    ){
      inc         <- inc[rowSums(inc) >= j, , drop = FALSE] # Keep only those members with j or more positions #! We need to subset with something different from base functions - something that always returns a
      inc         <- inc[, colSums(inc) >= k, drop = FALSE]  # Keep only those affiliations with more than k members
    }
    inc
  }

  inc <- level.up(inc, j = min.row, k = min.col)

  if(identical(silent, FALSE)) cat("\n", "Rows: ", nrow(inc), "(", nrow(inc)-nrow(incidence),")", "  Cols: ", ncol(inc), "(", ncol(inc)-ncol(incidence),")")

  inc
}

#' Convert simple triplet matrix to sparse matrix
#'
#' This function converts a simple triplet matrix in a sparse format to a sparse matrix object. It is implemented in other 'topicspace' functions to convert the 'docterm' document-term matrix to a sparse matrix object before passing it to the LDA model fitting process.
#'
#' @param simple_triplet_matrix_sparse A simple triplet matrix in a sparse format.
#'
#' @return A sparse matrix object created from the provided simple triplet matrix.
#' @export
#'
#' @examples
#' # Create a simple triplet matrix
#' stm <- list(
#'   i = c(1, 2, 3),
#'   j = c(1, 2, 3),
#'   v = c(4, 5, 6),
#'   nrow = 3,
#'   ncol = 3
#' )
#'
#' # Convert the simple triplet matrix to a sparse matrix
#' sparse_mat <- as.sparseMatrix(stm)
#' sparse_mat
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{topics.to.mca}}.
as.sparseMatrix <- function(simple_triplet_matrix_sparse) {

  Matrix::sparseMatrix(
    i = simple_triplet_matrix_sparse$i,
    j = simple_triplet_matrix_sparse$j,
    x = simple_triplet_matrix_sparse$v,
    dims = c(
      simple_triplet_matrix_sparse$nrow,
      simple_triplet_matrix_sparse$ncol
    ),
    dimnames = dimnames(simple_triplet_matrix_sparse)
  )

}


#' Step 1. Generate topic models from documents
#'
#' This function generates topic models from a dataset of abstracts. The dataset format is a slightly reformatted version of the export format from Zotero, see \link{abstracts}.
#'
#' @param abstracts A dataframe including the minimum required variables: ID, documents, title, start.year.
#' @param label An optional label for the output. Default is an empty string.
#' @param lda.k.seq A sequence of integers specifying the number of topics for each model.
#' @param lda.prune.row Minimum number of unique terms required in a document for it to be included in the model. Default is 25.
#' @param lda.prune.col Minimum number of unique documents required for a term to be included in the model. Default is 5.
#' @param filter.terms An optional list of terms to filter the models.
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{l.topics}{A list of topic models generated for different values of the lda.k.seq argument. Each element of the list corresponds to a different topic model, where the key is the number of topics used for that particular model.}
#'   \item{tidyterms}{A tibble containing cleaned and formatted tidy terms. See also 'topic.model.ready.tidy.terms' function.}
#'   \item{docterm}{A document-term matrix. Each row represents a document, each column represents a term, and values are the frequency of the corresponding term in the respective document. The document-term matrix is pruned to remove terms and documents that do not meet certain criteria before being used in the topic modeling process. See also 'prune.incidence' function.}
#' }
#' @export
#'
#' @examples
#' # Load example abstracts data
#' data(abstracts)
#' # Set minimum requirements for abstract data
#' example(abstracts)
#' # Run function with default settings
#' documents.to.topics(abstracts)
#' @seealso
#' \code{\link{prune.incidence}}, \code{\link{topic.model.ready.tidy.terms}}, \code{\link{docterm.to.lda}}.
documents.to.topics <- function(abstracts, label = "", lda.k.seq = seq(from = 50, to = 200, by = 50), lda.prune.row = 25, lda.prune.col = 5, filter.terms = NULL){

  tidyterms <- topic.model.ready.tidy.terms(abstracts)
  #tidyterms <- tidyterms %>% filter()
  docterm   <- tidyterms[tidyterms$tf_idf < 1 & tidyterms$tf_idf > 0, ] %>% tidytext::cast_dtm(ID, stem, n)

  # Ensure unique column names
  colnames(docterm) <- make.names(colnames(docterm), unique = TRUE)

  # Pruning
  # Terms and documents need to have at least 25 unique terms and a term has to be in at least 5 unique documents

  doc.set <- (as.sparseMatrix(docterm) > 0) %>% prune.incidence(., min.row = lda.prune.row, min.col = lda.prune.col)
  docterm <- docterm[rownames(doc.set) , colnames(doc.set)]

  # Fit model
  # List of topic models - with varying K
  k                      <- lda.k.seq
  l.topics               <- map(k, docterm.to.lda, docterm = docterm)
  names(l.topics)        <- k

  # Save output
  list(l.topics, tidyterms, docterm)
}

#' Make a list of MCA from a list of topic models
#'
#' This function performs multiple correspondence analysis (MCA) from a list of topic models as returned from the 'documents.to.topics' function.
#'
#' @param l.topics A list of topics. Refers to the first element of the list of three as returned from the 'documents.to.topics' function.
#' @param label An optional label for the output. Default is an empty string.
#' @param mca.cut.off.seq A sequence of values specifying the cut-off for MCA analysis.
#' @param mca.k.set A set of integers specifying the number of topics for the MCA.
#'
#' @return A list of MCA results for each specified combination of cut-off and number of topics considered.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' # Run function with default settings
#' topics.to.mca(tm_list[[1]])
#' @seealso
#' \code{\link{documents.to.topics}}.
topics.to.mca    <- function(l.topics, label = "", mca.cut.off.seq = seq(from = 0.01, to = 0.05, by = 0.01), mca.k.set = c(25, 50, 100, 150, 200)){

  s.seq            <- mca.cut.off.seq %>% set_names(.,.)
  set              <- mca.k.set
  set              <- which(as.numeric(names(l.topics)) %in% set)

  l.mca.seqs       <- purrr::map(l.topics[set], function(x, s.seq)  map(s.seq, mca.from.topics, topics = x), s.seq = s.seq)
  l.mca.seqs
}



# Section 2: Topic Modeling Functions
# ------------------------------------

#' Step 2. Perform LDA topic modeling from a defined number of topics
#'
#' This function performs LDA topic modeling from a defined number of topics.
#'
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function, weighted by the .
#' @param n_topics A numeric value for a k number of topics. Default is 100. For comparing between different solutions, please refer to the 'topic.pseudo.loglikelihood' function.
#'
#' @return A tibble containing the distribution of topics for each document.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' # Run function with default settings
#' docterm.to.lda(tm_list[[3]])
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{topic.pseudo.loglikelihood}}, \code{\link{topics.to.mca}}, \href{https://cran.r-project.org/web/packages/text2vec/index.html}{text2vec package}.
docterm.to.lda   <- function(docterm, n_topics = 100){

  lda_model       <- text2vec::LDA$new(n_topics = n_topics)
  doc_topic_distr <- lda_model$fit_transform(as.sparseMatrix(docterm), n_iter = 1000)

  topic.names               <- lda_model$get_top_words(lambda = 1) %>% as_tibble(.name_repair = "minimal") %>% map(head, 3) %>% map(paste, collapse = ".") %>% unlist()
  colnames(doc_topic_distr) <- topic.names

  topics                               <- doc_topic_distr %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  attr(topics, "pseudo_loglikelihood") <- lda_model$.__enclos_env__$private$calc_pseudo_loglikelihood()
  topics
}


# Section 3: Topics to space functions
# ------------------------------------

#' Step 3. Make MCA from topics
#'
#' This function performs an MCA based on the topics as returned from the 'docterm.to.lda' function.
#'
#' @param topics A dataframe containing the distribution of topics for each document as returned from the 'docterm.to.lda' function.
#' @param s A numeric value specifying the cut-off threshold for MCA analysis. Default is 0.01.
#' @param small Logical. If FALSE (the default), the function will retain additional matrices in the output for potential further analysis. If TRUE, additional matrices are excluded to save memory. This is useful when only the primary MCA results are needed and memory conservation is important.
#' @return An object of class 'soc.ca' containing the MCA results. For more information on the different elements of the 'soc.ca' object, please refer to the documentation of the 'soc.ca' package.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' # Run function with default settings
#' mca.from.topics(topics)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{map.topics}}, \code{\link{map.documents}}.
mca.from.topics  <- function(topics, s = 0.01, small = FALSE){
  d                <- select(topics, -ID) >= s
  rownames(d)      <- topics$ID
  d                <- prune.incidence(d, min.row = 4)
  d                <- data.frame(d)
  if(any(dim(d) == 0)) return(NA)

  res              <- soc.ca::soc.mca(active = d, identifier = rownames(d))

  if(identical(small, TRUE)) res$indicator.matrix.active  <- NULL  # We do this to save memory
  if(identical(small, TRUE)) res$indicator.matrix.passive <- NULL
  res$cor.ind                  <- res$cor.ind[, 1:5]
  res$ctr.ind                  <- res$ctr.ind[, 1:5]
  res$coord.ind                <- res$coord.ind[, 1:5]
  res
}


#' Sample Distinctive Terms
#'
#' This function returns the sample of the distinctive terms.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param s A numeric value indicating the threshold for selecting distinctive terms. Default is 0.5.
#'
#' @return A tibble containing the sample of distinctive terms including frequency and coordinates.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' # Run function with default settings
#' sample.distinctive.terms(result, tm_list[[3]])
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{mca.from.topics}}, , \code{\link{supterms}}.
sample.distinctive.terms <- function(result, docterm, s = 0.5){
  sup.terms      <- supterms(docterm, result)
  ls             <- purrr::map(sup.terms %>% dplyr::select(-ID), average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
  md             <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)
  md.z           <- purrr::map(sup.terms %>% dplyr::select(-ID), average.coord, object = result, dim = c(1,3)) %>% map(.x = ., as_tibble, .name_repair = "minimal") %>% bind_rows(.id = "word") %>% filter(label == 1)
  md$Z           <- md.z$Y
  mds <- md %>% mutate(x = sqrt(X^2), y = sqrt(Y^2), z = sqrt(Z^2)) %>% filter((x > s | y > s | z > s))
  mds
}



#' Get coordinates for supplementary terms and risk ratios for supplementary variables
#'
#' This functions calculates the coordinates for supplementary terms and the risk ratios for supplementary variables. Supplementary variables are part of the abstracts dataframe or the authors dataframe.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param vars A dataframe with two columns, namely ID and the supplementary variable.
#'
#' @return A tibble containing the coordinates and risk ratios for the supplementary variables.
#'
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' # Set vars argument
#' vars <- abstracts %>% mutate(var = Language %in% c("en")) %>% select(ID, var)
#' # vars <- authors %>% mutate(var = grepl("Bourdieu, P", Author)) %>% select(ID, var)
#' # Run function with default settings
#' supterms.coord.and.risk.ratios(result, tm_list[[3]], vars)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{mca.from.topics}}, \code{\link{map.supvar}}
supterms.coord.and.risk.ratios <- function(result, docterm, vars){
# It should have dim arguments - or supplementary.categories should give out Z coordinates by default.
  sup.terms   <- supterms(docterm, result)
  term.coords <- soc.ca::supplementary.categories(result, sup.terms %>% purrr::map_df(as.factor) %>% select(-ID))
  term.coords <- term.coords %>% mutate(term = Variable)

  vars        <- left_join(tibble(ID = result$names.ind), vars, by = "ID")
  sup.terms   <- left_join(tibble(ID = vars$ID), sup.terms, by = "ID")

  f.rr <- function(x, y){
    rr  <- epitools::riskratio(x, y)
    rr$measure[2, 1]
  }

  # This part could be made into a map function.
  coords     <- purrr::map(vars %>% select(-ID), ~purrr::map(sup.terms[, -1], f.rr, y = .x) %>% purrr::map(as_tibble, .name_repair = "minimal") %>% bind_rows(.id = "term") %>% left_join(term.coords, ., by = "term"))
  coords     <- coords %>% bind_rows(.id = "type")

  md             <- coords
  md$deviation   <- md$value
  md$deviation[md$value < 1] <- 1/md$deviation[md$value < 1]
  md             <- md[order(md$deviation, decreasing = TRUE),]

  md$RR          <- cut(md$value, breaks = c(0, 0.5, 0.75, 1, 1/0.75, 1/0.5, Inf), include.lowest = T) %>% factor()
  md
}

#' Define supplementary terms
#'
#' This function defines supplementary terms based on a minimum number of cases for a term to be considered significant.
#'
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param minimum.cases The minimum number of cases for a term to be considered significant.
#'
#' @return A tibble with the supplementary terms associated to each document.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' # Run function with default settings
#' supterms(tm_list[[3]], result)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{mca.from.topics}}, \code{\link{map.supterms}}
supterms <- function(docterm, result, minimum.cases = result$n.ind * 0.05) {

  row.set   <- rownames(docterm) %in% result$names.ind
  docterm   <- docterm[row.set,]
  cs        <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= minimum.cases)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set), by = "ID")
  sup.words
}




# Section 4: Plotting topic spaces
# --------------------------------

#' Plot the space of topics and the core features of the 'topicspace' package
#'
#' This function displays the core features of the 'topicspace' package. It visualizes the space of topics derived from MCA results including the distribution of topics, documents, supplementary terms, and year added as supplementary variable.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param label An optional label for the output. Default is an empty string.
#' @param one.plus Label for the positive direction in the first dimension.
#' @param one.minus Label for the negative direction in the first dimension.
#' @param two.plus Label for the positive direction in the second dimension.
#' @param two.minus Label for the negative direction in the second dimension.
#' @param three.plus Label for the positive direction in the third dimension.
#' @param three.minus Label for the negative direction in the third dimension.
#' @param guess.labels Logical. If TRUE, the function will attempt to guess labels for each direction based on the most contributing categories.
#' @param repel.text.size Size of text to be repelled in the plots.
#' @param browse Logical. If TRUE, the browser is invoked for debugging purposes.
#'
#' @return A list containing MCA results, coordinates and risk ratio for supplementary variable, together with a series of plots to visualize the distribution of topics and abstracts, supplementary terms, and year added as supplementary variable.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' docterm <- tm_list[[3]]
#'
#' # Run function with default settings
#' from.mca.to.results(result, docterm)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{mca.from.topics}}.
from.mca.to.results <- function(result, docterm, label = "",
                                one.plus = "+ unknown", one.minus = "- unknown", two.plus = "+ unknown", two.minus = "- unknown", three.plus = "+ unknown", three.minus = "- unknown",
                                guess.labels = TRUE, repel.text.size = 3, browse = FALSE){

  if(identical(browse, TRUE)) browser()

  # Data and abstracts
  #file <- paste0("saved/", label, "_", "1_data.Rda")
  #load(file)
  #file <- paste0("saved/", label, "_", "2_topics.Rda")
  #load(file)

  cs <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= 20)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set))
  #sup.words <- tibble("ID" = result$names.ind) %>% left_join(., doc.set)

  if(guess.labels == TRUE){
    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(3,1)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words.3   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    one.plus       <- md.sup.words %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    one.minus      <- md.sup.words %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("-", .)
    two.plus       <- md.sup.words %>% arrange(-Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    two.minus      <- md.sup.words %>% arrange(Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.plus     <- md.sup.words.3 %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.minus    <- md.sup.words.3 %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
  }

  dim1.label     <- paste("Dim 1:", one.plus, "/", one.minus)
  dim2.label     <- paste("Dim 2:", two.plus, "/", two.minus)
  dim3.label     <- paste("Dim 3:", three.plus, "/", three.minus)

  fix.logical.in.extract <- function(md){
    md$label     <- sub(pattern = ": TRUE", replacement = "", md$Modality) %>% sub("..", ", ",. , fixed = T) %>% sub("..", ", ",. , fixed = T)
    md$label     <- sub(pattern = ": FALSE", replacement = "", md$label) %>% sub("..", ", ",. , fixed = T) %>% sub("..", ", ",. , fixed = T)
    md$is.true   <- !grepl(pattern = ": FALSE", md$Modality)
    md
  }

  active.abstracts <- abstracts %>% filter(ID %in% result$names.ind)


  md           <- soc.ca::extract_mod(result) %>% fix.logical.in.extract()
  md           <- md[md$ctr.set,]


  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, color = is.true, label = label, family = "sans"), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Space of topics", subtitle = paste0("The ", nrow(md), " most contributing categories"))
  p.map.active.terms.12 <- p + ggplot2::scale_color_manual(values = c("red", "black"), name = "Has topic") + coord_fixed()
  p.map.active.terms.12


  md           <- soc.ca::extract_mod(result, c(2,3)) %>% fix.logical.in.extract()
  md           <- md[md$ctr.set,]

  p              <- soc.ca::map.ca.base(right = two.plus, left = two.minus, up = three.plus, down = three.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, color = is.true, label = label, family = "sans"), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Space of topics", subtitle = paste0("The ", nrow(md), " most contributing topics"))
  p.map.active.terms.23 <- p + ggplot2::scale_color_manual(values = c("red", "black"), name = "Has topic") + coord_fixed()
  p.map.active.terms.23

  # Most contributing abstracts -----
  md.ind       <- soc.ca::extract_ind(result, dim = c(1,2))
  md.ind.3     <- soc.ca::extract_ind(result, dim = c(1,3))
  md.ind$Z     <- md.ind.3$Y
  md.ind$quadrant12 <- soc.ca::create.quadrant(result)
  md.ind$quadrant23 <- soc.ca::create.quadrant(result, dim = 2:3)

  md.ind$ID <- md.ind$Individual
  md.join    <- left_join(md.ind, active.abstracts) %>% as_tibble(.name_repair = "minimal")

  md.join %>% select(X, title) %>% arrange(-X) %>% head(5)
  md.join %>% select(X, title) %>% arrange(X) %>% head(5)
  md.join %>% select(Y, title) %>% arrange(-Y) %>% head(5)
  md.join %>% select(Y, title) %>% arrange(Y) %>% head(5)

  # Space of abstracts -----

  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = X, y = Y), size = 0.5, color = "grey60") + geom_density_2d(data = md.join, aes(x = X, y = Y), color = "darkred")
  p.map.ind.12   <- p + coord_fixed() + ggtitle("")  + ggtitle("Cloud of documents", subtitle = "1st and 2nd dim.")

  p              <- soc.ca::map.ca.base(right = two.plus, left = two.minus, up = three.plus, down = three.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = Y, y = Z), size = 0.5, color = "grey60") + geom_density_2d(data = md.join, aes(x = Y, y = Z), color = "darkred")
  p.map.ind.23   <- p + coord_fixed() + ggtitle("Cloud of documents", subtitle = "2nd and 3rd dim.")

  # Space of abstracts with time ----
  # Ideally we project the years into the space - but not now
  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = X, y = Y, color = start.year), size = 1)
  p              <- p + ggplot2::scale_color_gradient2(low = "black", mid = "papayawhip", high = "darkred", midpoint = median(md.join$start.year), guide = "legend", name = "Year")
  p.map.ind.year.12   <- p + coord_fixed() + labs(title = "Cloud of documents by year", subtitle = "1st and 2nd dim.", caption = paste("Median year: ", md.join$start.year %>% median))

  p              <- soc.ca::map.ca.base(right = two.plus, left = two.minus, up = three.plus, down = three.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = Y, y = Z, color = start.year), size = 1)
  p              <- p + ggplot2::scale_color_gradient2(low = "black", mid = "papayawhip", high = "darkred", midpoint = median(md.join$start.year), guide = "legend", name = "Year")
  p.map.ind.year.23   <- p + coord_fixed() + labs(title = "Cloud of documents by year", subtitle = "2nd and 3rd dim.", caption = paste("Median year: ", md.join$start.year %>% median))

  # Space of supplementary terms -----

  # Term frequency
  ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
  md             <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, label = word, family = "sans", color = log(Freq)), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Frequent terms in the space of topics", subtitle = "1st and 2nd dim.", caption = "Terms have occurences in at least 25 abstracts")
  p.map.words.12 <- p + coord_fixed() + ggplot2::scale_color_gradient(high = "red", low = "black", guide = "legend")
  p.map.words.12

  ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(2,3)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
  md             <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)
  p              <- soc.ca::map.ca.base(right = two.plus, left = two.minus, up = three.plus, down = three.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, label = word, family = "sans", color = log(Freq)), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Frequent terms in the space of topics", subtitle = "2nd and 3rd dim.", caption = "Terms have occurences in at least 25 abstracts")
  p.map.words.23 <- p + coord_fixed() + ggplot2::scale_color_gradient(high = "red", low = "black", guide = "legend")
  p.map.words.23

  # Term development with RR

  md.join  <- md.join %>% mutate(added = start.year <= median(start.year))

  vars <- md.join %>% select(ID, added)
  md  <- supterms.coord.and.risk.ratios(result, docterm, vars)

  md             <- md %>% filter(label == "TRUE")
  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + geom_text(data = md, aes(x = X, y = Y, label = term, family = "sans", color = RR), check_overlap = T, size = 3)
  p              <- p + labs(title = "Date added for terms in the space of topics", subtitle = paste0("RR 1 is the median year (", median(md.join$start.year), ")"))
  p              <- p + scale_color_brewer(type = "div", palette = "RdYlBu")
  p.term.add.rr  <- p



  list(result, md.join,
       p.map.active.terms.12, p.map.active.terms.23,
       p.map.words.12, p.map.words.23,
       p.map.ind.12, p.map.ind.23,
       p.map.ind.year.12, p.map.ind.year.23,
       p.term.add.rr
  )

}


#' Map the space of topics
#'
#' This is a plotting function to map the space of topics.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param dim The dimensions in the order they are to be plotted. The first number defines the horizontal axis and the second number defines the vertical axis.
#' @param label An optional label for the output. Default is an empty string.
#' @param one.plus Label for the positive direction in the first dimension.
#' @param one.minus Label for the negative direction in the first dimension.
#' @param two.plus Label for the positive direction in the second dimension.
#' @param two.minus Label for the negative direction in the second dimension.
#' @param three.plus Label for the positive direction in the third dimension.
#' @param three.minus Label for the negative direction in the third dimension.
#' @param guess.labels Logical. If TRUE, the function will attempt to guess labels for each direction based on the most contributing categories.
#' @param repel.text.size Size of text to be repelled in the plots.
#' @param browse Logical. If TRUE, the browser is invoked for debugging purposes.
#'
#' @return A plot of the space of topics.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' docterm <- tm_list[[3]]
#' # Run function with default settings
#' map.topics(result, docterm)
map.topics <- function(result, docterm, dim = c(1,2), label = "",
                       one.plus = "+ unknown", one.minus = "- unknown", two.plus = "+ unknown", two.minus = "- unknown", three.plus = "+ unknown", three.minus = "- unknown",
                       guess.labels = TRUE, repel.text.size = 3, browse = FALSE){

  if(identical(browse, TRUE)) browser()

  cs <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= 20)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set))

  if(guess.labels == TRUE){
    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(3,1)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words.3   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    one.plus       <- md.sup.words %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    one.minus      <- md.sup.words %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("-", .)
    two.plus       <- md.sup.words %>% arrange(-Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    two.minus      <- md.sup.words %>% arrange(Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.plus     <- md.sup.words.3 %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.minus    <- md.sup.words.3 %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
  }

  dim1.label     <- paste("Dim 1:", one.plus, "/", one.minus)
  dim2.label     <- paste("Dim 2:", two.plus, "/", two.minus)
  dim3.label     <- paste("Dim 3:", three.plus, "/", three.minus)

  fix.logical.in.extract <- function(md){
    md$label     <- sub(pattern = ": TRUE", replacement = "", md$Modality) %>% sub("..", ", ",. , fixed = T) %>% sub("..", ", ",. , fixed = T)
    md$label     <- sub(pattern = ": FALSE", replacement = "", md$label) %>% sub("..", ", ",. , fixed = T) %>% sub("..", ", ",. , fixed = T)
    md$is.true   <- !grepl(pattern = ": FALSE", md$Modality)
    md
  }

  active.abstracts <- abstracts %>% filter(ID %in% result$names.ind)


  md           <- soc.ca::extract_mod(result, dim) %>% fix.logical.in.extract()
  md           <- md[md$ctr.set,]


  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, color = is.true, label = label, family = "sans"), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Space of topics", subtitle = paste0("The ", nrow(md), " most contributing categories"))
  p.map.active.terms.12 <- p + ggplot2::scale_color_manual(values = c("red", "black"), name = "Has topic") + coord_fixed()
  p.map.active.terms.12
}

#' Map supplementary terms in the space of topics
#'
#' This is a plotting function to map supplementary terms in the space of topics.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param dim The dimensions in the order they are to be plotted. The first number defines the horizontal axis and the second number defines the vertical axis.
#' @param label An optional label for the output. Default is an empty string.
#' @param one.plus Label for the positive direction in the first dimension.
#' @param one.minus Label for the negative direction in the first dimension.
#' @param two.plus Label for the positive direction in the second dimension.
#' @param two.minus Label for the negative direction in the second dimension.
#' @param three.plus Label for the positive direction in the third dimension.
#' @param three.minus Label for the negative direction in the third dimension.
#' @param guess.labels Logical. If TRUE, the function will attempt to guess labels for each direction based on the most contributing categories.
#' @param repel.text.size Size of text to be repelled in the plots.
#' @param browse Logical. If TRUE, the browser is invoked for debugging purposes.
#'
#' @return A plot of the supplementary terms.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' docterm <- tm_list[[3]]
#' # Run function with default settings
#' map.supterms(result, docterm)
map.supterms <- function(result, docterm, dim = c(1,2), label = "",
                         one.plus = "+ unknown", one.minus = "- unknown", two.plus = "+ unknown", two.minus = "- unknown", three.plus = "+ unknown", three.minus = "- unknown",
                         guess.labels = TRUE, repel.text.size = 3, browse = FALSE){

  if(identical(browse, TRUE)) browser()

  cs <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= 20)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set))

  if(guess.labels == TRUE){
    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(3,1)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words.3   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    one.plus       <- md.sup.words %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    one.minus      <- md.sup.words %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("-", .)
    two.plus       <- md.sup.words %>% arrange(-Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    two.minus      <- md.sup.words %>% arrange(Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.plus     <- md.sup.words.3 %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.minus    <- md.sup.words.3 %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
  }

  dim1.label     <- paste("Dim 1:", one.plus, "/", one.minus)
  dim2.label     <- paste("Dim 2:", two.plus, "/", two.minus)
  dim3.label     <- paste("Dim 3:", three.plus, "/", three.minus)


  ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim) %>% map(.x = ., as_tibble, .name_repair = "minimal")
  md             <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_text(data = md, aes(x = X, y = Y, label = word, family = "sans", color = log(Freq)), check_overlap = T, size = repel.text.size)
  p              <- p + labs(title = "Frequent terms in the space of topics", subtitle = "1st and 2nd dim.", caption = "Terms have occurences in at least 25 abstracts")
  p.map.words.12 <- p + coord_fixed() + ggplot2::scale_color_gradient(high = "red", low = "black", guide = "legend")
  p.map.words.12

}


#' Map supplementary variables in the space of topics
#'
#' This is a plotting function to map supplementary variables in the space of topics. The function calculates and maps the risk ratios of terms according to a supplementary variable from the 'abstracts' or the 'authors' dataframe.
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param vars A dataframe with two columns, namely ID and the supplementary variable with two categories (TRUE/FALSE).
#' @param label An optional label for the output. Default is an empty string.
#' @param one.plus Label for the positive direction in the first dimension.
#' @param one.minus Label for the negative direction in the first dimension.
#' @param two.plus Label for the positive direction in the second dimension.
#' @param two.minus Label for the negative direction in the second dimension.
#' @param three.plus Label for the positive direction in the third dimension.
#' @param three.minus Label for the negative direction in the third dimension.
#' @param guess.labels Logical. If TRUE, the function will attempt to guess labels for each direction based on the most contributing categories.
#' @param repel.text.size Size of text to be repelled in the plots.
#' @param browse Logical. If TRUE, the browser is invoked for debugging purposes.
#'
#' @return A plot displaying the risk ratios of terms according to a supplementary variable.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' docterm <- tm_list[[3]]
#' vars <- abstracts %>% mutate(var = Language %in% c("en")) %>% select(ID, var)
#' # vars <- authors %>% mutate(var = grepl("Bourdieu, P", Author)) %>% select(ID, var)
#' # Run function with default settings
#' map.supvar(result, docterm, vars)
map.supvar <- function(result, docterm, vars, label = "",
                       one.plus = "+ unknown", one.minus = "- unknown", two.plus = "+ unknown", two.minus = "- unknown", three.plus = "+ unknown", three.minus = "- unknown",
                       guess.labels = TRUE, repel.text.size = 3, browse = FALSE){

  #if(identical(browse, TRUE)) browser()

  cs <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= 20)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set))

  if(guess.labels == TRUE){
    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(3,1)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words.3   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    one.plus       <- md.sup.words %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    one.minus      <- md.sup.words %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("-", .)
    two.plus       <- md.sup.words %>% arrange(-Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    two.minus      <- md.sup.words %>% arrange(Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.plus     <- md.sup.words.3 %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.minus    <- md.sup.words.3 %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
  }

  dim1.label     <- paste("Dim 1:", one.plus, "/", one.minus)
  dim2.label     <- paste("Dim 2:", two.plus, "/", two.minus)
  dim3.label     <- paste("Dim 3:", three.plus, "/", three.minus)


  md  <- supterms.coord.and.risk.ratios(result, docterm, vars)

  md             <- md %>% filter(label == "TRUE")
  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_text(data = md, ggplot2::aes(x = X, y = Y, label = term, family = "sans", color = RR), check_overlap = TRUE, size = repel.text.size)

  # Check the number of unique colors (categories) in RR
  if(length(unique(md$RR)) == 2) {
    # Use red and blue for two categories of RR
    p <- p + ggplot2::scale_color_manual(values = c("#D73027", "#4575B4"))
  } else {
    # Use default RdYlBu palette for more than two categories
    p <- p + ggplot2::scale_color_brewer(type = "div", palette = "RdYlBu")
  }

  p.term.add.rr  <- p
  p.term.add.rr
}


#' Map the distribution of documents in the space of topics
#'
#' This is a plotting function to map the distribution of documents
#'
#' @param result An object of class 'soc.ca' containing the MCA results as returned from the 'mca.from.topics' function.
#' @param docterm A document-term matrix as returned from the 'documents.to.topics' function.
#' @param label An optional label for the output. Default is an empty string.
#' @param one.plus Label for the positive direction in the first dimension.
#' @param one.minus Label for the negative direction in the first dimension.
#' @param two.plus Label for the positive direction in the second dimension.
#' @param two.minus Label for the negative direction in the second dimension.
#' @param three.plus Label for the positive direction in the third dimension.
#' @param three.minus Label for the negative direction in the third dimension.
#' @param guess.labels Logical. If TRUE, the function will attempt to guess labels for each direction based on the most contributing categories.
#' @param repel.text.size Size of text to be repelled in the plots.
#' @param browse Logical. If TRUE, the browser is invoked for debugging purposes.
#'
#' @return A list with two plots: the distribution of documents on axes 1 and 2, and the distribution of documents on axes 2 and 3.
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' topics <- docterm.to.lda(tm_list[[3]])
#' result <- mca.from.topics(topics)
#' docterm <- tm_list[[3]]
#' # Run function with default settings
#' map.documents(result, docterm)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{docterm.to.lda}}, \code{\link{mca.from.topics}}.
map.documents <- function(result, docterm, label = "",
                                one.plus = "+ unknown", one.minus = "- unknown", two.plus = "+ unknown", two.minus = "- unknown", three.plus = "+ unknown", three.minus = "- unknown",
                                guess.labels = TRUE, repel.text.size = 3, browse = FALSE){

  if(identical(browse, TRUE)) browser()

  cs <- docterm %>% {as.sparseMatrix(.) > 0} %>% Matrix::colSums()
  set       <- which(cs >= 20)
  doc.set   <- docterm[, set] %>% {as.matrix(.) > 0} %>% as_tibble(rownames = "ID", .name_repair = "minimal")
  sup.words <- tibble("ID" = result$names.ind) %>% left_join(., as.data.frame(doc.set))
  #sup.words <- tibble("ID" = result$names.ind) %>% left_join(., doc.set)

  if(guess.labels == TRUE){
    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(1,2)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    ls             <- map(sup.words %>% select(-ID), soc.ca::average.coord, object = result, dim = c(3,1)) %>% map(.x = ., as_tibble, .name_repair = "minimal")
    md.sup.words.3   <- ls %>% bind_rows(.id = "word") %>% filter(label == 1)

    one.plus       <- md.sup.words %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    one.minus      <- md.sup.words %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("-", .)
    two.plus       <- md.sup.words %>% arrange(-Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    two.minus      <- md.sup.words %>% arrange(Y) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.plus     <- md.sup.words.3 %>% arrange(-X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
    three.minus    <- md.sup.words.3 %>% arrange(X) %>% head(3) %>% .[["word"]] %>% paste(collapse = ", ") %>% paste("+", .)
  }

  dim1.label     <- paste("Dim 1:", one.plus, "/", one.minus)
  dim2.label     <- paste("Dim 2:", two.plus, "/", two.minus)
  dim3.label     <- paste("Dim 3:", three.plus, "/", three.minus)

  active.abstracts <- abstracts %>% filter(ID %in% result$names.ind)

  # Most contributing abstracts -----
  md.ind       <- soc.ca::extract_ind(result, dim = c(1,2))
  md.ind.3     <- soc.ca::extract_ind(result, dim = c(1,3))
  md.ind$Z     <- md.ind.3$Y
  md.ind$quadrant12 <- soc.ca::create.quadrant(result)
  md.ind$quadrant23 <- soc.ca::create.quadrant(result, dim = 2:3)

  md.ind$ID <- md.ind$Individual
  md.join    <- left_join(md.ind, active.abstracts) %>% as_tibble(.name_repair = "minimal")

  md.join %>% select(X, title) %>% arrange(-X) %>% head(5)
  md.join %>% select(X, title) %>% arrange(X) %>% head(5)
  md.join %>% select(Y, title) %>% arrange(-Y) %>% head(5)
  md.join %>% select(Y, title) %>% arrange(Y) %>% head(5)

  # Space of abstracts -----

  p              <- soc.ca::map.ca.base(right = one.plus, left = one.minus, up = two.plus, down = two.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = X, y = Y), size = 0.5, color = "grey60") + geom_density_2d(data = md.join, aes(x = X, y = Y), color = "darkred")
  p.map.ind.12   <- p + coord_fixed() + ggtitle("Cloud of documents", subtitle = "1st and 2nd dim.")

  p              <- soc.ca::map.ca.base(right = two.plus, left = two.minus, up = three.plus, down = three.minus)
  p              <- p + ggplot2::geom_point(data = md.join, aes(x = Y, y = Z), size = 0.5, color = "grey60") + geom_density_2d(data = md.join, aes(x = Y, y = Z), color = "darkred")
  p.map.ind.23   <- p + coord_fixed() + ggtitle("Cloud of documents", subtitle = "2nd and 3rd dim.")

  list(
       p.map.ind.12, p.map.ind.23
  )

}




# Section 5: Supplementary functions
# ----------------------------------

#' Topic pseudo-loglikelihood
#'
#' This function provides a plot of pseudo-loglikelihood for different k number of topics.
#'
#' @param l.topics A list of topic models as returned from 'documents.to.topics' function. Topics are generated for different values of the lda.k.seq argument in the 'documents.to.topics' function.
#' @param lda.k.seq A sequence of integers specifying the number of topics for each model. Must be the same as defined in the 'documents.to.topics' function.
#' @param n_topics A numeric value specifying the k number of topics intended for LDA. Default is 100.
#'
#' @return A plot of pseudo-loglikelihood for different numbers of topics
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' # Run function with default settings
#' topic.pseudo.loglikelihood(tm_list[[1]])
#' @seealso
#' \code{\link{documents.to.topics}}, \href{https://cran.r-project.org/web/packages/text2vec/index.html}{text2vec package}.
topic.pseudo.loglikelihood <- function(l.topics, lda.k.seq = seq(from = 50, to = 200, by = 50), n_topics = 100) {
  ld <- l.topics %>% map(attr, "pseudo_loglikelihood") %>% unlist() %>% tibble::enframe()
  ld$name <- as.numeric(ld$name)
  ld$set <- ld$name %in% lda.k.seq
  ld$set <- factor(ld$set + (ld$name == n_topics)) %>% forcats::fct_recode("Testing" = "1", "LDA" = "2") %>% forcats::fct_rev()
  ld$set.labels <- as.character(ld$name)
  ld$set.labels[ld$set == "Not selected"] <- ""
  p <- ggplot(ld, aes(x = name, y = value)) + geom_line(size = 0.3) + geom_point(aes(fill = set, size = set), shape = 21)
  p <- p + geom_text(mapping = aes(label = set.labels), vjust = -1, hjust = 1, size = 3.5, family = "serif")
  p <- p + scale_fill_manual(values = c("darkred", "black", "white"), name = "Selected for:") + scale_size_manual(values = c(4, 2), name = "Selected for:")
  p <- p + ggthemes::theme_tufte() + xlab("Number of topics") + ylab("Pseudo-Loglikelihood") + labs(title = "Selecting the number of topics based on Pseudo Loglikelihood")
  p
}


#' Plot MCA variance
#'
#' This function provides a plot of MCA explained variance for axes 1 and 2 according to different mca cut-off and number of topics considered.
#'
#' @param l.mca.seqs A list of MCA results for each specified combination of cut-off and number of topics considered, as returned from the 'topics.to.mca' function.
#' @param n_topics A numeric value specifying the k number of topics intended for LDA. Default is 100.
#'
#' @return A ggplot of with the explained variance per cut off and number of topics
#' @export
#'
#' @examples
#' # Requirements
#' example(abstracts)
#' tm_list <- documents.to.topics(abstracts)
#' l.mca.seq <- topics.to.mca(tm_list[[1]])
#' # Run function with default settings
#' mca.variance.plot(l.mca.seq)
#' @seealso
#' \code{\link{documents.to.topics}}, \code{\link{topics.to.mca}}
mca.variance.plot <- function(l.mca.seqs, n_topics = 100) {
  # Helper function to clean the list
  f.clean <- function(x) {
    x[is.na(x)] <- NULL
    x
  }

  # Helper function to get the variance
  f.get.var <- function(x) {
    x$adj.inertia[1:2, "Adj.Var"] %>% sum() %>% unlist()
  }

  # Process the l.mca.seqs list
  ld.mca.var <- map_depth(l.mca.seqs, f.clean, .depth = 1) %>%
    map_depth(., f.get.var, .depth = 2) %>%
    map(unlist) %>%
    map(tibble::enframe, name = "Cutoff", value = "Variance") %>%
    bind_rows(.id = "Topics")

  # Prepare data for plotting
  ld.mca.var <- ld.mca.var %>%
    group_by(Topics) %>%
    mutate(topic.max = (Variance == max(Variance)) * 1)

  ld.mca.var$Topics <- as.numeric(ld.mca.var$Topics) %>% forcats::as_factor()
  ld.mca.var$label <- ld.mca.var$Cutoff
  ld.mca.var$label[ld.mca.var$topic.max == 0] <- ""
  ld.mca.var$selected <- as.numeric(ld.mca.var$Topics == n_topics) %>% forcats::as_factor()

  colors <- c("white", "black")
  line.colors <- paletteer::paletteer_d(palette = "awtools::a_palette")

  # Create the plot
  p <- ggplot(ld.mca.var, aes(x = as.numeric(Cutoff), y = Variance / 100, color = Topics, fill = selected, group = Topics, label = label)) +
    geom_line(linewidth = 1) +
    geom_point(aes(size = topic.max), shape = 21, color = "black") +
    geom_text(color = "black", hjust = 1, vjust = -1, family = "serif", size = 3.5) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1), labels = scales::percent, limits = c(0, 1)) +
    scale_fill_manual(values = colors, name = "Number of topics", guide = "none") +
    scale_color_manual(values = line.colors, name = "Number of topics") +
    scale_size(range = c(1.5, 3), guide = "none") +
    ggthemes::geom_rangeframe(color = "black") +
    xlab("Topic prediction cutoff") +
    ylab("% explained variance") +
    ggthemes::theme_tufte() +
    labs(title = "Selecting topic fit cutoff with explained variance on dimension 1 and 2",
         subtitle = "The selected MCA's are highlighted in black for the main analysis and in white for the tests")

  p
}
