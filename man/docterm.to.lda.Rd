% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{docterm.to.lda}
\alias{docterm.to.lda}
\title{Step 2. Perform LDA topic modeling from a defined number of topics}
\usage{
docterm.to.lda(docterm, n_topics = 100)
}
\arguments{
\item{docterm}{A document-term matrix as returned from the 'documents.to.topics' function, weighted by the .}

\item{n_topics}{A numeric value for a k number of topics. Default is 100. For comparing between different solutions, please refer to the 'topic.pseudo.loglikelihood' function.}
}
\value{
A tibble containing the distribution of topics for each document.
}
\description{
This function performs LDA topic modeling from a defined number of topics.
}
\examples{
# Requirements
example(abstracts)
tm_list <- documents.to.topics(abstracts)
# Run function with default settings
docterm.to.lda(tm_list[[3]])
}
\seealso{
\code{\link{documents.to.topics}}, \code{\link{topic.pseudo.loglikelihood}}, \code{\link{topics.to.mca}}, \href{https://cran.r-project.org/web/packages/text2vec/index.html}{text2vec package}.
}
