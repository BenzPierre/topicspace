#' topicspace: Mapping Topic Spaces Combining Latent Dirichlet Allocation with Multiple Correspondence Analysis
#'
#' @description
#' The `topicspace` package is designed for social scientists interested in exploring the relationships between discursive and social structures. It provides functions for pre-processing textual data, performing topic modeling, and mapping topic spaces by combining Latent Dirichlet Allocation (LDA) with Multiple Correspondence Analysis (MCA).
#'
#' To learn more about topicspace, please explore the vignette: browseVignettes(package = "topicspace")
#'
#' @details
#' This approach aims to enhance our understanding of the hierarchies among topics providing measures of the extent to which topics vary in their contribution to structuring topic spaces along several dimensions. Following the operational principles of MCA (Le Roux and Rouanet, 2004; Hjellbrekke, 2018), topics form the 'cloud of active variables,' while documents are represented in the 'cloud of individuals.' This allows for considering sufficiently prevalent terms as supplementary variables, referred to as 'supplementary terms,' for which we can calculate coordinates within the topic space. A distinctive feature of the package is its capability to calculate the risk ratio of these prevalent terms to be associated with both documents and author-related specific variables (Kropp and Larsen, 2023; Rossier et al., 2023; Benz et al., 2024).
#'
#' @section Package Structure:
#' The `topicspace` package is structured into four parameterizable blocks of functions, referred to as 'steps' to reflect their sequential order of execution:
#'
#' **Step 1. From Text to Terms**
#' This step transforms raw text data into a structured format, preparing it for topic modeling. It includes creating term frequency (TF), inverse document frequency (IDF), and TF-IDF values, and constructing the document-term matrix (DTM). Additionally, this step provides elements for evaluating different models, depending on the number of topics (k) and various MCA cut-offs.
#'
#' **Step 2. From Terms to Topics**
#' This step applies Latent Dirichlet Allocation (LDA) based on a specified number of topics (k). It uses the `text2vec` package by Selivanov et al. (2020) to implement LDA.
#'
#' **Step 3. From Topics to Spaces**
#' This step performs Multiple Correspondence Analysis (MCA) on the document-by-topic matrix generated in Step 2. It categorizes topics for each document based on set thresholds, resulting in a binary assignment (presence or absence of a topic) rather than probabilities. MCA is implemented with the `soc.ca` package by Larsen et al. (2021).
#'
#' **Step 4. Mapping Spaces**
#' This final step includes functions to map topic spaces, documents, supplementary terms, and risk ratios based on document or author-related variables.
#'
#' @author
#' Pierre Benz & Anton G. Larsen
#'
#' @references
#'## Some works using the `topicspace` package:
#'
#' Kropp, K., & Larsen, A. G. (2023). Changing the topics: the social sciences in EU-funded research projects.
#' Comparative European Politics, 21(2), 176–207. https://doi.org/10.1057/s41295-022-00313-5.
#'
#' Rossier, T., Benz, P., Grau Larsen, A., & Kropp, K. (2023). The Space of Research Topics in Economics:
#' Scientific Position-Takings and Individual Positions in Swiss Economic Science. Œconomia. History, Methodology, Philosophy, (13-2), 427–473.
#' https://doi.org/10.4000/oeconomia.15359.
#'
#' Benz, P., Kropp, K., Nobel, T. C., & Rossier, T. (2024). Homologies in fields of cultural production.
#' Evidence from the European scientific field. Poetics, 107, 101945. https://doi.org/10.1016/j.poetic.2024.101945.
#'
#'
#'## Some seminal references about LDA and MCA:
#'
#' Blei, D. M., Ng, A. Y., & Jordan, M. I. (2003). Latent Dirichlet Allocation. Journal of machine Learning research, 3, 993–1022.
#'
#' Hjellbrekke, J. (2018). Multiple Correspondence Analysis for the Social Sciences. Routledge.
#'
#' Le Roux, B., & Rouanet, H. (2004). Geometric Data Analysis: From Correspondence Analysis to Structured Data Analysis. Springer Science & Business Media.
#'
#' @seealso
#' Selivanov, D., Bickel, M., & Wang, Q. (2020). Package ‘text2vec’: Modern Text Mining Framework for R, 1-11.
#' Available at: <https://cran.r-project.org/web/packages/text2vec/index.html>
#'
#' Larsen, A. G., et al. (2021) Package ‘soc.ca’: Specific Correspondence Analysis for the Social Sciences.
#' Available at: <https://cran.r-project.org/web/packages/soc.ca/index.html>
#'
#' @note
#' This package is currently under development; functions and documentation may change in future versions.
#'
"_PACKAGE"
