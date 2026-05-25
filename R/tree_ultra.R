#' Phylogenetic tree (ultrametric): tree_ultra
#'
#' An ultrametric phylogenetic tree representing evolutionary distances among
#' 33 species. Species names have been anonymised using capital letters (A, B,
#' C, etc.) for illustrative purpose. This tree can be used for examples
#' involving phylogenetic models, such as \code{\link[phyr]{pglmm}} as used in \code{\link{dredge_pglmm}}.
#'
#' @format An object of class \code{phylo}, with 33 tips and branch lengths in
#' millions of years.
#' @details The tree was originally generated from the TimeTree database and
#' then made ultrametric using the \code{\link[phytools]{force.ultrametric}}.
#' Species names were replaced by generic labels to match those in the \code{\link{senescence}} dataset.
#'
#' @source TimeTree database (\url{https://www.timetree.org})
"tree_ultra"
