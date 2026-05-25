#' Small helper to write title boxes
#'
#' Internal helper function. This function allows to write title boxes.
#'
#' @param text A \code{character} string.
#' @param width A \code{numeric} setting the width of the title box.
#' @param border A \code{character} string to set the border type.
#'
#' @noRd
.title_box <- function(text, width = 80, border = "#") {

  middle_text <- paste0("### ", text, " ###")

  padding <- max(0, width - nchar(middle_text))
  padding_left <- floor(padding / 2)
  padding_right <- ceiling(padding / 2)

  middle_line <- paste0(strrep(" ", padding_left), middle_text, strrep(" ", padding_right))

  border_line <- strrep(border, width)

  cat("\n\n", border_line, "\n", middle_line, "\n", border_line, "\n\n", sep = "")
}

#' Small helper for separation lines
#'
#' Internal helper function. This function allows to add separation lines in print methods.
#'
#' @param line A \code{character} string to set line type.
#' @param width A \code{numeric} setting the width of the title box.
#'
#' @noRd
.sep_line <- function(line_type = "-", width = 80) {
  line <- strrep(line_type, width)
  cat(line, "\n\n")
}

#' Small helper to capitalize first letter of a string
#'
#' Internal helper function. This function allows to capitalize first letter of a string.
#'
#' @param text A \code{character} string.
#'
#' @noRd
.first_upp <- function(text) {
  substr(text, 1, 1) <- toupper(substr(text, 1, 1))
  return(text)
}
