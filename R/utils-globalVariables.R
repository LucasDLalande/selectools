# Global variables used in NSE code
# (e.g. MuMIn::subset.model.selection)
#
# Prevents R CMD check notes:
# "no visible binding for global variable"
#
# @noRd
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("df", "delta"))
}
