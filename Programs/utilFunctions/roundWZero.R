roundWZero = function (x, digits = 0)
{
  xstring = paste0(round(x, digits))
  for (j in 1:length(x)) {
    if (xstring[j] != "NA") {
      ln = str_length(str_split(xstring[j], "\\.")[[1]][2])
      if (!is.na(ln)) {
        if (ln < digits) {
          for (i in 1:(digits - ln)) {
            xstring[j] = paste0(xstring[j], "0")
          }
        }
      }
      else {
        if (digits > 0) {
          xstring[j] = paste0(xstring[j], ".")
          for (i in 1:(digits)) {
            xstring[j] = paste0(xstring[j], "0")
          }
        }
      }
    }
  }
  return(xstring)
}
