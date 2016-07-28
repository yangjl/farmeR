#' \code{get_file2tab}
#' Scan txt files and extract information into a data.frame.
#'
#' @param files Full names of the files to scan [vector].
#' @param features Unique features for grep to search [vector].
#'
#' @return return Extracted values [data.frame].
#'
#' @examples
#' files <- list.files(path = dir, pattern = fileptn, full.names = TRUE)
#' features <- c("C methylated in CHH context")
#' get_file2tab(files, features, replace=T )
#'
#' @export
get_file2tab <- function(files, features, replace=T ){

  out <- as.data.frame(matrix(0, nrow=length(files), ncol=length(features)))
  row.names(out) <- files
  names(out) <- features
  for(i in 1:length(files)){
    text <- readLines(files[i])
    for(j in 1:length(features)){
      val1 <- grep(features[j], text, value=TRUE)
      if(replace){val1 <- gsub(features[j], "", val1)}
      if(length(val1)==0) {
        out[i,j] <- NA
      }else{
        out[i,j] <- val1
      }
    }
  }
  return(out)
}
