select_newest <- function(path, file.pattern, by = NULL) {
  # some error dealing
  if (length(path) > 1) {
    stop("'path' needs to be a single string.")
  } else if(length(path) == 0L){
    stop("'path' is empty.")
  }
  if (length(file.pattern) > 1) {
    stop("'file.pattern' needs to be a single string.")
  } else if(length(file.pattern) == 0L){
    stop("'file.pattern' is empty.")
  }
  
  # read in files in directory
  files <- list.files(path, pattern = file.pattern)
  
  if(length(files) == 0L){
    stop("No file that matches pattern exists.")
  }
  
  if(!is.null(by)){
    # do we have several files per object? -> take newest version
    for (h in 1:length(by)) {
      versions <- files[grepl(by[h], files)]
      versions.sans <-tools::file_path_sans_ext(files)
      if (length(versions) > 1) {
        find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
          col.names <- paste0("V", seq(length = length(x)))
          as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
        }))
        dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
          is.date <- try(as.Date(x), silent = T)
          if(inherits(is.date, "try-error")){
            x <- NULL
          }
          return(x)
        }))))
        
        newest <- max(dates)
        
        files <-
          c(versions[grepl(newest, versions)], files[!grepl(by[h], files)])
      }
    }
  } else {
    versions <- files
    versions.sans <-tools::file_path_sans_ext(files)
    find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
      col.names <- paste0("V", seq(length = length(x)))
      as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
    }))
    dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
      is.date <- try(as.Date(x), silent = T)
      if(inherits(is.date, "try-error")){
        x <- NULL
      }
      return(x)
    }))))
    
    newest <- max(dates, na.rm =T)
    files <- versions[grepl(newest, versions)]
  }
  return(paste(path, files, sep = '/'))
}