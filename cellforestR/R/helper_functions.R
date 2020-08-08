get_prefix_from_path <- function(path) {
  last_slash_pos <- tail(which(strsplit(path, "")[[1]] == "/"), n = 1)
  prefix <- substr(path, start = 1, stop = last_slash_pos)

  return(prefix)
}

#' @title Processes that happen between two paths
#'
#' @param path_upstream String path of file upstream (e.g., normalize's RDS)
#' @param path_downstream String path of file downstream (e.g., reduce's meta.tsv)
#'
#' @importFrom stringr str_remove
processes_between_paths <- function(path_upstream, path_downstream) {
  folder_up <- get_prefix_from_path(path_upstream)
  folder_down <- get_prefix_from_path(path_downstream)
  path_diff <- str_remove(string = folder_down, pattern = folder_up)
  if (path_diff == "") {
    return()
  }

  processes_and_hashes <- strsplit(path_diff, "/")[[1]]
  processes <- c()
  # get process name and skip hashes
  for (i in 1:length(processes_and_hashes)) {
    if ((i %% 2) == 1) {
      processes <- append(processes, processes_and_hashes[i])
    }
  }

  return(processes)
}
