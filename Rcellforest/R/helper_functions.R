get_prefix_from_path <- function(path) {
  last_slash_pos <- tail(which(strsplit(path, "")[[1]] == "/"), n = 1)
  prefix <- substr(path, start = 1, stop = last_slash_pos)

  return(prefix)
}
