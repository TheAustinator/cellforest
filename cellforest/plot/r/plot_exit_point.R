dpi <- 150
width_in <- plot_width_px / dpi
height_in <- plot_height_px / dpi

if (is.null(plot_filepath) | plot_filepath == "None") {
  plot_filepath <- "/tmp/plot.png"
  print(paste0("using temp filepath: ", plot_filepath))
}
print(paste0("plot_filepath: ", plot_filepath))
ggsave(plot_filepath, width = width_in, height = height_in, units = "in", dpi = dpi, device = "png")
