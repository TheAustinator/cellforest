dpi = 150
width_in = plot_width_px / dpi
height_in = plot_height_px / dpi

ggsave(plot_filepath, width = width_in, height = height_in, units = "in", dpi = dpi)
