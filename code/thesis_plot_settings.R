#### Plot settings for printing in the dissertation (116 mm text width) ####
# Figures are sized to the full dissertation text width of 116 mm and
# at most ~18 cm high (leaving room for the caption).
# Font sizes are in pt and chosen to match ~10-11 pt body text after
# printing at 100 % scale (figures are not rescaled in LaTeX/Word).

library(ggplot2)

thesis_dpi = 600       # print resolution (all figures)
fig_width_cm = 11.6    # full text width (116 mm)
fig_max_height_cm = 18 # max. figure height incl. room for caption

title_size = 10  # panel titles
lab_size = 9     # axis titles
ax_size = 8      # axis tick labels
legend_size = 8  # legend text
annot_size = 8   # in-plot annotations (divided by .pt where used)
panel_label_size = 11 # panel labels (A, B, ...)

theme_set(theme_grey(base_size = 9))

# legend anchored inside the top-right corner of the panel, so it grows
# inwards (and never past the panel edge) when font sizes change
legend_top_right = theme(legend.position = "inside",
                         legend.position.inside = c(0.99, 0.98),
                         legend.justification = c(1, 1))

# arguments for panel labels of egg::ggarrange
panel_label_args = list(gp = grid::gpar(fontsize = panel_label_size, fontface = "bold"))

save_thesis_fig = function(filename, plot, height_cm, width_cm = fig_width_cm){
  #' save figure with thesis dimensions (116 mm text width)
  stopifnot(height_cm <= fig_max_height_cm)
  ggsave(filename = filename,
         plot = plot,
         width = width_cm,
         height = height_cm,
         units = "cm",
         dpi = thesis_dpi,
         bg = "white")
}
