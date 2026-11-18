#### Plot settings for printing in a B5 thesis ####
# B5 paper is 176 x 250 mm. With typical margins the text block is
# roughly 12.5 cm wide and 19 cm high, so figures are sized to the full
# text width and at most ~18 cm high (leaving room for the caption).
# Font sizes are in pt and chosen to match ~10-11 pt body text after
# printing at 100 % scale (figures are not rescaled in LaTeX/Word).

library(ggplot2)

dpi = 600              # print resolution
fig_width_cm = 12.5    # full text width on B5
fig_max_height_cm = 18 # max. figure height incl. room for caption

title_size = 9   # panel titles
lab_size = 8     # axis titles
ax_size = 7      # axis tick labels
legend_size = 7  # legend text
annot_size = 7   # in-plot annotations (divided by .pt where used)
panel_label_size = 10 # panel labels (A, B, ...)

theme_set(theme_grey(base_size = 8))

# arguments for panel labels of egg::ggarrange
panel_label_args = list(gp = grid::gpar(fontsize = panel_label_size, fontface = "bold"))

save_thesis_fig = function(filename, plot, height_cm, width_cm = fig_width_cm){
  #' save figure with B5 thesis dimensions
  stopifnot(height_cm <= fig_max_height_cm)
  ggsave(filename = filename,
         plot = plot,
         width = width_cm,
         height = height_cm,
         units = "cm",
         dpi = dpi,
         bg = "white")
}
