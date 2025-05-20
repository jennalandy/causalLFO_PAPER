.libPaths("~/apps/R_4.1.0/")
library(ggthemes)
library(ggplot2)

simulation_parameters <- readRDS("../application/simulation_parameters_BA.rds")
P <- simulation_parameters$P
ATE <- simulation_parameters$ATE_C

custom_colors <- c(
  "#000000",  # Neutral color light (black)
  "#999999",  # Neutral color dark (grey)

  "#56B4E9",  # Cool color light (blue)
  "#0072B2",  # Cool color dark (blue)

  "#E69F00",  # Warm color light (orange)
  "#D55E00",  # Warm color dark (orange)

  "#CC79A7"   # Warm color 2 (pink)
)


simulation_sig_order <- str_remove(colnames(P), "SBS") %>%
    as.numeric() %>%
    sort() %>%
    paste("SBS", ., sep = "")

truth_data <-  data.frame(
    sig = factor(names(ATE), levels = simulation_sig_order),
    ATE = ATE
)

truth_boxplot <- function(...) {
    geom_boxplot(
        aes(y = ATE, fill = 'black', alpha = 1, group = sig, color = 'black'),
        data = truth_data,
        show.legend = FALSE,
        ...
    )
}

method_order <- list(
    "oracle" = "Oracle",
    "observed_C" = "Observed Outcome",
    "all_data" = "All Data",
    "random_split_and_nmf" = "Random Split",
    "impute" = "Impute",
    "stabilize" = "Stabilize",
    "impute_and_stabilize" = "Impute And Stabilize",
    "impute_and_stabilize_and_decompose" = "Impute And Stabilize, Decompose"
)

set_colors <- function(name = "Okabe-Ito", subset = names(method_order)) {
    if (name == 'custom') {
        colors <- custom_colors
    } else {
        colors <- palette(name)
    }
    method_order_subset = method_order[subset]
    values <- setNames(colors[seq_along(names(method_order))], method_order)
    values <- values[unlist(method_order_subset)]
    scale_color_manual(
        values = values,
        breaks = method_order_subset
    )
}

set_fill <- function(name = "Okabe-Ito") {
    if (name == 'custom') {
        colors <- custom_colors
    } else {
        colors <- palette(name)
    }
    values <- setNames(colors[seq_along(names(method_order))], method_order)
    values_fill <- values
    values_fill['transparent'] <- 'transparent'
    method_order_fill <- method_order
    method_order_fill['transparent'] <- "pi = 0.8"
    scale_fill_manual(
        values = values_fill,
        breaks = method_order_fill
    )
}

my_theme <- function() {
    theme_bw() +
    theme(
        strip.background = element_rect(fill = 'white'),
        strip.placement = "outside",
        legend.position = 'top',
        legend.justification = 'left',
        text = element_text(size = 20),
        legend.box = "vertical",
        legend.box.just = "left"
    )
}

my_color_legend <- function(order = 1) {
    guides(color = guide_legend(nrow = 2, ncol = 4, byrow = FALSE, order = order))
}

my_fill_legend <- function(order = 2) {
    guides(fill = guide_legend(nrow = 2, ncol = 4, byrow = FALSE, order = order))
}

my_save <- function(
    file, width = 10, height = 6,
    width_scale = 1, height_scale = 1
) {
    ggsave(file, width = width*width_scale, height = height*height_scale)
}