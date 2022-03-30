

# GWAS/PheWAS -------------------------------------------------------------

#' Dummy data frame for `plot_manhattan()`
#'
#' @param type Either 'gwas' or 'phewas
#'
#' @return A data frame
#' @export
#'
#' @examples
#' dummy_manhattan_df()
dummy_manhattan_df <- function(type = "gwas") {

  match.arg(type,
            choices = c("gwas", "phewas"))

  result <- tibble::tribble(
    ~CHR, ~SNP, ~BETA,  ~SE,      ~P,
    "1",   "a",   0.1, 0.01,        0.05,
    "1",   "b",     0, 0.02,       1e-04,
    "1",   "c",  -0.1, 0.03,       1e-09,
    "2",   "d",   0.1, 0.01,        0.01,
    "2",   "e",     0, 0.02,        0.02,
    "2",   "f",  -0.1, 0.03,        0.03,
    "3",   "g",   0.1, 0.01,        0.04,
    "3",   "h",     0, 0.02,        0.05,
    "3",   "i",  -0.1, 0.03,        0.06,
    "11",  "j",   0.1, 0.01,        0.07,
    "11",  "k",     0, 0.02,        0.08,
    "11",  "l",  -0.1, 0.03,       1e-04,
    "11",  "m",   0.1, 0.01,        0.07,
    "11",  "n",     0, 0.02,        0.08,
    "11",  "o",  -0.1, 0.03,       1e-04,
  )

  if (type == "phewas") {
    result <- result %>%
      dplyr::mutate(
        "CHR" = dplyr::recode(.data[["CHR"]],
                              "1" = "Endocrine",
                              "2" = "CVD",
                              "3" = "Eye",
                              "11" = "Skin")
      )
  }

  return(result)
}

#' Take a sample of rows from non-significant results
#'
#' May be useful when plotting Manhattan or Q-Q plots for GWAS results to speed
#' up plotting.
#'
#' @param df A data frame.
#' @param pval_col Name of column in `df` containing p values.
#' @param sig Significance threshold. P values equal to or greater than this
#'   will be sampled.
#' @inheritParams dplyr::slice_sample
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' dummy_manhattan_df() %>%
#'   sample_nonsig_results(pval_col = "P",
#'                         sig = 0.001,
#'                         n = 2)
sample_nonsig_results <- function(df,
                                  pval_col = "P",
                                  sig = 5e-08,
                                  n,
                                  prop) {
  # validate
  if (!rlang::is_missing(n) & !rlang::is_missing(prop)) {
    stop("Supply exactly only one of `n` and `prop` arguments")
  }

  # split into significant and non-significant dfs
  df_split <- split(df,
                    df[[pval_col]] >= sig)

  # sample non-significant results
  if (!rlang::is_missing(prop)) {
    df_split[["TRUE"]] <-
      dplyr::slice_sample(df_split[["TRUE"]],
                          prop = prop,
                          replace = FALSE)
  } else if (!rlang::is_missing(n)) {
    df_split[["TRUE"]] <-
      dplyr::slice_sample(df_split[["TRUE"]],
                          n = n,
                          replace = FALSE)
  }

  # return result
  dplyr::bind_rows(df_split)
}

#' Basic Manhattan plotting function
#'
#' Chromosome names are labelled on the x axis in the centre of each chromosome
#' group.
#'
#' @param df A data frame.
#' @param chr Name of chromosome column. Should be type character or factor.
#' @param log10_p Name of -log10(p value) column. Should be on -log10 scale.
#' @param text Name of column with text annotation.
#' @param order_idx Name of column with order to arrange SNPs. Should be type
#'   integer.
#' @param geom_point_args A named list of arguments to be passed to
#'   [ggplot2::geom_point()].
#' @param labs_args List of arguments to be passed to [ggplot2::labs()].
#'
#' @return A ggplot plot object
plot_manhattan_basic <- function(df,
                                 chr = "CHR",
                                 log10_p = "P",
                                 text = "SNP",
                                 order_idx = "ORDER_IDX",
                                 geom_point_args = list(alpha = 1,
                                                        size = 2),
                                 labs_args = list(x = "Chromosome",
                                                  y = "-log10(p)")) {
  # validate args
  assertthat::assert_that(all(1:nrow(df) == sort(df[[order_idx]])),
                          msg = "Numbers in `order_idx` column should include all values between 1 and `nrow(df)`")

  # prepare x axis for plotting
  x_axis_df <- df %>%
    dplyr::group_by(.data[[chr]]) %>%
    dplyr::summarize("center" = (max(.data[[order_idx]]) + min(.data[[order_idx]])) / 2)

  # plot
  ggplot2::ggplot(df, ggplot2::aes(
    x = .data[[order_idx]],
    y = .data[[log10_p]],
    text = .data[[text]]
  )) +

    # Show all points
    ggplot2::geom_point(ggplot2::aes(colour = .data[[chr]]),
                        size = geom_point_args$size,
                        alpha = geom_point_args$alpha) +

    # Label X axis. Labels positioned at centre of each chromosome
    ggplot2::scale_x_continuous(label = x_axis_df[[chr]],
                                breaks = x_axis_df[["center"]]) +

    # Theme: remove background, legend, grid and panel borders
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +

    # Reduce space between plot area and x axis
    ggplot2::scale_y_continuous(expand = c(0, 0.5)) +

    # axis labels
    ggplot2::labs(x = labs_args$x,
                  y = labs_args$y)
}

plot_manhattan_gwas <- function(df,
                                chr = "CHR",
                                log10_p = "P",
                                text = "SNP",
                                order_idx = "ORDER_IDX",
                                geom_point_args = list(alpha = 1,
                                                       size = 2),
                                labs_args = list(x = "Chromosome",
                                                 y = "-log10(p value)"),
                                col_manual = c("blue4", "orange3")) {

  # validate args
  if (!is.character(df[[chr]])) {
    df[[chr]] <- as.factor(as.integer(df[[chr]]))
  }

  assertthat::assert_that(all(df[[chr]] %in% 1:22),
                          msg = "Values in `df[[chr]]` should be integers between 1 and 22")

  plot_manhattan_basic(df = df,
                       chr = chr,
                       log10_p = log10_p,
                       text = text,
                       order_idx = order_idx,
                       geom_point_args = geom_point_args,
                       labs_args = labs_args) +
    ggplot2::scale_color_manual(values = rep(col_manual,
                                             dplyr::n_distinct(df[[chr]])))

}


# prepare dummy data for basic Manhattan plot
# dummy_manhattan_df(type = "phewas") %>%
#   arrange(CHR,
#           P) %>%
#   mutate("ORDER_IDX" = 1:n()) %>%
#   mutate(P = -log10(P)) %>%
#   mutate(CHR = as.character(CHR)) %>%
#   plot_manhattan_basic() +
#   theme(axis.text.x = element_text(angle = 325, hjust=0, vjust=0))
