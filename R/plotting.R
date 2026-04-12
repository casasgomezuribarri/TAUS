#' Histogram of times to event, coloured by exit reason and faceted by a variable of interest
#'
#' Plots a faceted histogram of times to event (death or censoring), coloured by exit reason.
#' Useful for visually inspecting the distribution of survival times and censoring
#' across subgroups.
#'
#' @param data A data frame containing the survival data.
#' @param time_var Name of the time variable in \code{data} (string).
#' @param event_var Name of the event variable in \code{data} (string). Expected values are \code{"Dead"} and
#'   \code{"Censored"}.
#' @param var Name of the variable to facet by (string).
#' @param nbins Number of histogram bins. Defaults to \code{max(data[[time_var]])}.
#' @param ncol Number of columns in the facet layout. Defaults to 1.
#'
#' @return A \code{ggplot2} object, returned invisibly. The plot is also
#'   printed as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f_t_census(data = my_data, time_var = "time", event_var = "status", var = "group")
#' }
#'
f_t_census <- function(data, # the dataset
                       time_var, # name of the time variable (string)
                       event_var, # name of the event variable (string) - values will show in the legend
                       var, # variable to facet by (string)
                       nbins = max(data[[time_var]]), # max(t) by default
                       ncol = 1) {
    ft <- ggplot(data, aes(x = !!sym(time_var), fill = !!sym(event_var))) +
        geom_histogram(position = "dodge", bins = nbins, alpha = 1) +
        scale_fill_manual(values = c("Dead" = "red", "Censored" = "green")) +
        labs(
            title = "Histogram of survival times",
            x = "Time (days)",
            y = "Frequency",
            fill = "Exit reason"
        ) +
        facet_wrap(as.formula(paste("~", var)), ncol = ncol) + # parse var
        # theme_minimal() +
        theme(
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 13),
            strip.text = element_text(size = 16),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 13)
        )
    print(ft) # display the plot
    invisible(ft) # return without printing - so it can be saved in a variable
}


#' Plot conditional survival estimates
#'
#' Visualises the output of \code{cond_surv_mat()} in one of four modes
#' depending on the combination of \code{tau_value} and \code{collapse_time}:
#'
#' \itemize{
#'   \item \code{tau_value} specified, \code{collapse_time = FALSE}: plots
#'     \eqn{P(T > \tau | T > t)} against \eqn{t} for a fixed \eqn{\tau}.
#'   \item \code{tau_value} specified, \code{collapse_time = TRUE}: plots
#'     \eqn{P(T > \tau)} of a randomly-sampled individual (i.e. age unknown) for each group.
#'   \item \code{tau_value = NULL}, \code{collapse_time = FALSE}: plots a
#'     heatmap of \eqn{P(T > \tau | T > t)} over all \eqn{t} and \eqn{\tau}.
#'   \item \code{tau_value = NULL}, \code{collapse_time = TRUE}: plots
#'     \eqn{P(T > \tau)} of a randomly-sampled individual (i.e. age unknown) for all values
#'     of \eqn{\tau} for each group. Note that this is NOT a survival curve, despite the initial resemblance.
#' }
#'
#' @param cond_surv Output of \code{cond_surv_mat()}.
#' @param var_value Optional. Filter to a single group level (string). Must be
#'   one of the values of the grouping variable used in \code{cond_surv_mat()}.
#' @param tau_value Optional. Target age(s) to plot. Either a single numeric
#'   value (i.e. same \eqn{\tau} fro all groups), or a vector. If vector, it must
#'   be of the same length as there are levels in the grouping variable used in
#'   \code{cond_surv_mat()}}
#' @param color_var Optional. Name of the variable to colour lines or points
#'   by (string). Defaults to the grouping variable if \code{NULL}.
#' @param linetype_var Optional. Name of the variable to map to line type
#'   (string). If \code{NULL}, all lines are solid.
#' @param markershape_var Optional. Name of the variable to map to point shape
#'   (string). Only used when \code{collapse_time = TRUE}. If \code{NULL},
#'   no shape mapping is applied.
#' @param ncol_hm Number of columns in the facet grid of the heatmap. Only used when
#'   \code{tau_value = NULL} and \code{collapse_time = FALSE} (i.e. heatmap mode).
#'   Defaults to 1.
#' @param collapse_time If \code{TRUE}, plots \eqn{P(T > \tau)} marginalised
#'   over \eqn{t} using \eqn{A(t)} as weights, instead of the full
#'   time-dependent conditional survival curve. Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot2} object, returned invisibly. The plot is also
#'   printed as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Heatmap of all t and tau combinations
#' cond_surv_plot(cond_surv)
#'
#' # Conditional survival curve for a fixed tau (conditional on age, plotted in the x-axis)
#' cond_surv_plot(cond_surv, tau_value = 15)
#'
#' # Target-Age Unified Survival probability for a fixed tau, by group
#' cond_surv_plot(cond_surv, tau_value = 15, collapse_time = TRUE)
#'
#' # Target-Age Unified Survival curve across all tau values
#' cond_surv_plot(cond_surv, collapse_time = TRUE)
#' }
cond_surv_plot <- function(cond_surv, # the output of cond_surv_mat
                           var_value = NULL, # optional filter - must be one of the values of var
                           tau_value = NULL, # optional filter - must be a valid value within the range taus
                           color_var = NULL, # variable to color by (if NULL, use var)
                           linetype_var = NULL, # variable to linetype by (if NULL, no linetype)
                           markershape_var = NULL, # variable to use for marker shape (if NULL, no shape)
                           ncol_hm = 1, # number of columns if doing a heatmap (no filters)
                           collapse_time = FALSE # if TRUE, plot P(T>tau) instead of P(T>tau|T>t) for all t - collapsing t using A(t) as weights
) {
    # aggregate if appropriate
    if (cond_surv$args$aggregate_by_cats) {
        # compute maximum death observed in each group (levels in var)
        max_values <- cond_surv$args$data %>%
            filter(!!sym(cond_surv$args$event_var) == 1) %>% # only deaths
            group_by(!!sym(cond_surv$args$var), !!!syms(cond_surv$args$cats)) %>% # for each group of interest
            summarise(max_t = max(!!sym(cond_surv$args$time_var), na.rm = TRUE)) %>% # max_t is the latest death observed
            tidyr::unite("unique_label", all_of(c(cond_surv$args$var, cond_surv$args$cats)),
                sep = "_", remove = FALSE # in the future this could be made customisable for better legends
            )
    } else {
        # compute maximum death observed in each group (levels combinations in var and cats)
        max_values <- cond_surv$args$data %>%
            filter(!!(sym(cond_surv$args$event_var)) == 1) %>% # only deaths
            group_by(!!sym(cond_surv$args$var)) %>% # for each group of interest
            summarise(max_t = max(!!sym(cond_surv$args$time_var), na.rm = TRUE)) %>% # max_t is the latest death observed
            mutate(
                unique_label = paste0(!!sym(cond_surv$args$var)) # in the future this could be made customisable for better legends
            )
    }
    # useful for later
    max_time <- cond_surv$args$data %>%
        summarise(max_value = max(!!sym(cond_surv$args$time_var[[1]]), na.rm = TRUE))

    # if either filter is specified
    if (!is.null(tau_value) && collapse_time == FALSE) { # if tau_value is specified and collapse_time is FALSE the plot is P(T>tau|T>t) against t

        filtered_cond_surv <- cond_surv$cond_surv %>%
            filter(tau == tau_value)

        # for later
        title <- bquote("Time-dependant probability of outliving " ~ .(tau_value))

        if (!is.null(var_value)) { # if var_value is specified, filter by it too
            filtered_cond_surv <- filtered_cond_surv %>%
                filter(unique_label == var_value)

            # update
            title <- bquote("Time-dependant probability of outliving " ~ .(tau_value))
        }

        if (is.null(color_var)) { # if color_var isn't specified, use unique_var
            color_var <- "unique_label"
        }

        # plot the curve. If linetype_var is NULL:
        if (is.null(linetype_var)) { # if linetype_var isn't specified, all are solid
            plot <- ggplot(filtered_cond_surv, aes(
                x = time,
                y = Stau_St,
                group = unique_label,
                color = !!sym(color_var)
            ))
        } else { # but if it is, then use it
            plot <- ggplot(filtered_cond_surv, aes(
                x = time,
                y = Stau_St,
                group = unique_label,
                color = !!sym(color_var),
                linetype = !!sym(linetype_var)
            ))
        }

        # plot the curve
        plot <- plot +
            geom_line() +
            geom_ribbon(aes(ymin = Stau_St_lo, ymax = Stau_St_up, group = unique_label, fill = !!sym(color_var)), alpha = 0.2) +
            labs(x = "Time", y = bquote(O[tau == .(tau_value)](t)), title = title) +
            # scale_x_continuous(breaks = 1:30) +
            ylim(0, 1) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 18),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 13),
                strip.text = element_text(size = 16),
                legend.text = element_text(size = 13),
                legend.title = element_blank(),
                # legend.position = "bottom"
            )

        print(plot) # display the plot
        invisible(plot) # return without printing - so it can be saved in a variable
    } else if (!is.null(tau_value) && collapse_time == TRUE) { # if tau_value is specified and collapse_time is TRUE the plot is P(T>tau) for each var level

        # if var_values is NULL, tau_value is a number and is the same for all values of the variable of interest
        if (var_value %>% is.null()) {
            # tau_value must be a number
            if (!is.numeric(tau_value) || length(tau_value) != 1) {
                stop("tau_value must be a single numeric value or a vector of the same length as var_value")
            }
            # and this number can be in the plot title
            title <- bquote("Probability of outliving " ~ .(tau_value) ~ " when age is unknown")

            # filter cond_surv
            filtered_cond_surv <- cond_surv$cond_surv %>%
                filter(tau == tau_value)

            if (is.null(color_var)) { # if color_var isn't specified, use unique_var
                color_var <- "unique_label"
            }

            if (is.null(markershape_var)) { # if markershape_var isn't specified, no shape
                shape_mapping <- NULL
            } else {
                shape_mapping <- aes(shape = !!(sym(markershape_var)))
            }

            # plot it
            plot <- ggplot(filtered_cond_surv, aes(x = unique_label, y = O_tau, color = !!(sym(color_var)))) +
                geom_point(shape_mapping, position = position_dodge(width = 0.5), size = 3) +
                geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
                labs(x = "Group", y = bquote(O[tau == .(tau_value)]), title = title) +
                ylim(0, 1) +
                theme_minimal() +
                theme(
                    plot.title = element_text(size = 18),
                    axis.title = element_text(size = 16),
                    axis.text.x = element_text(angle = 13, hjust = 1),
                    axis.text = element_text(size = 13),
                    strip.text = element_text(size = 16),
                    strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
                    legend.text = element_text(size = 13),
                    legend.title = element_blank(),
                    # legend.position = "bottom"
                )

            print(plot) # display the plot
            invisible(plot) # return without printing - so it can be saved in a variable

            # if var_value has been specified, tau_values should be a vector of taus for each level, in the same order
        } else {
            # tau_value must be the same length
            if (length(tau_value) != length(var_value)) {
                stop("tau_value must be a single numeric value or a vector of the same length as var_value")
            }
            # and the title is a bit more generic
            title <- bquote("Probability of outliving " * tau * " when age is unknown")

            # first make a tibble with the pairs of tau and var
            pairs <- tibble(
                tau = tau_value,
                unique_label = var_value # rlang's !!sym()
            )

            # now filter cond_surv
            filtered_cond_surv <- cond_surv$cond_surv %>%
                semi_join(pairs, by = c("tau", "unique_label")) # nice and neat

            # create a variable that joiins the group and the tau choice for labelling, with proper handling of characters for downstream printing of the greek symbol tau
            filtered_cond_surv <- filtered_cond_surv %>%
                mutate(
                    label = paste0("plain('", unique_label, "') ~ (tau==", tau, ")")
                )

            if (is.null(color_var)) { # if color_var isn't specified, use unique_var
                color_var <- "unique_label"
            }

            if (is.null(markershape_var)) { # if markershape_var isn't specified, no shape
                shape_mapping <- NULL
            } else {
                shape_mapping <- aes(shape = !!(sym(markershape_var)))
            }

            # plot it
            plot <- ggplot(filtered_cond_surv, aes(x = label, y = O_tau, color = !!(sym(color_var)))) +
                geom_point(shape_mapping, position = position_dodge(width = 0.5), size = 3) +
                geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
                labs(x = "Group", y = bquote(O[tau]), title = title) +
                scale_x_discrete(labels = function(x) parse(text = x)) + # to print the tau symbol, not the "tau" string
                ylim(0, 1) +
                theme_minimal() +
                theme(
                    plot.title = element_text(size = 18),
                    axis.title = element_text(size = 16),
                    axis.text.x = element_text(angle = 13, hjust = 1),
                    axis.text = element_text(size = 13),
                    strip.text = element_text(size = 16),
                    strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
                    legend.text = element_text(size = 13),
                    legend.title = element_blank(),
                    # legend.position = "bottom"
                )

            print(plot) # display the plot
            invisible(plot) # return without printing - so it can be saved in a variable
        }
    } else if (is.null(tau_value) && collapse_time == FALSE) { # if tau_value is not specified and collapse_time is FALSE, the plot is a heatmap.

        filtered_cond_surv <- cond_surv$cond_surv
        title <- bquote("Conditional survival per group")

        if (!is.null(var_value)) { # if var_value is specified, filter by it too
            filtered_cond_surv <- filtered_cond_surv %>%
                filter(unique_label == var_value)

            # update title
            title <- bquote("Conditional survival. Group: " ~ filtered_cond_surv$unique_label)
        }
        n_groups <- length(unique(filtered_cond_surv$unique_label)) # number of groups to plot (useful later)
        n_rows <- ceiling(n_groups / ncol_hm) # number of rows in the grid, useful later

        # heatmap
        plot <- ggplot(filtered_cond_surv, aes(x = time, y = tau, fill = Stau_St)) +
            geom_tile(color = NA) +
            geom_vline(data = max_values, aes(xintercept = max_t), color = "#00f891", linetype = "solid", linewidth = 1) + # adds a vertical line at the time of the last death observed
            scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), na.value = "grey") +
            labs(
                title = title,
                x = "Age (t)",
                y = bquote("Age to outlive (" * tau * ")"),
                fill = bquote("P(T > " * tau * " | T > t)")
            ) +
            facet_wrap(~unique_label, scales = "free", ncol = ncol_hm) +
            theme_minimal() +
            theme(
                aspect.ratio = 1,
                plot.title = element_blank(),
                axis.title = element_text(size = max(6, 11 - 1.5 * n_rows) + 1), # depends on number of groups
                axis.text = element_text(size = max(6, 11 - 1.5 * n_rows)), # depends on number of groups
                strip.text = element_text(size = max(6, 11 - 1.5 * n_rows)), # depends on number of groups
                legend.text = element_text(size = 10),
                legend.title = element_blank(),
            ) +
            guides(
                fill = guide_colorbar(
                    barwidth = 0.7,
                    barheight = n_rows * 4, # depends on number of rows
                    title.hjust = 0.5
                )
            )

        print(plot) # display the plot
        invisible(plot) # return without printing - so it can be saved in a variable
    } else if (is.null(tau_value) && collapse_time == TRUE) { # if tau_value is not specified and collapse_time is TRUE, the plot is P(T>tau) against tau

        filtered_cond_surv <- cond_surv$cond_surv
        title <- bquote("Probability of outliving " * tau * " when age is unknown")

        if (!is.null(var_value)) { # if var_value is specified, filter by it too
            filtered_cond_surv <- filtered_cond_surv %>%
                filter(unique_label == var_value)

            # update title
            title <- bquote("P(T>" ~ tau ~ ") when age is unknown. Group: " ~ filtered_cond_surv$unique_label)
        }

        if (is.null(color_var)) { # if color_var isn't specified, use unique_var
            color_var <- "unique_label"
        }

        # plot the curve. If linetype_var is NULL:
        if (is.null(linetype_var)) { # if linetype_var isn't specified, all are solid
            plot <- ggplot(filtered_cond_surv, aes(
                x = tau,
                y = O_tau,
                group = unique_label,
                color = !!sym(color_var)
            ))
        } else { # but if it is, then use it
            plot <- ggplot(filtered_cond_surv, aes(
                x = tau,
                y = O_tau,
                group = unique_label,
                color = !!sym(color_var),
                linetype = !!sym(linetype_var)
            ))
        }

        # plot the curve
        plot <- plot +
            geom_line() +
            geom_ribbon(aes(ymin = O_tau_lo, ymax = O_tau_up, group = unique_label, fill = !!sym(color_var)), alpha = 0.2) +
            labs(x = expression(tau), y = expression(O(tau)), title = title) +
            xlim(0, max_time[[1]]) +
            ylim(0, 1) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 18),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 13),
                strip.text = element_text(size = 16),
                strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
                legend.text = element_text(size = 13),
                legend.title = element_blank(),
                # legend.position = "bottom"
            )

        print(plot) # display the plot
        invisible(plot) # return without printing - so it can be saved in a variable
    }
}

#' Kaplan-Meier survival curve
#'
#' Fits and plots a Kaplan-Meier survival curve using \code{survfit2()} and
#' \code{ggsurvfit()}. Supports grouping, confidence intervals, custom colours
#' and linetypes, and optional median survival annotations.
#'
#' @param data A data frame containing the survival data.
#' @param title Plot title (string).
#' @param colors Optional. Named or unnamed character vector of colours for
#'   each group level. Required if \code{median_surv = TRUE}.
#' @param linetypes Optional. Named or unnamed character vector of linetypes
#'   for each group level. If \code{NULL}, all lines are solid.
#' @param title_size Font size for the plot title. Defaults to 20.
#' @param conf_int If \code{TRUE}, adds a shaded confidence interval around
#'   each survival curve. Defaults to \code{TRUE}.
#' @param group Name of the grouping variable (string). Use \code{NA} for an
#'   overall survival curve with no grouping.
#' @param time_var Name of the time variable (string).
#' @param event_var Name of the event variable (string). Expected to be a
#'   binary variable where 1 indicates the event (death) and 0 censoring.
#' @param xmax Maximum value for the x-axis.
#' @param legend Legend position. Passed to \code{ggplot2::theme()}.
#'   Common values are \code{"none"}, \code{"bottom"}, \code{"right"}.
#'   Defaults to \code{"none"}.
#' @param median_surv If \code{TRUE}, adds median survival annotations and
#'   reference lines to the plot. Requires \code{colors} to be specified.
#'   Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot2} object, returned invisibly. The plot is also
#'   printed as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Overall survival curve, no grouping
#' create_survival_curve(
#'     data = my_data, title = "Overall survival",
#'     group = NA, time_var = "time", event_var = "status",
#'     xmax = 30, median_surv = FALSE
#' )
#'
#' # Grouped survival curve with median annotations
#' create_survival_curve(
#'     data = my_data, title = "Survival by group",
#'     group = "treatment", time_var = "time", event_var = "status",
#'     colors = c("red", "blue"), xmax = 30
#' )
#' }
create_survival_curve <- function(data,
                                  title,
                                  colors = NULL,
                                  linetypes = NULL,
                                  title_size = 20,
                                  conf_int = TRUE,
                                  group,
                                  time_var,
                                  event_var,
                                  xmax,
                                  legend = "none",
                                  median_surv = TRUE) {
    if (!(is.na(group))) {
        # convert group to a factor, define formula, fit KM model
        data[[group]] <- as.factor(data[[group]])
        surv_formula <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", group))
        fit <- survfit(surv_formula, data = data)

        # extract median survival for each level in data[[group]]
        surv_data <- summary(fit)$table
        annotations <- data.frame(
            group = rownames(surv_data),
            x = surv_data[, "median"]
        )
        labels <- levels(data[[group]]) # labels for the annotation
    } else {
        # convert group to a factor, define formula, fit KM model
        surv_formula <- as.formula("data$event ~ 1")
        fit <- survfit2(surv_formula)

        # extract median survival for each level in data[[group]]
        surv_data <- summary(fit)$table
        annotations <- data.frame(
            group = c("All"),
            x = surv_data["median"]
        )
        labels <- c("All")
    }

    use_linetype_aes <- !is.null(linetypes)

    # Plot the survival curves
    p <- survfit2(surv_formula, data = data) %>%
        ggsurvfit(type = "survival", linetype_aes = use_linetype_aes) +
        labs(
            x = "Time (days)",
            y = "Proportion alive"
        ) +
        xlim(0, xmax) + # Make sure all plots have the same axis limits
        ylim(0, 1) + # Make sure all plots have the same axis limits
        (if (!is.null(colors)) scale_fill_manual(values = colors) else NULL) + # Set fill colors based on the palette
        (if (!is.null(colors)) scale_color_manual(values = colors) else NULL) + # Set line colors based on the palette
        (if (!is.null(linetypes)) scale_linetype_manual(values = linetypes) else NULL) +
        theme(
            axis.text = element_text(size = 13), # Font size for axis ticks
            axis.title = element_text(size = 16), # Adjust font of labels
            legend.position = legend, # Try "bottom" or "none"
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 16),
            plot.margin = margin(8, 8, 8, 8), # Plot margins (t, r, b, l)
            plot.title = element_text(size = title_size, hjust = 0.5) # Title settings
        ) +
        ggtitle(title)

    if (conf_int == TRUE) p <- p + add_confidence_interval(alpha = 0.1)


    if (median_surv == TRUE) {
        if (is.na(colors)) {
            stop("Median survival annotations need color specifications")
        }

        # Add a single "Median survival" annotation
        p <- p + annotate("text",
            x = 0, y = 0.25,
            label = "Median survival:",
            size = 8, color = "black", fontface = "bold", hjust = 0
        )

        # add the median survival for each level, and the lines
        for (i in 1:nrow(annotations)) {
            # add horizontal line
            p <- p + geom_segment(
                x = 0, xend = annotations$x[i], y = 0.5, yend = 0.5,
                color = colors[i], linetype = "solid"
            )

            # add vertical line
            p <- p + geom_segment(
                x = annotations$x[i], xend = annotations$x[i], y = 0.5, yend = 0,
                color = colors[i], linetype = "solid"
            )

            # Add annotation text
            p <- p + annotate("text",
                x = 0, y = 0.25 - i * 0.05, # Adjust y position for each label
                label = paste(labels[i], ": ", round(annotations$x[i], 1)),
                size = 8, color = colors[i], hjust = 0
            )
        }
    }
    print(p)
    invisible(p)
}
