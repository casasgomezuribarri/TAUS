# i love this

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
