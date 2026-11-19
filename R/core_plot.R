# Core layer: plotting.
#
# Every plot function takes the ggtheme jamovi passes in (so the base look follows the user's
# Plots-tab theme) and adds its own colour scales last, because jamovi's theme object carries
# a colour scale of its own. `palette` is "color" (colour-blind-safe Okabe-Ito) or "gray"
# (shades, shapes and line types, for black-and-white printing).

dc_method_colors <- c(
    "TWFE"               = "#D55E00",
    "Callaway-Sant'Anna" = "#0072B2",
    "Sun-Abraham"        = "#009E73",
    "Borusyak-Jaravel-Spiess" = "#CC79A7",
    "Gardner (did2s)"    = "#E69F00",
    "Stacked DiD"        = "#555555",
    "Triple difference (DDD)" = "#D55E00",
    "DiD, eligible only"      = "#0072B2",
    "DiD, ineligible only"    = "#009E73")

dc_gray_colors  <- c("#000000", "#4D4D4D", "#808080", "#999999", "#262626", "#666666", "#000000", "#4D4D4D", "#808080")
dc_gray_shapes  <- c(16, 17, 15, 18, 4, 8, 16, 17, 15)
dc_gray_lines   <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "solid", "dashed", "dotted")

# colours, shapes and line types per estimator for a palette
dc_pal <- function(palette = c("color", "gray")) {
    palette <- match.arg(palette)
    nm <- names(dc_method_colors)
    gray <- palette == "gray"
    list(colors = if (gray) stats::setNames(dc_gray_colors, nm) else dc_method_colors,
         shapes = stats::setNames(if (gray) dc_gray_shapes else rep(16, length(nm)), nm),
         lines  = stats::setNames(if (gray) dc_gray_lines else rep("solid", length(nm)), nm),
         gray = gray)
}

# `dyn`: long data frame with columns method, e, est, lo, hi (lo/hi NA for the reference period).
# `ci`: "bars" (error bars, dodged by `dodge`) or "ribbon" (lines with shaded bands).
dc_plot_event <- function(dyn, ggtheme = NULL, palette = "color", ci = c("bars", "ribbon"), dodge = 0.5) {
    ci <- match.arg(ci)
    pal <- dc_pal(palette)
    dyn$method <- factor(dyn$method, levels = names(dc_method_colors))
    p <- ggplot2::ggplot(dyn, ggplot2::aes(x = e, y = est, colour = method, shape = method)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = -0.5, colour = "grey60", linewidth = 0.4, linetype = "dashed")
    if (ci == "ribbon") {
        dyn$lo[is.na(dyn$lo)] <- dyn$est[is.na(dyn$lo)]     # keep the band continuous through e = -1
        dyn$hi[is.na(dyn$hi)] <- dyn$est[is.na(dyn$hi)]
        p <- ggplot2::ggplot(dyn, ggplot2::aes(x = e, y = est, colour = method, shape = method)) +
            ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
            ggplot2::geom_vline(xintercept = -0.5, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = method), alpha = 0.15,
                                 colour = NA, na.rm = TRUE, show.legend = FALSE) +
            ggplot2::geom_line(ggplot2::aes(linetype = method), linewidth = 0.7, na.rm = TRUE) +
            ggplot2::geom_point(size = 2)
    } else {
        pd <- ggplot2::position_dodge(width = dodge)
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.25,
                                        position = pd, na.rm = TRUE) +
            ggplot2::geom_point(position = pd, size = 2)
    }
    p <- p + ggplot2::scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]))) +
        ggplot2::labs(x = "Event time (periods relative to first treatment)",
                      y = "Estimated effect", caption = "Reference period: e = -1")
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p <- p + ggplot2::scale_colour_manual(values = pal$colors, drop = TRUE, name = NULL) +
        ggplot2::scale_shape_manual(values = pal$shapes, drop = TRUE, name = NULL) +
        ggplot2::theme(legend.position = "bottom")
    if (ci == "ribbon")
        p <- p + ggplot2::scale_fill_manual(values = pal$colors, guide = "none") +
            ggplot2::scale_linetype_manual(values = pal$lines, drop = TRUE, name = NULL)
    p
}

dc_bacon_colors <- c(
    "Treated vs Untreated"     = "#009E73",
    "Earlier vs Later Treated" = "#0072B2",
    "Later vs Earlier Treated" = "#D55E00")
dc_bacon_grays <- c(
    "Treated vs Untreated"     = "#000000",
    "Earlier vs Later Treated" = "#808080",
    "Later vs Earlier Treated" = "#BFBFBF")

# `detail`: bacon 2x2 comparisons (weight, estimate, type); `twfe`: static TWFE estimate
dc_plot_bacon <- function(detail, twfe, ggtheme = NULL, palette = "color") {
    cols <- if (identical(palette, "gray")) dc_bacon_grays else dc_bacon_colors
    detail$type <- factor(detail$type, levels = names(dc_bacon_colors))
    p <- ggplot2::ggplot(detail, ggplot2::aes(x = weight, y = estimate, colour = type, shape = type)) +
        ggplot2::geom_hline(yintercept = twfe, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(x = "Weight in the TWFE estimate", y = "2x2 DiD estimate",
                      caption = "Dashed line: static TWFE estimate")
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p + ggplot2::scale_colour_manual(values = cols, drop = TRUE, name = NULL) +
        ggplot2::scale_shape_manual(values = c(16, 17, 15), drop = TRUE, name = NULL) +
        ggplot2::theme(legend.position = "bottom")
}

# `df`: columns label (factor order = plot order), lb, ub, original (logical)
dc_plot_honest <- function(df, ggtheme = NULL, palette = "color") {
    cols <- if (identical(palette, "gray")) c(`TRUE` = "#808080", `FALSE` = "#000000")
            else c(`TRUE` = "#555555", `FALSE` = "#0072B2")
    df$label <- factor(df$label, levels = unique(df$label))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = label, ymin = lb, ymax = ub, colour = original)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
        ggplot2::geom_errorbar(width = 0.25, linewidth = 0.8) +
        ggplot2::labs(x = NULL, y = "Robust confidence interval",
                      caption = "Mbar: post-treatment violations at most Mbar times the largest pre-treatment violation")
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p + ggplot2::scale_colour_manual(values = cols, guide = "none")
}

# `df`: columns method, est, lo, hi
dc_plot_forest <- function(df, ggtheme = NULL, palette = "color") {
    pal <- dc_pal(palette)
    df$method <- factor(df$method, levels = rev(names(dc_method_colors)[names(dc_method_colors) %in% df$method]))
    p <- ggplot2::ggplot(df, ggplot2::aes(y = method, x = est, xmin = lo, xmax = hi, colour = method, shape = method)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
        ggplot2::geom_errorbar(orientation = "y", width = 0.25, linewidth = 0.8) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(x = "Overall ATT with confidence interval", y = NULL)
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p + ggplot2::scale_colour_manual(values = pal$colors, guide = "none") +
        ggplot2::scale_shape_manual(values = pal$shapes, guide = "none")
}

# `gt`: columns g (cohort), t (period), att, se -- Callaway & Sant'Anna group-time effects
dc_plot_gt <- function(gt, ggtheme = NULL, palette = "color") {
    gray <- identical(palette, "gray")
    gt <- gt[is.finite(gt$att), ]
    gt$g <- factor(gt$g, levels = sort(unique(gt$g)))
    lim <- max(abs(gt$att), na.rm = TRUE)
    gt$fillv <- if (gray) abs(gt$att) else gt$att        # gray: shade = size, the sign stays in the printed number
    p <- ggplot2::ggplot(gt, ggplot2::aes(x = t, y = g, fill = fillv)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]))) +
        ggplot2::labs(x = "Period", y = "Cohort (first treated period)",
                      fill = if (gray) "|ATT(g,t)|" else "ATT(g,t)")
    if (nrow(gt) <= 150)
        p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", att),
                                                 colour = if (gray) fillv > 0.6 * lim else FALSE), size = 2.8, show.legend = FALSE) +
            ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"))
    if (!is.null(ggtheme))
        p <- p + ggtheme
    fillscale <- if (gray) ggplot2::scale_fill_gradient(low = "white", high = "#262626", limits = c(0, lim))
                 else ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim))
    p + fillscale + ggplot2::theme(legend.position = "right")
}

# `means`: group (0/1), post (0/1), mean, n from dc_means22(); `labels`: c(before, after)
dc_plot_22 <- function(means, labels, ggtheme = NULL, palette = "color") {
    gray <- identical(palette, "gray")
    means$period <- factor(labels[means$post + 1], levels = labels)
    means$grp <- factor(ifelse(means$group == 1, "Treated", "Comparison"), levels = c("Comparison", "Treated"))
    c0 <- means[means$group == 0, ]; t1 <- means[means$group == 1, ]
    change_c <- c0$mean[c0$post == 1] - c0$mean[c0$post == 0]
    cf <- data.frame(period = factor(labels, levels = labels),
                     mean = c(t1$mean[t1$post == 0], t1$mean[t1$post == 0] + change_c))
    p <- ggplot2::ggplot(means, ggplot2::aes(x = period, y = mean, colour = grp, group = grp, shape = grp)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_line(data = cf, ggplot2::aes(x = period, y = mean, group = 1), inherit.aes = FALSE,
                           colour = if (gray) "#A6A6A6" else "#CC79A7", linetype = "dashed", linewidth = 0.8) +
        ggplot2::labs(x = NULL, y = "Mean outcome",
                      caption = "Dashed: treated group's path if it had changed like the comparison group")
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p + ggplot2::scale_colour_manual(values = if (gray) c(Comparison = "#808080", Treated = "#000000")
                                                  else c(Comparison = "#0072B2", Treated = "#D55E00"), name = NULL) +
        ggplot2::scale_shape_manual(values = c(Comparison = 17, Treated = 16), name = NULL) +
        ggplot2::theme(legend.position = "bottom")
}

# `st`: columns row (unit position, sorted by cohort), t, status
dc_plot_status <- function(st, ggtheme = NULL, palette = "color") {
    st$status <- factor(st$status, levels = c("Untreated", "Treated", "Not observed"))
    p <- ggplot2::ggplot(st, ggplot2::aes(x = t, y = row, fill = status)) +
        ggplot2::geom_raster() +
        ggplot2::scale_y_reverse(expand = c(0, 0)) +
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = function(l) seq(ceiling(l[1]), floor(l[2]))) +
        ggplot2::labs(x = "Period", y = "Units (sorted by cohort)", fill = NULL)
    if (!is.null(ggtheme))
        p <- p + ggtheme
    cols <- if (identical(palette, "gray")) c(Untreated = "#D9D9D9", Treated = "#333333", `Not observed` = "#FFFFFF")
            else c(Untreated = "#C6DBEF", Treated = "#D55E00", `Not observed` = "#FFFFFF")
    p + ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
        ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
}

# `curve`: dose, att, att_se, acrt, acrt_se; `z`: normal critical value
dc_plot_dose <- function(curve, z, ggtheme = NULL, palette = "color") {
    gray <- identical(palette, "gray")
    long <- rbind(
        data.frame(panel = "Effect of dose d: ATT(d)", dose = curve$dose, est = curve$att,
                   lo = curve$att - z * curve$att_se, hi = curve$att + z * curve$att_se),
        data.frame(panel = "Marginal effect of dose: ACRT(d)", dose = curve$dose, est = curve$acrt,
                   lo = curve$acrt - z * curve$acrt_se, hi = curve$acrt + z * curve$acrt_se))
    long$panel <- factor(long$panel, levels = unique(long$panel))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = dose, y = est)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = if (gray) "#808080" else "#0072B2", alpha = 0.2) +
        ggplot2::geom_line(colour = if (gray) "#000000" else "#0072B2", linewidth = 0.9) +
        ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) +
        ggplot2::labs(x = "Dose", y = "Estimated effect")
    if (!is.null(ggtheme))
        p <- p + ggtheme
    p
}
