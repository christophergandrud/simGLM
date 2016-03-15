#' Simulate quantities of interest from generalised linear models
#'
#'
#' @examples
#' # Normal Linear Model example
#' library(car) # Contains data
#' m1 <- lm(prestige ~ education + type,
#'          data = Prestige)
#'
#' fitted_prestige <- expand.grid(education = 6:16, typewc = 1)
#'
#' sim_glm(obj = m1, newdata = fitted_prestige, x_var = 'education')
#'
#' fitted_prestige <- expand.grid(education = 6:16,
#'                                typewc = 0:1)
#' sim_glm(obj = m1, newdata = fitted_prestige, x_var = 'education',
#'         group_var = 'typewc')
#'
#' # Logistic Model example
#' URL <- 'http://www.ats.ucla.edu/stat/data/binary.csv'
#' Admission <- read.csv(URL)
#' Admission$rank <- as.factor(Admission$rank)
#'
#' m2 <- glm(admit ~ gre + gpa + rank,
#'           data = Admission, family = 'binomial')
#'
#' fitted_admit_1 <- with(Admission,
#'                        expand.grid(gre = seq(220, 800, by = 10),
#'                                    gpa = mean(gpa),
#'                                    rank2 = 0:1))
#'
#' sim_glm(obj = m2, newdata = fitted_admit_1, x_var = 'gre',
#'         group_var = 'rank2', model = 'logit')
#'
#' fitted_admit_2 <- expand.grid(gre = seq(220, 800, by = 10), gpa = c(2, 3.5),
#'                               rank4 = 1)
#' sim_glm(obj = m2, newdata = fitted_admit_2, model = 'logit',
#'         x_var = 'gre', group_var = 'gpa')
#'
#' @import ggplot2
#' @importFrom MASS mvrnorm
#' @importFrom dplyr %>% group_by summarise
#'
#' @export


sim_glm <- function(obj, newdata, x_var, group_var, n = 1000, model = 'lm',
                    col_pal) {

    require(MASS)


    if (!(model %in% c('lm', 'logit'))) stop('model must be either "lm" or "logit.',
                                             call. = FALSE)

    if (missing(x_var)) {
        if (missing(x_var) & ncol(newdata) == 1) {
            x_var <- names(newdata)
            message(sprintf('%s set for x_var.', names(newdata)))
        }
        if (missing(x_var)) stop("x_var must be specified to determine the simulation plot's x-axis.",
                                 call. = FALSE)
    }

    if (!missing(group_var) & length(unique(newdata[, group_var])) == 1) stop(
        'Your group_var only has one value, so there is no need to set group_var.',
        call. = FALSE
    )

    # Find point estimates and variance-covariance matrix
    coef_est <- coef(obj) %>% matrix
    vcov_est <- vcov(obj)

    # Draw simulations
    drawn <- mvrnorm(n = n, mu = coef_est, Sigma = vcov_est) %>% data.frame

    # Ensure fitted variables match
    drawn_names <- names(drawn)
    if (any(!(names(newdata) %in% drawn_names))) stop('All column names in newdata must match *coefficient names* in the fitted model.',
                                                 call. = FALSE)
    if (!(x_var %in% drawn_names)) stop('x_var must be a *coefficient name* in the fitted model',
                                        call. = FALSE)
    if (!missing(group_var)){
        if (!(group_var %in% drawn_names)) stop('group_var must be a *coefficient name* in the fitted model',
                                            call. = FALSE)
    }

    # Mark fitted values as distinct from the point estimates
    names(newdata) <- sprintf('%s_fitted_', names(newdata))

    # Merge simulated point estimates and fitted values
    drawn_fitted <- merge(drawn, newdata)

    # Find quantity of interest
    names_fitted <- names(drawn_fitted)
    drawn_fitted$qi_ <- 0
    for (i in drawn_names[-1]) {
        fitted_i <- sprintf('%s_fitted_', i)

        # Check if fitted variable exists
        if (fitted_i %in% names_fitted) {
            temp_comb <- drawn_fitted[, i] * drawn_fitted[, fitted_i]
        }
        else if (!(fitted_i %in% names_fitted)) {
            temp_comb <- 0
            message(sprintf('%s fitted at 0.', i))
        }

        drawn_fitted$qi_ <- drawn_fitted$qi_ + temp_comb
    }
    drawn_fitted$qi_ <- drawn_fitted[, 1] + drawn_fitted$qi_ # add to intercept

    if (model == 'lm') {
        qi_name <- 'Predicted Value of y\n'
    }
    if (model == 'logit') {
        # Find probabilities of y = 1
        drawn_fitted$qi_ <- exp(drawn_fitted$qi_) / (1 - exp(drawn_fitted$qi_))
        # Drop simulations outside of [0, 1]
        drawn_fitted <- subset(drawn_fitted, qi_ < 1)
        drawn_fitted <- subset(drawn_fitted, qi_ > 0)
        qi_name <- 'Pr(y = 1)\n'
    }

    # Find original distribution of x_var for rug plot --------
    original_data <- model.frame(obj)
    rug <- original_data[, x_var]
    if (missing(group_var)) {
        rug <- data.frame(xvar__ = rug, mean_sim = 1)
        names(rug) <- c('xvar__', 'mean_sim')
    }
    else if (!missing(group_var)) {
        rug <- data.frame(xvar__ = rug, mean_sim = 1,
                          group_var__ = as.factor(min(newdata[, sprintf('%s_fitted_',
                                                              group_var)])))
        names(rug) <- c('xvar__', 'mean_sim', 'group_var__')
    }

    # Find central intervals
    xvar_position <- match(sprintf('%s_fitted_', x_var), names_fitted)
    names(drawn_fitted)[xvar_position] <- 'xvar__'

    # Plot with no groups ----------
    if (missing(group_var)) {
        central <- drawn_fitted %>% group_by(xvar__) %>%
            summarise(mean_sim = mean(qi_),
                      lower_95 = quantile(qi_, probs = 0.025),
                      upper_95 = quantile(qi_, probs = 0.975),
                      lower_90 = quantile(qi_, probs = 0.05),
                      upper_90 = quantile(qi_, probs = 0.95)
            )

        # Plot
        out_plot <- ggplot(central, aes(xvar__, mean_sim)) +
            geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.05) +
            geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.05) +
            geom_line()
    }

    # Plot with groups ----------
    else if (!missing(group_var)) {
        group_position <- match(sprintf('%s_fitted_', group_var), names_fitted)
        if (is.na(group_position)) stop('group_var must be a variable in newdata.',
                                        call. = FALSE)
        names(drawn_fitted)[group_position] <- 'group_var__'

        central <- drawn_fitted %>% group_by(xvar__, group_var__) %>%
            summarise(mean_sim = mean(qi_),
                      lower_95 = quantile(qi_, probs = 0.025),
                      upper_95 = quantile(qi_, probs = 0.975),
                      lower_90 = quantile(qi_, probs = 0.05),
                      upper_90 = quantile(qi_, probs = 0.95)
            )
        central$group_var__ <- as.factor(central$group_var__)

        # Organise colour palette
        ## Modified Darjeeling palette from wesanderson package
        if (missing(col_pal)) col_pal <- c("#F98400", "#5BBCD6",
                                           "#FF0000", "#00A08A", "#F2AD00")
        n_levels <- length(unique(as.factor(central$group_var__)))
        if(length(col_pal) > length(n_levels)) {
            col_pal <- col_pal[1:n_levels]
        }
        else if (length(col_pal) < length(n_levels)) sprintf(
            'The colour palette needs %s values, but there are only %s',
            length(n_levels), length(col_pal)
        )

        # Plot
        out_plot <- ggplot(central, aes(xvar__, mean_sim,
                            group = group_var__,
                            fill = group_var__)) +
            geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.05) +
            geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.05) +
            geom_line(aes(colour = group_var__)) +
            scale_colour_manual(values = col_pal, name = group_var) +
            scale_fill_manual(values = col_pal, name = group_var)

    }
    out_plot + xlab(sprintf('\n%s', x_var)) + ylab(qi_name) +
                geom_rug(data = rug, aes(x = xvar__, y = mean_sim),
                         sides = 'b', alpha = 0.2) +
                theme_bw()
}



