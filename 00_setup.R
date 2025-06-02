suppressMessages({
    library(here)
    library(rlang)
    library(birdie)
    library(tidyverse)
    library(scales)
    library(wacolors)
    library(patchwork)
    library(geomtextpath)
})

races = c(white="White", black="Black", hisp="Hispanic", asian="Asian",
          aian="Native", other="Other")

PAL_R = wacolors$rainier
names(PAL_R) = NULL

theme_paper = function() {
    theme_bw(base_family="Times", base_size=10) +
        theme(plot.background=element_blank())
}



tidy_thresh = function(r_probs, x, nm=deparse(substitute(x)), prefix="pr_", se_boot=0) {
    x = as.factor(x)
    N = length(x)

    r_names =  substring(colnames(r_probs), nchar(prefix)+1L)
    stopifnot(nrow(r_probs) == N)

    r_est = predict(r_probs)
    m = prop.table(table(x, r_est), 2)

    out = tibble(X = rep(levels(x), ncol(m)),
                 race = rep(r_names, each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = nm

    if (se_boot > 0) {
        ests = matrix(nrow=length(m), ncol=se_boot)
        for (i in seq_len(se_boot)) {
            wts = birdie:::rdirichlet(rep(1, N), N)
            ests[, i] = c(prop.table(xtabs(wts ~ x + r_est), 2))
        }
        out$se = sqrt(diag(cov(t(ests))))
    }

    out
}

tidy_ols = function(r_probs, x, grp, nm=deparse(substitute(x)), prefix="pr_", se_boot=0) {
    x = as.factor(x)
    N = length(x)

    r_names = substring(colnames(r_probs), nchar(prefix)+1L)
    r_probs = as.matrix(select(r_probs, starts_with(prefix)))
    stopifnot(nrow(r_probs) == N)

    grp_locs = vctrs::vec_group_loc(grp)$loc
    grp_ids = vctrs::vec_group_id(grp)
    p_gr = prop.table(xtabs(r_probs ~ grp_ids), 2)
    ffill = function(x) ifelse(is.na(x), 0, x)
    m = map(levels(x), function(l) {
        grp_m = imap(grp_locs, function(idx, i) {
            # .lm.fit(r_probs[idx, , drop=FALSE], x[idx] == l)$coefficients
            m = lm(x[idx] == l ~ 0 + r_probs[idx, , drop=FALSE])
            cbind(ffill(coef(m)),
                  ffill(sqrt(diag(vcov(m)))))
        })
        grp_ests = do.call(rbind, lapply(grp_m, function(x) x[, 1]))
        grp_var = do.call(rbind, lapply(grp_m, function(x) x[, 2]^2))
        cbind(
            est = colSums(grp_ests * p_gr),
            var = colSums(grp_var * p_gr^2)
        )
    })

    m_var = do.call(rbind, lapply(m, function(x) x[, 2]))
    m = do.call(rbind, lapply(m, function(x) x[, 1]))

    out = tibble(
        X = rep(levels(x), ncol(m)),
        race = rep(r_names, each=nrow(m)),
        estimate = as.numeric(m),
        se = sqrt(as.numeric(m_var)),
    )
    names(out)[1] = nm

    if (se_boot > 0) {
        ests = matrix(nrow=length(m), ncol=se_boot)
        for (i in seq_len(se_boot)) {
            wts = birdie:::rdirichlet(rep(1, N), N)
            ests[, i] = map(levels(x), function(l) {
                grp_ests = imap(grp_locs, function(idx, i) {
                    lm.wfit(r_probs[idx, , drop=FALSE], x[idx] == l, wts[idx])$coefficients
                }) |>
                    do.call(rbind, args=_)
                colSums(grp_ests * p_gr)
            }) |>
                do.call(rbind, args=_) |>
                c()
        }
        out$se = sqrt(diag(cov(t(ests))))
    }

    out
}

tidy_true = function(race, x, nm=deparse(substitute(x))) {
    x = as.factor(x)
    N = length(x)

    r_names = levels(race)
    stopifnot(length(race) == N)

    m = prop.table(table(x, race), 2)

    out = tibble(X = rep(levels(x), ncol(m)),
                 race = rep(r_names, each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = nm
    out
}

plot_calib = function(r_probs, race, group="white", bins=16) {
    y_fit = r_probs[[str_c("pr_", group)]]
    tibble(fitted = cut(y_fit, bins),
           y_res = race) %>%
        group_by(fitted) %>%
        summarize(mean = mean(y_res == group),
                  n = n()) %>%
        mutate(fitted = parse_number(as.character(fitted))) %>%
        ggplot(aes(x=fitted, y=mean, size=n)) +
        geom_point() +
        labs(title=str_c("Race: ", group)) +
        scale_size_area(guide="none")  +
        theme_bw()
}

# Helpful function adapted from: https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
#' Calculate AUC
#'
#' @param x the predictor
#' @param y a binary indicator
#'
#' @returns the scalar AUC
fastAUC <- function(x, y) {
    x1 = x[y==1]; n1 = length(x1);
    x2 = x[y==0]; n2 = length(x2);
    r = rank(c(x1, x2))
    exp(log(sum(r[1:n1]) - n1*(n1 + 1)/2) - log(n1) - log(n2))
}

calc_roc <- function(d_pr) {
    imap(p_r, function(x, r) {
        fastAUC(d_pr[[str_c("pr_", r)]], d$race == r)
    }) |>
        as_tibble()
}

log_score <- function(x) {
    pr_act = as.matrix(x)[cbind(1:nrow(x), as.integer(d$race))]
    pr_act[pr_act == 0] = 1e-6
    mean(log(pr_act))
}
