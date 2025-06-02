# OLS: Finite sample performance, weak-instrument test

# OLS by size ----------
X = d$party
X_name = "party"

d_true = tidy_true(d$race, X, X_name)

ols_cty = function(idx) {
    tidy_ols(r_probs$county[idx, ], X[idx], d$GEOID_county[idx], X_name) |>
        suppressWarnings() |>
        rename(est = estimate) |>
        left_join(d_true, by=c("party", "race")) |>
        rename(est_true = estimate) |>
        mutate(level = "county", method = length(idx)) |>
        eval_fit_tv(p_r)
}

n = 1000 * round(10^seq(0, 3, length.out=10))
cells = tibble(
    n = n,
    n_per_cty = map_dbl(n, ~ .x / n_distinct(d$GEOID_county[1:.x]))
)
res = map(n, ~ ols_cty(1:.x), .progress = TRUE) |>
    bind_rows() |>
    rename(n = method) |>
    left_join(cells, by="n") |>
    mutate(race = fct_inorder(race))



race_labs = c(overall="Overall", races)
d_lab = filter(res, n==1e6) |>
    arrange(desc(tv)) |>
    mutate(lab = race_labs[as.character(race)],
           y = tv + c(0, 0.01, -0.005, 0.006, -0.002, -0.007, 0))

ggplot(res, aes(n_per_cty, tv, group=race)) +
    geom_line(aes(color=(race=="overall"), lwd=(race=="overall"))) +
    geom_point(size=1.0) +
    geom_text(aes(label=lab, y=y), data = d_lab, hjust=-0.2, size=3.0, family="Times") +
    scale_linewidth_manual(values=c(0.5, 1.0), guide="none") +
    scale_color_manual(values=c("#506064", "black"), guide="none") +
    scale_x_log10("Average observations per county", labels=scales::comma,
                  expand=expansion(mult=c(0.05, 0.15))) +
    scale_y_log10("Total variation distance") +
    theme_paper()
ggsave(here("paper/figures/nc_ols_n.pdf"), width=7, height=4.0)


# Variance inflation ------------
counties = vctrs::vec_group_loc(d$GEOID_county)
vifs = map(counties$loc, function(idx) {
    R = as.matrix(r_probs$county[idx, ])
    diag(solve(crossprod(R))) * colSums(R)
}) |>
    do.call(rbind, args=_)
p_gr = map(counties$loc, function(idx) {
    R = as.matrix(r_probs$county[idx, ])
    colSums(R)
}) |>
    do.call(rbind, args=_)

colSums(vifs * prop.table(p_gr, 2)^2) / colSums(prop.table(p_gr, 2)^2)

d_vif = as_tibble(vifs) |>
    mutate(county = counties$key) |>
    pivot_longer(-county, names_prefix="pr_", names_to="race", values_to="vif")
d_gr = as_tibble(prop.table(p_gr, 2)) |>
    mutate(county = counties$key) |>
    pivot_longer(-county, names_prefix="pr_", names_to="race", values_to="p_gr")
breaks = 10^seq(0, 3.6, length.out=19)
p_indiv = left_join(d_vif, d_gr, by=c("county", "race")) |>
    ggplot(aes(vif, weight=p_gr)) +
    facet_grid(~ fct_inorder(races[race])) +
    geom_histogram(breaks=breaks, fill="#777777", color="white", lwd=0.2) +
    scale_x_log10("Variance inflation due to BISG") +
    scale_y_continuous("Fraction of individuals", labels=scales::percent,
                       expand=expansion(mult=c(0.0, 0.02))) +
    theme_paper()
p_cty = ggplot(d_vif, aes(vif)) +
    facet_grid(~ fct_inorder(races[race])) +
    geom_histogram(breaks=breaks, fill="#777777", color="white", lwd=0.2) +
    scale_x_log10("Variance inflation due to BISG") +
    scale_y_continuous("Number of counties", labels=scales::comma,
                       expand=expansion(mult=c(0.0, 0.02))) +
    theme_paper()
(p_indiv / p_cty) & theme(plot.margin = margin(l=1, t=2, b=4))
ggsave(here("paper/figures/nc_vif.pdf"), width=7, height=4.0)
