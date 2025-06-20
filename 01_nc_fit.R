# Setup --------
voters = read_rds(here("data-raw/nc_voters.rds"))
d_cens = read_rds(here("data-raw/nc_block_vap_race.rds")) %>%
    drop_na()
d_county = group_by(d_cens, county) %>%
    summarize(across(vap:vap_other, sum), .groups="drop") %>%
    mutate(tract = NA, block = NA)

# voters = select(voters, last_name, zip:county, race, party, lic, turnout=voted_2020_11)

# make a R|GZ table for predict_race_sgz
make_p_rgx = function(voters, level=c("block", "tract", "county", "zip"), counts=FALSE) {
    if (level == "block") {
        voters_geo = distinct(voters, county, tract, block)
        d = bind_rows(d_cens, d_county) %>%
            right_join(voters_geo, by=c("county", "tract", "block")) %>%
            mutate(GEOID = ifelse(is.na(block), as.character(county),
                                  str_c(county, tract, block)),
                   .before="county") %>%
            select(-county, -tract, -block)
    } else if (level == "tract") {
        voters_geo = distinct(voters, county, tract)
        d = group_by(d_cens, county, tract) %>%
            summarize(across(vap:vap_other, function(x) sum(x, na.rm=TRUE)), .groups="drop") %>%
            bind_rows(d_county) %>%
            right_join(voters_geo, by=c("county", "tract")) %>%
            mutate(GEOID = ifelse(is.na(tract), as.character(county),
                                  str_c(county, tract)),
                   .before="county") %>%
            select(-county, -tract, -block)
    } else if (level == "zip") {
        voters_geo = voters %>%
            mutate(zip = na_if(proc_zip(zip), "<none>")) %>%
            distinct(county, zip) %>%
            mutate(GEOID = if_else(is.na(zip), str_c("cty", county), zip)) %>%
            distinct(GEOID)
        d_county = mutate(d_county, GEOID = str_c("cty", county)) %>%
            select(-county, -tract, -block)
        d = census_race_geo_table("zcta", year=2010, counts=TRUE) %>%
            rename_with(~ str_c("vap_", .), white:other) %>%
            mutate(vap = rowSums(across(vap_white:vap_other), na.rm=TRUE)) %>%
            bind_rows(d_county) %>%
            right_join(voters_geo, by="GEOID")
    } else {
        d = rename(d_county, GEOID = county) %>%
            select(-tract, -block)
    }

    if (isFALSE(counts)) {
        d = mutate(d, across(vap_white:vap_other, ~ (0.01 + .) / (vap + 0.06)))
    }

    d %>%
        select(-vap) %>%
        rename_with(~ str_sub(., 5), vap_white:vap_other) %>%
        mutate(GEOID = as.character(GEOID))
}

# Subsample a table and format geography and outcome variables
set.seed(5118)

d = voters |>
    semi_join(bind_rows(d_cens, d_county), by=c("county", "tract", "block")) |>
    mutate(GEOID_county = as.character(county),
           GEOID_tract = if_else(is.na(tract), GEOID_county, str_c(county, tract)),
           GEOID_block = if_else(is.na(block), GEOID_county, str_c(county, tract, block)),
           GEOID_zip = if_else(is.na(zip), str_c("cty", county), as.character(zip)),
           voted_20 = c("no", "yes")[voted_2020_11 + 1],
           voted_21 = c("no", "yes")[voted_2021_11 + 1],
           party=coalesce(party, "ind")) |>
    slice_sample(n=1e6)
print(as.character(head(d$last_name))) # ensure seed is working


# Do BISG --------
p_r = prop.table(table(d$race))

geo_levels = c("county", "zip", "tract", "block")

r_probs = map(geo_levels, function(level) {
    cat("Doing BISG at the", level, "level\n")
    bisg(~ nm(last_name) + GEOID,
         data=rename(d, GEOID=str_c("GEOID_", level)),
         p_rgx=make_p_rgx(d, level), p_r=p_r)
}) %>%
    set_names(geo_levels)

# BISG quality
log_score_baseline = mean(log(p_r[as.integer(d$race)]))
log_scores = map_dbl(r_probs, log_score)
acc_thresh = map_dbl(r_probs, ~ mean(max.col(.) == as.integer(d$race)))

if (!file.exists(path <- here("paper/data/nc_bisg.rds"))) {
    list(base_score = log_score_baseline,
         score = log_scores,
         acc = acc_thresh) |>
        write_rds(path)
}

if (FALSE) {
    d_acc = imap_dfr(r_probs, function(tbl, lbl) {
        d_auc = tbl |>
            summarize(across(everything(), function(pr) {
                fastAUC(pr, d$race == str_sub(cur_column(), 4))
            })) |>
            rename_with(~ str_sub(., 4)) |>
            t()
        tibble(level = lbl,
               acc = mean(predict(tbl) == d$race),
               race = rownames(d_auc),
               auc = d_auc[, 1])
    })

    ggplot(d_acc, aes(fct_inorder(level), auc)) +
        facet_wrap(~ fct_inorder(race)) +
        geom_col() +
        scale_y_continuous(expand=c(0, 0), limits=c(0, 1))

    d_acc |>
        filter(level=="county") |>
    ggplot(aes(fct_inorder(toupper(race)), auc)) +
        geom_col(fill="#2c3f4a") +
        scale_y_continuous(name="AUROC", limits=c(0, 1),
                           expand=c(0, 0), labels=percent) +
        labs(x="Race") +
        theme_bw(base_family="IBM Plex Sans", base_size=24)
        # theme_paper()

    ggplot(d, aes(fct_relabel(race, toupper), after_stat(count)/nrow(d))) +
        geom_bar(fill="#2c3f4a") +
        scale_y_continuous(name="Share of voter file", limits=c(0, 1),
                           expand=c(0, 0), labels=percent) +
        labs(x="Race") +
        theme_bw(base_family="IBM Plex Sans", base_size=24)
}


# Party ID -----------
ctrl = birdie.ctrl(abstol=1e-5)
prior_mm = list(
    scale_int = 5,
    scale_beta = 2 * p_r,
    scale_sigma = 0.5
)
prior_sat_gibbs = list(alpha = matrix(0.5, nrow = 4, ncol = 6))

fits_party = imap(r_probs, function(d_pr, level) {
    cat(level, "\n")

    if (level == "block") level = "tract"
    # algo = if (level == "county") "em_boot" else "em"

    link_vec = "GEOID"
    names(link_vec) = str_c("GEOID_", level)
    # fill in pct white and black for mixed model
    d_mm = left_join(d, make_p_rgx(d, level), by=link_vec) |>
        rename(GEOID=names(link_vec)) |>
        mutate(across(white:other, ~ coalesce(., p_r[cur_column()]))) # fill in unlinked

    suppressWarnings(list(
        sat = birdie(d_pr, party ~ GEOID, data=d_mm,
                     algorithm="gibbs", iter=500, warmup=25,
                     # algorithm=algo, iter=100,
                     family=cat_dir(), prior=prior_sat_gibbs, ctrl=ctrl),
        mmm = birdie(d_pr, party ~ white + black + (1 | GEOID), data=d_mm,
                     family=cat_mixed(), prior=prior_mm, ctrl=ctrl)
    ))
})

# for (level in names(r_probs)) {
# # for (level in "county") {
#     cat(level, "\n")
#     nm_level = level
#     if (level == "block") level = "tract"
#
#     link_vec = "GEOID"
#     names(link_vec) = str_c("GEOID_", level)
#     # fill in pct white and black for mixed model
#     d_mm = left_join(d, make_p_rgx(d, level), by=link_vec) |>
#         rename(GEOID=names(link_vec)) |>
#         mutate(across(white:other, ~ coalesce(., p_r[cur_column()]))) # fill in unlinked
#
#     fits_party[[nm_level]]$sat_gibbs = birdie(d_pr, party ~ GEOID, data=d_mm,
#                                            algorithm="gibbs", iter=100, warmup=25,
#                                            family=cat_dir(), prior=prior_sat_gibbs, ctrl=ctrl)
# }

print(max(map_dbl(geo_levels, ~ fits_party[[.]]$sat$algo$runtime)))
print(max(map_dbl(geo_levels, ~ fits_party[[.]]$mmm$algo$runtime)))

rm(make_p_rgx, voters, d_cens)
