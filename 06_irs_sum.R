library(easycensus)

# IRS text output ------------------------
morg_int_brk = c(1, 3785.9, 5681, 7204, 8734, 10390, 12346, 14618, 17298, 22301, 750000)

res_txt = "
        white       black        hisp       asian        aian      other
1  0.89394365 0.934647100 0.953928413 0.878645558 0.869804696 0.82222368
2  0.01063682 0.011650768 0.003147550 0.005412941 0.017885690 0.01714904
3  0.01056057 0.010476624 0.004052581 0.006082662 0.015768905 0.01929709
4  0.01013059 0.009146048 0.004552202 0.006604051 0.018265482 0.01978764
5  0.01047017 0.008145958 0.005276701 0.007706196 0.016134144 0.02019548
6  0.01011557 0.007111155 0.005526864 0.009067576 0.014205106 0.02039068
7  0.01035289 0.005784415 0.005576239 0.011114815 0.012223831 0.01700871
8  0.01078163 0.004272466 0.005389839 0.013982239 0.009505774 0.01853845
9  0.01069592 0.003655878 0.005197625 0.017011964 0.009813602 0.01657461
10 0.01095254 0.002905778 0.004598053 0.021007179 0.008118151 0.01665905
11 0.01135965 0.002203808 0.002753933 0.023364818 0.008274620 0.01217556
"

res_amt = c(1428.9228, 633.5371, 597.4768, 2096.2345, 1482.6316, 2215.9076)


# process ------------------------------
fmt_brk = label_dollar(accuracy=0.1, scale=1e-3, suffix="k")
morg_int_lbl = str_c("[", fmt_brk(lag(morg_int_brk, 1, default=0)),
                     ", ", fmt_brk(morg_int_brk), ")")
morg_int_lbl[1] = "$0"
morg_int_lbl[2] = str_c("[$1, ", fmt_brk(morg_int_brk[2]), ")")
morg_int_lbl[11] = str_c("[", fmt_brk(morg_int_brk[10]), ", $750k]")

res_m = str_split_1(res_txt, "\n") |>
    str_split("\\s+") |>
    keep(~ length(.) > 1) |>
    map(~ .[-1]) # remove row numbers

d_est = res_m[-1] |>
    map(function(row) {
        tibble::new_tibble(list(
            race = races[res_m[[1]]],
            est = as.numeric(row)
        ))
    }) |>
    bind_rows(.id = "morg_int") |>
    mutate(race = fct_inorder(race),
           morg_int = fct_inorder(morg_int_lbl[as.integer(morg_int)]))

d_amt = tibble(race = fct_inorder(races),
               est = res_amt / (1 - d_est$est[1:6]))

res_m[-1] |>
    map(as.numeric) |>
    do.call(rbind, args=_) |>
    `colnames<-`(res_m[[1]]) |>
    `rownames<-`(morg_int_lbl) |>
    write_rds(here("paper/data/irs_pct.rds"), compress="gz")
(res_amt / (1 - as.numeric(res_m[[2]]))) |>
    setNames(res_m[[1]]) |>
    write_rds(here("paper/data/irs_amt.rds"), compress="gz")


# Census benchmarks --------------------------

# cens_find_dec(tenure)
d_mortg = cens_get_dec(tables_sf1$H4, geo="zcta")
idx_own = which(str_detect(d_mortg$tenure[1:3], "free and clear"))
idx_mort = which(str_detect(d_mortg$tenure[1:3], "mortgage"))
d_mortg = d_mortg |>
    group_by(GEOID) |>
    summarize(pct_mortg = value[idx_mort] / (value[idx_own] + value[idx_mort]))

# Check: can we use racial comp. to predict fraction of owners w/ a mortgage?
d_ei = census_race_geo_table("zcta") |>
    left_join(d_mortg, by="GEOID") |>
    mutate(total = white + black + hisp + asian + aian + other,
           across(white:other, ~ . / total))
ggplot(d_ei, aes(white, pct_mortg, size=total)) +
    geom_point(alpha=0.5) +
    scale_size_area(max_size=2) +
    geom_smooth()
# Answer: no. So probably OK to assume this fraction is const. across races
lm(pct_mortg ~ white + black + hisp + asian + aian + other - 1, data=d_ei) |>
    summary()


tables_sf2 = cens_parse_tables("dec/sf1", 2010)
# cens_find(tables_sf2, tenure, race)
d_ten_race = cens_get_dec(tables_sf2$HCT1, geo="zcta")
d_ten_race = d_ten_race |>
    transmute(GEOID = GEOID,  tenure = tenure, value = value,
              hisp = tidy_ethnicity(hispanic_or_latino_origin_of_householder),
              race = tidy_race(str_remove(race_of_householder, "^householder who is "))) |>
    filter(hisp != "total", !is.na(hisp),
           race != "total", !is.na(race)) |>
    mutate(race = case_when(
        hisp == "hisp" ~ "hisp",
        race %in% c("nhpi", "other", "two") ~ "other",
        TRUE ~ race
    )) |>
    count(GEOID, tenure, race, wt=value)

# join, impute with average of nearby ZIPs, and tally to nation
d_mortg_race = d_ten_race |>
    left_join(d_mortg, by="GEOID") |>
    group_by(zip_g3 = str_sub(GEOID, 1, 3)) |>
    mutate(pct_mortg = coalesce(pct_mortg, mean(pct_mortg, na.rm=TRUE))) |>
    group_by(zip_g2 = str_sub(GEOID, 1, 2)) |>
    mutate(pct_mortg = coalesce(pct_mortg, mean(pct_mortg, na.rm=TRUE))) |>
    group_by(race) |>
    summarize(total = sum(n),
              pct_owners = sum(n[tenure == "owner occupied"]) / total,
              pct_mortg = sum((n * pct_mortg)[tenure == "owner occupied"]) / total) |>
    mutate(race = factor(race, levels=names(races), labels=races)) |>
    arrange(race)


d_mortg_race |>
    mutate(across(pct_owners:pct_mortg, ~ percent(., accuracy=0.1))) |>
    transmute(`Race`=race, `Total households`=comma(total),
           `Fraction owner-occupied`=pct_owners,
           `Fraction with mortgage`=pct_mortg) |>
    write_rds(here("paper/data/race_mortg.rds"), compress="gz")


# plot ---------------------------------

# expected rate if all groups had the same prop. having mortgages
d_mort_adj = d_est |>
    filter(as.integer(morg_int) == 1, # only non-users
           !race %in% c("Native", "Other")) |>
    left_join(d_mortg_race, by="race") |>
    group_by(morg_int) |>
    arrange(race) |>
    mutate(est = 1 - est, # pct that use
           est_nodisp = est[1] * pct_mortg / pct_mortg[1]) |>
           # est_nodisp = weighted.mean(est, total) * pct_mortg /
           #     weighted.mean(pct_mortg, total)) |>
    ungroup() |>
    select(-morg_int, -total, -pct_owners)

p1 = d_est |>
    filter(as.integer(morg_int) > 1, # remove non-users
           !race %in% c("Native", "Other")) |>
    mutate(race = fct_reorder(race, -est)) |>
ggplot(aes(race, est, fill=morg_int)) +
    geom_col(color="white", linewidth=0.2) +
    geom_text(aes(label=morg_int, group=morg_int,
                  color=as.integer(morg_int) > 6,
                  alpha=as.integer(morg_int) %% 3 == 2),
              position=position_stack(vjust = 0.5), size=1.9, family="Times") +
    geom_segment(aes(x = as.integer(race) - 0.5, xend = as.integer(race) + 0.5,
                     y = est_nodisp, yend = est_nodisp),
                 data=d_mort_adj, inherit.aes=FALSE,
                 linewidth=0.8, linetype="11") +
    scale_fill_wa_d("volcano", name="Home mortgage\ninterest deduction", reverse=TRUE) +
    scale_color_manual(values=c("black", "white"), guide="none") +
    scale_alpha_manual(values=c(0, 1), guide="none") +
    scale_y_continuous("Estimated fraction of filers", labels=percent,
                       minor_breaks=seq(0, 1, 0.01),
                       expand=expansion(mult=c(0, 0.05))) +
    labs(x=NULL) +
    theme_paper() +
    theme(panel.grid.major.x=element_blank(),
          legend.position="left",
          legend.margin=margin())

p2 = d_amt |>
    filter(!race %in% c("Native", "Other")) |>
    left_join(filter(d_est, morg_int == "$0"),
              by="race", suffix=c("", "0")) |>
ggplot(aes(reorder(race, est0), est)) +
    geom_col(fill=wacolors$volcano[3]) +
    scale_y_continuous("Estimated average HMID amount", labels=dollar,
                       minor_breaks=seq(0, 20000, 1000),
                       expand=expansion(mult=c(0, 0.05))) +
    labs(x=NULL) +
    theme_paper() +
    theme(panel.grid.major.x=element_blank(),
          plot.title.position="plot")

p = (p1 + labs(title="(a) Deduction distribution")) +
    (p2 + labs(title="       (b) Average deduction when taken")) +
    plot_layout(widths=c(2, 1.5)) &
    theme(legend.margin=margin(),
          plot.margin=margin(l=4, r=2))
ggsave(here("paper/figures/irs_hmid.pdf"), plot=p, width=8, height=4.25)
