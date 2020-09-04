library(gt); library(tidyverse)
data("gtcars")

gtcars_8 <-
  gtcars %>%
  dplyr::group_by(ctry_origin) %>%
  dplyr::top_n(2) %>%
  dplyr::ungroup() %>%
  dplyr::filter(ctry_origin != "United Kingdom")

gtcars_8 %>% 
  group_by(ctry_origin) %>% 
  gt() %>% 
  gtsave('test.rtf')

gt::render_gt(gtcars_8 %>% 
                group_by(ctry_origin) %>% 
                gt())


res_rd %>% 
  filter(param %in% c('rdark_30','Q10')) %>% 
  select(param, Treat, `50%`, sd) %>% 
  group_by(Treat) %>% 
  gt()

stan_summary <- function(stanfit, probs=c(0.1,0.5,0.9), treatment){
  library(tidyverse); library(rstan); 
  tmp <- summary(stanfit, probs=probs)$summary
  par_names <- rownames(tmp)
  out <- tmp %>% as_tibble() %>%  
    mutate(Parameter = par_names, Treatment=treatment)
  # out$Parameter <- par_names
  return(out)
}
v_cc <- stan_summary(vc400_c.c, treatment="Control to Control")
v_cw <- stan_summary(vc400_c.w, treatment="Control to Treatment")
v_wc <- stan_summary(vc400_w.c, treatment="Treatment to Control")
v_ww <- stan_summary(vc400_w.w, treatment="Treatment to Treatment")
tmp <- bind_rows(v_cc, v_cw, v_wc, v_ww)

j_cc <- stan_summary(jm400_c.c, treatment="Control -> Control")
j_cw <- stan_summary(jm400_c.w, treatment="Control -> Treatment")
j_wc <- stan_summary(jm400_w.c, treatment="Treatment -> Control")
j_ww <- stan_summary(jm400_w.w, treatment="Treatment -> Treatment")
tmp2 <- bind_rows(j_cc, j_cw, j_wc, j_ww)

# tmp3 <- bind_rows(tmp %>% mutate(Response="Vcmax"), tmp2 %>% mutate(Response="Jmax"))

#! this one works for Vcmax
out_vcmax <- tmp %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  mutate(Parameter = case_when(
    Parameter == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
    Parameter == "Vcmax_Topt" ~ "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
    Parameter == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
    Parameter == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
    Parameter == "sigma" ~ "V<sub>Cmax</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>% 
  gt(groupname_col = 'Parameter') %>% 
  row_group_order(groups = c("T<sub>opt</sub> (&deg;C)",
                             "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
                             "V<sub>Cmax</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>%
  tab_header(title = html("V<sub>Cmax<sub>")) %>%
  fmt_number(columns = vars(`10%`,`50%`,`90%`), decimals=2)
out_vcmax

#! 
out_jmax <- tmp2 %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  mutate(Parameter = case_when(
    Parameter == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
    Parameter == "Jmax_Topt" ~ "J<sub>max</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
    Parameter == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
    Parameter == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
    Parameter == "sigma" ~ "J<sub>max</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>% 
  gt(groupname_col = 'Parameter') %>% 
  row_group_order(groups = c("T<sub>opt</sub> (&deg;C)",
                             "J<sub>max</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
                             "J<sub>max</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>%
  tab_header(title = html("J<sub>max<sub>")) %>%
  fmt_number(columns = vars(`10%`,`50%`,`90%`), decimals=2)

out_vcmax_L <- out_vcmax %>% as_latex()  
out_vcmax_L[[1]] %>% writeLines(con = "doc/Table_Vcmax_params.tex")

out_jmax %>% as_latex() %>% as.character() %>% writeLines(con = "doc/Table_Jmax_params.tex")

out_vcmax %>% as_raw_html() %>% as.character() %>% writeLines(con = "doc/Table_Vcmax_params.html")


tmp %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  mutate(Parameter = case_when(
      Parameter == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
      Parameter == "Vcmax_Topt" ~ "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
      Parameter == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
      Parameter == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
      Parameter == "sigma" ~ "V<sub>Cmax</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>% 
  # arrange(factor(Treatment)) %>% 
  # group_by(Parameter) %>% 
  gt(groupname_col = 'Parameter') %>% 
  # gt(
  #   # groupname_col = ("Treatment")
  #   groupname_col = "Parameter"
  #   # rowname_col='Parameter'
  # ) %>% 
  row_group_order(groups = c("T<sub>opt</sub> (&deg;C)",
                             "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
                             "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
                             "V<sub>Cmax</sub> &sigma; (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)")) %>%
  tab_header(title = html("V<sub>Cmax<sub>")) %>%
  fmt_number(columns = vars(`10%`,`50%`,`90%`), decimals=2) #%>%
  # text_transform(
  #   locations = cells_stub(        # cells_data
  #     # columns = vars(Parameter)),
  #     rows = vars(Parameter),
  #   fn = function(x) {
  #     dplyr::case_when(
  #       x == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
  #       x == "Vcmax_Topt" ~ "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
  #       x == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
  #       x == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
  #       x == "sigma" ~ "&sigma;")
  #   }
  # ))


tmp %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  # arrange(factor(Treatment)) %>% 
  # group_by(Treatment) %>%
  gt(
    # groupname_col = ("Treatment")
    groupname_col = "Parameter"
    # rowname_col='Parameter'
     ) %>% 
  row_group_order(groups = c("Vcmax_Topt","Topt","Ha","kopt","sigma")) %>%
  tab_header(title = html("V<sub>Cmax<sub>")) %>%
  fmt_number(columns = vars(`10%`,`50%`,`90%`), decimals=2) %>%
  text_transform(
    locations = cells_data(        # cells_data
      columns = vars(Parameter)),
    fn = function(x) {
        dplyr::case_when(
          x == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
          x == "Vcmax_Topt" ~ "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
          x == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
          x == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)",
          x == "sigma" ~ "&sigma;")
    }
  )

tmp %>% filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  arrange(factor(Treatment)) %>% 
  group_by(Treatment) %>%
  gt(
    # groupname_col = ("Treatment")
    # rowname_col='Parameter'
  ) %>% 
  # row_group_order(c("Vcmax_Topt","Topt","Ha","sigma")) %>% 
  tab_header(title = html("V<sub>Cmax<sub>")) %>% 
  fmt_number(columns = vars(`10%`,`50%`,`90%`), decimals=2) %>% 
  text_transform(
    locations = cells_data(        # cells_data
      columns = vars(Parameter)),
    fn = function(x) {
      dplyr::case_when(
        x == "Topt"  ~ "T<sub>opt</sub> (&deg;C)",
        x == "Vcmax_Topt" ~ "V<sub>Cmax</sub> at T<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)",
        x == "kopt"   ~ "k<sub>opt</sub> (&#956;mol m<sup>-2</sup>s<sup>-1</sup>)", 
        x == "Ha" ~ "H<sub>a</sub> (kJ mol<sup>-1</sup>)", 
        x == "sigma" ~ "&sigma;")
    }
  ) %>% 
  tab_spanner(label='Interval', columns=vars(`10%`,`50%`,`90%`))




%>%
  # tab_row_group(group=md("**blah**"), rows = Parameter == "kopt")
  # tab_row_group(group=html("Temp,<br>&deg;F"), rows = Parameter == "kopt")
  cols_label(kopt = html("Kopt,<br>umol"))


tmp %>% filter(Parameter != "lp__") %>% 
  select(Treatment,Parameter, `10%`, `50%`, `90%`) %>% 
  group_by(Parameter) %>% 
  gt(groupname_col = Parameter) %>% 
  tab_header(title = html("V<sub>Cmax<sub>"))




tmp %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment, Parameter, `10%`, `50%`, `90%`) %>% 
  # gather(interval, estimate, -Treatment, -Parameter) %>% 
  # group_by(Treatment) %>% 
  gt(rowname_col="Parameter", groupname_col = "Treatment", 
     rownames_to_stub = F)

tmp %>% 
  filter(Parameter != "lp__") %>% 
  select(Treatment, Parameter, `10%`, `50%`, `90%`) %>% 
  gather(Level,Estimate, -Treatment,-Parameter) %>%
  spread(Level,Estimate)
  # group_by(Treatment) %>% 
  # gt(groupname_col = Parameter)
  # spread(Parameter, estimate, -Treatment)


library(dplyr)
stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)
stocksm <- stocks %>% gather(stock, price, -time)
stocksm %>% spread(stock, price)
stocksm %>% spread(time, price)


sza %>%
  dplyr::filter(latitude == 20) %>%
  dplyr::select(-latitude) %>%
  dplyr::filter(!is.na(sza)) %>%
  tidyr::spread(key = "tst", value = sza) %>%
  gt(rowname_col = "month") %>%
  fmt_missing(
    columns = TRUE,
    missing_text = ""
  ) %>%
  tab_stubhead_label(label = html("month<br>(20&deg;N)")) %>%
  tab_header(title = html("&#x2600; Solar Zenith Angles &#x2600;")) %>%
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    row.padding = px(3)
  )

order_countries <- c("Germany", "Italy", "United States", "Japan")

# Use dplyr functions to get the car with the best city gas mileage;
# this will be used to target the correct cell for a footnote
best_gas_mileage_city <- 
  gtcars %>% 
  dplyr::arrange(desc(mpg_c)) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(car = paste(mfr, model)) %>%
  dplyr::pull(car)

# Use dplyr functions to get the car with the highest horsepower
# this will be used to target the correct cell for a footnote
highest_horsepower <- 
  gtcars %>% 
  dplyr::arrange(desc(hp)) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(car = paste(mfr, model)) %>%
  dplyr::pull(car)
# Create a display table with `gtcars`, using all of the previous
# statements piped together + additional `tab_footnote()` stmts
tab <-
  gtcars %>%
  dplyr::arrange(
    factor(ctry_origin, levels = order_countries),
    mfr, desc(msrp)
  ) %>%
  dplyr::mutate(car = paste(mfr, model)) %>%
  dplyr::select(-mfr, -model) %>%
  dplyr::group_by(ctry_origin) %>%
  gt(rowname_col = "car") %>%
  cols_hide(columns = vars(drivetrain, bdy_style)) %>%
  cols_move(
    columns = vars(trsmn, mpg_c, mpg_h),
    after = vars(trim)
  ) %>%
  tab_spanner(
    label = md("*Performance*"),
    columns = vars(mpg_c, mpg_h, hp, hp_rpm, trq, trq_rpm)
  ) %>%
  cols_merge(
    col_1 = vars(mpg_c),
    col_2 = vars(mpg_h),
    pattern = "{1}c<br>{2}h"
  ) %>%
  cols_merge(
    col_1 = vars(hp),
    col_2 = vars(hp_rpm),
    pattern = "{1}<br>@{2}rpm"
  ) %>%
  cols_merge(
    col_1 = vars(trq),
    col_2 = vars(trq_rpm),
    pattern = "{1}<br>@{2}rpm"
  ) %>%
  cols_label(
    mpg_c = "MPG",
    hp = "HP",
    trq = "Torque",
    year = "Year",
    trim = "Trim",
    trsmn = "Transmission",
    msrp = "MSRP"
  ) %>%
  fmt_currency(
    columns = vars(msrp),
    currency = "USD",
    decimals = 0
  ) %>%
  cols_align(
    align = "center",
    columns = vars(mpg_c, hp, trq)
  ) %>%
  tab_style(
    style = cells_styles(text_size = px(12)),
    locations = cells_data(columns = vars(trim, trsmn, mpg_c, hp, trq))
  ) %>%
  text_transform(
    locations = cells_data(columns = vars(trsmn)),
    fn = function(x) {
      
      speed <- substr(x, 1, 1)
      
      type <-
        dplyr::case_when(
          substr(x, 2, 3) == "am" ~ "Automatic/Manual",
          substr(x, 2, 2) == "m" ~ "Manual",
          substr(x, 2, 2) == "a" ~ "Automatic",
          substr(x, 2, 3) == "dd" ~ "Direct Drive"
        )
      
      paste(speed, " Speed<br><em>", type, "</em>")
    }
  ) %>%
  tab_header(
    title = md("The Cars of **gtcars**"),
    subtitle = "These are some fine automobiles"
  ) %>%
  tab_source_note(
    source_note = md(
      "Source: Various pages within [edmunds.com](https://www.edmunds.com).")
  ) %>%
  tab_footnote(
    footnote = md("Best gas mileage (city) of all the **gtcars**."),
    locations = cells_data(
      columns = vars(mpg_c),
      rows = best_gas_mileage_city)
  ) %>%
  tab_footnote(
    footnote = md("The highest horsepower of all the **gtcars**."),
    locations = cells_data(
      columns = vars(hp),
      rows = highest_horsepower)
  ) %>%
  tab_footnote(
    footnote = "All prices in U.S. dollars (USD).",
    locations = cells_column_labels(
      columns = vars(msrp))
  )
tab
