# plan = drake_plan(
# 
#   # data used to make qq plot
#   pval_data = runif(1000),  # in reality: n > 10^9
# 
#   # other data dependencies of the Rmd
#   somedata = runif(10),
# 
#   # make qq plot as png file
#   qq = make_qq(v_pval = pval_data, fig_path = file_out("qq.png")),
# 
#   # generate report, which is to include the png
#   report = rmarkdown::render(
#     knitr_in("EDA.Rmd"),
#     output_file = file_out(here::here("EDA_report.html")),
#     quiet = TRUE
#   )
# )


# with explicit trigger

plan = drake_plan(

  # data used to make qq plot
  pval_data = runif(10),  # in reality: n > 10^9

  # other data dependencies of the Rmd
  somedata = runif(10),

  # make qq plot as png file
  qq = make_qq(v_pval = pval_data, fig_path = file_out("qq.png")),

  # generate report
  report = target(
    command = rmarkdown::render(
      knitr_in("EDA.Rmd"),
      output_file = file_out(here::here("EDA_report.html")),
      quiet = TRUE
    ),
    trigger = trigger(change = digest(file.info(file_in("qq.png"))))  # problem: seems to interfere with automatic dependency detection in Rmd
  )
)
