helptab_ui <- fluidPage(
  id = "tab1_content",
  titlePanel("User Manual"),
  fluidRow(
    column(6,
           div(
            h4(
              "Manual for Flow Cytometry in SCPline: ",
              a("Manual", href = "data/help_doc/Flow Cytometry.pdf", target = "_blank")
            ),
            h4(
              "Manual for Single Cell Multiomics in SCPline: ",
              a("Manual", href = "data/help_doc/scMultimics.pdf", target = "_blank")
            ),
            h4(
              "Manual for Single MASS Workflow in SCPline: ",
              a("Manual", href = "data/help_doc/MASS.pdf", target = "_blank")
            )
           )
    ),
    column(6,
           img(src = "data/overview_img/homepage.png", width = "500px", height = "500px")
    ),
  )
)
