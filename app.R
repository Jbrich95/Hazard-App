# app.R

library(shiny)
library(rhandsontable)
library(splines2)

#------------------------------------------------------------
# Default data
#------------------------------------------------------------
a <- c(0, 6, 12, 24, 48, 96, 120)
b <- c(6, 12, 24, 48, 96, 120, 144)


default_data <- data.frame(
  time_period = c(
    "0-6 months",
    "6-12 months",
    "1-2 years",
    "2-4 years",
    "4-8 years",
    "8-10 years",
    "10-12 years"
  ),
  Probability  = c(0.20, 0.10, 0.10, 0.30, 0.1, 0.1, 0.1),
  Days_lost_average     = c(70, 50, 50, 75, 30, 30, 30),
  Prob_multiplier = c(2, 3, 1, 1, 1, 1, 1),
  Days_lost_multiplier = c(1, 3, 1, 1, 1, 1, 1)
)

#------------------------------------------------------------
# Model fitting function
#------------------------------------------------------------

fit_hazard_model <- function(dat)
{
  
  
  cond.prob <- dat$Probability
  Y <- dat[["Days_lost_average"]]  
  knots <- unique(a[-1])
  knots <- knots[knots < max(b)]
  
  if(length(knots) < 1)
    knots <- median(a)
  
  boundary_max <- max(b)
  basis_fun <- function(t) {
    iSpline(
      t,
      degree = 3,
      knots = knots,
      Boundary.knots = c(0, boundary_max),
      intercept = TRUE
    )
  }
  
  h_basis_fun <- function(t) {
    mSpline(
      t,
      degree = 3,
      knots = knots,
      Boundary.knots = c(0, boundary_max),
      intercept = TRUE
    )
  }
  
  obj <- function(par)
  {
    H_fun <- function(t)
      basis_fun(t) %*% par
    
    pred.cum.hazard <- sapply(
      seq_along(a),
      function(i)
        H_fun(b[i]) - H_fun(a[i])
    )
    
    pred.prob <- 1 - exp(-pred.cum.hazard)
    
    sum((pred.prob - cond.prob)^2)
  }
  
  p0 <- rep(0.01, ncol(basis_fun(1)))
  
  fit <- nlminb(
    p0,
    obj,
    lower = 0,
    control = list(iter.max = 2000)
  )
  
  par_hat <- fit$par
  
  H_hat <- function(t)
    as.vector(basis_fun(t) %*% par_hat)
  
  h_hat <- function(t)
    as.vector(h_basis_fun(t) %*% par_hat)
  
  S_hat <- function(t)
    exp(-H_hat(t))
  
  f_hat <- function(t)
    S_hat(t) * h_hat(t)
  
  fitted_interval <- sapply(
    seq_along(a),
    function(i) {
      1 - exp(-(H_hat(b[i]) - H_hat(a[i])))
    }
  )
  
  probs_interval <- sapply(
    seq_along(a),
    function(i) {
      S_hat(a[i]) - S_hat(b[i])
    }
  )
  
  #E_unconditional <- sum(Y * probs_interval)
  
  Y_func = function(t){
  #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y[1],Y, Y[length(Y)]), xout = t)$y
  }
  f_hat <- function(t)
    S_hat(t) * h_hat(t)
  Risk_hat <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func(x))
  E_unconditional <- integrate(function(t)
    Risk_hat(t), 0, max(b))$value
  
  list(
    h_hat = h_hat,
    H_hat = H_hat,
    S_hat = S_hat,
    f_hat = f_hat,
    a = a,
    b = b,
    Y = Y,
    Y_func= Y_func,
    Risk_hat = Risk_hat,
    fitted_interval = fitted_interval,
    E_unconditional = E_unconditional,
    fit = fit
  )
}

#------------------------------------------------------------
# Model fitting function
#------------------------------------------------------------

fit_hazard_model_adjusted <- function(dat)
{
  
  mult = dat$Prob_multiplier
  cond.prob <- 1-(1-dat$Probability)^mult
  Y <- dat[["Days_lost_average"]]  * dat[["Days_lost_multiplier"]] 
  knots <- unique(a[-1])
  knots <- knots[knots < max(b)]
  
  if(length(knots) < 1)
    knots <- median(a)
  
  boundary_max <- max(b)
  basis_fun <- function(t) {
    iSpline(
      t,
      degree = 3,
      knots = knots,
      Boundary.knots = c(0, boundary_max),
      intercept = TRUE
    )
  }
  
  h_basis_fun <- function(t) {
    mSpline(
      t,
      degree = 3,
      knots = knots,
      Boundary.knots = c(0, boundary_max),
      intercept = TRUE
    )
  }
  
  obj <- function(par)
  {
    H_fun <- function(t)
      basis_fun(t) %*% par
    
    pred.cum.hazard <- sapply(
      seq_along(a),
      function(i)
        H_fun(b[i]) - H_fun(a[i])
    )
    
    pred.prob <- 1 - exp(-pred.cum.hazard)
    
    sum((pred.prob - cond.prob)^2)
  }
  
  p0 <- rep(0.01, ncol(basis_fun(1)))
  
  fit <- nlminb(
    p0,
    obj,
    lower = 0,
    control = list(iter.max = 2000)
  )
  
  par_hat <- fit$par
  
  H_hat <- function(t)
    as.vector(basis_fun(t) %*% par_hat)
  
  h_hat <- function(t)
    as.vector(h_basis_fun(t) %*% par_hat)
  
  S_hat <- function(t)
    exp(-H_hat(t))
  
  f_hat <- function(t)
    S_hat(t) * h_hat(t)
  
  fitted_interval <- sapply(
    seq_along(a),
    function(i) {
      1 - exp(-(H_hat(b[i]) - H_hat(a[i])))
    }
  )
  
  probs_interval <- sapply(
    seq_along(a),
    function(i) {
      S_hat(a[i]) - S_hat(b[i])
    }
  )
  

  Y_func = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y[1],Y, Y[length(Y)]), xout = t)$y
  }
 # E_unconditional <- sum(Y * probs_interval)
 
  f_hat <- function(t)
    S_hat(t) * h_hat(t)
  Risk_hat <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func(x))
  E_unconditional <- integrate(function(t)
    Risk_hat(t), 0, max(b))$value

  
  list(
    h_hat = h_hat,
    H_hat = H_hat,
    S_hat = S_hat,
    f_hat = f_hat,
    a = a,
    b = b,
    Y = Y,
    Y_func = Y_func,
    Risk_hat = Risk_hat,
    fitted_interval = fitted_interval,
    E_unconditional = E_unconditional,
    fit = fit
  )
}
#------------------------------------------------------------
# UI
#------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Injury Risk / Time-Loss Tool"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Interval Inputs"),
      
      rHandsontableOutput("tbl"),
      
      br(),
      
      
      width = 4
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Baseline Hazard",
          plotOutput("hazard_plot", height = 500)
        ),
    
        tabPanel(
          "Baseline Density",
          plotOutput("density_plot", height = 500)
        ),
        tabPanel(
          "Baseline Risk",
          plotOutput("risk_plot", height = 500)
        ),
        tabPanel(
          "Compounded Hazard",
          plotOutput("compare_hazard_plot", height = 500)
        ),
        tabPanel(
          "Compounded Density",
          plotOutput("compare_density_plot", height = 500)
        ),
        tabPanel(
          "Compounded Risk",
          plotOutput("compare_risk_plot", height = 500)
        ),
        
        
        
        tabPanel(
          "Expected days lost",
          tableOutput("results_tbl")
        )
      )
    )
  )
)

#------------------------------------------------------------
# Server
#------------------------------------------------------------

server <- function(input, output, session)
{
  
  rv <- reactiveValues(
    dat = default_data
  )
  
  output$tbl <- renderRHandsontable({
    
    rhandsontable(
      rv$dat,
      rowHeaders = F  ) |>
      hot_col("time_period", readOnly = TRUE,   renderer = "
        function(instance, td, row, col, prop, value, cellProperties) {
          Handsontable.renderers.TextRenderer.apply(this, arguments);
          td.style.background = '#f0f0f0';
          td.style.color = '#666';
        }
      "
    ) |>
      hot_row("Probability") |>
      hot_row("Days_lost_average") |>
      hot_row("Prob_multiplier")|>
      hot_row("Days_lost_multiplier")
  })
  
  observeEvent(input$tbl$changes$changes, {
    tbl <- hot_to_r(input$tbl)
    
    if (!is.null(tbl)) {
      rv$dat <- tbl
    }
  })
  
  
  
  fit <- reactive({
    
    req(nrow(rv$dat) > 0)
    
    validate(
      need(all(rv$dat$Probability >= 0 & rv$dat$Probability <= 1),
           "Probabilities must be between 0 and 1"),
      need(all(rv$dat$Prob_multiplier >= 0),
           "Multipliers must be larger than 0"),
      need(all(rv$dat$Days_lost_multiplier >= 0),
           "Multipliers must be larger than 0"),
      need(all(rv$dat$Days_lost_average >= 0),
           "Cannot lose negative days"),
      need(!any(is.na(rv$dat)),
           "Missing values are not allowed")
    )
    
    list(fit_hazard_model(rv$dat), fit_hazard_model_adjusted(rv$dat))
  })
  
  output$hazard_plot <- renderPlot({
    
    mod <- fit()[[1]]
    
    tt <- seq(
      1e-6,
      max(b),
      length.out = 1000
    )
    
    plot(
      tt,
      mod$h_hat(tt),
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Hazard",
      main = "Continuous Hazard"
    )
    
    abline(
      v = unique(c(
        a,
        b
      )),
      lty = 2
    )
    
  })
  
  output$compare_hazard_plot <- renderPlot({
    
    mod <- fit()[[1]]
    mod2 <- fit()[[2]]
    
    tt <- seq(
      1e-6,
      max(b),
      length.out = 1000
    )
    
    h_hat1 = mod$h_hat(tt)
    h_hat2 = mod2$h_hat(tt)
    
    par(mar = c(6, 4, 4, 2) + 0.1)
    
    
    plot(
      tt,
      mod$h_hat(tt),
      type = "l",
      ylim = range(h_hat1,h_hat2),
      lwd = 3,
      xlab = "Months",
      ylab = "Hazard",
      main = "Continuous Hazard"
    )
    points( tt,
            mod2$h_hat(tt), col = "red", type = "l", lwd= 3, lty = 2)
    
    
    abline(
      v = unique(c(
        a,
        b
      )),
      lty = 2
    )
    legend(
      "topright",
      legend = c("Baseline hazard", "Compounded hazard"),
      col = c("black", "red"),
      lty = c(1, 2),
      lwd = 3,
      bty = "n"
    )
  })
  output$risk_plot <- renderPlot({
    
    mod <- fit()[[1]]
    
    tt <- seq(
      0,
      max(b),
      length.out = 1000
    )
    
    plot(
      tt,
      mod$Risk_hat(tt),
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Risk",
      main = "Risk Curve"
    )
    
  })
  
  output$compare_risk_plot <- renderPlot({
    
    mod <- fit()[[1]]
    mod2 <- fit()[[2]]
    
    tt <- seq(
      0,
      max(b),
      length.out = 1000
    )
    R_Hat1 =  mod$Risk_hat(tt)
    R_Hat2 =  mod2$Risk_hat(tt)
    
    plot(
      tt,
      R_Hat1,
      ylim = range(R_Hat1, R_Hat2),
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Risk",
      main = "Risk Curve"
    )
    points(
      tt,
      R_Hat2,
      type = "l",
      lwd = 3,
      lty = 2,
      col="red"
    )
    abline(
      v = unique(c(
        a,
        b
      )),
      lty = 2
    )
    
    legend(
      "topright",
      legend = c("Baseline risk", "Compounded risk"),
      col = c("black", "red"),
      lty = c(1, 2),
      lwd = 3,
      bty = "n"
    )
  })
  
  output$density_plot <- renderPlot({
    
    
    
    mod <- fit()[[1]]
    
    
    
    tt <- seq(
      1e-6,
      max(b),
      length.out = 1000
    )
    
    f_dens =  mod$f_hat(tt)
    
    Y_to_h <- function(y) y * (max(f_dens) / max(mod$Y))
    plot(
      tt,
      f_dens,
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Density",
      main = "Injury Density"
    )
    
    # for (i in seq_along(a)) {
    #   if(!is.infinite(b[i])){
    #     segments(a[i], Y_to_h(mod$Y[i]), mod$b[i], Y_to_h(mod$Y[i]),
    #              lwd = 6, col = "blue", lty = 2)
    #   }else{
    #     segments(a[i], Y_to_h(mod$Y[i]), 120, Y_to_h(mod$Y[i]),
    #              lwd = 6, col = "blue", lty = 2)
    #     
    #   }
    #   
    # }
     points(
       tt,
       Y_to_h(mod$Y_func(tt)),
       col="blue",
       lty=1,
       type = "l",
       lwd = 3
     )
    axis(side = 4,
         at = Y_to_h(pretty(c(0, 100))),
         labels = pretty(c(0, 100)),
         col.axis = "blue",
         col = "blue", tick=F)
    mtext(side = 4, "Consequence", cex = 2, col="blue")
    abline(v=c(a,b),lty = 2)
  })
  
  output$compare_density_plot <- renderPlot({
    
    
    
    mod <- fit()[[1]]
    mod2 <- fit()[[2]]
    
    
    
    tt <- seq(
      1e-6,
      max(b),
      length.out = 1000
    )
    
    f_dens1 =  mod$f_hat(tt)
    f_dens2 =  mod2$f_hat(tt)
    
    Y_to_h <- function(y) y * (max(f_dens2, f_dens1) / max(mod2$Y, mod$Y))
    plot(
      tt,
      f_dens1,
      ylim=range(f_dens1, f_dens2),
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Density",
      main = "Injury Density"
    )
    points(
      tt,
      f_dens2,
     col="red",
     lty=2,
      type = "l",
      lwd = 3
    )
    
    points(
      tt,
      Y_to_h(mod$Y_func(tt)),
      col="blue",
      lty=1,
      type = "l",
      lwd = 3
    )
    points(
      tt,
      Y_to_h(mod2$Y_func(tt)),
      col="blue",
      lty=2,
      type = "l",
      lwd = 3
    )
   #  # 
   #  # 
   #  # for (i in seq_along(a)) {
   #  #   if(!is.infinite(b[i])){
   #  #     segments(a[i], Y_to_h(mod$Y[i]), mod$b[i], Y_to_h(mod$Y[i]),
   #  #              lwd = 6, col = "blue", lty = 2)
   #  #   }else{
   #  #     segments(a[i], Y_to_h(mod$Y[i]), 120, Y_to_h(mod$Y[i]),
   #  #              lwd = 6, col = "blue", lty = 2)
   #  #     
   #  #   }
   #  #   
   # # }
    
    axis(side = 4,
         at = Y_to_h(pretty(c(0, max(mod$Y, mod2$Y)))),
         labels = pretty(c(0, max(mod$Y, mod2$Y))),
         col.axis = "blue",
         col = "blue", tick=F)
    mtext(side = 4, "Days lost", cex = 2, col="blue")
    abline(v=c(a,b),lty = 2)
    
    legend(
      "topright",
      legend = c("Baseline density", "Compounded density", "Baseline days lost", "Compounded days lost"),
      col = c("black", "red", "blue", "blue"),
      lty = c(1, 2, 1 , 2),
      lwd = 3,
      bty = "n"
    )
  })
  
  
  output$results_tbl <- renderTable({
    
    mod <- fit()[[1]]
    mod2 <- fit()[[2]]
    
    data.frame(
      Metric = c(
        "Expected days lost",
        "Compounded expected days lost"
      ),
      Value = round(
        c(
          mod$E_unconditional,
          mod2$E_unconditional
        ),
        3
      )
    )
    
  })
  
}

shinyApp(ui, server,  options = list(launch.browser = TRUE))

