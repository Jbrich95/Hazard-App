# app.R

library(shiny)
library(rhandsontable)
library(splines2)

#------------------------------------------------------------
# Default data
#------------------------------------------------------------
a <- c(0, 6, 12, 24, 48, 96)
b <- c(6, 12, 24, 48, 96, 144)


default_data <- data.frame(
  time_period = c(
    "0-6 months",
    "6-12 months",
    "1-2 years",
    "2-4 years",
    "4-8 years",
    "8-12 years"
  ),
  Probability  = c(0.30, 0.10, 0.20, 0.50, 0.75, 0.9),
  Days_lost_average     = c(21, 21, 54, 60, 90, 90),
  Days_lost_lower = c(7, 7, 7, 7, 7, 7),
  Days_lost_upper = c(180, 180, 270, 365, 540, 540),
  Prob_multiplier = c(1, 1, 1, 1, 1, 1),
  Days_lost_multiplier = c(1, 1, 1, 1, 1, 1)
)

knots <- unique(a[-1])
knots <- c(6, 12, 24)
#------------------------------------------------------------
# Model fitting function
#------------------------------------------------------------

fit_hazard_model <- function(dat)
{
  
  
  cond.prob <- dat$Probability
  Y <- dat[["Days_lost_average"]]  
  Y.upper <- dat[["Days_lost_upper"]]  
  Y.lower <- dat[["Days_lost_lower"]]  
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
  
  
  f_hat <- function(t)
    S_hat(t) * h_hat(t)
  
  Y_func = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y[1],Y, Y[length(Y)]), xout = t)$y
  }
  Risk_hat <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func(x))
  E_unconditional <- integrate(function(t)
    Risk_hat(t), 0, max(b))$value
  
  Y_func.upper = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y.upper[1],Y.upper, Y.upper[length(Y.upper)]), xout = t)$y
  }
  Risk_hat.upper <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func.upper(x))
  E_unconditional.upper <- integrate(function(t)
    Risk_hat.upper(t), 0, max(b))$value
  
  Y_func.lower = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y.lower[1],Y.lower, Y.lower[length(Y.lower)]), 
           xout = t)$y
  }
  Risk_hat.lower <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func.lower(x))
  E_unconditional.lower <- integrate(function(t)
    Risk_hat.lower(t), 0, max(b))$value
  
  E_1yr <- integrate(function(t)
    Risk_hat(t), 0,12)$value
  E_2yr <- integrate(function(t)
    Risk_hat(t), 0,24)$value
  E_5yr <- integrate(function(t)
    Risk_hat(t), 0,120)$value
  
  E_1yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,12)$value
  E_2yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,24)$value
  E_5yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,120)$value
  
  E_1yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,12)$value
  E_2yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,24)$value
  E_5yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,120)$value
  
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
    Y_func.upper = Y_func.upper,
    Risk_hat.upper = Risk_hat.upper,
    Y_func.lower= Y_func.lower,
    Risk_hat.lower = Risk_hat.lower,
    fitted_interval = fitted_interval,
    E_unconditional = E_unconditional,
    E_1yr = E_1yr,
    E_2yr = E_2yr,
    E_5yr = E_5yr,
    E_unconditional.upper = E_unconditional.upper,
    E_1yr.upper = E_1yr.upper,
    E_2yr.upper = E_2yr.upper,
    E_5yr.upper = E_5yr.upper,
    E_unconditional.lower = E_unconditional.lower,
    E_1yr.lower = E_1yr.lower,
    E_2yr.lower = E_2yr.lower,
    E_5yr.lower = E_5yr.lower,
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
  Y.upper <- dat[["Days_lost_upper"]]   * dat[["Days_lost_multiplier"]] 
  Y.lower <- dat[["Days_lost_lower"]]  * dat[["Days_lost_multiplier"]] 
  
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
  
  Y_func = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y[1],Y, Y[length(Y)]), xout = t)$y
  }
  Risk_hat <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func(x))
  E_unconditional <- integrate(function(t)
    Risk_hat(t), 0, max(b))$value
  
  Y_func.upper = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y.upper[1],Y.upper, Y.upper[length(Y.upper)]), xout = t)$y
  }
  Risk_hat.upper <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func.upper(x))
  E_unconditional.upper <- integrate(function(t)
    Risk_hat.upper(t), 0, max(b))$value
  
  Y_func.lower = function(t){
    #  Y[max(which(a <= t))]
    approx(c(0,(a+b)/2,max(b)), c(Y.lower[1],Y.lower, Y.lower[length(Y.lower)]), 
           xout = t)$y
  }
  Risk_hat.lower <- function(t)
    f_hat(t)*apply(as.matrix(t), 1, function(x) Y_func.lower(x))
  E_unconditional.lower <- integrate(function(t)
    Risk_hat.lower(t), 0, max(b))$value
  
  E_1yr <- integrate(function(t)
    Risk_hat(t), 0,12)$value
  E_2yr <- integrate(function(t)
    Risk_hat(t), 0,24)$value
  E_5yr <- integrate(function(t)
    Risk_hat(t), 0,120)$value
  
  E_1yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,12)$value
  E_2yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,24)$value
  E_5yr.upper <- integrate(function(t)
    Risk_hat.upper(t), 0,120)$value
  
  E_1yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,12)$value
  E_2yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,24)$value
  E_5yr.lower <- integrate(function(t)
    Risk_hat.lower(t), 0,120)$value
  
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
    Y_func.upper = Y_func.upper,
    Risk_hat.upper = Risk_hat.upper,
    Y_func.lower= Y_func.lower,
    Risk_hat.lower = Risk_hat.lower,
    fitted_interval = fitted_interval,
    E_unconditional = E_unconditional,
    E_1yr = E_1yr,
    E_2yr = E_2yr,
    E_5yr = E_5yr,
    E_unconditional.upper = E_unconditional.upper,
    E_1yr.upper = E_1yr.upper,
    E_2yr.upper = E_2yr.upper,
    E_5yr.upper = E_5yr.upper,
    E_unconditional.lower = E_unconditional.lower,
    E_1yr.lower = E_1yr.lower,
    E_2yr.lower = E_2yr.lower,
    E_5yr.lower = E_5yr.lower,
    fit = fit
  )
}
#------------------------------------------------------------
# UI
#------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Injury Risk (Lateral Meniscectomy)"),
  
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
      hot_row("Days_lost_lower") |>
      hot_row("Days_lost_upper") |>
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
      need(all(rv$dat$Days_lost_upper >= 0),
           "Cannot lose negative days"),
      need(all(rv$dat$Days_lost_lower >= 0),
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
    
    risk      <- mod$Risk_hat(tt)
    risk.upper <- mod$Risk_hat.upper(tt)
    risk.lower <- mod$Risk_hat.lower(tt)
    
    plot(
      tt,
    log = "y",
      risk,
      type = "l",
      ylim = c(min(risk.lower[risk.lower>0]), max(risk.upper)),
      lwd = 3,
      xlab = "Months",
      ylab = "Risk",
      main = "Risk Curve"
    )
    
    polygon(
      c(tt, rev(tt)),
      c(risk.upper, rev(risk.lower)),
      col = adjustcolor("steelblue", alpha.f = 0.25),
      border = NA
    )
    
    lines(
      tt,
      risk,
      lwd = 3,
      col = "black"
    )
    
    lines(
      tt,
      risk.upper,
      lty = 2,
      col = "steelblue"
    )
    
    lines(
      tt,
      risk.lower,
      lty = 2,
      col = "steelblue"
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
      legend = c("Average risk", "Risk range"),
      col = c("black", "steelblue"),
      lty = c(1, 1),
      lwd = 3,
      bty = "n"
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
    
    risk1      <- mod$Risk_hat(tt)
    risk.upper1 <- mod$Risk_hat.upper(tt)
    risk.lower1 <- mod$Risk_hat.lower(tt)
    risk2      <- mod2$Risk_hat(tt)
    risk.upper2 <- mod2$Risk_hat.upper(tt)
    risk.lower2 <- mod2$Risk_hat.lower(tt)
    
    plot(
      tt,
      risk1,
      log = "y",
      ylim = c(min(risk.lower1[risk.lower1 > 0], risk.lower2[risk.lower2 >
                                                               0]),
               max(risk.upper2, risk.upper1)),    
      type = "l",
      lwd = 3,
      xlab = "Months",
      ylab = "Risk",
      main = "Risk Curve"
    )
    points(
      tt,
      risk2,
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
    
    polygon(
      c(tt, rev(tt)),
      c(risk.upper1, rev(risk.lower1)),
      col = adjustcolor("steelblue", alpha.f = 0.25),
      border = NA
    )
    polygon(
      c(tt, rev(tt)),
      c(risk.upper2, rev(risk.lower2)),
      col = adjustcolor("red", alpha.f = 0.25),
      border = NA
    )
    

    lines(
      tt,
      risk1,
      lwd = 3,
      col = "black"
    )
    
    lines(
      tt,
      risk.upper1,
      lty = 2,
      col = "steelblue"
    )
    lines(
      tt,
      risk.lower1,
      lty = 2,
      col = "steelblue"
    )

    
    lines(
      tt,
      risk2,
      lwd = 3,
      lty = 2,
      col = "red"
    )
    
    lines(
      tt,
      risk.upper2,
      lty = 2,
      col = "red"
    )
    lines(
      tt,
      risk.lower2,
      lty = 2,
      col = "red"
    )
    
    
    legend(
      "topright",
      legend = c("Baseline average risk", "Compounded average risk",
                 "Baseline risk range", "Compounded risk range"),
      col = c("black", "red", "steelblue", adjustcolor("red", alpha.f = 0.25)),
      lty = c(1, 2, 1, 1),
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
    mtext(side = 4, "Days lost", cex = 2, col="blue")
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
    mod  <- fit()[[1]]
    mod2 <- fit()[[2]]
    
    
    mean_vals <- c(
      mod$E_1yr,
      mod2$E_1yr,
      NA,
      mod$E_2yr,
      mod2$E_2yr,
      NA,
      mod$E_5yr,
      mod2$E_5yr,
      NA,
      mod$E_unconditional,
      mod2$E_unconditional
    )
    
    lower_vals <- c(
      mod$E_1yr.lower,
      mod2$E_1yr.lower,
      NA,
      mod$E_2yr.lower,
      mod2$E_2yr.lower,
      NA,
      mod$E_5yr.lower,
      mod2$E_5yr.lower,
      NA,
      mod$E_unconditional.lower,
      mod2$E_unconditional.lower
    )
    
    upper_vals <- c(
      mod$E_1yr.upper,
      mod2$E_1yr.upper,
      NA,
      mod$E_2yr.upper,
      mod2$E_2yr.upper,
      NA,
      mod$E_5yr.upper,
      mod2$E_5yr.upper,
      NA,
      mod$E_unconditional.upper,
      mod2$E_unconditional.upper
    )
    
    data.frame(
      Metric = c(
        "Expected days lost (1 year)",
        "Compounded expected days lost (1 year)",
        "",
        "Expected days lost (2 years)",
        "Compounded expected days lost (2 years)",
        "",
        "Expected days lost (5 years)",
        "Compounded expected days lost (5 years)",
        "",
        "Expected days lost (12 years)",
        "Compounded expected days lost (12 years)"
      ),
      
      Value = ifelse(
        is.na(mean_vals),
        "",
        sprintf("%.1f (%.1f - %.1f)", mean_vals, lower_vals, upper_vals)
      ),
      check.names = FALSE
      
    )
    
  })
}
shinyApp(ui, server,  options = list(launch.browser = TRUE))

