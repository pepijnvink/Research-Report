sim_scm <- function(timepoints = 5,
                    burnin = 45,
                    N = 500,
                    ndat = 1,
                    phi,
                    betac,
                    betac2 = NULL,
                    time_beta_change = NULL,
                    psi,
                    intercepts = c(0,0),
                    meanc = c(0,0),
                    varC = c(1,1),
                    seed = NULL
){
  if (!is.null(betac2) & is.null(time_beta_change)){
    cli::cli_abort(c("The function does not know when to change the effect of the confounders.",
                     "i" = "Please provide a time for the new effects."))
  }
  if (is.null(betac2) & is.null(time_beta_change)){
    betac2 <- betac
    time_beta_change <- 1 # time_beta_change will be at t=1
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  tot_timepoints <- burnin+timepoints
  variable.names <- purrr::map(.x = c(-(burnin-1):(timepoints)), 
                               function(x) purrr::map2(.x = c("x", "y"),
                                                       .y = x,
                                                       .f = paste0)) %>% 
    unlist() %>%
    c("ID", ., "C1", "C2")
  datas <- tibble(data = rep(NA, ndat))
  for (i in 1:ndat){
    # dataset i
    tab <- matrix(nrow=N,ncol=2*tot_timepoints)
    
    # values on confounders
    c <- mvrnorm(n = N, mu = meanc, Sigma = diag(varC))
    
    for (j in 1:N){
      # datamatrix for person j
      mat <- matrix(NA,
                    nrow = 2,
                    ncol = tot_timepoints)
      
      # residuals
      error <- mvrnorm(n = tot_timepoints,
                       mu = rep(0, 2),
                       Sigma = diag(psi))
      mat[,1] <- intercepts + betac%*%c[j,] + error[1,]
      # rest of burnin period + first beta
      for (k in 2:(burnin+time_beta_change-1)){
        mat[, k] <- intercepts + phi%*%mat[,k-1] + betac%*%c[j,] + error[k,]
      }
      # second beta
      for (l in (burnin+time_beta_change):tot_timepoints){
        mat[, l] <- intercepts + phi%*%mat[,l-1] + betac2%*%c[j,] + error[l,]
      }
      tab[j, ] <- c(mat)
    }
    # add ID and confounder values
    tab <- cbind(1:N, tab, c) %>%
      as_tibble()
    # nest dataset
    datas[i, ] <- nest(tab)
  }
  # set variable names
  datas$data <- purrr::map(datas$data,
                           purrr::set_names,
                           variable.names)
  return(datas)
}