lav_lavaan_step16_rotation <- function(lavoptions, lavmodel, lavpartable, lavh1, lavdata,
                                     x, lavvcov, VCOV, lavcache, lavimplied, lavsamplestats) {
  # # # # # # # # # # # 
  # #  16. rotation # # 
  # # # # # # # # # # # 
  # if lavmodel@nefa > 0L and lavoptions$rotation not "none"
  #   store unrotated  solution in partable (column est.unrotated)
  #   rotate lavmodel via lav_model_efa_rotate and overwrite column est in partable 
  #   if lavoptions$se not in none, bootstrap, external, twostep 
  #     if lavoptions$rotation.se == "delta"
  #       re-compute vcov with delta rule (*)
  #       re-compute SE and store them in lavpartable (*)
  #     else if lavoptions$rotation.se == "bordered"
  #       create 'new' partable where the user = 7/77 parameters are free (*)
  # 
  # (*) code too complicated to summarize here
  if ((.hasSlot(lavmodel, "nefa")) && (lavmodel@nefa > 0L) &&
      (lavoptions$rotation != "none")) {
    
    # store unrotated solution in partable
    lavpartable$est.unrotated <- lavpartable$est
    
    # rotate, and create new lavmodel
    if (lavoptions$verbose) {
      cat("rotating EFA factors using rotation method =",
          toupper(lavoptions$rotation), "...")
    }
    x.unrotated <- as.numeric(x)
    lavmodel.unrot <- lavmodel
    lavmodel <- lav_model_efa_rotate(lavmodel   = lavmodel,
                                     lavoptions = lavoptions)
    # overwrite parameters in @ParTable$est
    lavpartable$est <- lav_model_get_parameters(lavmodel = lavmodel,
                                                type = "user", extra = TRUE)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
    
    # VCOV rotated parameters
    if (!lavoptions$se %in% c("none", "bootstrap", "external", "two.step")) {
      if (lavoptions$verbose) {
        cat("computing VCOV for se =", lavoptions$se,
            "and rotation.se =", lavoptions$rotation.se, "...")
      }
      
      # use delta rule to recompute vcov
      if (lavoptions$rotation.se == "delta") {
        # Jacobian
        JAC <- numDeriv::jacobian(func = lav_model_efa_rotate_x,
                                  x = x.unrotated, lavmodel = lavmodel.unrot,
                                  init.rot = lavmodel@H, lavoptions = lavoptions,
                                  type = "user", extra = FALSE,
                                  method.args = list(eps = 0.0050),
                                  method = "simple") # important!
        
        # force VCOV to be pd, before we transform (not very elegant)
        VCOV.in <- lav_matrix_symmetric_force_pd(lavvcov$vcov,
                                                 tol = 1e-10)
        #VCOV.in <- as.matrix(Matrix:::nearPD(x = lavvcov$vcov)$mat)
        
        # apply Delta rule
        VCOV.user <- JAC %*% VCOV.in %*% t(JAC)
        
        # re-compute SE and store them in lavpartable
        tmp <- diag(VCOV.user)
        min.idx <- which(tmp < 0)
        if (length(min.idx) > 0L) {
          tmp[min.idx] <- as.numeric(NA)
        }
        tmp <- sqrt(tmp)
        # catch near-zero SEs  (was ^(1/2) < 0.6)
        zero.idx <- which(tmp < .Machine$double.eps^(1 / 3))
        if (length(zero.idx) > 0L) {
          tmp[zero.idx] <- 0.0
        }
        lavpartable$se <- tmp
      } else if (lavoptions$rotation.se == "bordered") {
        # create 'new' partable where the user = 7/77 parameters are free
        PT.new <- lavpartable
        
        user7.idx <- which(PT.new$user == 7L |
                             PT.new$user == 77L)
        PT.new$free[user7.idx] <- 1L
        PT.new$free[PT.new$free > 0L] <-
          seq_len(sum(PT.new$free > 0L))
        
        # create 'new' lavmodel (where user7/77 parameters are free)
        lavmodel.new <- lav_model(lavpartable = PT.new,
                                  lavoptions = lavoptions,
                                  th.idx     = lavmodel@th.idx)
        lavmodel.new@GLIST    <- lavmodel@GLIST
        lavmodel.new@H        <- lavmodel@H
        lavmodel.new@lv.order <- lavmodel@lv.order
        
        # create 'border' for augmented information matrix
        x.rot <- lav_model_get_parameters(lavmodel.new)
        JAC <- numDeriv::jacobian(func = lav_model_efa_rotate_border_x,
                                  x = x.rot, lavmodel = lavmodel.new,
                                  lavoptions = lavoptions,
                                  lavpartable = lavpartable,
                                  #method.args = list(eps = 0.0005),
                                  #method = "simple")
                                  method = "Richardson")
        # store JAC
        lavmodel.new@ceq.efa.JAC <- JAC
        
        # no other constraints
        if (length(lavmodel@ceq.linear.idx) == 0L &&
            length(lavmodel@ceq.nonlinear.idx) == 0L &&
            length(lavmodel@cin.linear.idx)    == 0L &&
            length(lavmodel@cin.nonlinear.idx) == 0L) {
          lavmodel.new@con.jac <- JAC
          attr(lavmodel.new@con.jac, "inactive.idx") <- integer(0L)
          attr(lavmodel.new@con.jac, "ceq.idx") <- seq_len(nrow(JAC))
          attr(lavmodel.new@con.jac, "cin.idx") <- integer(0L)
          lavmodel.new@con.lambda <- rep(0, nrow(JAC))
          
          # other constraints
        } else {
          inactive.idx <- attr(lavmodel@con.jac, "inactive.idx")
          ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
          cin.idx <- attr(lavmodel@con.jac, "cin.idx")
          lambda <- lavmodel@con.lambda
          nbord <- nrow(JAC)
          
          # recompute con.jac
          if (is.null(body(lavmodel.new@ceq.function))) {
            ceq <- function(x, ...) {
              return(numeric(0))
            }
          } else {
            ceq <- lavmodel.new@ceq.function
          }
          if (is.null(body(lavmodel.new@cin.function))) {
            cin <- function(x, ...) {
              return(numeric(0))
            }
          } else {
            cin <- lavmodel.new@cin.function
          }
          CON.JAC <-
            rbind(JAC,
                  numDeriv::jacobian(ceq, x = x.rot),
                  numDeriv::jacobian(cin, x = x.rot))
          attr(CON.JAC, "cin.idx") <- cin.idx + nbord
          attr(CON.JAC, "ceq.idx") <- c(1:nbord, ceq.idx + nbord)
          attr(CON.JAC, "inactive.idx") <- inactive.idx + nbord
          
          lavmodel.new@con.jac <- CON.JAC
          lavmodel.new@con.lambda <- c(rep(0, nbord), lambda)
        }
        
        # overwrite lavpartable/lavmodel with rotated version
        #lavmodel    <- lavmodel.new
        #lavpartable <- PT.new
        
        # compute VCOV, taking 'rotation constraints' into account
        VCOV <- lav_model_vcov(lavmodel        = lavmodel.new,
                               lavsamplestats  = lavsamplestats,
                               lavoptions      = lavoptions,
                               lavdata         = lavdata,
                               lavpartable     = PT.new,
                               lavcache        = lavcache,
                               lavimplied      = lavimplied,
                               lavh1           = lavh1)
        
        # compute SE and store them in lavpartable
        tmp <- lav_model_vcov_se(lavmodel = lavmodel.new,
                                 lavpartable = PT.new, VCOV = VCOV)
        lavpartable$se <- tmp
        
        # store rotated VCOV
        #tmp.attr <- attributes(VCOV)
        #VCOV1 <- VCOV
        #attributes(VCOV1) <- tmp.attr["dim"]
        #lavvcov <- list(se = tmp,
        #                information = lavoptions$information,
        #                vcov = VCOV1)
      }
      if (lavoptions$verbose) {
        cat(" done.\n")
      }
    } # vcov
  } # efa
  return(list(
    lavpartable = lavpartable,
    lavmodel = lavmodel
  ))
}

