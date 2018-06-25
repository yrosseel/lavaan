lav_samplestats_step2 <- function(UNI               = NULL,
                                  ov.names          = NULL, # error message only

                                  # polychoric and empty cells
                                  zero.add          = c(0.5, 0.0),
                                  zero.keep.margins = TRUE,
                                  zero.cell.warn    = FALSE,

                                  # keep track of tables with zero cells?
                                  zero.cell.tables  = TRUE,

                                  optim.method      = "nlminb") {

    nvar <- length(UNI)
    COR <- diag(nvar)

    if(zero.cell.tables) {
        zero.var1 <- character(0L)
        zero.var2 <- character(0L)
    }

    # one-by-one (for now)
    for(j in seq_len(nvar-1L)) {
        for(i in (j+1L):nvar) {
            #if(verbose) { cat(" i = ", i, " j = ", j,
            #                  "[",ov.names[i], "-", ov.names[j], "] ",
            #                  "(",ov.types[i], "-", ov.types[j], ")\n") }
            #pstar.idx <- PSTAR[i,j]
            #COR.NAMES[pstar.idx] <- paste(ov.names[i],"~~",ov.names[j],sep="")
            if(class(UNI[[i]]) == "lavOLS" && class(UNI[[j]]) == "lavOLS") {
                if(UNI[[i]]$nexo > 0L) {
                    Y1 <- UNI[[i]]$y - UNI[[i]]$yhat
                    Y2 <- UNI[[j]]$y - UNI[[j]]$yhat
                } else {
                    Y1 <- UNI[[i]]$y; Y2 <- UNI[[j]]$y
                }
                COR[i,j] <- COR[j,i] <- cor(Y1, Y2, use="pairwise.complete.obs")
            } else if(class(UNI[[i]]) == "lavOLS" &&
                      class(UNI[[j]]) == "lavProbit") {
                # polyserial
                out <- ps_cor_TS(fit.y1=UNI[[i]], fit.y2=UNI[[j]])
                COR[i,j] <- COR[j,i] <- out
            } else if(class(UNI[[j]]) == "lavOLS" &&
                      class(UNI[[i]]) == "lavProbit") {
                # polyserial
                out <- ps_cor_TS(fit.y1=UNI[[j]], fit.y2=UNI[[i]])
                COR[i,j] <- COR[j,i] <- out
            } else if(class(UNI[[i]]) == "lavProbit" &&
                      class(UNI[[j]]) == "lavProbit") {
                # polychoric correlation
                out <- pc_cor_TS(fit.y1=UNI[[i]], fit.y2=UNI[[j]],
                                 method = optim.method,
                                 zero.add = zero.add,
                                 zero.keep.margins = zero.keep.margins,
                                 zero.cell.warn = zero.cell.warn,
                                 zero.cell.flag = zero.cell.tables,
                                 Y1.name = ov.names[i],
                                 Y2.name = ov.names[j])
                if(zero.cell.tables) {
                    if(attr(out, "zero.cell.flag")) {
                        zero.var1 <- c(zero.var1, ov.names[j])
                        zero.var2 <- c(zero.var2, ov.names[i])
                    }
                    attr(out, "zero.cell.flag") <- NULL
                }
                COR[i,j] <- COR[j,i] <- out
            }
            # check for near 1.0 correlations
            if(abs(COR[i,j]) > 0.99) {
                warning("lavaan WARNING: correlation between variables ", ov.names[i], " and ", ov.names[j], " is (nearly) 1.0")
            }
        }
    }

    # keep track of tables with zero cells
    if(zero.cell.tables) {
        zero.cell.tables <- cbind(zero.var1, zero.var2)
        attr(COR, "zero.cell.tables") <- zero.cell.tables
    }

    COR
}
