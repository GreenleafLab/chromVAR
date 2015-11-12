
setMethod("show",
          signature="deviationResult",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Deviations for ", object@tf, " for ",
                length(object@normalized_deviations), " cells/samples. \n", sep = "")
            invisible(NULL)
          })

setMethod("variability", "deviationResult",
          function(object){object@variability})


setGeneric("deviationFDR", function(object, ...) standardGeneric("deviationFDR"))

setMethod("deviationFDR", "deviationResult",
          function(object, returns = c("object","result")){
            returns = match.arg(returns)
            p = pnorm(object@normalized_deviations)
            p = ifelse(p > 0.5, (1-p)*2, p*2)
            fdr = p.adjust(p, method="BH")
            if (returns =="result"){
              return(fdr)
            } else{
              object@fdrs = fdr
              return(object)
            }
            })


setGeneric("differentialSamples", function(object, fdr = 0.1) standardGeneric("differentialSamples"))

setMethod("differentialSamples", "deviationResult",
          function(object, fdr = 0.1){
            fdrs = deviationFDR(object)
            sig = which(fdrs < fdr)
            up = intersect(sig, which(object@normalized_deviations > 0))
            down = intersect(sig, which(object@normalized_deviations < 0))
            return(list(up = up, down = down))
          })








