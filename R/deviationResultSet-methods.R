

setMethod("show",
          signature="deviationResultSet",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("deviationResults for ", length(object@results), " annotation sets. \n ", sep = "")
            invisible(NULL)
          })

setGeneric("variability", function(object) standardGeneric("variability"))

setMethod("variability", "deviationResultSet",
          function(object){sapply(object@results, function(x) x@variability)})



setGeneric("variability_error", function(object, lower = T) standardGeneric("variability_error"))

setMethod("variability_error", "deviationResultSet",
          function(object, lower = T){
            if (lower){
              out = sapply(object@results, function(x) x@variability_error[1])
            } else{
              out = sapply(object@results, function(x) x@variability_error[2])
            }
            return(out)
          })


setGeneric("plot_variability", function(object, top = NA, ...) standardGeneric("plot_variability"))

setMethod("plot_variability", "deviationResultSet",
          function(object, top = 3, xlab = "Sorted TFs"){

            res_df = data.frame(var = variability(object),
                                #min = variability(object) - variability_error(object, lower = T),
                                #max = variability(object) + variability_error(object, lower = F),
                                min = variability_error(object, lower = T),
                                max = variability_error(object, lower = F),
                                tf = object@tfs,
                                #type = "Data",
                                ranks = rank(-1 * variability(object)))

            out = ggplot2::ggplot(res_df, ggplot2::aes(x = ranks, y = var, min = min, max = max)) + ggplot2::geom_point()+ ggplot2::geom_errorbar() +
              ggplot2::xlab(xlab) + ggplot2::ylab("Normalized Variability") + ggplot2::scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))

            if (!is.na(top)){
              top_df = subset(res_df, ranks <= top)
              out = out + ggplot2::geom_text(data = top_df, ggplot2::aes(x = ranks, y = var, label = tf), size = 2, hjust=-0.15,col = "Black")
            }

            return(out)

          })



setGeneric("adjust_p_values", function(object) standardGeneric("adjust_p_values"))

setMethod("adjust_p_values","deviationResultSet",
          function(object){
            p_values = sapply(object@results, function(x) x@p_value)
            object@adjusted_p_values = p.adjust(p_values,method="BY")
            sim_p_values = sapply(object@results, function(x) x@sim_p_value)
            object@sim_adjusted_p_values = p.adjust(sim_p_values,method="BY")
            return(object)
          })
