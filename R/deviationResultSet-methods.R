

setMethod("show",
          signature="deviationResultSet",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("deviationResults for ", length(object@results), " TFs. \n ", sep = "")
            invisible(NULL)
          })

setGeneric("z_sd", function(object, permuted = F) standardGeneric("z_sd"))

setMethod("z_sd", "deviationResultSet",
          function(object, permuted = F){
            if (permuted){
              out = sapply(object@results, function(x) x@sim_z_sd)
            } else{
              out = sapply(object@results, function(x) x@z_sd)
            }
            return(out)
          })

setGeneric("normalized_variability", function(object, permuted = F) standardGeneric("normalized_variability"))

setMethod("normalized_variability", "deviationResultSet",
          function(object, permuted = F){
            if (permuted){
              out = sapply(object@results, function(x) x@sim_normalized_variability)
            } else{
              out = sapply(object@results, function(x) x@normalized_variability)
            }
            return(out)
          })



setGeneric("error_variability", function(object, permuted = F) standardGeneric("error_variability"))

setMethod("error_variability", "deviationResultSet",
          function(object, permuted = F){
            if (permuted){
              out = sapply(object@results, function(x) x@sim_error_variability)
            } else{
              out = sapply(object@results, function(x) x@error_variability)
            }
            return(out)
          })

setGeneric("error_z_sd", function(object, tail, permuted = F) standardGeneric("error_z_sd"))

setMethod("error_z_sd", "deviationResultSet",
          function(object, tail, permuted = F){
            if (tail=="lower"){
              index = 1
            } else if (tail=="upper"){
              index = 2
            } else{
              stop('argument "tail" must be equal to "lower" or "upper"')
            }
            if (permuted){
              out = sapply(object@results, function(x) x@sim_error_z_sd[index])
            } else{
              out = sapply(object@results, function(x) x@error_z_sd[index])
            }
            return(out)
          })


setGeneric("plot_z_scores", function(object, plot_extra = TRUE, top = NA, ...) standardGeneric("plot_z_scores"))

setMethod("plot_z_scores", "deviationResultSet",
          function(object, plot_extra = TRUE, top = 3){
            
            res_df = data.frame(z = z_sd(object),
                                min = error_z_sd(object, "lower"),
                                max = error_z_sd(object, "upper"),
                                tf = object@tfs,
                                type = "Data",
                                ranks = rank(-1 * z_sd(object)))
            
            if (!is.na(top)){
              top_df = subset(res_df, ranks <= top)
            }
            
            
            if (plot_extra){
              extra_df = data.frame(z = z_sd(object, permuted = T),
                                    min = error_z_sd(object, "lower", permuted = T),
                                    max =  error_z_sd(object, "upper", permuted = T),
                                    tf = object@tfs,
                                    type = "Permuted",
                                    ranks = rank(-1 * z_sd(object, permuted = T)))
              
              res_df = rbind(res_df, extra_df)
              
              
              out = ggplot(res_df, aes(x = ranks, y = z, col = type, ymin =min, ymax = max)) + geom_point() + geom_errorbar() + 
                xlab("Sorted TFs") + ylab("Standard Deviation of Z-scores") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05)) +
                theme_classic(base_size = 12) + theme(legend.title = element_blank(), legend.position = c(1,1), legend.justification = c(1,1))
              
            } else{
              out = ggplot(res_df, aes(x = ranks, y = z, ymin =min, ymax = max)) + geom_point()+  geom_errorbar() + 
                xlab("Sorted TFs") + ylab("Standard Deviation of Z-scores") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))
            }
            
            if (!is.na(top)){
              out = out + geom_text(data = top_df, aes(x = ranks, y = z, label = tf),size = 2, hjust=-0.15,col = "Black")
            }
            
            return(out)
            
          })


setGeneric("plot_variability", function(object, plot_extra = TRUE, top = NA, ...) standardGeneric("plot_variability"))

setMethod("plot_variability", "deviationResultSet",
          function(object, plot_extra = TRUE, top = 3){
            
            res_df = data.frame(var = normalized_variability(object),
                                min = normalized_variability(object) - error_variability(object),
                                max = normalized_variability(object) + error_variability(object),
                                tf = object@tfs,
                                type = "Data",
                                ranks = rank(-1 * normalized_variability(object)))
            
            if (!is.na(top)){
              top_df = subset(res_df, ranks <= top)
            }
            
            
            if (plot_extra){
              extra_df = data.frame(var = normalized_variability(object, permuted = T),
                                    min = normalized_variability(object, permuted = T) - error_variability(object, permuted = T),
                                    max = normalized_variability(object, permuted = T) + error_variability(object, permuted = T),
                                    tf = object@tfs,
                                    type = "Permuted",
                                    ranks = rank(-1 * normalized_variability(object, permuted = T)))
              
              res_df = rbind(res_df, extra_df)
              
              out = ggplot(res_df, aes(x = ranks, y = var, col = type, min = min, max = max)) + geom_point() +  geom_errorbar() + 
                xlab("Sorted TFs") + ylab("Normalized Variability") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))+
                theme_classic(base_size = 12) + theme(legend.title = element_blank(), legend.position = c(1,1), legend.justification = c(1,1))
              
            } else{
              out = ggplot(res_df, aes(x = ranks, y = var, min = min, max = max)) + geom_point()+ geom_errorbar() + 
                xlab("Sorted TFs") + ylab("Normalized Variability") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))
            }
            
            if (!is.na(top)){
              out = out + geom_text(data = top_df, aes(x = ranks, y = var, label = tf), size = 2, hjust=-0.15,col = "Black")
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
