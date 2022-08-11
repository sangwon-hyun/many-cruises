##' Print the alpha table, as a plot, or as a table.
##'
##' @param orig_cluster_labels Ordered clustered labels.
##'
##' @param bestres One flowmix object.
print_alpha <- function(bestres, type = c("plot", "table"), cluster_labels=NULL, orig_cluster_labels=NULL){
  type = match.arg(type)
  alpha = bestres$alpha %>% t()
  numclust = ncol(alpha)
  alpha = alpha[-(rownames(alpha) == "intp"), ]
  if(is.null(cluster_labels)){
    colnames(alpha) = paste0("clust-", (1:numclust) %>% sapply(toString))
  } else {
    colnames(alpha) = cluster_labels
    if(!is.null(orig_cluster_labels)){
      orig_cluster_labels = orig_cluster_labels[which(orig_cluster_labels %in% cluster_labels)]
      alpha = alpha[,orig_cluster_labels]
    }
  }
  max_abs = max(abs(alpha))
  cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(11)
  brk <- lattice::do.breaks(c(-max_abs, max_abs), 11)
  if(type == "plot"){
    p = drawmat_precise_ggplot(alpha, colours = cols) + ggtitle("Alpha coefficients") ##+ xlab("cluster-dimension")
    p = p + scale_fill_gradient2(low = ("red"),
                                 mid = "white",
                                 midpoint = 0,
                                 high = ("blue"))
    p = p + geom_hline(yintercept = 1:1000 - 0.5, col = rgb(0,0,0,0.3))
    p = p + geom_vline(xintercept = 1:1000 - 0.5, col = rgb(0,0,0,0.3))
    return(p)
    drawmat_precise_barebones(alpha, par.settings = list(fontsize = list(text = 20, points = 4)),
                              col.regions = cols,
                              at = brk,
                              colorkey = list(col = cols,
                                              at = brk),
                              xlab = "Cluster",
                              ylab = "Covariate",
                              main = list('Alpha coefficients', side = 1, line = 0.5))
  }
  if(type == "table"){
    tab = alpha %>% process_table(3) %>% kable(caption = "Alpha coefficients")
    return(tab)
  }
}


##' Plots Beta coefficient.
plot_beta <- function(bestres, cluster_labels = NULL, orig_cluster_labels = NULL,
                      cytogram_dimnames = NULL){

  ## Setup
  dimdat = ncol(bestres$beta[[1]])
  numclust = length(bestres$beta)

  ## More setup
  if(is.null(cytogram_dimnames)) cytogram_dimnames = 1:dimdat

  ## Even more setup
  if(is.null(cluster_labels)){
    ordered_cluster_labels = cluster_labels = 1:numclust
  } else {
    if(!is.null(orig_cluster_labels)){
      ordered_cluster_labels = orig_cluster_labels[which(orig_cluster_labels %in% cluster_labels)]
    }
  }

  ## Get the coefficients
  names(bestres$beta) = cluster_labels
  beta = bestres$beta %>% .[ordered_cluster_labels] %>% do.call(cbind,.)
  colnames(beta) = paste0(rep(ordered_cluster_labels, each = dimdat), "-", rep(c(cytogram_dimnames), numclust))
  beta = beta[-1, ]

  ## Formatting for the ggplot
  max_abs = max(abs(beta))
  brk <- lattice::do.breaks(c(-max_abs, max_abs), 11)
  cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(11)

  ## Make the ggplot
  p = drawmat_precise_ggplot(beta, colours = cols) + ggtitle("Beta coefficients") ##+ xlab("cluster-dimension")
  p = p + geom_vline(xintercept = dimdat * (1:numclust) + .5, col = rgb(0,0,0,0.3))
  p = p + geom_hline(yintercept = 1:1000 - 0.5, col = rgb(0,0,0,0.3))
  ## p = p + scale_fill_gradientn(colours = cols, guide="colorbar")
  p = p + scale_fill_gradient2(low = ("red"),
                               mid = "white",
                               midpoint = 0,
                               high = ("blue"))
  return(p)

  ## Plot the coefficients
  drawmat_precise_barebones(beta, par.settings = list(fontsize = list(text = 20, points = 4)),
                            col.regions = cols,
                            at = brk,
                            colorkey = list(col = cols,
                                            at = brk),
                            xlab = "Cluster-dimension",
                            ylab = "Covariate",
                            main = list('Beta coefficients', side = 1, line = 0.5))

}



##' @param bestres flowmix object
print_beta <- function(bestres, cluster_labels = NULL, orig_cluster_labels = NULL, cytogram_dimnames = NULL){

  ## Setup
  if(is.null(cytogram_dimnames)) cytogram_dimnames = 1:dimdat

  ## Even more setup
  if(is.null(cluster_labels)){
    ordered_cluster_labels = cluster_labels = 1:numclust
  } else {
    if(!is.null(orig_cluster_labels)){
      ordered_cluster_labels = orig_cluster_labels[which(orig_cluster_labels %in% cluster_labels)]
    }
  }

  ## Even more setp
  numclust = length(bestres$beta)
  dimdat = ncol(bestres$beta[[1]])
  stopifnot(dimdat == 3)

  ## Format beta coefficient table
  names(bestres$beta) = cluster_labels
  beta = bestres$beta %>% .[ordered_cluster_labels] %>% do.call(cbind,.)
  colnames(beta) = paste0(rep(ordered_cluster_labels, each = dimdat), "-", rep(c(cytogram_dimnames), numclust))

  ## beta = bestres$beta %>% do.call(cbind,.)
  ## colnames(beta) = paste0(rep(1:numclust, each = dimdat), "-", rep(c(cytogram_dimnames), numclust))

  tab = beta %>% process_table(3) %>% knitr::kable(format = "html", caption = "Beta coefficients") %>%
    kableExtra::kable_styling(font_size = 9) %>%
    kableExtra::column_spec(seq(from=1,to=ncol(beta), by=3), border_left = FALSE, border_right = TRUE)##, color=rgb(0,0,0,0.4))
  return(tab)

  ## ## OLD
  ## for(iclust in 1:numclust){
  ##   ## cat("####", paste0("Cluster ", iclust))
  ##   columns = (iclust-1) * 3 + (1:3)
  ##   beta %>% process_table(3) %>% .[,columns] %>% `colnames<-`(cytogram_dimnames) %>%
  ##     kable(caption = paste0("Beta coefficients, cluster ", iclust)) %>% print()
  ## }
  ## ## end of OLD
}


drawmat_precise_barebones <- function (mat, ...)
{
    if (is.null(colnames(mat))) {
        colnames(mat) <- paste(rep("col\n", ncol(mat)), c(1:ncol(mat)),
            sep = " ")
        ownames(mat) <- paste(rep("row", nrow(mat)), c(1:nrow(mat)),
            sep = " ")
    }
    lattice::levelplot(t(mat[c(nrow(mat):1), ]), las = 2, ...)
}



process_table <- function(mat, digit=3){
  apply(mat, 2, function(mycol){
    sapply(mycol, function(entry) if(entry==0) "." else round(entry, digit))
  })
}


plot_prob_ggplot <- function(bestres, cluster_labels = NULL, times = NULL){
  prob = bestres$prob
  numclust = ncol(prob)
  colnames(prob) = paste0("clust-", 1:numclust)
  prob = prob %>% as_tibble()
  if(!is.null(cluster_labels)) colnames(prob) = cluster_labels

  TT = nrow(prob)
  if(is.null(times)) times = 1:TT
  if(!is.null(times)){
    ## stopifnot(all(times == as_datetime(times)))
    times = times %>% lubridate::as_datetime()
  }
  prob = prob %>% add_column(times = times)
  p = prob %>% pivot_longer(-any_of("times")) %>%
    ggplot() +
    ## facet_wrap(~name) +
    ## facet_wrap(~factor(name, levels=paste0("clust-", 1:numclust))) +
    facet_wrap(~name)+
    geom_line(aes(x=times, y=value))
    ## ylim(c(0,1))

  ## p = p + scale_x_datetime(date_breaks = "6 hour", labels = scales::date_format("%b %d - %H:%M"))
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 1.0, hjust = 1.0))
  p = p + ggtitle("Cluster probabilities")
  return(p)

}




##' Takes estimated models from subsamples (these are themselves summaries from
##' cross-validation), and return the problem.
##'
##' @param cruise_id Cruise id e.g. "KM1508"
##' @param outputdir Defauts to \code{outputdir = "~/repos/flowmixapp/data/02-estimate/"}
##' @param outputdir_cytograms Defaults to \code{"~/repos/flowmixapp/data/01-cytograms"}
##' @param dimdat Dimension of data,d efaults to 3.
##' @param nsim Number of subsampling replicates; defaults to 100.
##' @param verbose Whether to print the progress. Defaults to \code{FALSE}.
##' @param save Whether to save results to file. Defaults to \code{FALSE}.
##'
##' @return Lists (two for everything, one for alpha, and another for beta)
##'   containing frequency tables, coefficient values, cluster means and
##'   probabilities, and scale-free statistics. Also returns model from
##'   original data, as well as those from the 100 subsampling experiments.
##'
subsample_summarize <- function(cruise_id,
                                cluster_labels,
                                outputdir = "~/repos/flowmixapp/data/02-estimate",  ## This is defined in the Rmd file.
                                outputdir_cytograms = "~/repos/flowmixapp/data/01-cytograms",  ## This is defined in the Rmd file.
                                dimdat = 3,
                                nsim = 100,
                                verbose = FALSE,
                                save = FALSE){


  ## Setup and basic checks
  summarydir = file.path(outputdir, cruise_id, "summaries")
  filename = paste0(cruise_id, "-grand-summary.RDS")
  if(file.exists(file.path(summarydir, filename)) & verbose){
    print(paste("Subsamples' grand summary for", cruise_id, "already done!"))
    return(NULL)
  } else {
    print(paste("Computing subsamples' grand summary for", cruise_id, " now!"))
  }

  ## Load model
  origres = readRDS(file.path(summarydir, "summary-orig.RDS"))$bestres
  X = origres$X

  ## Also load particle-level data.
  filename = paste0(cruise_id, "-particle-datobj.RDS") ## What is this???
  particledir = file.path(outputdir_cytograms, cruise_id, "data")
  datobj_particle = readRDS(file = file.path(particledir, filename))
  ylist_particle = datobj_particle$ylist %>% purrr::map(as.matrix)

  ######################
  ###  Frequencies #####
  ######################
  res_list = stat_list = beta_list = alpha_list = mn_list = prob_list = list()
  start.time = Sys.time()
  for(isim in 1:nsim){
    ## if(verbose) print_progress(isim, nsim, start.time = start.time)
    if(verbose) printprogress(isim, nsim, start.time = start.time)

    ## Load the data.
    prefix = "summary-subsample-"
    filename = paste0(prefix, isim, ".RDS")
    resfile = file.path(summarydir, filename)
    if(!file.exists(resfile)) next
    ## cvres = readRDS(file = resfile)
    cvres = tryCatch({ readRDS(file = resfile)}, error = function(e){  NULL })
    if(is.null(cvres)){ next }

    ## Take, model (estimated on subsampled data), recalculate on the full data.
    cnames = colnames(cvres$bestres$X)
    newres = predict(cvres$bestres, newx = X[,cnames])
    ## class(newres) = "flowmix"

    ## Reorder the new res.
    ## filename = paste0(cruise_id, "-datobj.RDS") ## What is this???
    ## particledir = file.path(outputdir_cytograms, cruise_id, "data")
    ## datobj = readRDS(file = file.path(particledir, filename))
    ## ylist = datobj$ylist
    newres = newres %>% reorder_kl(origres, ylist_particle, fac = 100, verbose = FALSE)

    ## ## Calculate a scale-free measure
    ## nvar = ncol(X)
    ## cnames = colnames(X)
    ## numclust = newres$numclust
    ## stats_array = array(NA, dim = c(numclust, dimdat, nvar))
    ## for(iclust in 1:numclust){
    ##   if(verbose) print_progress(iclust, numclust, "clusters")
    ##   for(idim in 1:dimdat){
    ##     stats = get_scalefree_stat(newres, iclust, idim)
    ##     stats = stats[cnames]
    ##     ## viz_scalefree_stats(stats)
    ##     stats_array[iclust, idim,] = stats
    ##   }
    ## }

    ## Store the beta and alpha coefficients
    alpha_list[[isim]] = newres$alpha
    beta_list[[isim]] = newres$beta
    mn_list[[isim]] = newres$mn
    prob_list[[isim]] = newres$prob
    res_list[[isim]] = newres

    ## ## Also store the new, scale-less measure.
    ## stat_list[[isim]] = stats_array
  }

  ## Remove empty elments, if any.
  alpha_list = alpha_list %>% purrr::compact()
  beta_list = beta_list %>% purrr::compact()
  mn_list = mn_list %>% purrr::compact()
  nsim = length(mn_list)
  nsim = length(alpha_list)

  threshold <- function(x, tol){abs(x) > tol}
  alpha_nonzero = lapply(alpha_list, threshold, 1E-8) %>% Reduce("+", .)

  ## Combine all betas into (K x (p+1) x d) arrays
  beta_arrays = lapply(beta_list, abind, along = 0)

  ## Calculate frequencies
  beta_freq = beta_arrays %>% purrr::map(. %>% threshold(1E-8)) %>% Reduce("+", .) %>% `/`(nsim)

  ## Combine all alphas
  alpha_freq = alpha_list %>% purrr::map(. %>% threshold(1E-8)) %>% Reduce("+", .) %>% t() %>% `/`(nsim)

  ## Save to file
  obj = list(alpha_freq = alpha_freq,
             beta_freq = beta_freq,
             alpha_list = alpha_list,
             beta_list = beta_list,
             res_list = res_list,
             mn_list = mn_list,
             prob_list = prob_list,
             ## stat_list = stat_list,
             numclust = numclust,
             dimdat = dimdat,
             origres = origres)

  ## Save to file
  if(save){
    filename = paste0(cruise_id, "-grand-summary.RDS")
    saveRDS(obj, file = file.path(summarydir, filename))
    if(verbose) cat("saved to ", file.path(summarydir, filename))
  }
  invisible(obj)
}

##' Helper function  to obtain scale-free statistic.
##' @param newres A flowmix object
##' @param iclust Cluster
##' @param idim Dimension
##'
##' @return A p-length vector of scale-free statistics, whose \code{names()}
##'   show the variables.
get_scalefree_stat <- function(newres, iclust, idim){
  X = newres$X
  bet = newres$beta %>% .[[iclust]] %>%
    .[-1,]  %>%  ## remove intercept
    .[,idim]  ## get one dimension
  full_mean = X %*% bet

  nvar = ncol(X)
  partial_means = sapply(1:nvar, function(ivar){
    onebet = newres$beta %>% .[[iclust]] %>%
      .[-1,] %>% ## remove intercept
      .[,idim] %>%  ## get one dimension
      .[ivar] ## Get one coefficient
    one_partial_mean = onebet * X[,ivar]
    one_partial_mean = one_partial_mean + rnorm(length(one_partial_mean), 0, 1E-8)
  })

  ## Useful plot
  ## cor(cbind(full_mean, partial_means)) %>% drawmat_precise_ggplot(col = c('red', 'white', 'blue'), limits = c(-1,1))

  ## Calculate all the variances between all of the partial means
  varlist = diag(cov(cbind(full_mean, partial_means)))
  stats = varlist[-1] / varlist[1]
  names(stats) = colnames(X)
  return(stats)
}



##' Takes the result from \code{get_scalefree_stat()}, and produce plot.
##'
##' @param stats A numeric vector whose \code{names()} are the variables.
##'
##' @return A ggplot.
viz_scalefree_stat <- function(stats){
  p = data.frame(stats) %>%
    tibble::rownames_to_column(var="variable")  %>% as_tibble() %>%
    ggplot() + geom_point(aes(x = variable, y = stats)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1.0, hjust = 1.0)) +
    ylim(c(0,1))
  return(p)
}




##' Helper function to reformat a numeric vector (into a character vector).
print_sparsify <- function(mycol, digit=3){
  sapply(mycol, function(entry) if(entry==0) "." else round(entry, digit))
}



##' Make a data frame with all the information.
##'
##' Let's design this now; the columnsare
##'
##' 1. Cluster (Pro, Syn, picoek1, ..)
##'
##' 2. Coefficient value.
##'
##' 3. Stability value.
##'
##' 4. type (beta, alpha)
##'
##' 5. dimension (alpha, beta)
##'
##' 6. Covariate name (Si, NO3, ...)
##'
aggregate_info <- function(){

  ## Read in the grand summary
  all_cruise_id = c("KM1508", "KM1512", "KM1513", "KOK1515")## KM1518
  cruise_id = all_cruise_id[1]

  ## Must have a reference model for matching.

  ## Collect all the grand summaries
  all_grand_summaries <- lapply(all_cruise_id, function(cruise_id){
    filename = paste0(cruise_id, "-grand-summary.RDS")
    summarydir = file.path(outputdir, cruise_id, "summaries")
    obj = readRDS(file = file.path(summarydir, filename))
    return(obj)
  })

  ## Reordered grand summary.
  grand_summary_reordered = Map(reorder_grand_summary,
                                all_grand_summaries,
                                cluster_orders)

  ## Aggregate into a single dataset.
  all_alpha_lists = grad_summary_reordered %>% purrr::map(.%>% pluck("alpha_list"))
  all_alpha_lists = all_alpha_lists %>% abind(3) ## Something like this


  return(agg)
}





##' It seems like it's a good idea to wrap the loading + reordering into one function.
##'
##' @param cruise_id Main cruise of interest.
##' @param cluster_labels (NOT USED NOW)Defaults to NULL, otherwise a vector of 10 strings.
##' @param outputdir output directory
##' j
##' @return  Returns the correct model
load_model <- function(cruise_id,
                       outputdir = "/home/sangwonh/repos/flowmixapp/data/02-estimate"){

  ## Load the cruise of interest
  summarydir = file.path(outputdir, cruise_id, "summaries")
  bestres = readRDS(file.path(summarydir, ## "no-wind",
                              "summary-orig.RDS"))$bestres ## temporary
  return(bestres)
}






##' 2d Plotting functionality using ggplot2 (only 2d data).
##'
##' @param datobj_2d A data matrix with X on grid; three columns are assumed;
##'   the first two are the coordinates of the 2d data grid (e.g. "diam" and
##'   "chl"); the thrid column is named "counts".
##'
##' @return A ggplot object.
##'
##' @export
##'
##' @import grDevices
bin_plot_2d <- function(datobj_2d, mn=NULL, sigma=NULL, prob=NULL, labels = NULL, fac = 20, colours = NULL, mn_colours = NULL){

  ## Basic checks
  if(is.null(colours)) colours = c("white", "blue")

  ## Get variable names
  varnames = datobj_2d %>% colnames()
  stopifnot(length(varnames) == 3) ## two data columns, one called "counts"
  stopifnot(varnames[3] == "counts")
  varname1 = varnames[1]
  varname2 = varnames[2]

  ## Information about clusters
  dt = data.frame(x = mn[1,], y = mn[2,], prob = prob * fac)
  if(!is.null(mn)){
    numclust = dim(mn)[2]
    if(is.null(labels)) labels = 1:numclust
    dt = cbind(dt, label = labels)
  }

  p = datobj_2d %>%
    ggplot() +
    theme_minimal() +
    geom_raster(aes(x = !!sym(varname1), y=!!sym(varname2), fill = counts)) +
    scale_fill_gradientn(colours = colours, guide="colorbar")+
    xlim(c(0,8)) + ylim(c(0, 8)) +
    theme(legend.position = "none")


  ## Add model.
  if(!is.null(mn)){

    ## Temporary
      mn_colours = labels %>% sapply(function(label){
        allcols = RColorBrewer::brewer.pal(6, "Set2")
        if(grepl("pro", label)){
          return(allcols[1])
        } else if (grepl("pico", label)){
          return(allcols[2])
        } else if (grepl("croco", label)){
          return(allcols[3])
        }  else if (grepl("syn", label)) {
          return(allcols[4])
        }  else if (grepl("bead", label)) {
          return(allcols[5])
        }  else if (grepl("small", label)) {
          return(allcols[6])
        } else {
          return(grDevices::rgb(0,0,0,0.4))
        }
      })

    for(iclust in 1:numclust){
      el = ellipse::ellipse(x = sigma[iclust,,], centre = mn[,iclust]) %>% as_tibble()
      ## p = p + geom_path(aes(x = x, y = y), data = el, colour = "red", lty = 2,
      ##                   lwd = pmin(prob[iclust] * 5, 0.5))
      p = p + geom_path(aes(x = x, y = y), data = el, colour = mn_colours[iclust], lty = 2,
                        lwd = pmin(prob[iclust] * 8, 0.8))
    }

    ## Add points
    ## p = p + geom_point(aes(x = x, y = y, size = (prob)),
    ##                    data = dt, colour = 'red') +

    p = p + geom_point(aes(x = x, y = y, size = (prob)),
                       data = dt, colour = mn_colours) +

    ## scale_size_identity()
      scale_size_area()


    ## Add labels
    ## p = p + geom_text(aes(x = x, y = y, label = label), size = 5, hjust = 0, vjust = 0, data = dt,
    ##                   fontface = "bold", col='black')

    ## Improve label style
    cex = ifelse(mn_colours==rgb(0,0,0,0.5), rel(3), rel(4))
    p = p + ggrepel::geom_text_repel(aes(x = x, y = y, label = label, point.size = sqrt(prob)),
                                     ## col = "black",
                                     col = mn_colours,
                                     cex = cex,
                                     bg.color = "white",
                                     bg.r = 0.1,
                                     fontface = "bold",
                                     ## point.size = NA,
                                     force_pull   = 5, # do not pull toward data points
                                     data = dt,
                                     seed = 1)
  }

  return(p)
}



##' Plot three 2d panels of data, optionally with a model (from \code{obj}).
##'
##' @param obj flowmix object. If NULL, only data is drawn.
##' @param ylist Data.
##' @param countslist Defaults to NULL.
##' @param obj A flowmix object.
##' @param tt time point of interest, out of 1 through \code{length(ylist)}.
##'
##' @return A grob object containing a 3-panel plot.
##'
##' @export
plot_2d_threepanels <- function(obj = NULL, ## Understandably, data (ylist) might not be in the object.
                                ylist,
                                countslist = NULL, ## The time point of interest, out of 1:TT
                                tt,
                                labels = NULL,
                                plist_return = FALSE,
                                colours = NULL,
                                cruise_id = NULL){

  ## Basic checks
  if(!is.null(obj)) stopifnot("flowmix" %in% class(obj))

  ## Setup
  TT = length(ylist)
  assertthat::assert_that(tt %in% 1:TT)
  mn = sigma = prob = NULL

  ## Scale the biomass (|countslist|) by the total biomass in that cytogram.
  counts_sum = sapply(countslist, sum)
  countslist = lapply(countslist, function(counts)counts/sum(counts))

  ###############################
  ## Make the three data plots ##
  ###############################
  dimslist = list(c(1:2), c(2:3), c(3,1))
  plist = lapply(dimslist, function(dims){

    ## Collapse to 2d
    y = ylist[[tt]]
    counts = countslist[[tt]]
    datobj_2d = collapse_3d_to_2d(y = y, counts = counts, dims = dims) %>% as_tibble()

    ## Extract model data
    ## obj = bestres
    if(!is.null(obj)){
      mn = obj$mn[tt,dims,]
      sigma = obj$sigma[,dims,dims]
      prob = obj$prob[tt,]
    }

    ## Make the heatmap
    p = bin_plot_2d(datobj_2d, mn, sigma, prob, labels, colours = colours)

    return(p)
  })

  if(plist_return) return(plist)

  ###############################
  ## Return a 3 x 1 data panel ##
  ###############################
  mytitle = mydatetime = names(ylist)[tt]
  if(!is.null(cruise_id)) mytitle = paste0(cruise_id,"  ", mydatetime)
  main_text_grob = grid::textGrob(mytitle,
                                  gp = grid::gpar(fontsize = 20, fontface = "plain"),
                                  x = 0, hjust = 0)
  gridExtra::arrangeGrob(plist[[1]], plist[[2]], plist[[3]],
                         ncol = 3,
                         top = main_text_grob)
}
