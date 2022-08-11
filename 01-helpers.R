##' Read from one vct file (which is one hours' worth of data, at a finer e.g. 3
##' minute level), and fully process this particle-level data. The processing
##' includes (1) taking the relevant columns and high-quality measurements from
##' the VCT files, (2) removing beads (3) removing the edge points, and (3)
##' relabeling dates.
##'
##' @param one_vct_filename File name of a single VCT file.
##'
##' @return One data frame with four columns named (diam, chl, pe, qc).
##' @export
read_from_one_vct_file <- function(one_vct_filename){

  ## Read from the hourly level data file
  df <- arrow::read_parquet(one_vct_filename)

  ## Retain the proper quality (q50) measurements
  df <- df %>% dplyr::filter(q50 == TRUE) %>% dplyr::select(-q50)

  ## Remove the beads
  df = df %>% dplyr::filter(pop_q50 != "beads")

  ## Obtain raw data
  y <- df %>% dplyr::select(diam = contains("diam_mid_q50"),
                            chl = contains("chl"),
                            pe = contains("pe"),
                            qc = contains("Qc_mid_q50")) ##%>% log()
  ## Log transform
  y <- y %>% dplyr::mutate(diam = log(diam),
                           chl = log(chl),
                           pe = log(pe),
                           qc = qc)

  ## Remove edge points
  y <- y %>% process_y(verbose = TRUE) %>% as_tibble()

  ## Get the single hour associated with this file
  date = df %>% dplyr::pull(date) %>% lubridate::ymd_hms() %>% unique() %>%
    lubridate::floor_date(unit = "hours") %>% unique()
  stopifnot(length(date) == 1)

  ## Return the results
  return(list(date = date, tab = y))
}


##' Processes a cytogram (matrix that has the columns chl, pe and diam) by
##' getting rid of points that lie on the boundaries of each axis.
##'
##' @param y Matrix whose rows represent cytogram particles.
##' @param ... Additional arguments.
##'
##' @return Processed matrix.
##'
##' @export
process_y <- function(y, verbose = FALSE){

  ## Basic checks
  if(nrow(y) < 10) return(y)
  stopifnot(all(c("chl", "pe", "diam") %in% colnames(y)))

  ## Setup
  rid.inds = list()

  ## Eliminate boundary points one dimension at a time.
  for(dimname in c("chl", "pe", "diam")){
    print(dimname)
    rid.inds[[dimname]] = identify_boundary_points_1d(y1d = y[,dimname, drop=TRUE],
                                                      qc = y[,"qc", drop=TRUE])
  }
  rid.inds = unique(unlist(rid.inds))

  if(verbose){
    cat(fill = TRUE)
    cat(length(rid.inds), "points on the edge out of", nrow(y), "are deleted", fill=TRUE)
  }
  return(as.matrix(y[-rid.inds,]))
}


##' Takes vector and identifies the ones that are (1) exactly on the min, and
##' (2) near the max.
##'
##' @param y1d vector of points.
##'
##' @return Indices of the points on the boundaries.
identify_boundary_points_1d_no_qc <- function(y1d, shave_prop = 0.005){
  y1d = unlist(y1d)
  rng = max(y1d) - min(y1d)
  rid.min.ind = which(y1d == min(y1d))
  ## rid.max.ind = which(y1d == max(y1d))
  rid.max.ind = which(y1d > max(y1d) - rng * shave_prop )
  rid.ind = c(rid.min.ind, rid.max.ind)
  return(unique(rid.ind))
}


##' Removes data from the top until the largest bin (of a 100-bin weighted
##' histogram) is not too unusual (smaller than mean + 1 std of the other bins).
##'
##' @param y1d 1 dimensional particles
##' @param qc accompanying qc for particles
##'
##' @return Indices to remove.
identify_boundary_points_1d <- function(y1d, qc){
  for(prop in seq(from=0.001, to=0.05, by = 0.001)){
    ind_to_delete = identify_boundary_points_1d_no_qc(y1d, prop)
    y1d_trimmed = y1d[-ind_to_delete]
    qc_trimmed = qc[-ind_to_delete]
    res = plotrix::weighted.hist(y1d_trimmed, qc_trimmed, breaks = 100, plot = FALSE)
    tail_count = res$density[length(res$density)]
    avg = res$density[-length(res$density)] %>% mean()
    std = res$density[-length(res$density)] %>% sd()
    if(tail_count < avg + std) break
  }
  print(paste(prop, "deleted!"))
  return(ind_to_delete)
}



##' Removes data from the top until the largest bin (of a 100-bin weighted
##' histogram) is not too unusual (smaller than mean + 1 std of the other bins).
##'
##' @param y1d 1 dimensional particles
##' @param qc accompanying qc for particles
##'
##' @return Indices to remove.
identify_boundary_points_1d <- function(y1d, qc){
  for(prop in seq(from=0.001, to=0.05, by = 0.001)){
    ind_to_delete = identify_boundary_points_1d_no_qc(y1d, prop)
    y1d_trimmed = y1d[-ind_to_delete]
    qc_trimmed = qc[-ind_to_delete]
    res = plotrix::weighted.hist(y1d_trimmed, qc_trimmed, breaks = 100, plot = FALSE)
    tail_count = res$density[length(res$density)]
    avg = res$density[-length(res$density)] %>% mean()
    std = res$density[-length(res$density)] %>% sd()
    if(tail_count < avg + std) break
  }
  print(paste(prop, "deleted!"))
  return(ind_to_delete)
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



##' 2d Plotting functionality using ggplot2 (only 2d data).
##'
##' @param datobj_2d A single three-column matrix; the first two are the
##'   coordinates of the 2d data grid (e.g. "diam" and "chl"); the thrid column
##'   is named "counts".
##' @param mn Array of all means.
##' @param sigma (numclust x dimdat x dimdat) array.
##' @param colours Colors plotting the binned data (e.g. \code{c("white", "blue")}).
##' @param mn_colours Colors for cluster means.
##' @param xlim X limits; defaults to NULL.
##' @param ylim Y limits; defaults to NULL.
##' @param labels Labels for the clusters.
##' @param fac Magnification factor for the mean points; defaults to 20.
##' @param prob Cluster probabilities.
##'
##' @return A ggplot object.
##'
##' @export
##'
##' @import grDevices
bin_plot_2d <- function(datobj_2d,
                        mn = NULL, sigma = NULL, prob = NULL, labels = NULL,
                        fac = 20, colours = NULL, mn_colours = NULL,
                        xlim = NULL, ylim = NULL){

  counts = x = y = label = NULL ## fixing check()

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
    geom_tile(aes(x = !!sym(varname1), y=!!sym(varname2), fill = counts)) +
    scale_fill_gradientn(colours = colours, guide="colorbar")+
    theme(legend.position = "none")


  if(!is.null(ylim))  p = p + ylim(ylim)
  if(!is.null(xlim))  p = p + xlim(xlim)

  ## Add model.
  if(!is.null(mn)){

    mn_colours =
    labels %>% sapply(function(label){
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
    ## mn_colours = rep(grDevices::rgb(0,0,0,0.4), numclust)

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
