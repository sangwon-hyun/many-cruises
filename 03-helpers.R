##' Plotting a long matrix with columns "lon", "lat", and "values".
##'
##' @param pred_long Long matrix
##' @param jitter_fac How much to jitter the points.
##'
##' @return a ggplot object.
plot_pipe <- function(pred_long, jitter_fac = 1, col_limits = NULL){

  ## Setup
  latrange = pred_long %>% pull(lat) %>% range(na.rm=TRUE) %>% diff()
  jitter = ggplot2::position_jitter(width = jitter_fac * .1 * latrange/15,
                                    height = jitter_fac * .1 * latrange/15)
  legend_title = "Value"
  world <- ggplot2::map_data("world")

  ## Make plot
  p =
    pred_long %>%
    add_column(size = 1) %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, col = value, size = size), alpha = 1,
               position = jitter) +
    scale_colour_gradientn(colours = c(RColorBrewer::brewer.pal(9, "YlGnBu"), "black"),
                           limits = col_limits,
                           ## oob = scales::censor,
                           oob = scales::squish,
                           na.value = "red") +
    ## scale_size_manual(guide = "none") +
    guides(col = guide_colourbar(title = legend_title)) +
    annotation_map(map_data("world"), fill = "antiquewhite", colour = "darkgrey") +
    ## coord_map("azequalarea", orientation = c(-36.92, 200, 0)) +
    coord_fixed() +
    theme_gray() +
    ylab("Latitude") +
    xlab("Longitude") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "bottom") +
    scale_size_continuous(guide = "none")
    ## scale_size(guide="none")

  ## if(continuous_size) p = p + scale_size_continuous(range = c(1E-3, 2))
  return(p)
}



##' Helper function for collapsing e.g. \code{pro1} and \code{pro1_other} into
##' \code{pro1}.
##' @param pred_long Long data matrix.
##' @param str_subset "pro1"
##' @param str_result "pro1"
##'
##' @return Combines all predictions (relative abundance and means)
##'
##' @export
combine_pop <- function(pred_long, str_subset, str_result){

  pred_long_new = pred_long %>%
    dplyr::mutate(population2 = ifelse(grepl(!!(str_subset), population),
                                      str_result, population))  %>%
    dplyr::select(-population) %>%
    group_by(across(c(-value))) %>%
    dplyr::summarize(sumvalue = sum(value), meanvalue = mean(value)) %>%
    ungroup() %>%
    rename(population = population2) %>%
    mutate(value = case_when(grepl("mu", type) ~ meanvalue,
                             TRUE ~ sumvalue)) %>%
    select(-sumvalue, -meanvalue)
}



## Table with metadata for the cruises
tab = do.call(rbind,
    list(c("FK180310-1", "yes", "no", "lagrangian"),
    c("FK180310-2", "yes", "no", "lagrangian"),
    c("KM1513", "yes", "no", "diel"),
    c("KM1709", "yes", "no", "mesoscale_eddies"),
    c("KM1713", "yes", "no", "gradients"),
    c("KM1508", "yes", "yes", "hot"),
    c("KM1512", "yes", "yes", "hot"),
    c("KM1518", "yes", "yes", "hot"),
    c("KM1601", "yes", "yes", "hot"),
    c("KM1602", "yes", "yes", "hot"),
    c("KM1708", "yes", "yes", "hot"),
    c("KM1717", "yes", "yes", "hot"),
    c("KM1802", "yes", "yes", "hot"),
    c("KOK1515", "yes", "yes", "hot"),
    c("KOK1604", "yes", "yes", "hot"),
    c("KOK1606", "yes", "yes", "gradients"),
    c("KOK1607", "yes", "yes", "hot"),
    c("KOK1608", "yes", "yes", "hot"),
    c("KOK1609", "yes", "yes", "hot"),
    c("KOK1801", "yes", "yes", "hot"),
    c("KOK1803", "yes", "yes", "hot"),
    c("MGL1704", "yes", "yes", "gradients")))
colnames(tab) = c("cruise_id", "main", "subsample", "cruise_type")
tab = as_tibble(tab) %>% select(cruise_id, cruise_type)



##' Dividing and averaging data.
##'
##' @param pred_long Long data matrix.
##' @param cutwidth The size of a lon/lat box to use to divide then average the
##'   data to.
##'
##' @return Data matrix with the resulting averages.
average_by_space <- function(pred_long, cutwidth = 0.1){

  ## Plot setup
  ## jitter = position_jitter(width = .1, height = .1) ## not used now

  ## Make plot
  pred_long_new =
    pred_long %>%
    na.omit() %>%
    mutate(loncut = ggplot2::cut_width(lon, width = cutwidth)) %>%
    mutate(latcut = ggplot2::cut_width(lat, width = cutwidth))  %>%
    ## group_by(population, cruise_id, type,
    ##          loncut, latcut) %>%
    group_by(across(-c(starts_with("value"),
                       starts_with("time"),
                       starts_with("mo"),
                       "lon",
                       "cruise_id",
                       "lat"
                       ))) %>%
    ## Summing and averaging in one swipe, for mu and pi
    dplyr::summarize(## sumvalue = sum(value),
                     value = mean(value),
                     time = mean(time),
                     lon = mean(lon),
                     lat = mean(lat))
  return(pred_long_new)
}



## ######################################################
## ###  cut+averaging data in some lat/lon regularity ###
## ######################################################
## cut2 <- function(x, breaks) {
##   r <- range(x)
##   b <- seq(r[1], r[2], length=2*breaks+1)
##   brk <- b[0:breaks*2+1]
##   mid <- b[1:breaks*2]
##   brk[1] <- brk[1]-0.01
##   k <- cut(x, breaks=brk, labels=FALSE)
##   mid[k]
## }
