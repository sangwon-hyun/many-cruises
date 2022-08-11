##' Median smooth a numeric vector whenever it is too small (smaller than 0) or
##' too large (larger than 2.5).
##'
##' @param vec A numeric vector
##'
##' @return Median smooth each point by the surrounding 10 points.
median_smooth <- function(vec){
  vec_new = vec
  for(ii in 1:length(vec)){
    if(is.na(vec[ii])) next
    if(vec[ii] > 2.5 | vec[ii] < 0){
      vec_new[ii] = median(vec[ii + c(-5,-4,-3,-2,-1,1,2,3,4,5)], na.rm=TRUE)
    }
  }
  return(vec_new)
}


##' Takes data frame |df| that contains |time| and |id| (latter must be a single
##' value), and aggregates it to the hourly level. Also optionally pads it so
##' that missing hours are present.
##'
##' @param df Data frame at a finer time scale than hourly (e.g. 3 minutes).
##'
##' @return Hourly data frame.
coarsen_to_hourly <- function(df, pad = TRUE){

  ## Basic checks
  if(is.null(df)) return(df)
  stopifnot("time" %in% colnames(df))
  stopifnot("id" %in% colnames(df))
  if(nrow(df)==0) return(df)
  stopifnot(length(id)==1)

  ## Aggregate to hourly
  df_hourly = df %>%
    mutate(time = lubridate::as_datetime(time)) %>%
    group_by(hour = cut(time, "1 hour")) %>%
    summarize(across(-c("time", "id"), ~ mean(.x, na.rm = TRUE)), id = unique(id))  %>%
    ungroup() %>%
    mutate(hour = lubridate::as_datetime(hour)) %>%
    rename(time = hour)
  if(pad){
    df_hourly = df_hourly %>% padr::pad() %>% fill(id)
  }
  return(df_hourly)
}




##' Lagging the PAR variable, creating four new variables \code{"lag_3"} through
##' \code{"lag_12"}. No other processing is done.
##'
##' @param df data frame with \code{par} variable.
##'
##' @return df with additional PAR variables.
lag_par <- function(df){

  ## Basic check: is the data hourly?
  per = df %>% pull(time) %>% xts::periodicity()
  stopifnot(per$units == "hours")

  ## Create the new data frames.
  lags = c(3,6,9,12)
  for(lag_amount in lags){
    df =  df %>% mutate( "par_{lag_amount}" :=  lag(par, lag_amount))
  }
  return(df)
}


##' Uses \code{lag_par} to create a data frame with four lagged PAR variables,
##' additionally completed (using nearby +-2 days, using \code{fill_in_par()})
##' so that there's complete PAR data in all the lagged PAR variables in the
##' original time points.
##'
##' @param df Data frame with \code{par} variable}
##'
add_all_lagged_par <- function(df){

  ## Lag the data by 3,6,9, and 12 hours.
  orig_times = df %>% na.omit(par) %>% arrange(time) %>% pull(time)
  df = df %>% arrange(time) %>% lag_par()

  ## Save the original times once.

  ## Use "fill_in_par()" to add lagged and completed data.
  varnames = c("par_3", "par_6", "par_9", "par_12")
  new_dfs = sapply(varnames, function(varname){
    df %>% select(id, time, par = !!(as.name(varname))) %>%
      fill_in_par(fill_all = TRUE)
  }, simplify = FALSE, USE.NAMES = TRUE)


  for(varname in varnames){
    df[,varname] = new_dfs[[varname]] %>% pull(par)
  }

  df = df %>% filter(time %in% orig_times)

  return(df)
}



##' Takes in a single PAR table from a single cruise. Then, fills in PAR by
##' using nearby, +-2 days.
##'
##' @param partable table with par, time, and id
##' @param id ID of cruise
##'
##' @return
fill_in_par <- function(partable, id = NULL,  fill_all = FALSE, fill_none = FALSE, makeplot = FALSE){

  ## Basic check
  stopifnot(all(sort(names(partable)) == c("id", "par", "time")))
  stopifnot(length(unique(partable$id)) == 1)

  ## Temporary
  ## partable = df %>% select(par, time, id)
  ## End of temporary

  ## partable = df %>% select(id, par = par_3, time)

  ## Find missing data proportions, by day
  partables = partable %>% group_by(day = cut(time, "1 day") ) %>% group_split()
  missing = partables %>% purrr::map(.%>% pull(par) %>% is.na() %>% mean()) %>% unlist()

  ## Get the days where PAR is missing, but not by too much
  missing_days = which(0 < missing & missing < 0.8)
  if(fill_all)  missing_days = 1:length(partables)
  if(fill_none) missing_days = c()
  partables_new = partables
  ddmax = length(partables)
  for(dd in missing_days){

    ## Surrounding +-2 days
    surrounding_dd = (dd + c(-2,-1,1,2)) %>% pmax(1) %>% pmin(ddmax) %>% unique()

    ## If the surrounding days aren't absurdly missing
    if(!all(missing[surrounding_dd] > 0.3)){
      replacement = partables[surrounding_dd] %>% bind_rows() %>% mutate(hr = lubridate::hour(time))%>%
        group_by(hr) %>% summarize(avg_par = mean(par, na.rm = TRUE)) %>% ungroup()
    } else {
      next
    }

    ## Format into a single table
    replacement_table = partables[[dd]] %>%
      mutate(hr = lubridate::hour(time)) %>%
      left_join(replacement, by=c("hr")) %>%
      mutate(par_complete = ifelse(is.na(par), avg_par, par)) %>%
      dplyr::select(time, par = par_complete, id, day)
    if(any(is.na(replacement_table$time))) browser()

    partables_new[[dd]] <- replacement_table
  }

  if(!is.null(id) & makeplot){
    g = bind_rows(partables %>% bind_rows() %>% add_column(type = "incomplete"),
              partables_new %>% bind_rows() %>% add_column(type='complete')) %>%
      ggplot()  +
      geom_line(aes(x=time, y=par), data = .%>% filter(type=="complete"), col = rgb(0,0,0,0.1)) +
      geom_point(aes(x=time, y=par), data = .%>% filter(type=="complete"), col='green', cex=.5) +
      geom_point(aes(x=time, y=par), data = .%>% filter(type=="incomplete"), cex=.5) +
      geom_vline(aes(xintercept = time),
                 data = . %>% mutate(hr = lubridate::hour(time)) %>% filter(hr == 0), col=rgb(0,0,0,0.5), lty=2) +
      xlab("") +
      ggtitle(id) +
      theme(plot.title=element_text(size=rel(.6), face="bold"))
    return(g)
    ## ggsave(paste0(id, "-complete-par.png"), height=5, width = 20)
  }

  return(partables_new %>% bind_rows() %>% select(-day))
}




##' Convenience function.
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))

  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}



##' List of datasets that each contain id.
##'
##' @param list_of_datasets List of data frames from cruises. Must have cruise
##'   ID as names (strings).
##'
##' @return ggplot of data completeness.
plot_frequency <- function(list_of_datasets){

  ## Basic check
  ids = names(list_of_datasets)
  assertthat::assert_that(!is.null(ids))
  assertthat::assert_that("character" %in% class(ids))

  ## Make frequency matrix
  freqmat = list_of_datasets %>%
    purrr::map(get_freq) %>%
    ## lapply(., function(freq){
    ##   freq[which(freq == 0)] = NA
    ##   freq
    ## }) %>%
    data.table::rbindlist(use.names = TRUE, fill = TRUE, idcol = "id")

  ## Plot the frequency matrix
  freqmat %>%  tibble::column_to_rownames('id') %>%
    ## select(sort(current_vars())) %>%
    drawmat_precise_ggplot(colours = coul) +
    scale_fill_gradientn(limits = c(0,1), colours = coul) +
    geom_vline(xintercept = seq(from=0, to=100, by=1) + 0.5, col = rgb(0,0,0,0.1)) +
    geom_hline(yintercept = seq(from=0, to=100, by=1) + 0.5, col = rgb(0,0,0,0.1))
}




##' Gets the frequencey of non-NA values.
##'
##' @param df Data frame.
##'
##' @return Completeness of each column.
##'
get_freq <- function(df){
  ## stopifnot("time" %in% colnames(df))
  ## stopifnot("lat" %in% colnames(df))
  ## stopifnot("lon" %in% colnames(df))
  if(nrow(df) == 0){
    df = df %>% select(-one_of("time", "lat", "lon", "id"))##-time,-lat,-lon, -id)
    df[1,] = NA
    return(df)
  } else {
    freq = df %>% select(-one_of("time", "lat", "lon", "id")) %>% summarize_all(is.na) %>% summarize_all(mean)
    freq[which(freq == 1)] = NA
    return(data.frame(1-freq))
  }
}
