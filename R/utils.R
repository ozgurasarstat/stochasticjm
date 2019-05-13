reduce_long <- function(x, frac = 0.1, id){

  idlist   <- x %>% select(id) %>% unique() %>% unlist() %>% as.numeric()
  nsubj    <- idlist %>% length()
  nrepeats <- with(x, tapply(id, id, function(i) table(i))) %>% unlist() %>% as.numeric()

  n_sel <- ceiling(nrepeats * frac)

  out <- c()

  for(i in 1:nsubj){

    x_i <- filter(x, id == idlist[i])
    x_i_sub <- x_i[c(1, sample(2:nrepeats[i], n_sel[i]) %>% sort()), ]

    out <- rbind(out, x_i_sub)

  }

  return(out)

}
