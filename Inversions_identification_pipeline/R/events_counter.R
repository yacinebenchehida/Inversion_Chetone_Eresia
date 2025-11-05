events_counter <- function(block_annotation){
  list=list()
  event = 1
  events <- ann_block$role
  if(!(any(events %in% c("inversion","translocation","inversion+translocation")))){
      if(any(events %in% c("colinear"))){
        list[[event]]="colinear"
        counted_events <- as.data.frame(do.call(rbind,list))
        colnames(counted_events) = c("event(s)")
        return(counted_events)
        break
      }else{
         return(NULL)
       break
      }
  }
  if(any(events %in% c("inversion"))){
    list[[event]]="inversion"
    event = event + 1
  }
  if(any(events %in% c("inversion+translocation"))){
    list[[event]]="inversion+translocation"
    event = event + 1
  }
  if(any(events %in% c("translocation"))){
    list[[event]]="translocation"
    event = event + 1
  }
  counted_events <- as.data.frame(do.call(rbind,list))
  colnames(counted_events) = c("event(s)")
  return(counted_events)
}
