# Function that summarizes which structural event types are present in a block
events_counter <- function(block_annotation){
  
  # initialize an empty list that will store event names
  list = list()
  
  # counter index for inserting into the list
  event = 1
  
  # extract the event types from the input data frame
  events <- block_annotation$role
  
  # --- CASE 1: No event of interest except possibly "colinear" ---
  # Check if none of the major events (inversion/translocation/etc.) are present
  if(!(any(events %in% c("inversion","translocation","inversion+translocation")))){
    
    # If only "colinear" is present, return a one-row result containing "colinear"
    if(any(events %in% c("colinear"))){
      list[[event]] = "colinear"
      counted_events <- as.data.frame(do.call(rbind, list))
      colnames(counted_events) = c("event(s)")
      return(counted_events)   # stop function here
    } else {
      # Otherwise, there is no event worth reporting → return NULL
      return(NULL)             # stop function here
    }
  }
  
  # --- CASE 2: At least one real structural event is present ---
  # Check each event type and record it if present
  
  if(any(events %in% c("inversion"))){
    list[[event]] = "inversion"
    event = event + 1
  }
  
  if(any(events %in% c("inversion+translocation"))){
    list[[event]] = "inversion+translocation"
    event = event + 1
  }
  
  if(any(events %in% c("translocation"))){
    list[[event]] = "translocation"
    event = event + 1
  }
  
  # Convert the collected events to a clean data frame for output
  counted_events <- as.data.frame(do.call(rbind, list))
  colnames(counted_events) = c("event(s)")
  
  # return the summary
  return(counted_events)
}
