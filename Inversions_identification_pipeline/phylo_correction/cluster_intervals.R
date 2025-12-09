cluster_intervals <- function(left_start, left_end, right_start, right_end) {  # cluster intervals that overlap on both left and right
    
  n <- length(left_start)                                                 # number of intervals
  adj <- matrix(FALSE, n, n)                                              # adjacency matrix for overlap graph
  
  for (i in 1:(n - 1)) {                                                  # loop over first interval index
    for (j in (i + 1):n) {                                                # loop over second interval index
      left_overlap  <- max(left_start[i],  left_start[j])  <              # check strict overlap on left interval
                       min(left_end[i],    left_end[j])
      right_overlap <- max(right_start[i], right_start[j]) <              # check strict overlap on right interval
                       min(right_end[i],   right_end[j])
      if (left_overlap && right_overlap) {                                # require overlap on both sides
        adj[i, j] <- TRUE                                                 # mark i and j as adjacent
        adj[j, i] <- TRUE                                                 # symmetric adjacency
      }
    }
  }
  
  cluster_id <- rep(NA_integer_, n)                                       # cluster id per interval
  current_cluster <- 0L                                                   # cluster counter
  
  for (i in 1:n) {                                                        # breadth first search over adjacency graph
    if (is.na(cluster_id[i])) {                                           # if interval not yet assigned
      current_cluster <- current_cluster + 1L                             # start a new cluster
      queue <- c(i)                                                       # initialize queue with current interval
      cluster_id[i] <- current_cluster                                    # assign cluster id
      while (length(queue) > 0) {                                         # process queue until empty
        k <- queue[1]                                                     # take first element in queue
        queue <- queue[-1]                                                # remove it from queue
        neighbors <- which(adj[k, ] & is.na(cluster_id))                  # unassigned neighbors of k
        if (length(neighbors) > 0) {                                      # if there are neighbors to assign
          cluster_id[neighbors] <- current_cluster                        # assign them to the same cluster
          queue <- c(queue, neighbors)                                    # add neighbors to queue
        }
      }
    }
  }
  
  return(cluster_id)                                                      # return cluster ids for all intervals
}
