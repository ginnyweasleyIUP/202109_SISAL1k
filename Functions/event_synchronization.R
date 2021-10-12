#event synchronization functions

#Extreme event are all the events that are above or below the q/2 (1-q/2) percentiles of a distribution. The value of the event is now unimportant. Only time t\*_x (t\*_y) is considered.

#Temporal threshold \tau defines within what timespan an event is connsidered as synchronous. 
#\tau = max(\delta t^x, min(\delta t\*_x, \delta t\*_y)/2)
#\delta t^x is the mean sampling rate of X, 
#\delta t\*_x are the interevent times in X or Y

#The co-occurence of an event in X and Y is counted and summed as C(X|Y) = \sum \sum J^{xy}_{lm},
#where the counter variable J = 1 (if 0<t\*_x-t\*_y<\tau), = 1/2 if t\*_x-t\*_y = 0, and = 0 otherwise

#The strength Q = (C(X|Y)+C(Y|X))/\sqrt(N_x*N_y), where N_i is the number of events in each TS
#The direction q = (C(X|Y)-C(Y|X))/\sqrt(N_x*N_y)

ident.event <- function(TS = list(time = list(), value = list()), q){
  #q in percent
  q.value_above = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(100-q)]
  q.value_below = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(q)]
  
  return_TS = list(time = TS$time[TS$value<=q.value_below | TS$value>=q.value_above], value = TS$value[TS$value<=q.value_below | TS$value>=q.value_above])
  
  return(return_TS)
}


ident.event_minmax <- function(TS = list(time = list(), value = list()), q){
  #q in percent
  q.value_above = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(100-q)]
  q.value_below = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(q)]
  
  time_above = TS$time[TS$value >= q.value_above]
  time_below = TS$time[TS$value <= q.value_below]
  
  return_TS = list(time_max = time_above, time_min = time_below)
  
  return(return_TS)
}

ident.event_max <- function(TS = list(time = list(), value = list()), q){
  #q in percent
  q.value_above = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(100-q)]
  q.value_below = quantile(TS$value, probs = seq(0,1,length.out = 100), na.rm = T)[floor(q)]
  
  time_above = TS$time[TS$value >= q.value_above]
  time_below = TS$time[TS$value <= q.value_below]
  
  return_TS = list(time = time_above, value = TS$value[TS$value>=q.value_above])
  
  return(return_TS)
}

threshold <- function(TS_time_1, e.time_1, e.time_2){
  return(max(mean(diff(TS_time_1, na.rm = T), na.rm = T), min(c(diff(e.time_1), diff(e.time_2)),na.rm = T)/2))
}
threshold.local <- function(e.time_1, e.time_2, ll, mm){
  if(ll>1 & mm>1){tau = min(abs(e.time_1[ll+1]-e.time_1[ll]), abs(e.time_1[ll]-e.time_1[ll-1]), abs(e.time_2[mm+1]-e.time_2[mm]), abs(e.time_2[mm]-e.time_2[mm-1]),na.rm = T)/2.}
  if(ll == 1 & mm>1){tau = min(abs(e.time_1[ll+1]-e.time_1[ll]), abs(e.time_2[mm+1]-e.time_2[mm]), abs(e.time_2[mm]-e.time_2[mm-1]),na.rm = T)/2.}
  if(ll > 1 & mm == 1){tau = min(abs(e.time_1[ll+1]-e.time_1[ll]), abs(e.time_1[ll]-e.time_1[ll-1]), abs(e.time_2[mm+1]-e.time_2[mm]),na.rm = T)/2.}
  if(ll == 1 & mm == 1){tau = min(abs(e.time_1[ll+1]-e.time_1[ll]), abs(e.time_2[mm+1]-e.time_2[mm]),na.rm = T)/2.}
  #hard threshold as mean uncertainty of dating...
  if(tau>50){tau = 50}
  return(tau)
}

co.occurence <- function(TS_1 = list(time = list(), value = list()), TS_2= list(time = list(), value = list()), q, TS_1_extreme = NULL, TS_2_extreme = NULL, volc = F){
  if(is.null(TS_1_extreme)){
    TS_1_extreme = ident.event(TS_1,q)
  }
  if(is.null(TS_2_extreme)){
    TS_2_extreme = ident.event(TS_2,q)
  }
  if(volc){
    TS_2_extreme = ident.event_max(TS_2,q)
  }
  tau = threshold(TS_1$time, TS_1_extreme$time, TS_2_extreme$time)
  C = 0
  for(ll in 1:length(TS_1_extreme$time)){
    for(mm in 1:length(TS_2_extreme$time)){
      J = 0
      tau = threshold.local(TS_1_extreme$time, TS_2_extreme$time, ll, mm)
      if(is.na(TS_1_extreme$time[ll]) | is.na(TS_2_extreme$time[mm])){next}
      if(0<(TS_1_extreme$time[ll]-TS_2_extreme$time[mm]) & tau > (TS_1_extreme$time[ll]-TS_2_extreme$time[mm])){J = 1}
      if((TS_1_extreme$time[ll]-TS_2_extreme$time[mm]) == 0){J=1/2}
      C = C + J
    }
  }
  return(C)
}

synchro_time <- function(TS_1 = list(time = list(), value = list()), TS_2= list(time = list(), value = list()), q, TS_1_extreme = NULL, TS_2_extreme = NULL, volc = F){
  
  if(volc){TS_2_extreme = ident.event_max(TS = TS_2, q = 5)}
  
  if(is.null(TS_1_extreme)){TS_1_extreme = ident.event(TS_1,q)}
  if(is.null(TS_2_extreme)){TS_2_extreme = ident.event(TS_2,q)}
  list_time = c()
  list_tau = c()
  #tau = threshold(TS_1$time, TS_1_extreme$time, TS_2_extreme$time)
  for(ll in 1:length(TS_1_extreme$time)){
    for(mm in 1:length(TS_2_extreme$time)){
      tau = threshold.local(TS_1_extreme$time, TS_2_extreme$time, ll, mm)
      #tau = 0.5
      if(is.na(TS_1_extreme$time[ll]) | is.na(TS_2_extreme$time[mm])){next}
      if((TS_1_extreme$time[ll]-TS_2_extreme$time[mm]) == 0){
        list_time = c(list_time, TS_1_extreme$time[ll])
        list_tau = c(list_tau, tau)
        next
        }
      if(tau >= abs(TS_1_extreme$time[ll]-TS_2_extreme$time[mm])){
        list_time = c(list_time, TS_1_extreme$time[ll])
        list_tau = c(list_tau, tau)
      }
      
    }
  }
  
  return_list = list(list_time = as.numeric(list_time), tau = as.numeric(list_tau))
  return(return_list)
}

synchro_time_random.cut <- function(TS_1 = list(time = list(), value = list()), TS_2= list(time = list(), value = list()), q, 
                                    TS_1_extreme = NULL, TS_2_extreme = NULL, volc = F){
  
  if(volc){
    random.cut =round(runif(1)*length(TS_1$value-1)+0.50)
    TS_1$value = c(TS_1$value[(random.cut+1):length(TS_1$value)],TS_1$value[1:random.cut])
    TS_2_extreme = ident.event_max(TS = TS_2, q = 5)
  }else{
    random.cut =round(runif(1)*length(TS_2$value-1)+0.50)
    TS_2$value = c(TS_2$value[(random.cut+1):length(TS_2$value)],TS_2$value[1:random.cut])
    
  }
  
  if(is.null(TS_1_extreme)){TS_1_extreme = ident.event(TS_1,q)}
  if(is.null(TS_2_extreme)){TS_2_extreme = ident.event(TS_2,q)}
  list_time = c()
  list_tau = c()
  for(ll in 1:length(TS_1_extreme$time)){
    for(mm in 1:length(TS_2_extreme$time)){
      tau = threshold.local(TS_1_extreme$time, TS_2_extreme$time, ll, mm)
      #tau = 0.5
      if(is.na(TS_1_extreme$time[ll]) | is.na(TS_2_extreme$time[mm])){next}
      if((TS_1_extreme$time[ll]-TS_2_extreme$time[mm]) == 0){
        list_time = c(list_time, TS_1_extreme$time[ll])
        list_tau = c(list_tau, tau)
        next
      }
      if(tau >= abs(TS_1_extreme$time[ll]-TS_2_extreme$time[mm])){
        list_time = c(list_time, TS_1_extreme$time[ll])
        list_tau = c(list_tau, tau)
      }
    }
  }
  return_list = list(list_time = as.numeric(list_time), tau = as.numeric(list_tau))
  return(return_list)
}

ES.strength <- function(TS_1 = list(time = list(), value = list()), 
                        TS_2= list(time = list(), value = list()), q, cf = F, volc = F){
  N_1 = length(ident.event(TS_1,q)$time)
  if(volc){
    N_2 = length(ident.event_max(TS_2,q)$time)
  }else{N_2 = length(ident.event(TS_2,q)$time)}
  

  strength = (co.occurence(TS_1, TS_2, q, volc = volc)+co.occurence(TS_2, TS_1,q, volc = volc))/sqrt(N_1*N_2)
  
  
  if(cf){
    bstrap <- c()
    for (i in 1:1000){
      TS_1_extreme = ident.event(TS_1,q)
      TS_1_extreme$time = sort(sample(TS_1_extreme$time, length(TS_1_extreme$time), replace = T))
      TS_2_extreme = ident.event(TS_2,q)
      TS_2_extreme$time = sort(sample(TS_2_extreme$time, length(TS_2_extreme$time), replace = T))
      bstrap <- c(bstrap,(co.occurence(TS_1, TS_2, q, TS_1_extreme = TS_1_extreme, TS_2_extreme = TS_2_extreme) + 
                            co.occurence(TS_2, TS_1,q, TS_1_extreme = TS_1_extreme, TS_2_extreme = TS_2_extreme))/sqrt(N_1*N_2))
    }
    conf_interval = c(quantile(bstrap,0.05), quantile(bstrap,0.95))
  } else{
    conf_interval = NULL
  }
  
  return(list(strength = strength, conf_interval = conf_interval))
}
ES.direction <- function(TS_1 = list(time = list(), value = list()), TS_2= list(time = list(), value = list()), q, cf = F){
  N_1 = length(ident.event(TS_1,q)$time)
  N_2 = length(ident.event(TS_2,q)$time)
  direction = (co.occurence(TS_1, TS_2, q)-co.occurence(TS_2, TS_1,q))/sqrt(N_1*N_2)
  
  if(cf){
    bstrap <- c()
    for (i in 1:1000){
      TS_1_extreme = ident.event(TS_1,q)
      TS_1_extreme$time = sort(sample(TS_1_extreme$time, length(TS_1_extreme$time), replace = T))
      TS_2_extreme = ident.event(TS_2,q)
      TS_2_extreme$time = sort(sample(TS_2_extreme$time, length(TS_2_extreme$time), replace = T))
      bstrap <- c(bstrap,(co.occurence(TS_1, TS_2, q, TS_1_extreme = TS_1_extreme, TS_2_extreme = TS_2_extreme) - 
                            co.occurence(TS_2, TS_1,q, TS_1_extreme = TS_1_extreme, TS_2_extreme = TS_2_extreme))/sqrt(N_1*N_2))
    }
    conf_interval = c(quantile(bstrap,0.05), quantile(bstrap,0.95))
  } else{
    conf_interval = NULL
  }
  
  return(list(direction = direction, conf_interval = conf_interval))
}

irr_detrend <- function(TS = list(time= list(), value = list())){
  lm_trend = lm(TS$value ~ TS$time)
  value_new = TS$value -lm_trend$coefficients[[2]]*TS$time
  return(list(time = TS$time, value = value_new))
}