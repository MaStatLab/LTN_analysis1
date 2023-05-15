
bfdr_d = function(d, pmap){
    # d: binary vector
    # pmap: vector of pmaps
    if (sum(d) == 0) {
        return(0)
    } else {
        return(sum(d * (1 - pmap)) / sum(d))
        }
}

bfdr_c = function(pmap, c){
    # decision rule: d(A) = 1 iff PMAP(A) >=c
    d = as.numeric(pmap >=c)
    return(bfdr_d(d, pmap))
}

bfdr_sort = function(pmap, c){
    if (sum(pmap >=c) == 0){
        return(0)
    }else{
    pmap = sort(pmap, decreasing = T)
    return(1-mean(pmap[1:sum(pmap>=c)]))
    }
}

# bfdr_control = function(pmap, level){
#     pmap = sort(pmap, decreasing = T)
#    rej = NULL
#    cum_bfdr = cumsum(1 - pmap) / (1:length(pmap))
#    return(sum(cum_bfdr < level))
#}


bfdr_control = function(pmap, level){
    # ties
    pmap = sort(unique(pmap), decreasing = T)
    rej = NULL
    cum_bfdr = cumsum(1 - pmap) / (1:length(pmap))
    return(pmap[sum(cum_bfdr <= level)])
}
