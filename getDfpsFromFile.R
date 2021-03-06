getDfpsFromFiles <- function(idx) {
    res <- list()
    for(i in idx) {
        fn <- paste("paramList", i, ".Rdata", sep = "")
        load(fn)
        res[[i]] <- list(
            skipFactor = paramList$skipFactor,
            zeta = paramList$zeta,
            piVal = paramList$piVal,
            overlapping = paramList$overlapping,
            dfps = paramList$dfps,
            numGenes = nrow(paramList$dfps)
        )
        print(paste(i, "done."))
    }
    return(res)
}

getDfpResults <- function(idx){
    df <- data.frame()
    browser()
    for(i in idx) {
        fn <- paste("paramList", i, ".Rdata", sep = "")
        load(fn)
        df = rbind(df,data.frame(paramList$skipFactor, paramList$zeta, paramList$piVal, paramList$overlapping, nrow(paramList$dfps),fn))
        print(paste(i, "done."))
    }
    colnames(df) <- c("skipFactor","zeta","piVal","overlapping","numGenes","fileName")
    return(df)
}
