

create.design.matrix.binary.trt <- function(x, pi.x, trt, method, reference.trt = NULL)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control


    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    # if not specified, set reference treatment
    # to be the last one
    if (is.null(reference.trt))
    {
        reference.trt <- unique.trts[1L]
    }

    which.reference <- which(unique.trts == reference.trt)

    if (n.trts != 2) stop("two trtment levels only for binary trt design matrix function")

    # construct modified design matrices
    # depending on what method is used
    if (method == "weighting")
    {
        # create 1 and -1 version of treatment vector
        trt2 <- 2 * (trt != reference.trt) - 1

        x.tilde <- trt2 * cbind(1, x)
    } else
    {   # A-learning method
        x.tilde <- (1 * (trt != reference.trt) - pi.x) * cbind(1, x)
    }
    x.tilde
}

create.weights.binary.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control

    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    if (n.trts != 2) stop("two trtment levels only for binary trt weighting function")

    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x * (trt == unique.trts[2L]) + (1 - pi.x) * (trt == unique.trts[1L]))
    } else
    {   # A-learning method
        wts     <- rep(1, length(pi.x))
    }
    wts
}






create.weights.mult.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as factor (actually not anymore!)

    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x)
    } else
    {   # A-learning method
        wts     <- rep(1, length(pi.x))
    }
    wts
}


create.design.matrix <- function(x, pi.x, trt, y, method, reference.trt = NULL)
{
    # check if multiple treatments or not
    if (is.factor(trt))
    {
        n.trts <- length(levels(trt))
    } else
    {
        n.trts <- length(unique(trt))
    }

    is.mult.trt <- n.trts > 2

    if (is.mult.trt)
    {
        # set to factor for multiple trtment trt vector if it isn't already
        if (!is.factor(trt)) trt <- as.factor(trt)

        return( create.design.matrix.mult.trt(x             = cbind(1, x),
                                              pi.x          = pi.x,
                                              trt           = trt,
                                              #y             = y,
                                              method        = method,
                                              reference.trt = reference.trt) )
    } else
    {
        return( create.design.matrix.binary.trt(x      = x,
                                                pi.x   = pi.x,
                                                trt    = trt,
                                                method = method,
                                                reference.trt = reference.trt) )
    }
}

create.weights <- function(pi.x, trt, method)
{
    # check if multiple treatments or not
    if (is.factor(trt))
    {
        n.trts <- length(levels(droplevels(trt)))
    } else
    {
        n.trts <- length(unique(trt))
    }

    is.mult.trt <- n.trts > 2

    if (is.mult.trt)
    {
        # set to factor for multiple trtment trt vector if it isn't already
        if (!is.factor(trt)) trt <- as.factor(trt)

        return( create.weights.mult.trt(pi.x   = pi.x,
                                        trt    = trt,
                                        method = method) )
    } else
    {
        return( create.weights.binary.trt(pi.x   = pi.x,
                                          trt    = trt,
                                          method = method) )
    }
}

