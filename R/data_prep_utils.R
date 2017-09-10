

create.design.matrix.binary.trt <- function(x, pi.x, trt, method)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control


    # construct modified design matrices
    # depending on what method is used
    if (method == "weighting")
    {
        # create 1 and -1 version of treatment vector
        trt2 <- 2 * trt - 1

        x.tilde <- trt2 * cbind(1, x)
    } else
    {   # A-learning method
        x.tilde <- (trt - pi.x) * cbind(1, x)
    }
    x.tilde
}

create.weights.binary.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control


    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x * (trt == 1) + (1 - pi.x) * (trt == 0))
    } else
    {   # A-learning method
        wts     <- rep(1, length(pi.x))
    }
    wts
}






create.weights.mult.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as factor

    stop("this doesn't work yet")

    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x * (trt == 1) + (1 - pi.x) * (trt == 0))
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
                                                method = method) )
    }
}

create.weights <- function(pi.x, trt, method)
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

