

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

