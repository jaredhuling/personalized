

create.design.matrix.split <- function(x, pi.x, trt, y, method, reference.trt = NULL)
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

        return( create.design.matrix.mult.trt.split(x             = cbind(1, x),
                                                    pi.x          = pi.x,
                                                    trt           = trt,
                                                    #y             = y,
                                                    method        = method,
                                                    reference.trt = reference.trt) )
    } else
    {
        return( create.design.matrix.binary.trt.split(x      = x,
                                                      pi.x   = pi.x,
                                                      trt    = trt,
                                                      method = method,
                                                      reference.trt = reference.trt) )
    }
}







create.design.matrix.binary.trt.split <- function(x, pi.x, trt, method, reference.trt = NULL)
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
        trt.multiplier <- 2 * (trt != reference.trt) - 1

        x.tilde <- cbind(1, x)
    } else
    {   # A-learning method
        trt.multiplier <- (1 * (trt != reference.trt) - pi.x)
        x.tilde <- cbind(1, x)
    }
    list(x = x.tilde,
         trt.multiplier = trt.multiplier)
}






