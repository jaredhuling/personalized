
# this creates a block matrix of contrasts where the reference
# treatment is treated as a "control".
# For example, if matrix x is supplied and treatment with three
# levels where the third level is the reference treatment,
# create.block.matrix.mult.trt() will return:
# |x_1   0  |
# |0     x_2|
# |-x_3 -x_3|
#
# where x = (x_1', x_2', x_3')' and x_i is the submatrix
# of x of observations with treatment status == i
#
# with binary treatments, this simplifies to the normal
# | x_1|
# |-x_2|
create.block.matrix.mult.trt <- function(x, trt, y, reference.trt = NULL)
{
    trt.levels      <- levels(trt)
    n.trts          <- length(trt.levels)
    trt.idx         <- vector(mode = "list", length = n.trts)
    sample.sizes    <- numeric(n.trts)

    # if not specified, set reference treatment
    # to be the last one
    if (is.null(reference.trt))
    {
        reference.trt <- trt.levels[n.trts]
    }
    which.reference <- which(trt.levels == reference.trt)

    # re-order treatments so that
    # the specified reference will be treated as the reference
    trt.idx.vec <- 1:n.trts
    trt.idx.vec <- c(trt.idx.vec[-which.reference], trt.idx.vec[which.reference])
    trt.levels  <- trt.levels[trt.idx.vec]

    for (t in 1:n.trts)
    {
        trt.idx[[t]]    <- which(levels(trt)[trt] == trt.levels[t])
        sample.sizes[t] <- length(trt.idx[[t]])
    }
    n.obs  <- NROW(x)
    n.vars <- NCOL(x)

    # if (sd(x[,1]) != 0) stop("x must have intercept as first column")
    stopifnot(n.obs == length(trt))
    stopifnot(n.obs == sum(sample.sizes))

    x.return <- array(0, dim = c(n.obs, n.vars * (n.trts - 1)))

    var.idx.list        <- vector(mode = "list", length = n.trts - 1)
    names(var.idx.list) <- trt.levels[-n.trts] # remove reference treatment
    n.vars.cumsum       <- c(0, cumsum(rep(n.vars, n.trts)))
    n.obs.cumsum        <- c(0, cumsum(sample.sizes))

    y.return <- numeric(n.obs)
    for (t in 1:(n.trts - 1))
    {
        idx.obs.cur  <- (n.obs.cumsum[t] + 1):n.obs.cumsum[t + 1]
        idx.vars.cur <- (n.vars.cumsum[t] + 1):n.vars.cumsum[t + 1]
        x.return[idx.obs.cur, idx.vars.cur] <- x[trt.idx[[t]],]
        y.return[idx.obs.cur] <- y[trt.idx[[t]]]
    }
    t <- n.trts

    for (r in 1:(n.trts - 1))
    {
        # keep observation index static
        idx.obs.cur  <- (n.obs.cumsum[t] + 1):n.obs.cumsum[t + 1]
        # replicate columns
        idx.vars.cur <- (n.vars.cumsum[r] + 1):n.vars.cumsum[r + 1]
        x.return[idx.obs.cur, idx.vars.cur] <- -x[trt.idx[[t]],]

        if (r == 1)
        {
            y.return[idx.obs.cur] <- y[trt.idx[[t]]]
        }
    }
    list(x            = x.return,     # design matrix
         y            = y.return,
         trt.levels   = trt.levels,   # treatment levels (re-ordered based off of reference treatment)
         coef.indices = var.idx.list) # list of indices for the coefficients specific
                                      # to specific treatment contrasts
}
