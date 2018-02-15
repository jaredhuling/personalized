

convert.cutpoint <- function(cutpoint, benefit.scores)
{
    if (is.character(cutpoint))
    {
        if (cutpoint == "median")
        {
            cutpoint <- median(benefit.scores, na.rm = TRUE)
        } else if (grepl("^quant[0-9]+$", cutpoint))
        {
            matches <- gregexpr('[0-9]+', cutpoint)
            qval    <- as.numeric(regmatches(cutpoint, matches)[[1]])
            if (qval >= 100 | qval < 1) stop("Invalid quantile value for cutpoint.")
            qval <- qval / 100
            cutpoint <- quantile(benefit.scores, probs = qval)
        } else
        {
            stop("Invalid cupoint supplied.")
        }
    } else if (is.numeric(cutpoint))
    {
        cutpoint <- cutpoint[1]
    } else
    {
        stop("Invalid value specified for cutpoint.")
    }
    cutpoint
}
