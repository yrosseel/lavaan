# generic statistical model -- YR 10 july 2012


# super class -- virtual statistical model
lavRefModel <- setRefClass("lavRefModel",

# fields
fields = list(
    npar            = "integer",     # number of free model parameters
    theta           = "numeric",     # the model parameters (free only)
    theta.labels    = "character"    # parameter names (if any)
),

# methods
methods = list(

show = function(header=TRUE) {
    if(header)
        cat(class(.self), "model parameters (theta):\n")
    out <- theta # avoid changing theta by giving names
    if(length(theta.labels) > 0L)
        names(out) <- theta.labels
    print(out)
}

))

