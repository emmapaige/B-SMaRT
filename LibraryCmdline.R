cmdline.integer <- function(...)
# return the value part of a command line option as an integer
{
    s <- cmdline.option(...)
    if (!is.null(s)) return (as.integer(s))
}

cmdline.flag <- function(arg)
# return TRUE iff command line flag was present
# Note: --arg=value are not considered flags
{
    0 != length(grep(paste("^--",arg,"$", sep=""), commandArgs()))
}

cmdline.logical <- function(...)
# return value part of command line option as numeric
{
    s <- cmdline.option(...)
    if (!is.null(s)) return (as.logical(s))
}

cmdline.numeric <- function(...)
# return value part of command line option as numeric
{
    s <- cmdline.option(...)
    if (!is.null(s)) return (as.numeric(s))
}

cmdline.has.option <- function(...)
# return true if option was specified
{
    !is.null(cmdline.option(..., allow.omit=TRUE))
}

cmdline.option <- function(option, default=NULL,
        die.on.fail=TRUE, # deprecated -- use stop.on.fail instead
        stop.on.fail=die.on.fail,
        allow.omit=!stop.on.fail,
        allowed.values=NULL)
# return the value part of a command line option
{
    ca <- grep("=", grep("^--", commandArgs(), value=T), value=T)

    keys <- sub(pattern="=.*", replacement="", ca)
    keys <- sub(keys, pattern="--", replacement="")
    i <- match(option, keys)
    if (is.na(i))
    {
        if (is.null(default))
        {
            if (!allow.omit) stop("Could not find option ", option, "\n")
            return (NULL)
        }
        else
        {
            return (default)
        }
    }
    values <- sub(pattern=".*=", replacement="", ca)

    if (!is.null(allowed.values))
    {
        ok <- values[i] %in% allowed.values
        if (!all(ok))
        {
            stop("Illegal values for option ", option, ": ",
                    paste(sep=" ,", values[i][!ok]), "\n")
        }
    }

    return (values[i])
}

cmdline.strings <- function(option, ...)
# return comma separated values
{
    out <- cmdline.option(option, ...)
    if (is.null(out)) return (NULL)
    unlist(strsplit(out, split=",", perl=TRUE))
}
