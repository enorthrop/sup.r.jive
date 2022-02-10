

# Always use package name: r.jive::jive()

# never use par(),
#     Instead use withr::with_par(mfrow=c(1,1), plot)
#     Or in the first line of the function,
#            write withr::local_par(mfrow=c(1,1))
#            so it's exists locally within a function


# Write EVERYTHING in this .R files as a function
#     never do "x <- 1" because of compiling errors


# Setting a seed within a function
#myfunction <- function () {
#  old <- .Random.seed
#  on.exit( { .Random.seed <<- old } )
#  set.seed(2)
#  runif(1)
#}
