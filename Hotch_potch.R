# This is a file of old R code that I didn't end up using but has some interesting stuff that may be useful


# Grabs the diagonal of a HiC file and rotates it so the diagonal is horizontal, also just takes upper half
HiC.diagonalize <- function(x, width){
	x.horiz <- diag(x)
	x.width <- length(x[1,])
	x.height <- length(x[,1])
	for (i in 1:width){
		new.row <- c(rep(0,i), x[row(x) == col(x) + i])
		x.horiz <- rbind(new.row, x.horiz)
	}
	return(x.horiz)
}