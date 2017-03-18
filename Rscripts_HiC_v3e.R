#Takes a vector of distances and fits a model for scaling s^-1 (fractal globule) or s^-3/2 (equilibrium globule)
HiC.distance.fit <- function(distances, n){
	x <- distances[distances > 5000 & distances < 1e07] #reduce our space to something reasonable, within 1 Mb
	x.h <- hist(x,breaks=seq(0,1.1e07, 1e04)) #histogram function with 10 kb bins2:(n+1)
	plot(x.h$breaks[2:(n+1)], x.h$counts[1:n], type="l", lwd=3, lty=1, ylab="Counts", xlab="Distance")
	data.m <- data.frame(x.h$breaks[2:(n+1)], x.h$counts[1:n], x.h$breaks[2:(n+1)]^-1, x.h$breaks[2:(n+1)]^-1.5) #build a data frame with s^-1 and s^-3/2 as columns
	colnames(data.m) <- c('X1','X2','X3','X4')
	mod1 <- lm(X2 ~ X3, data = data.m)
	points(x.h$breaks[2:(n+1)],mod1$coefficients[[2]] * (x.h$breaks[2:(n+1)])^-1 + mod1$coefficients[[1]],type="l", col=rgb(0,0,1,0.7), lwd=3,lty=2)
	mod32 <- lm(X2 ~ X4, data = data.m)
	points(x.h$breaks[2:(n+1)],mod32$coefficients[[2]] * (x.h$breaks[2:(n+1)])^-1.5 + mod32$coefficients[[1]],type="l", col=rgb(1,0.5,0,0.7), lwd=3,lty=2)
	text.mod1 <- paste0("p(c) ~ s^-1   (R^2 = ", round(summary(mod1)$r.squared,4), ")")
	text.mod32 <- paste0("p(c) ~ s^-3/2   (R^2 = ", round(summary(mod32)$r.squared,4), ")")
	legend("topright", c(text.mod1, text.mod32), fill = c(rgb(0,0,1,.7), rgb(1,0.5,0,0.7)))
	return(list(mod1, mod32))	
}

# This is my goto. Basic heatmap that I'm using for all Hi-C matrices, implemented using heatmap.2 and a variety of color schemes. 
heatmap.natural <- function(x){
	require(gplots)
	require(RColorBrewer)
	
	yellow.orange.red <- c("white",brewer.pal(9, "YlOrRd"))
	orange.and.blue <- c(brewer.pal(9, "Oranges")[7:1], brewer.pal(9, "Blues")[1:9])
	yellow.blue <- c(brewer.pal(9, "YlOrBr")[9:1], brewer.pal(9, "Blues")[1:9])
	blue.yellow <- yellow.blue[18:1]
	blues <- brewer.pal(9, "Blues")
	
	heatmap.2(x,dendrogram='none', Rowv=FALSE, Colv=FALSE,symm=TRUE,key=FALSE,keysize=0.5,key.title=NA,key.xlab=NA,key.ylab=NA,trace='none',scale='none',labRow=NA,labCol=NA, col=colorRampPalette(blue.yellow)(10000))
}


#Implementation of Lieberman-Aiden's basic coverage normalizaiton for a Hi-C matrix. Each position is simply divided by the sum of its column and the sum of its row. First gets sums of all columns and rows, then adds pseudocounts (0.5) for all zeros, then does normalization.
HiC.matrix.vanilla.normalize <- function(x){
	num.rows <- length(x[,1])
	num.cols <- length(x[1,])
	row.sums <- numeric()
	col.sums <- numeric()
	for (a in 1:num.rows){
		row.sums <- c(row.sums, sum(x[a,]))
	}
	
	for (b in 1:num.cols){
		col.sums <- c(col.sums, sum(x[,b]))
	}
	
	#add pseudocounts
	x[x==0] <- 0.5
			
	x1 <- matrix(nrow = num.rows, ncol = num.cols)
	for (i in 1:num.rows){
		for (j in 1:num.cols){
			R.sum <- row.sums[i]
			C.sum <- col.sums[j]
			x1[i,j] <- x[i,j] * (1/R.sum) * (1/C.sum) * 1E10
		}
	}
	rownames(x1) <- rownames(x)
	colnames(x1) <- colnames(x)
	return(x1)
}

# Take the log of every position of a matrix. Also replaces all zero entries with minimum value (could just use 1, log value is 0)
HiC.matrix.scale.log <- function(x){
	x1 <- as.matrix(x)
	sorted <- sort(x1)
	#find smallest non-zero entry
	min.log = 0
	for (i in 1:length(sorted)){
		if (sorted[i] > 0){
			min.log = log(sorted[i])
			break
		}
	}

	for (i in 1:length(x[1,])){
		for (j in 1:length(x[,1])){
			if(x1[j,i] > 0){
				x1[j,i] <- log(x[j,i])
			}
			else{
				x1[j,i] <- min.log
			}
		}
	}
	return (x1)
}

# "compresses" a matrix by collapsing all values below the supplied bottom percentile to the value at the percentile,
# and the same for the top. Contrast enhancing.
HiC.matrix.compress <- function(x, top.percentile=1, bottom.percentile=0.60){
	top.compress <- quantile(x, top.percentile)
	bottom.compress <- quantile(x, bottom.percentile)
	x1 <- x
	x1[x1 > top.compress] <- top.compress
	x1[x1 < bottom.compress] <- bottom.compress
	return(x1)
}

# Wrapper for basic processing of a binned Hi-C file. Reads in, does vanilla norm, log, compression, hist. equalizing.
HiC.matrix.processfromfile <- function(filename,top,bottom){
	x <- read.matrix(filename)
	x <- HiC.matrix.process.standard(x)
	return(x)
}

# Wrapper for processing a Hi-C matrix that is already loaded. Does vanilla norm, log, compression, hist. equalizing.
HiC.matrix.process.standard <- function(x,top,bottom){	
	x <- HiC.matrix.vanilla.normalize(x)
	x <- HiC.matrix.scale.log(x)
	x <- HiC.matrix.compress(x, top, bottom)
	x <- histogram.equalize(x)
	return(x)
}

# Processes and creates heatmap from a Hi-C binned file. Uses proessing routine above.
HiC.heatmap.plotfromfile <- function(filename, top=1, bottom=0.50, outfile){
	x <- HiC.matrix.processfromfile(filename,top,bottom)
	#jpeg(paste(gsub('.txt','',filename),'_bot', bottom,'.jpeg',sep=''),4000,4000)
	#jpeg(paste(gsub('.txt','',filename),'.jpeg',sep=''),4000,4000)
	jpeg(outfile ,4000,4000)
	heatmap.natural(x)
	dev.off()
}

# One-off code for trying a series of lower compression values for a Hi-C binned folder.
HiC.heatmap.compression.series <- function(filename){
	x <- read.matrix(filename)
	x <- HiC.matrix.vanilla.normalize(x)
	x <- HiC.matrix.scale.log(x)
	for (i in seq(0.25, 0.75, 0.05)){
		x1 <- HiC.matrix.compress(x, 1, i)
		x1 <- histogram.equalize(x1)
		jpeg(paste(gsub('.txt','',filename),'_bot', i,'.jpeg',sep=''),4000,4000)
		heatmap.natural(x1)
		dev.off()
	}
}



# Performs histogram equalization on a matrix.
histogram.equalize <- function(x){
	x1 <- round(x, digits=2)
	sorted <- sort(as.matrix(x1))
	counts <- table(sorted)
	cur.total <- 0
	cdf <- list()
	for (i in 1:length(counts)){
		cur.total <- cur.total + counts[i]
		cdf <- c(cdf,cur.total)
		names(cdf)[i] <- names(counts)[i]
		#if (i %% 100000 == 0){
			#print(i)
		#}
	}
	#return(cdf)
	numrows <- length(x1[,1])
	numcols <- length(x1[1,])
	denom <- (numrows * numcols) - 1
	cdf.min <- cdf[[1]]
	x2 <- matrix(nrow = numrows, ncol=numcols)
	for (j in 1:numrows){
		for (k in 1:numcols){
			#CDF function: CDF - CDF.min / NxM -1 * 255
			x2[j,k] <- (((cdf[[as.character(x1[j,k])]] - cdf.min) / denom) * 255)
		}
	}
	rownames(x2) <- rownames(x)
	colnames(x2) <- colnames(x)
	return(x2)
}


# A bit janky, but lets you click on a HiC heat map of all Drosophila chromosomes adn returns the 
# coordinates. Very hard-coded, only works for Drosophila and for data prepared exactly as I prep it.
# Just a tool for convenience. Left room to build it out for single chromosome arms.
HiC.click.coord <- function(chr.choice='all'){
	x <- locator()
	#if (chr == 'all')
	for (i in x[[1]]){
		#Heatmap isn't plotted 0 to 1, but rather something random to something random, depending on the dimensions of the quartz window. The correction below is for plotting using dev.new(width=9,height=6.5). If using a different sized window, use locator() function to identify the x coordinate of the right and left edge and replace corrections.
		left.correction <- 0.02378131
		right.correction <- 0.9480718
		index <- (i - left.correction) * (1/(right.correction - left.correction))
		index1 <- 0
		index2 <- 0
		chr <- ''
		if (chr.choice == 'all'){
			index <- index * 120381547
			if (index >= 0 & index < 22422828){	
				index1 <- index + 50000
				index2 <- index - 50000
				chr <- 'X'
			}
			if (index >= 22422828 & index < 45434372){
				index1 <- index - 22422828 + 50000
				index2 <- index - 22422828 - 50000
				chr <- '2L'
			}
			if (index >= 45434372 & index < 66581080){
				index1 <- index - 45434372 + 50000
				index2 <- index - 45434372 - 50000
				chr <- '2R'
			}
			if (index >= 66581080 & index < 91124637){
				index1 <- index - 66581080 + 50000
				index2 <- index - 66581080 - 50000
				chr <- '3L'
			}
			if (index >= 91124637 & index < 119029690){
				index1 <- index - 91124637 + 50000
				index2 <- index - 91124637 - 50000
				chr <- '3R'
			}
			if (index >= 119029690 & index < 120381547){
				index1 <- index - 119029690 + 50000
				index2 <- index - 119029690 - 50000
				chr <- '4'
			}
		}
		if (chr.choice == '3R'){
			index <- index * 27905053
			index1 <- index + 50000
			index2 <- index - 50000
			chr <- '3R'
		}
		window <- paste("chr",chr, ':', round(index2, digits=0), '-', round(index1, digits=0), sep = "")
		exact <- paste("chr",chr, ':', round(index, digits=0), sep = "")
		print(paste(exact, window, sep = '     '))
	}
}

#Zeros out everything but some given number of diagonal rows. So a width of 3 would leave the diagonal and the +1 and +2 diagonals, make everything else zero
HiC.zeroOffDiag <- function(x, width){
	x.rows <- length(x[,1])
	x.cols <- length(x[1,])
	x1 <- x
	cols.to.zero.above <- width + 2 # leave two diagonals above the diagonal
	cols.to.zero.below <- 1 # leave two diagonals below the diagonal
	for (i in 1:x.rows){
		if (cols.to.zero.above <= x.cols){
			x1[i,cols.to.zero.above:x.cols] <- 0
		}
		cols.to.zero.above <- cols.to.zero.above + 1
		
		if (i > width + 1){
			x1[i,1:cols.to.zero.below] <- 0
			cols.to.zero.below <- cols.to.zero.below + 1
		}		
	}
	return(x1)
}

#Zeros the diagonal and additional (+1, +2, etc. diagonal) diagonals as specified by width; leaves everything else.
HiC.zeroDiag <- function(x, width){
	x.rows <- length(x[,1])
	x.cols <- length(x[1,])
	x1 <- x
	for (i in 1:x.rows){
		start <- max(1, i - width)
		stop <- min(x.cols, i + width)
		x1[i, start:stop] <- 0
	}	
	return(x1)
}

# Extracts from a binned Hi-C matrix just a single chromosome or arm (intxs with self only)
grab.chr <- function(x, chr){
	hits <- grep(chr,rownames(x))
	return(x[hits,hits])
}

# Extracts single chromosome from a "linear" dataset, i.e. a one column dataframe like DNase data. Requires row name handling because of the stupid one-column issue (subsetting one column dataframe turns it into a vector with no names).
grab.chr.linear <- function(x, chr){
	hits <- grep(chr,rownames(x))
	x1 <- as.data.frame(x[hits,])
	rownames(x1) <- rownames(x)[hits]
	return(x1)
}

# extracts a matrix containing a certain number of columns around the diagonal. So for width=20, each row will be the 20 columns to the left of the diagonal(col=row), the diagonal value, and 20 cols to the right of the diagonal. Sort of straightens the matrix, allows looking at patterns of decay within certain distances.
HiC.matrix.extractMiddle <- function(x, width){
	x1 <- data.frame()
	for(i in 1:nrow(x)){
		if(i > width & i <= (nrow(x) - width)){
			new.row <- x[i,(i - width):(i + width)]
			#print(new.row)
			x1 <- rbind(x1, new.row)
		}
	}
	rownames(x1) <- rownames(x)[(width+1) : (nrow(x) - width)]
	#colnames(x1) <- rownames(x1)
	return(as.matrix(x1))
}

# This is for taking a one-column dataframe, e.g. DNase data, and grabbing the rows to match a Hi-C matrix from which the middles were extracted using HiC.matrix.extractMiddle 
HiC.vector.matchMiddles <- function(x,width){
	x.range <- (width + 1):(nrow(x) - width)
	return(as.data.frame(x[x.range,],row.names = row.names(x)[x.range]))
}

# Normalization for looking at diagonal values, mostly. Simply divides each value by the sum of its row
HiC.matrix.normbyrowonly <- function(x){
	x1 <- x
	for (i in 1:nrow(x)){
		sum <- sum(x[i,])
		if (sum > 0){
			x1[i,] <- x[i,] / sum
		}
		else{
			x1[i,] <- 0	
		}		
	}
	return(as.matrix(x1))
}


# Makes a uniformly decaying Hi-C matrix with a decay profile equal to experimentally supplied data in X. Calculates the average values at various distances from teh diagonal in the experimental data, then assigns those values to every row around the diagonal.
HiC.make.uniform <- function(x,width,size){
	x.middles <- HiC.matrix.extractMiddle(x, width)
	#x.decay <- x.middles[,51:101]
	x.means <- apply(x.middles, MARGIN=2, mean)
	x.new <- matrix(nrow=size,ncol=size)
	for (i in 1:size){
		start <- max((width + 2) - i, 1)
		end <-  min((width + 1), size - i + 1) + width
		start.zeros <- i + start - width - 2
		end.zeros <- size - (i + end - width) + 1
		new.row <- c(rep(0, start.zeros), x.means[start:end], rep(0, end.zeros))
		x.new[i,] <- new.row
	}
	return(x.new)
}

# Makes little blocks of 'doped' values along the diagonal. Makes boxes of size n by n (n defined by size variable) that have all counts increased by multiplying by supplied factor. Boxes alternate, with an altered box followed by an unaltered box.
diagonal.dope <- function(x, size, factor){
	x1 <- x
	for (i in seq(1, nrow(x), 2*size)){
		for (j in 0:(size - 1)){
			for (k in 0:(size - 1)){
				x1[i + j, i + k] <- factor * x[i + j, i + k]
			}
		}
	}
	return(x1)
}

# Makes heatmap plots for all .txt files in a supplied folder. Designed for "local" maps. Data can be logged or histogram equalized depending on preference.
HiC.localPlot.from.folder <- function(folder, outfolder, size, bottom.compress=0, LOG=FALSE,HISTEQ=FALSE){
	files <- list.files(folder)
	for (file in files){
		
		file.stem <- gsub('.txt','',file)
		outfile <- paste(outfolder, '/', file.stem, sep='')
		file.path <- paste(folder, file, sep='/')
		x <- read.matrix(file.path)
		center <- round(nrow(x) / 2, 0)
		range <- (center - size):(center + size)
		x <- x[range, range]
		x <- HiC.matrix.compress(x, 1, bottom.compress)
		if(LOG){
			x <- HiC.matrix.scale.log(x)
			outfile <- paste(outfile, '_log', sep='')
		}
		if (HISTEQ){
			x <- histogram.equalize(x)
			outfile <- paste(outfile, '_histeq', sep='')
		}
		median <- quantile(x, 0.9)[1]
		outfile <- paste(outfile, '.jpeg', sep='')
		jpeg(outfile, width = 4000, height = 4000)
		heatmap.natural(x)
		abline(v=seq(0.1,1,0.1),lty=2, lwd=1.2, col="white") #v is 0.1 to 0.97, h is 0 to 0.9
		abline(h=seq(0,0.9,0.1),lty=2,lwd=1.2, col="white")
		points(0.535,0.46,pch=19,col="green",cex=3)
		dev.off()
	}
}

# For binary ChIP matrix, sorts sequentially on columns and makes a heatmap
chip.binary.sortedheatmap <- function(x, method.choice='euclidean'){
	x <- as.matrix(x)
	x <- x[order(x[,1],x[,2],x[,3],x[,4],x[,5],x[,6], decreasing=TRUE),]
	#plot_ly(z=x, type="heatmap", x = colnames(x))
	chip.binary.heatmap(x)
}

# For binary ChIP data, which is a bunch of rows representing boundary elements and columns representing presence/absece (1/0) of a ChIP factor. Makes 6 heatmaps, one sorted on each column. Currently hard-coded for 6 insulator proteins.
chip.binary.singlysorted <- function(x, stem){
	blues <- brewer.pal(9, "Blues")
	ordering <- c(1:6,1:6)
	for (i in 1:ncol(x)){
		x1 <- x[,ordering[i:(i + 5)]]
		x1 <- x1[order(x1[,1], decreasing=TRUE),]
		jpeg(paste(stem, i,'.jpeg',sep=''),height=3000, width=3000)
		heatmap.2(x1,Colv=NA,Rowv=NA,trace="none", dendrogram="none",labRow=NA,key=FALSE, cexCol=8, margins=c(50,50), col=blues)
		dev.off()
	}
}

# I think this is deprecated...
chip.binary.heatmap <- function(x){
	blues <- brewer.pal(9, "Blues")
	heatmap.2(x,symm=TRUE,Colv=FALSE,Rowv=FALSE, trace="none",dendrogram="none",labRow=NA,key=FALSE,colsep=0:(ncol(x)+1),sepcolor="black",sepwidth=c(0.005,0.005), col=blues,margins=c(50,50),cexCol=8)
}

# Takes a binary matrix of factor occupancy, calculates for every pairwise combination of columns (every factor)
# combination, the probability of finding B bound given A is bound, the probability of finding B bound given
# A is not bound, and the ratio (enrichment)
ChIP.conditional.occupancies <- function(x){
	out <- c(1,1,1,1,1)
	for (i in 1:ncol(x)){
		for (j in 1:ncol(x)){
			if (i != j){
				compare <- x[,i] == x[,j]
				p.j.given.i <- signif(nrow(x[x[,i] == 1 & x[,j] == 1,]) / nrow(x[x[,i] == 1,]), 3)
				p.j.given.not.i <- signif(nrow(x[x[,i] == 0 & x[,j] == 1,]) / nrow(x[x[,i] == 0,]), 3)
				enrichment <- signif(p.j.given.i / p.j.given.not.i, 3)
				out <- rbind(out, c(colnames(x)[i], colnames(x)[j], p.j.given.i, p.j.given.not.i, enrichment))
			}
		}
	}
	return(out)
}

# Calculates the frequency of each unique combination of factor binding (actually counts occurrences of unique rows of any matrix)
chip.factor.combofreq <- function(x){
	library("plyr")
	counts <- count(x)
	return(counts[order(counts$freq, decreasing=TRUE),])
}

# Plots the cumulative profiles for ChIP data around a position. Takes as input a folder that has multiple .txt files that contain two columns: the position and the cumulative count. Makes plots for all files, puts them on a single output file (but comments can be made to return to a separate output per file). Prints as PDF.
chip.metaprofile.multiple.fromfiles <- function(folder, outfile){
	files <- list.files(folder)
	pdf(outfile,100,100)
	rows <- round((length(files[grep('.txt',files)]) / 4),0)
	par(mfcol=c(rows,4)) #comment for single
	par(mai=c(1.5,1.5,1.5,1.5)) #comment for single
	for (file in files){
		if(grepl('.txt',file)){
			path <- paste(folder, file,sep='/')
			x <- read.table(path)
			gene.name <- gsub('.txt','',file)
			plot(x[,1],x[,2],type="l",lwd=10, main=gene.name,xlab="",ylab="",cex.main=15, cex.axis=9,xaxt="n")
			abline(v=0, lty=2,col=2)
		}
	}
	dev.off()
}

#This appears to be the same as above...no idea why there are separate functions.
chip.metaprofile.multiple.singlePlot <- function(folder){
	files <- list.files(folder)
	pdf(outfile,100,100)
	rows <- round((length(files[grep('.txt',files)]) / 4),0)
	par(mfcol=c(rows,4)) #comment for single
	par(mai=c(1.5,1.5,1.5,1.5)) #comment for single
	for (file in files){
		if(grepl('.txt',file)){
			path <- paste(folder, file,sep='/')
			x <- read.table(path)
			gene.name <- gsub('.txt','',file)
			#outfile <- paste(folder, '/',gene.name,'.pdf',sep='') #uncomment for single
			#jpeg(outfile, 1000,100)
			#pdf(outfile) #uncomment for single
			#plot(x[,1],x[,2],type="l",lwd=2.5, main=gene.name,xlab="Distance from Boundary",ylab="Cumulative Enrichment") #uncomment for single
			plot(x[,1],x[,2],type="l",lwd=10, main=gene.name,xlab="",ylab="",cex.main=15)
			abline(v=0, lty=2,col=2)
			#dev.off()
			#return()
		}
	}
	dev.off()
}

# Takes a folder full of text files representing the binned ChIP values around a series of positions (one pos per row) and a width. Makes a heatmap for each txt file, printing a JPEG. Width determines number of bins around middle (width=40 gives 81 columns). Resulting plots are sorted becaues they are processed with chip.heatmaps.process()
chip.heatmaps.folder.all <- function(folder, width){
	files <- list.files(folder)
	for (file in files){
		if(grepl('.txt',file)){
			path <- paste(folder, file,sep='/')
			x <- read.table(path)
			x <- chip.heatmaps.process(x, width,40,40)
			gene.name <- gsub('.txt','',file)
			#pdf(paste(folder,'/',gene.name,'.pdf',sep=''),15,15)
			jpeg(paste(folder,'/',gene.name,'.jpeg',sep=''),4000,4000)
			chip.heatmap(x, width, gene.name)
			dev.off()
		}
	}
}

# Hard-codey script that takes a folder full of text files representing the binned ChIP values around a series of positions (one pos per row), grabs the ones corresponding to the list in genes. Can actually put anything there. Script processes each individual heatmap by compressing top 10%, making all negative values 0, and scaling by normalizing to sum. Combines these individual matrices into a larger matrix, each experiment separated by 50 columns of 0's. Then sorts on each column independently, prints a heatmap for sorting on each column. Pretty slow because these processes are slow. Just added feature that makes a separate file with the sorted  directionality Hi-C data, printed in red
chip.heatmaps.insulators.fromfile <- function(folder, width=80){
	genes <- c('BEAF-32','CP190','CTCF','GAF','mod(mdg4)','Su(Hw)','DNase_S5_rep1')
	#genes <- c('BEAF-32','DNase_S5_rep1')
	profiles <- list()
	big.profile <- data.frame()
	for (i in 1:length(genes)){
		gene <- genes[i]
		path <- paste(folder, '/',gene,'.txt',sep='')
		x <- read.table(path,row.names=1)
		profiles[[i]] <- x
		#process these guys a little
		middle <- round(ncol(x) / 2,0)
		x[x < 0] <- 0
		top.compress <- quantile(x[,middle],0.9)
		x[x > top.compress] <- top.compress
		x <- (100000 * x) / sum(x) # just a scaling function since different datasets have different scales
		if(length(big.profile) > 0){
			big.profile <- data.frame(big.profile,matrix(rep(0,nrow(x) * 50),ncol=50),x)
		}
		else{
			big.profile <- x
		}	
	}
	directionality.path <- paste(folder, '/','nc14_HiC_directionality','.txt',sep='')
	directional <- read.table(directionality.path, row.names=1)
	middle <- round(ncol(directional) / 2,0)
	directional <- directional[,(middle - width):(middle + width)]
	top.compress <- quantile(directional[,middle],0.9)
	bottom.compress <- quantile(directional[,middle],0.1)
	directional[directional > top.compress] <- top.compress
	directional[directional < bottom.compress] <- bottom.compress
	directional <- as.matrix(directional)
	#x <- x + 0.25
	#x[x < 0] <- 0
	#x <- (100000 * x) / sum(x) 
	#big.profile <- data.frame(big.profile,matrix(rep(0,nrow(x) * 50),ncol=50),x)
	#return(big.profile)
	
	for (i in 1:length(profiles)){
		middle <- round(ncol(profiles[[i]]) / 2,0)
		row.sums <- apply(profiles[[i]][,(middle - 40):(middle + 40)], MARGIN=1, sum)
		ordering <- order(row.sums,decreasing=TRUE)
		x1 <- as.matrix(big.profile[ordering,])
		gene <- genes[i]
		jpeg(paste(folder,'/sortedby_', gene, ".jpeg", sep=''),4000,4000)
		chip.heatmap(x1, '')
		dev.off()
		
		jpeg(paste(folder,'/directionality_sortedby_', gene, ".jpeg", sep=''),round(4000 / length(profiles),0),4000,)
		chip.heatmap(directional[ordering,],'',"red")
		dev.off()
		return()
	}
}

# subfunction that takes a single ChIP individual value matrix, compresses top 10%, makes all negative values 0, sorts by the middle 81 columns.
chip.heatmaps.process <- function(x, width, order.pos1, order.pos2){
	middle <- round(ncol(x) / 2,0)
	#print(middle)
	top.compress <- quantile(x[,middle],0.9)
	x1 <- x[,(middle - width):(middle + width)]
	x1 <- as.matrix(x1)
	
	row.sums <- apply(x[,(middle - order.pos1):(middle + order.pos2)], MARGIN=1, sum)
	#row.sums <- apply(x[,(middle - 10):(middle + 10)], MARGIN=1, sum)
	x1 <- x1[order(row.sums,decreasing=TRUE),]
	x1[x1 > top.compress] <- top.compress
	x1[x1 < 0] <- 0
	return(x1)
}

# uses heatmap.2 to print out a heatmap for ChIP values. Variety of colors available.
chip.heatmap <- function(x, title, colour="blue"){
	require(gplots)
	require(RColorBrewer)

	yellow.orange.red <- c("white",brewer.pal(9, "YlOrRd"))
	orange.and.blue <- c(brewer.pal(9, "Oranges")[7:1], brewer.pal(9, "Blues")[1:9])
	yellow.blue <- c(brewer.pal(9, "YlOrBr")[9:1], brewer.pal(9, "Blues")[1:9])
	blue.yellow <- yellow.blue[18:1]
	blues <- brewer.pal(9, "Blues")
	reds <- brewer.pal(9, "Reds")
	
	if (colour[1] == "blue"){ colour <- blues}
	if (colour[1] == "yellow.orange.red"){ colour <- yellow.orange.red}
	if (colour[1] == "red"){ colour <- reds}
	#heatmap.2(x1)
	heatmap.2(x,dendrogram='none', main=title, Rowv=FALSE, Colv=FALSE,symm=TRUE,key=FALSE,keysize=0.5,key.title=NA,key.xlab=NA,key.ylab=NA,trace='none',scale='none',labRow=NA,labCol=NA, col=colorRampPalette(colour)(100))

}

# master script for making decay figure for paper
decay.figure.make <- function(nc14.filename, nc12.filename, width, bin.size,  outfile){
	x <- read.matrix(nc14.filename)
	nc14 <- HiC.matrix.extractMiddles.allchr(x, width) #not sure what width I was using for figures originally...
	x <- read.matrix(nc12.filename)
	nc12 <- HiC.matrix.extractMiddles.allchr(x, width)
	
	nc12 <- HiC.matrix.normbyrowonly(nc12)
	nc14 <- HiC.matrix.normbyrowonly(nc14)
	pdf(outfile)
	HiC.plot.decay(nc14, bin.size, "n", TRUE, rgb(0,0,1,0.1,0.7), -4, 0.3)
	HiC.plot.decay.addMeans(nc14, bin.size, "blue")
	HiC.plot.decay.addMeans(nc12, bin.size, rgb(1,0.4,0,0.7))
	#4.8 is a good spacing for "o" style plotting, on outside of default balls. 3.2 is good for "l"
	HiC.plot.decay(nc14, bin.size, "p", FALSE, rgb(0,0,1,0.025), -4.6, 1)
	HiC.plot.decay(nc12, bin.size, "p", FALSE,rgb(1,0.4,0,0.025), 4.6, 1)
	dev.off()
}

# For making overall decay plots. It's important not to have the ends of chromosomes in the decay data, so this uses the extract middle routine which cuts off the ends that don't have sufficient columns to left or right to match width, so a width of 40 will take only the middle of the chromosomes with the 40 bins on the ends left off.
HiC.matrix.extractMiddles.allchr <- function(x, width){
	x1 <- rbind(
	HiC.matrix.extractMiddle(grab.chr(x,'X'),width), 
	HiC.matrix.extractMiddle(grab.chr(x,'2L'),width), 
	HiC.matrix.extractMiddle(grab.chr(x,'2R'),width), 
	HiC.matrix.extractMiddle(grab.chr(x,'3L'),width), 
	HiC.matrix.extractMiddle(grab.chr(x,'3R'),width)
	)
	return(x1)
}

# Takes a matrix of values (design is for them to be centered on diagonal values) and plots each row independently. Recommended colors rgb(0,0,1,0.025) and rgb(1,0.6,0,0.025). Bin size is in kilobases.
HiC.plot.decay <- function(x, bin.size=100, plot.type="n",new=TRUE,color, offset=0, expansion=0.75){
	max <- max(x)
	width <- ncol(x)
	height <- nrow(x)
	x.vals <- seq((-1 * bin.size) *((width - 1) / 2), bin.size *((width - 1) / 2), bin.size)
	#print(x.vals)
	x.vals <- x.vals + offset
	#return(x.vals)
	if(new){
		plot(x.vals,x.vals, type="n", ylim = c(0,0.7), xlab = "", ylab = "",cex.axis=2)
	}
	for (i in 1:height){
		points(x.vals, x[i,], type = plot.type,pch=16,cex=expansion, lwd=1,col=color)
	}
	
}

# Overlay means on the decay plots, as a line currently
HiC.plot.decay.addMeans <- function(x, bin.size, color){
	width <- ncol(x)
	height <- nrow(x)
	x.vals <- seq((-1 * bin.size) *((width - 1) / 2), bin.size *((width - 1) / 2), bin.size)
	means <- apply(x, MARGIN=2, mean)
	points(x.vals, means, col=color, type="l", lty=1,lwd=2.5)
}

# Shortcut to get my "cool region" at the 5' end of 3R that I'm using as a showcase. Input file is a binned matrix file.
cool.region.get <- function(filename){
	x <- read.matrix(filename)
	x <- grab.chr(x, '3R')
	x <- HiC.matrix.vanilla.normalize(x)
	x <- HiC.matrix.scale.log(x)
	x <- x[55:390,55:390]
	x <- histogram.equalize(x)
	return(x)
}
# Shortcut to make a heatmap of the cool region. Compression of 0.2 generally looks good.
cool.region.make.heatmap <- function(filename, outfile, compress=0.2){
	x <- cool.region.get(filename)
	x <- HiC.matrix.compress(x, 1, compress)
	jpeg(outfile, 4000, 4000)
	heatmap.natural(x)
	dev.off()
}

# Makes a plot of an "epigenomic" dataset (e.g., ChIP or DNase) that matches a region of Hi-C data. Hi-C matri supplied in hic and epigenomic file in epifile must have row names in the format of "chr_bin". The epigenomic file is normalized (subtract minimum, divide by max) and can optionally be expanded or contracted for visual purposes by raising to exponent. This probably shouldn't be used but some of the tracks are more "even" than others and it makes viewing them together difficult. As long as this change is noted, I think it's ok. Just changes the visuals, and for this purposes we're just trying to show where a signal is high or low so I think it's ok.
epigenomic.match.plot <- function(hic, epifile, outfile, exponent=1){
	epi <- read.table(epifile, row.names=1)
	epi <- epi[row.names(hic),1]
	#epi <- epi / mean(epi)
	epi <- epi - min(epi)
	epi <- epi  / max(epi)
	epi <- epi^exponent
	pdf(outfile,25,6)
	par(yaxp=c(0,10,10))
	plot(epi, type="l", bty="n",xaxt="n",yaxt="n",xlab="",ylab="",lwd=10)
	axis(side=2, col="black", at=c(0,1),labels=c(0,1))
	dev.off()
}

# Hard-codey routine to make a distance decay plot for a Hi-C matrix with the bins sorted by DNase (or other epigenomic) data. Currently written as plotting lower, middle and upper thirds separately but could easily be tweaked.
decay.plot.sort.dnase <- function(hic, epifile, bin.size, width, outfile){
	hic <- HiC.matrix.extractMiddles.allchr(hic, width)
	epi <- read.table(epifile, row.names=1)
	namelogical <- rownames(epi) %in% rownames(hic)
	epi <- data.frame(epi[namelogical,], row.names=rownames(epi)[namelogical])
	percentile66 <- quantile(epi[,1], 0.6666)
	percentile33 <- quantile(epi[,1], 0.33333)
	upper.third <- rownames(epi)[epi[,1] > percentile66]
	middle.third <- rownames(epi)[epi[,1] < percentile66 & epi[,1] > percentile33]
	bottom.third <- rownames(epi)[epi[,1] < percentile33]
	#return(hic[upper.third,]) # issue is that some of the middle third names aren't in  hic
	
	upper <- HiC.matrix.normbyrowonly(hic[upper.third,])
	middle <- HiC.matrix.normbyrowonly(hic[middle.third,])
	bottom <- HiC.matrix.normbyrowonly(hic[bottom.third,])
	pdf(outfile)
	HiC.plot.decay(upper, bin.size, "n", TRUE, rgb(0,0,1,0.1,0.7), -4, 0.3)
	
	#4.8 is a good spacing for "o" style plotting, on outside of default balls. 3.2 is good for "l"
	HiC.plot.decay(upper, bin.size, "p", FALSE, rgb(0,0,1,0.025), 4.6, 1)
	HiC.plot.decay(upper, bin.size, "p", FALSE, rgb(0.9,0.9,0,0.025), 0, 1)
	HiC.plot.decay(bottom, bin.size, "p", FALSE,rgb(1,0.4,0,0.025), -4.6, 1)
	HiC.plot.decay.addMeans(upper, bin.size, "blue")
	HiC.plot.decay.addMeans(middle, bin.size, rgb(1,1,0,0.7))
	HiC.plot.decay.addMeans(bottom, bin.size, rgb(1,0.4,0,0.7))
	dev.off()
}

# Makes shading for "compartments" from Hi-C data. Takes as input a file of O/E matrix dumped from juicebox and converted to a proper matrix in python. Shading can be matched to matrix heatmap in Illustrator. Instructions: 1) Dump the OE matrix from juicebox, making sure you're in upper left so it starts with 0. Convert to matrix with python script. Needs to be 25 kb to work for same same. 
compartment.shading <- function(OE.matrix.file, outfile){
	x <- read.matrix(OE.matrix.file)
	x[x > 5] <- 5 #gets rid of weird outliers
	km <- kmeans(x, 2)
	ids <- km$cluster[55:390]
	ids <- c(ids, 0.5) #just to make the last one print...really dumb I know
	pdf(outfile,25,6)
	plot(ids,type="n", bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	
	curr.start <- 0
	curr.id <- 0
	for (i in 1:length(ids)){
		id <- ids[i]
		if (id != curr.id){
			if(curr.id != 0){
				if (curr.id == 1){ colour = "gray"}
				if (curr.id == 2){ colour = "light gray"}
				rect(curr.start- 0.5, 0, i-0.5, 2, col=colour, border = NA)
			}
			curr.id <- id
			curr.start <- i
			
		}
	}
	dev.off()
}

HiC.heatmap.difference <- function(filename1, filename2, name1, name2){
	process <- function(filename){
		x <- read.matrix(filename)
		x <- grab.chr(x, '3R')
		x <- HiC.matrix.vanilla.normalize(x)
		x <- HiC.matrix.scale.log(x)
		x <- x[55:390,55:390]
		return(x)
	}
	plot.subtraction <- function(y1, y2, y1.name, y2.name){
		date.string <- gsub('-','',Sys.Date())
		y3 <- y1 - y2
		#y3 <- histogram.equalize(y3)
		jpeg(paste(date.string, '_', y1.name, '_minus_',y2.name,'.jpeg',sep=''),4000,4000)
		heatmap.natural(y3)
		dev.off()
	}
	x1 <- process(filename1)
	x2 <- process(filename2)
	plot.subtraction(x1, x2, name1, name2)	
	plot.subtraction(x2, x1, name2, name1)	
	#x.1min2 <- 
}
########################################################################
# HELPER FUNCTIONS
########################################################################

# Reads delimited file in as a matrix instead of dataframe
read.matrix <- function(x){
	return(as.matrix(read.table(x)))
}

# Standard formatted write.table with quotes off and tab sep
write.fileout <- function(x, file){
	write.table(x, file, sep = '\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# returns reverse complement
rev.comp <- function(x){
	x <- toupper(x)
	x <- chartr('GATC','CTAG', x)
	x <- sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
	return(x)
}