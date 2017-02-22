
# Standard distance histogram with 1 kb steps
HiC.dist.hist <- function(x, title){
	hist(x$V1,breaks=c(seq(0,100000,1000),100000000),xlim=c(0,50000), col = 4, ylab="Counts",xlab="Distance", main=title)
	abline(v=c(1000,5000,10000),lty=2,col=2)
	pct.1kb <- length(x[x[,1] > 10000,1]) / length(x[,1]) * 100
	pct.1kb <- paste(as.character(round(pct.1kb, digits=2)),'% > 1kb', sep='')
	pct.10kb <- length(x[x[,1] > 100000,1]) / length(x[,1]) * 100
	pct.10kb <- paste(as.character(round(pct.10kb, digits=2)),'% > 10kb', sep='')
	legend("topright", c(pct.1kb, pct.10kb))
}

pct.greaterthan <- function(x,num){
	x1 <- length(x[x[,1] > num,]) / length(x[,1])
	print(x1)
}

# Custom histogram that represents fraction of total reads at 100 kb intervals, also tosses out things <10 kb

HiC.dist.hist2 <- function(x, title){
	x <- x[x[,1] > 10000,]
	#return(x)
	x.hist <- hist(x,breaks=c(seq(0,10000000,100000),100000000),xlim=c(0,10000000), col = 4, ylab="Counts",xlab="Distance", main=title)
	
}

HiC.distance.distr.1 <- function(x, color){
	x <- x[x[,1] > 2000,] #filter away cis reads (crude)
	x.hist <- hist(x,breaks=seq(0,100000000,100000),plot=FALSE)
	points(x.hist$counts/sum(x.hist$counts),type="l",col=color)
	
}

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
	return(x1)
}

# Scales a matrix to the median of the diagonal
HiC.matrix.scale.diagonalMedian <- function(x){
	x1 <- x
	diagonal <- numeric()
	for (i in 1:length(x1[1,])){
		diagonal <- c(diagonal, x[i,i])
	}
	m <- median(diagonal)
	return (x1 / m)
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
# 0.1 to 1
temp <- function(x){
	x.vals <- seq(0.1, 1, 0.9 / nrow(x))
	
	plot(x.vals,x[,1],type="l",lwd=2)
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

#Extracts from a binned Hi-C matrix just a single chromosome or arm (intxs with self only)
grab.chr <- function(x, chr){
	hits <- grep(chr,rownames(x))
	return(x[hits,hits])
}

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

# Takes a matrix of values (design is for them to be centered on diagonal values) and plots each row independently. Recommended colors rgb(0,0,1,0.025) and rgb(1,0.6,0,0.025)
HiC.plot.decay <- function(x,new=TRUE,color){
	max <- max(x)
	width <- ncol(x)
	height <- nrow(x)
	x.vals <- seq(-100 *((width - 1) / 2), 100 *((width - 1) / 2),100)
	#return(x.vals)
	if(new){
		plot(x.vals,x.vals, type="n", ylim = c(0,max), xlab = "Distance (kb)", ylab = "Fraction of Total Reads")
	}
	for (i in 1:height){
		points(x.vals, x[i,], type = "l",lwd=10,col=color)
	}
	
}

HiC.plot.decay.addMeans <- function(x, color){
	width <- ncol(x)
	height <- nrow(x)
	x.vals <- seq(-100 *((width - 1) / 2), 100 *((width - 1) / 2),100)
	means <- apply(x, MARGIN=2, mean)
	points(x.vals, means, col=color, type="l", lty=2,lwd=0.8)
}

HiC.make.uniform <- function(x,width,size){
	x.middles <- HiC.matrix.extractMiddle(x, width)
	#x.decay <- x.middles[,51:101]
	x.means <- apply(x.middles, MARGIN=2, mean)
	x.new <- matrix(nrow=size,ncol=size)
	for (i in 1:size){
		start <- max((width + 2) - i, 1)
		end <-  min((width + 1), size - i + 1) + width
		#print(paste(start,end))
		start.zeros <- i + start - width - 2
		end.zeros <- size - (i + end - width) + 1
		new.row <- c(rep(0, start.zeros), x.means[start:end], rep(0, end.zeros))
		x.new[i,] <- new.row
		#print(paste(start, end, start.zeros, end.zeros))
	}
	return(x.new)
}

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

HiC.localPlot.from.folder <- function(folder, outfolder, LOG=FALSE,HISTEQ=FALSE){
	files <- list.files(folder)
	for (file in files){
		
		file.stem <- gsub('.txt','',file)
		outfile <- paste(outfolder, '/', file.stem, sep='')
		file.path <- paste(folder, file, sep='/')
		x <- read.matrix(file.path)
		if(LOG){
			x <- HiC.matrix.scale.log(x)
			outfile <- paste(outfile, '_log', sep='')
		}
		if (HISTEQ){
			x <- histogram.equalize(x)
			outfile <- paste(outfile, '_histeq', sep='')
		}
		median <- quantile(x, 0.9)[1]
		#x[x < median] <- median
		outfile <- paste(outfile, '.jpeg', sep='')
		jpeg(outfile, width = 4000, height = 4000)
		heatmap.natural(x)
		abline(v=seq(0.1,1,0.1),lty=2, lwd=1.2, col="white") #v is 0.1 to 0.97, h is 0 to 0.9
		abline(h=seq(0,0.9,0.1),lty=2,lwd=1.2, col="white")
		points(0.535,0.46,pch=19,col="green",cex=3)
		dev.off()
	}
}

# searches a DNA string for supplied dinucleotide (or its reverse complement) and converts into a digital string,
# 1 where the dinuc is present and 0 where it isn't. For downstream frequency analysis.
dinucleotide.digitalString <- function(seq, dinuc){
	len <- nchar(seq)
	seq <- toupper(seq)
	dinuc <- toupper(dinuc)
	#sdinuc.rc <- rev.comp(dinuc) #comment out for no RC
	matches.f <- gregexpr(dinuc,seq)
	#matches.r <- gregexpr(dinuc.rc,seq) #comment out for no RC
	hits <- rep(0,len - 1)
	if(attributes(matches.f[[1]])$match.length[1] != -1){
		hits[as.numeric(matches.f[[1]])] <- 1
	}
	#if(attributes(matches.r[[1]])$match.length[1] != -1){ #comment out for no RC
	#	hits[as.numeric(matches.r[[1]])] <- 1			   #comment out for no RC
	#}                                                     #comment out for no RC
	return(hits)
}

# Takes a digital vector of dinucleotide positions, plots either autocorrelation or fourier transform.
dinucleotide.plot.acf <- function(seq, dinuc, funct='ACF'){
	if(funct == 'ACF'){
		ac <- acf(dinucleotide.digitalString(seq, dinuc),plot=FALSE)
		plot(ac[2:length(ac[[1]])], main=dinuc)
	}
	if(funct == 'FFT'){
		fft <- fft(dinucleotide.digitalString(seq, dinuc))
		plot.frequency.spectrum(fft, dinuc,c(0,2.5 * max(Mod(fft)[2:50])),c(0,30))
	}
}

# Plots the autocorrelation function for all dinucleotides in the supplied sequence. Takes as input a file of
# just sequences, no carrots
dinucleotide.plotall.acf <- function(file, funct='ACF'){
	seq <- seq.asSingleString.fromFile(file)
	par(mfcol=c(4,4))
	nts <- c('G','A','T','C')
	for (nt1 in nts){
		for (nt2 in nts){
			dinuc <- paste(nt1, nt2, sep="")
			dinucleotide.plot.acf(seq, dinuc, funct)			
		}
	}
}




########################################################################
# HELPER FUNCTIONS
########################################################################

# Reads delimited file in as a matrix instead of dataframe
read.matrix <- function(x){
	return(as.matrix(read.table(x)))
}

write.fileout <- function(x, file){
	write.table(x, file, sep = '\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}

rev.comp <- function(x){
	x <- toupper(x)
	x <- chartr('GATC','CTAG', x)
	x <- sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
	return(x)
}

plot.frequency.spectrum <- function(X.k, title, ylimits = c(0,range(Mod(X.k))[2]),xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))

  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main=title, 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=ylimits)
}

seq.asSingleString.fromFile <- function(file){
	seq <- read.table(file)
	seq <- seq[grep('>',seq[,1], invert=TRUE),] #matches all rows without '>'
	s <- ''
	for (i in 1:length(seq)){
		s <- paste(s, seq[i],sep='')
	}
	seq <- s

}