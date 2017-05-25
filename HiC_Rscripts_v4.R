########################################################################
# Figure making functions
########################################################################

Fig.allchr.heat <- function(infilename, outfilename){
	a <- HiC.matrix.processfromfile(infilename, 1, 0.7)
	jpeg(outfilename, 4000,4000)
	heatmap.natural(a)
	dev.off()
}

########################################################################
# Main functions
########################################################################

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
	x <- HiC.matrix.process.standard(x, top, bottom)
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


# Makes a plot of an "epigenomic" dataset (e.g., ChIP or DNase) that matches a region of Hi-C data. Hi-C matri supplied in hic and epigenomic file in epifile must have row names in the format of "chr_bin". The epigenomic file is normalized (subtract minimum, divide by max) and can optionally be expanded or contracted for visual purposes by raising to exponent. This probably shouldn't be used but some of the tracks are more "even" than others and it makes viewing them together difficult. As long as this change is noted, I think it's ok. Just changes the visuals, and for this purposes we're just trying to show where a signal is high or low so I think it's ok.
epigenomic.match.plot <- function(hic, epifile, outfile, colour, exponent=1){
	epi <- read.table(epifile, row.names=1)
	epi <- epi[row.names(hic),1]
	#epi <- epi / mean(epi)
	epi <- epi - min(epi)
	epi <- epi  / max(epi)
	epi <- epi^exponent
	pdf(outfile,25,6)
	par(yaxp=c(0,10,10))
	plot(epi, type="l", bty="n",xaxt="n",yaxt="n",xlab="",ylab="",lwd=6.5, col=colour)
	axis(side=2, col="black", at=c(0,1),labels=c('',''))
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

# Takes two files of Hi-C binned matrices, produces a simple difference map for the "cool region" at left end of 3R. I also set it up to do autonaming for hte output file with the date string at the beginning. Not sure why.
cool.region.difference.heatmap <- function(filename1, filename2, name1, name2){
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
}


# Takes two Hi-C binned matrix files, converts each into a single numeric vector, calculates teh correlation and makes a scatter plot. The plots contain too many points because the matrices are large, so I used sampling to just plot some of the points. The correlation is on the full dataset.
HiC.matrix.correlation.plot <- function(file1, file2, outfile, sample.size=10000, return=FALSE){
	x <- unlist(read.table(file1))
	y <- unlist(read.table(file2))
	x <- x / sum(x) * 1e8
	y <- y / sum(y) * 1e8
	R.corr <- cor(x,y)
	print(R.corr)
	sample <- sample.int(length(x), sample.size)
	pdf(outfile)
	plot(x[sample], y[sample], log="xy", xlab="", ylab="", pch=19, col=rgb(0,0,0,0.2))
	legend("topleft", as.character(R.corr))
	dev.off()
	if(return){
		return(list(x, y, sample))
	}
}


# prints a simple key for the heatmap natural color scheme.
heatmap.key.make <- function(outfilename){
	x <- as.matrix(rbind(1:200,1:200))
	jpeg(outfilename,200,200)
	heatmap.natural(x)
	dev.off()
}


#numbers are for 1300 bp bins
polytene.plot.heat <- function(filename, outname, topcompress=0.995, bottomcompress=0.75, startbin=565, endbin=645, LOG=TRUE){
	x <- read.matrix(filename)
	x <- x[startbin:endbin,startbin:endbin]
	x[is.na(x)] <- 0
	x <- HiC.matrix.compress(x,topcompress, bottomcompress)
	if (LOG){
		x <- HiC.matrix.scale.log(x)
	}
	jpeg(outname, 4000,4000)
	heatmap.natural(x)
	dev.off()
}

boundaries.call <- function(x, left.thresh, right.thresh){
	boundaries.merge <- function(x1, merge.dist=2000){
		#return(x1)
		merged <- data.frame(chr=character(), Lmost = integer(), Rmost = integer(), name = character(), boundaryScore = numeric())
		current.start <- 1
		current.end <- 1
		current.max.score <- x1[1, 6]
		
		
		for(i in 2:nrow(x1)){
			if (x1[i, 2] == x1[current.end, 3] + 1){
				current.end <- i
				if (x1[i,6] > current.max.score){
					current.max.score <- x1[i,6]
				}
				#if (x1[i,6] > current.line[,6]){current.line <- x1[i,]} #comment this to select left-most boundary, leave for highest scoring boundary
			}
			else{
				midpoint <- round((current.start + current.end) / 2, 0)
				newline <- cbind(x1[midpoint,1:3], 'boundary', current.max.score)
				
				#print(newline)
				merged <- rbind(merged, newline)
				#print(merged)
				current.start <- i
				current.end <- i
				current.max.score <- x1[i,6]
			}
		}
		return(merged)
	}
	left.thresh <- -1 * left.thresh
	called <- x[x[,4] <= left.thresh & x[,5] >= right.thresh,]
	called <- cbind(called, called[,5] - called[,4])
	return(boundaries.merge(called))
}

local.heatmap.plot.series <- function(folder){
	files <- list.files(folder)
	for (file in files){
		if(! grepl('location', file)){
			date.string <- gsub('-','',Sys.Date())
			#outstem <- gsub('.*/', '', filename, perl = TRUE)
			#outstem <- gsub('.txt','', outstem)
			if (substr(folder, nchar(folder), nchar(folder)) != '/'){ folder <- paste(folder, '/', sep='')}
			outstem <- gsub('.txt','', file)
			
			x <- read.matrix(paste(folder, file, sep=''))
			x[is.na(x)] <- 0
			for (top in c(1, 0.99, 0.98, 0.97)){
				for (bottom in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)){
					outfile <- paste(date.string, '_', outstem, '_', top, '_', bottom, '.jpg' ,sep='')
					x1 <- HiC.matrix.compress(x, top, bottom)
					diag(x1) <- max(x1)
					jpeg(outfile, 4000,4000)
					heatmap.natural(x1)
					dev.off()
				}
			}
		}
	}
}

chip.heatmaps.control <- function(folder, outfolder, number){
	files <- list.files(folder)
	x <- read.table(paste(folder, files[1], sep="/"), row.names=1)
	sortby <- rownames(x)[sample(1:nrow(x), number)]
	chip.heatmaps.folder.sortbyX(folder, sortby, outfolder)
}

chip.heatmaps.polycomb <- function(folder, outfolder){
	files <- list.files(folder)
	for (file in files){
		if(grepl('.txt',file)){
			path <- paste(folder, file,sep='/')
			x <- read.table(path, row.names=1)
			left.means <- apply(x[,1:99], MARGIN=1, mean)
			right.means <- apply(x[,101:200], MARGIN=1, mean)
			thresh <- quantile(apply(x, MARGIN=1, mean), 0.6)
			over.under <- rownames(x)[left.means > thresh & right.means < thresh]
			under.over <- rownames(x)[left.means < thresh & right.means > thresh]
			over.over <- rownames(x)[left.means > thresh & right.means > thresh]
			under.under <- rownames(x)[left.means < thresh & right.means < thresh]
			sortby <- c(over.under, under.over, over.over, under.under)
			x <- x[sortby, ]
			x <- as.matrix(x)
			x <- chip.heatmap.compress(x,0.99,0.7)
			gene.name <- gsub('.txt','',file)
			#pdf(paste(folder,'/',gene.name,'.pdf',sep=''),15,15)
			jpeg(paste(outfolder,'/',gene.name,'.jpeg',sep=''),4000,4000)
			chip.heatmap(x, 'test', 'red')
			dev.off()
		}
	}
}

chip.heatmaps.folder.sortbyBoundaryStrength <- function(folder, boundaryfile, outfolder){
	boundaries <- read.table(boundaryfile)
	boundaries <- boundaries[order(boundaries[,5], decreasing=TRUE),]
	sortby <- paste(boundaries[,1], ':', boundaries[,2], '-', boundaries[,3], sep="")
	chip.heatmaps.folder.sortbyX(folder, sortby, outfolder)
}

chip.heatmaps.folder.sortbyDNase <- function(folder, DNase.file, outfolder){
	dnase <- read.table(DNase.file, row.names=1)
	window <- 5
	middle <- round(ncol(dnase) / 2, 0)
	means <- apply(dnase[,(middle - window):(middle + window)], MARGIN=1, mean)
	dnase <- dnase[order(means, decreasing=TRUE),]
	sortby <- rownames(dnase)
	chip.heatmaps.folder.sortbyX(folder, sortby, outfolder)
}

chip.heatmaps.folder.sortbyX <- function(folder, sortby, outfolder){
	files <- list.files(folder)
	for (file in files){
		if(grepl('.txt',file)){
			path <- paste(folder, file,sep='/')
			x <- read.table(path, row.names=1)
			x <- x[sortby, ]
			x <- as.matrix(x)
			x <- chip.heatmap.compress(x,0.99,0.7)
			gene.name <- gsub('.txt','',file)
			#pdf(paste(folder,'/',gene.name,'.pdf',sep=''),15,15)
			jpeg(paste(outfolder,'/',gene.name,'.jpeg',sep=''),4000,4000)
			chip.heatmap(x, 'test', 'blue')
			dev.off()
		}
	}
}

chip.heatmap.compress <- function(x, top=0.95, bottom=0.05){
	max.val <- quantile(unlist(x), top, na.rm=TRUE)
	min.val <- quantile(x, bottom, na.rm=TRUE)
	x[x > max.val] <- max.val
	x[x < min.val] <- min.val
	return(x)
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
	heatmap.2(x,dendrogram='none', main=title, Rowv=FALSE, Colv=FALSE,symm=TRUE,key=FALSE,keysize=0.5,key.title=NA,key.xlab=NA,key.ylab=NA,trace='none',scale='none',labRow=NA,labCol=NA, col=colorRampPalette(colour)(1000))
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

