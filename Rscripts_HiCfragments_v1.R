
# Retrieves the elements in the enhancer element list tied to the queried gene
Get.enhancer.elements <- function(gene){
	return(enhancer.element.list[enhancer.element.list[,2] == gene,])
}

#Loads RF contact tables with correct settings
RF.contacts.load <- function(file){
	return(read.table(file, as.is=TRUE, fill=TRUE))
}

# Produces a Hi-C contact plot for a supplied element
Plot.element.contacts <- function(element, hiC.data, xmin, xmax, plot.elements=FALSE){
	element.loc <- enhancer.element.list[element,4]
	element.info <- split.chr.loc.diffseps(element.loc, ":", "\\.\\.")
	xmin2 <- as.numeric(element.info[2]) + xmin
	xmax2 <- as.numeric(element.info[3]) + xmax
	chr <- element.info[1]

	for (i in 1:length(hiC.data[,1])){
		fragment.info <- split.chr.loc.sameseps(hiC.data[i,1], "_")
		if (fragment.info[1] == chr){ #if same chromosome
			#print(hiC.data[i,1])
			if((element.info[2] < fragment.info[2] & element.info[3] < fragment.info[2])
			| (element.info[2] > fragment.info[3] & element.info[3] > fragment.info[3])
			){ #if both ends are to left of Lmost, or both ends are to R of Rmost, these don't overlap
				#Do nothing
			}
			else{
				#dev.new()
				#HiC.plot.fragmentrow(hiC.data[i,],fragment.info[1], xmin2, xmax2) #second part is chromosome
				HiC.table <- convert.fragrow.toTable(hiC.data[i,],fragment.info[1])
				HiC.table <- add.zeros.toHiCTable(HiC.table)
				HiC.plot1(HiC.table, xmin2, xmax2, element.loc)
				#return() #eventually add support for launching multiple frags
				if (plot.elements == TRUE){ Add.element.boxes(chr, xmin2, xmax2) }
				return()
			}
		}
	}
}

Plotall.element.contacts <- function(elem, xmin, xmax){
	par(mfrow=c(2,2))
	Plot.element.contacts(elem,fragment.data.hic[['c14.a.hic13']],xmin, xmax)
	Plot.element.contacts(elem,fragment.data.hic[['c14.p.hic14']],xmin, xmax)
	Plot.element.contacts(elem,fragment.data.hic[['c14.a.hic15']],xmin, xmax)
	Plot.element.contacts(elem,fragment.data.hic[['c14.p.hic16']],xmin, xmax)
	par(mfrow=c(1,1))
}

# Splits chromosome location blocks that use two separators, e.g. "2R:13450_14450", returns vector: chr, LeftMost, RightMost
split.chr.loc.diffseps <- function(combined, splitter1, splitter2){
	x <- strsplit(combined, splitter1)
	chr <- x[[1]][[1]]
	y <- strsplit(x[[1]][[2]], splitter2)
	Lmost <- as.numeric(y[[1]][[1]])
	Rmost <- as.numeric(y[[1]][[2]])
	return (c(chr, Lmost, Rmost))
}

# Splits chromosome location blocks that use identical separators, e.g. "2R_13450_14450", returns vector: chr, LeftMost, RightMost
split.chr.loc.sameseps <- function(combined, splitter1){
	x <- strsplit(combined, splitter1)
	#print(x[[1]])
	chr <- x[[1]][[1]]
	Lmost <- as.numeric(x[[1]][[2]])
	Rmost <- as.numeric(x[[1]][[3]])
	return (c(chr, Lmost, Rmost))
}


# Slightly cumbersome function that converts a row of fragment contacts to a table
convert.fragrow.toTable <- function(row, chr){
	output <- c(0,0,0)
	
	for (i in 2:length(row)){
		if(nchar(row[,i]) > 0){
			x <- strsplit(as.character(row[,i]), ':')
			#print(x)
			if (!is.na(x[[1]][1])){
				count <- as.numeric(x[[1]][[2]])
				y <- strsplit(x[[1]][[1]], "_")
				chr2 = y[[1]][[1]]
				Lmost = as.numeric(y[[1]][[2]])
				Rmost = as.numeric(y[[1]][[3]])
				if (chr == chr2){
					output <- rbind(output, c(Lmost, Rmost, count))
				}
			}
		}
	}
	output <- output[2:length(output[,1]),] #takes off blank first row
	output <- output[order(output[,1]),] #sorts by leftmost
	return(output)
}


# Fills in Hi-C contact table with 0s for fragments not reported
add.zeros.toHiCTable <- function(table.nozeros){
	table.withzeros <- table.nozeros[1,]
	prev.Rmost <- table.nozeros[1,2]
	output <- table.nozeros[1,]
	for (i in 2:length(table.nozeros[,1])){
		curr.line <- table.nozeros[i,]
		curr.Lmost <- curr.line[1]
		if (curr.Lmost == prev.Rmost){ #no gap
			output <- rbind(output, curr.line)
		}
		else{
			gap.line <- c(prev.Rmost, curr.Lmost, 0)
			output <- rbind(output, gap.line)
			output <- rbind(output, curr.line)
		}
		prev.Rmost <- curr.line[2]
	}
	return(output)
}

Add.element.boxes <- function(chr, min, max){
	for (i in 1:length(enhancer.element.list[,1])){
		element.loc <- enhancer.element.list[i,4]
		element.name <- row.names(enhancer.element.list)[i]
		element.info <- split.chr.loc.diffseps(element.loc, ":", "\\.\\.")
		if (chr == element.info[1]){
			L <- as.numeric(element.info[2])
			R <- as.numeric(element.info[3])
			print(c(L,R))
			if (R > min & L < max ){
				polygon(c(L,R,R,L),c(-10,-10,100,100),col=rgb(1,1,0.8,0.3))
				plot.coords <- par("usr")
				text.vert.pos <- (plot.coords[4] - plot.coords[3]) * 0.85
				text((R+L)/2,text.vert.pos,element.name,srt=90)
			}
		}
	}
	
}



# Plots Hi-C contact data from a table. Draws boxes for restriction fragments.
HiC.plot1 <- function(x,xmin, xmax, title){
	#print(x)
	#return()
	mids <- (x[,1] + x[,2]) / 2
	plot(mids,x[,3],type="o",xlim = c(xmin, xmax),lwd=2,ylab="",xlab="", main=title)
	#This part draws boxes
	#HiC.plot1.drawElementBoxes(x, xmin, xmax)
}

HiC.plot1.drawElementBoxes <- function(x, xmin, xmax){
	for (i in 1:length(x[,1])){
		window.height <- max(x[,3])
		box.height <- window.height * 0.01
		color = "blue"
		if (i %% 2 == 0){
			color = "lightblue"
		}
		#rect(x[i,1],box.height,x[i,2],-1 * box.height,col=color)
		xLeft <- x[i,1]
		xRight <- x[i,2]
		yBot <- -1 * box.height
		yTop <- box.height
		if (xRight > xmin & xLeft < xmax){
			polygon(c(xLeft, xRight, xRight, xLeft),c(yBot, yBot, yTop, yTop), col = color)
		}
	}
}

enhancer.shade <- function(floor,L, R, color){
	rect(L,floor,R,1E6,col=color,lwd=0.0)
}

shade.eve.enhancers <- function(){
	enhancer.shade(1,5865267,5865750,rgb(1,1,0,0.25))
	enhancer.shade(1,5863006,5863516,rgb(1,1,0,0.25))
	enhancer.shade(1,5873440,5874240,rgb(1,1,0,0.25))
	enhancer.shade(1,5874230,5875033,rgb(1,1,0,0.25))
	enhancer.shade(1,5871404,5872203,rgb(1,1,0,0.25))
}
