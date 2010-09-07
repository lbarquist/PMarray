setMethod("initialize", "PMset", function(.Object, files, idField) {
	.Object@files <- files
	.Object@idField <- idField
	orgList = list()
	plateList = list()
	tempList = list()
	for(file in files){
		handle = file(description=file, open="r")
		fileData = readLines(con = handle, n= -1)
		close(handle)
		i = 0
		
		while (i < length(fileData) - 1){
			i = i +1
			sub("^M","",fileData[i]);
			if (regexpr("Data\\sFile", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				dataFile = gsub("\\s","", splat[[1]][2], perl=TRUE)
			} else {
				stop(paste("Data File field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Set\\sup\\sTime", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				setupTime = splat[[1]][2]
			} else {
				stop(paste("Set up Time field not found, file ", file, ",line ", i,"\n"))
			} 
			i = i + 1 
			if (regexpr("Position", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				pos = splat[[1]][2]
			} else {
				stop(paste("Position field not found, file ", file, ",line ", i,"\n"))
			}	
			i = i + 1 
			if (regexpr("Plate Type", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				plateType = gsub("\\s","",splat[[1]][2], perl=TRUE)
				if (! (plateType %in% plateList) ){
					plateList[length(plateList) + 1] = plateType
				}
			} else {
				stop(paste("Plate Type field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Strain Type", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				if (regexpr("NOT APPLICABLE",splat[2], perl=TRUE) == 1)
					strainType = NA
				else
					strainType = splat[[1]][2]
			} else {
				stop(paste("Strain Type field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Sample Number", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				sampleNumber = splat[[1]][2]
				if (regexpr(idField, "sampleNumber", fixed=TRUE) == 1 && !(sampleNumber %in% orgList)){
					orgList[length(orgList)+1] <- sampleNumber
				}	 
			} else {
				stop(paste("Sample Number field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Strain Name", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				strainName = splat[[1]][2]
				if (regexpr(idField, "strainName", fixed=TRUE) == 1 && !(strainName %in% orgList)){
					orgList[length(orgList)+1] <- strainName
				}
			} else {
				stop(paste("Strain Name field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Strain Number", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				strainNumber = splat[[1]][2]
				if (regexpr(idField, "strainNumber", fixed=TRUE) == 1 && !(strainNumber %in% orgList)){
					orgList[length(orgList)+1] <- strainNumber
				}
			} else {
				stop(paste("Strain Number field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 1 
			if (regexpr("Other", fileData[i],perl=TRUE) == 1){
				splat = strsplit(fileData[i], ",", perl=TRUE)
				other = splat[[1]][2]
			} else {
				stop(paste("Other field not found, file ", file, ",line ", i,"\n"))
			}
			i = i + 2
			dataList = list()
			if (regexpr("Hour", fileData[i]) == 1){
				colNames <- strsplit(fileData[i],",",perl=TRUE)
				colNames = gsub(" ", "", colNames[[1]])
			} else {
				stop(paste("incorrect file format, file ", file,",line ", i,"\n"))
			}
			i = i + 1
			rowNames = list(c())
			localPos = 1
			fileData[i]
			while(!(regexpr("^$", fileData[i], perl=TRUE) == 1)){
				temp <- strsplit(fileData[i],",",perl=TRUE)
				for (j in 1:length(temp[[1]])){
					if (j != 1){
						dataList[[length(rowNames[[1]])]][j-1] = temp[[1]][j]
					} else {
						rowNames[[1]][length(rowNames[[1]]) + 1] = temp[[1]][1]
						dataList = c(dataList, list(c()))
					}
				}
				localPos = localPos + 1
				i = i + 1
			}
			.Object@organisms <- orgList
			.Object@plateTypes <- plateList
			.Object@timeStep <- as.numeric(rowNames[[1]][2])
			.Object@time <- as.numeric(rowNames[[1]][length(rowNames[[1]])]) + as.numeric(rowNames[[1]][2])
			data <- data.frame(matrix(do.call("rbind", lapply(dataList, as.numeric)), nrow = length(rowNames[[1]]), dimnames=list(rowNames[[1]], colNames[2:length(colNames)])))
			.Object@wellNames <- colNames[2:length(colNames)]
			i = i + 1
			tempList[length(tempList) +1] <- new("PlateFrame", dataFile=dataFile, setupTime=setupTime, pos=pos, plateType=plateType, strainType=strainType, sampleNumber=sampleNumber, strainName=strainName, strainNumber=strainNumber, other=other, data=data)
		}		 		
	}
	repMat = matrix(data=0, nrow=length(plateList), ncol=length(orgList))
	for ( i in tempList ) {
			row <- grep(i@plateType,plateList)
			col <- grep(paste("^",slot(i, idField),"$",sep=""), orgList)
			repMat[row,col] <- repMat[row,col] + 1
			#print(paste(paste("^",slot(i, idField),"$",sep=""),i@strainNumber, i@plateType, repMat[row,col]))
	}
	warning("Assuming ", repMat[1,1], " replicates\n")
	#todo: Need sanity check on repVec
	.Object@replicates <- repMat[1,1]
	#repList = list()
	#sortList = list()
	#for (i in 1:length(plateList)){
	#	thisList=list()
	#	for (j in 1:length(orgList)){
	#		thisList = c(thisList, list(c()))
	#	}
	#	sortList = c(sortList, list(c()))
	#}
	#for (i in orgList){
	#	for (j in 1:repVec[1]){
	#		repList[length(repList) + 1] <- paste(i,".",j, sep="")
	#	}
	#}
	#topList = vector("list",length(plateList))
	#containerFrame <- data.frame(row.names = plateList)
	#for( i in repList ){
	#	containerFrame[[i]] <- vector(mode="list", length=length(plateList))
	#}
	#row.names(containerFrame) <- plateList 
	#containerFrame = data.frame(matrix(data=NA,nrow=length(plateList), ncol=length(repList)), row.names=plateList);
	#names(containerFrame) = repList
	#print(repList)
	for ( i in 1:length(tempList)) {
		#colList = list()
		row <- grep(tempList[[i]]@plateType,plateList)
		col <- grep(paste("^",slot(tempList[[i]],idField),"$",sep=""),orgList)
		tempList[[i]]@replicateTag <- paste(slot(tempList[[i]], idField),"^",repMat[row,col], sep="")
		#print(paste("Row: ",row," Column: ", col," plateType: ", paste(slot(tempList[[i]], idField),"^",repMat[row,col], sep=""),"\n", sep=""))
		repMat[row,col] = repMat[row,col] - 1
	#	containerFrame[[paste(slot(i, idField),".",repVec[row], sep="")]][row] <- i
		#colList <- sortList[column]
		#colList[[row]] <- i
		#sortList[column] <- colList # INDIVIDUAL LISTS FOR EACH ROW --- data.frame eg. data.frame(252.1 = list()..., row.names)
	}
	#print(sortList)
	#.Object@plates <- data.frame(matrix(do.call("rbind", sortList), nrow = length(plateList), dimnames=list(plateList,repList)))
	.Object@plates <- tempList
	
	.Object
})


"getReplicate" <- function(PMset="PMset", plate="character", replicateTag="character"){
	for (i in PMset@plates){
		if(regexpr(replicateTag, i@replicateTag, fixed=TRUE) == 1 ){
			if(regexpr(plate, i@plateType, fixed=TRUE) == 1){
				return(i)
			}
		}
	}
	
	stop(paste("No matching entry found:", replicateTag, plate))	
}

"modelWell" <- function(PMset="PMset", plate = "character", replicateTag="character", well="character"){
	A = getReplicate(PMset, plate, replicateTag)
	cont = grofit.control(suppress.messages=TRUE)
	return(gcFitSpline(row.names(A@data), A@data[[well]], control=cont))
}


"signalValue" <- function(PMset="PMset", plate = "character", replicateTag="character", well="character", zero=FALSE){
	A = getReplicate(PMset, plate, replicateTag)
	normCount = as.integer(2/PMset@timeStep)
	sum = 0
	for ( i in c(1:normCount)){
		sum = sum + A@data[[well]][i]
	}
	normVal = sum / normCount
	sum = 0
	max = 0
	for ( i in c(1:length(A@data[[well]]))){
		sum = sum + A@data[[well]][i]
	#	print(A@data[[well]][i])
		if ( A@data[[well]][i] > max){
			max = A@data[[well]][i]
		}
	}
	ave = sum / length(A@data[[well]])
#	print (paste(ave, max, normVal))
	val = (ave + max)/2 - normVal
	if(zero){
		val = val - signalValue(PMset, plate, replicateTag, "A01")
	}
	return(val)
}


"reportOnWells" <- function (PMset="PMset", plateVec="vector", orgVec="vector", repVec="vector", file="character", zero=FALSE){
	noWells = length(PMset@wellNames) + 2
	thisMat = matrix(data = NA,nrow = length(orgVec), ncol = (noWells*(length(plateVec))))
	rownames(thisMat) = orgVec
	distVec <- c()
	for(i in plateVec){
		for (j in orgVec){
			for (k in repVec){
				distVec <- c(distVec, getSignalVector(PMset, plate = i, replicateTag=paste(j,k,sep='^'), lg=TRUE, zero=zero))
			}
		}
	}
	truehist(distVec,nbins="FD")
	ok = FALSE
	value <- readline(prompt="Choose an on/off threshold: ")
	offVec <- c()
	onVec <- c()
	for(i in distVec){
		if(i < as.numeric(value)){
			offVec <- c(offVec, i)
			#print(paste("off", i))
		}else{
			onVec <- c(onVec, i)
		}
	}
	offDist = fitdistr(offVec, "normal")

	truehist(offVec, nbins="FD")
	x = seq(min(offVec),max(offVec), length=200)
	curve(dnorm(x,offDist$estimate[1], offDist$estimate[2]), add=TRUE,type="l",col="red")
	readline(prompt="press enter")
	#print(offVec)
	onDist = fitdistr(onVec, "normal")
	x = seq(min(onVec),max(onVec), length=200)
	truehist(onVec,nbins="FD")
	curve(dnorm(x,onDist$estimate[1], onDist$estimate[2]), add=TRUE,type="l",col="red")
	
	for(i in 1:length(orgVec)){
		for(k in 1:length(plateVec)){
				for(l in 1:(noWells-2)){
					probOn <- c()
					probOff <- c()
					for(j in 1:length(repVec)){
						value <- log(signalValue(PMset, plate = plateVec[k], replicateTag=paste(orgVec[i],repVec[j],sep="^"), well=PMset@wellNames[l]), base=2)
						if(!is.finite(value)){
							value = 0;
						}
						#value <- signalValue(PMset, plate = k, replicateTag=paste(i,j,sep="^"), well=l)
						probOn <- c(probOn, pnorm(value,onDist$estimate[1], onDist$estimate[2], lower.tail=TRUE))
						probOff <- c(probOff,  pnorm(value,offDist$estimate[1], offDist$estimate[2], lower.tail=FALSE))
						#print(paste(i,k,l,j,value, probOn, probOff))
					}
					sumLogOdds = 0
					#probProdOn = probOn[1]
					#probProdOff = probOff[1]
					#print(length(probOn))
					count = 0
					for(m in 1:length(probOn)){
						#probProdOn = probProdOn * probOn[m]
						#print(probProdOn)
						#probProdOff = probProdOff * probOff[m]
						#sumLogOdds = sumLogOdds + log((probOn[m] * (1 - probOff[m])) / (probOff[m] * (1 - probOn[m])), base=2)
						#print(probProdOff)
						if (log((probOn[m] * (1 - probOff[m])) / (probOff[m] * (1 - probOn[m])), base=2) > 2){
							count = count + 1;
						}
					}
					#print(paste(probProdOn, probProdOff))
					
					#odds = (probProdOn * (1 - probProdOff)) / (probProdOff * (1 - probProdOn))
					if(count == length(repVec)){
						thisMat[i,(k-1)*noWells + l] <- 1
						cat(orgVec[i],plateVec[k], PMset@wellNames[l],"\n", file=file, append=TRUE,sep=",")
					}else{
						thisMat[i,(k-1)*noWells + l] <- 0
					}
				}	
				thisMat[i,(k-1)*noWells + l + 1] <- NA
				thisMat[i,(k-1)*noWells + l + 2] <- NA
		}
			
	} 
	pairs.breaks = seq(0,1, by=0.5)
	myPalette = colorpanel(n=2, low="black", high="yellow")
	#print(myPalette)
	heatmap.2(thisMat, Rowv=TRUE, Colv=FALSE, key=FALSE, ylab="Biotypes", col=myPalette, trace="none", density.info=c("none"), dendrogram=c("row"), margins=c(10,8), breaks=pairs.breaks, labCol="")
	#return(thisMat)
}

"plotWell" <- function (PMset="PMset", plate = "character", well="character", repVec, colVec="vector"){
	cont = grofit.control(suppress.messages=TRUE,interactive=FALSE)
	time = seq(from = 0.00, to = (PMset@time - PMset@timeStep), by = PMset@timeStep)
	#timeMat = matrix(replicate((length(PMset@organisms) * PMset@replicates), time), nrow=(length(PMset@organisms) * PMset@replicates), byrow=TRUE)
	timeMat = matrix(replicate(length(repVec),time), nrow=(length(repVec)), byrow=TRUE)
	#print(timeMat)
	df <- NULL
	chVec = vector("numeric", length= (length(PMset@organisms)*PMset@replicates))
	#colVec = vector("numeric", length= (length(PMset@organisms)*PMset@replicates))
	x = 1
	#for ( i in 1:length(PMset@organisms)){
	#	for(j in 1:PMset@replicates){
	for (z in repVec){
			#replicateTag = paste(PMset@organisms[i],"^",j, sep="")
			A = getReplicate(PMset, plate, z)
			#print(paste(z, plate, well))
			#print(A@data[[well]])
			#rowVec = rbind(c("", paste(PMset@organisms[i],"replicate",j), "", A@data[[well]]))
			rowVec = rbind(c("", paste(z), "", A@data[[well]]))
			#print(rowVec)
			df = rbind(df, rowVec)
			#print(df)
			chVec[x] = x
			#colVec[x] = i
			x = x + 1
	}
	#	}
	#}
	print ("Fitting Data")
	df=data.frame(df)
	#df[,4:] <- as.numeric(df[,4:])
	#print(df)
	#return(df)
	#print(min(timeMat), na.rm=TRUE)
	#for (i in 4:(length(PMset@wellNames) - 1)){
	#	print(df[[i]])
		
		#df[[i]] <- as.numeric(df[[i]])
	#}
	model = gcFit(time=timeMat, data=df, control = cont)
	#for (i in 1:length(model$gcFittedModels)){
	#	model$gcFittedModels$raw.data = as.numeric(model$gcFittedModels$raw.data)
	#}
	print("plotting data")
	#return(model)
	#print(model)

	plot.gcFit(model, opt="s", pch = chVec, cex = 1, colModel=1, colData=colVec)
}
	

"reportSig" <- function(PMset="PMset", plateVec="vector", reference="character", orgVec="vector",file="character", pval){
	eset <- getESet(PMset, plateVec, append(reference, orgVec), TRUE, FALSE)
	nameVec <- c()
	refName = make.names(reference)
	noWells <- length(PMset@wellNames)+2
	thisMat = matrix(data = NA,nrow = length(orgVec), ncol = (noWells*(length(plateVec))))
	sequence <- vector(mode="numeric")
	for(i in 1:(length(orgVec)+1)){
		sequence <- append(sequence, seq(from=i,to=i, length=PMset@replicates))
		if(i < length(orgVec) + 1){
			nameVec <- append(nameVec, paste(make.names(orgVec[i]),"-",refName,sep=""))
		}
	}
	#print(nameVec)
	design <- model.matrix(~ -1+factor(sequence))
	colnames(design) <- make.names(c(reference, orgVec))
	#print(design)
	fit <- lmFit(eset, design)
	contrast.matrix <- makeContrasts(contrasts=nameVec, levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	
	results <- decideTests(fit2, p.value=pval)
	#print(results)
	#print(fit2$coefficients)
	for(i in 1:length(orgVec)){
		
		pvals <- p.adjust(as.vector(fit2$p.value[,i]), method="BH")
		#print(pvals)
		for(k in 1:length(plateVec)){
			for(l in 1:(noWells-2)){
				if(results[(k-1)*length(PMset@wellNames) + l, i] != 0){
					#print(fit2$coefficients[(k-1)*l + l,i])
					thisMat[i, (k-1)*noWells + l] <- fit2$coefficients[(k-1)*length(PMset@wellNames) + l,i]
					cat(orgVec[i], PMset@wellNames[l], pvals[(k-1)*length(PMset@wellNames) + l], fit2$coefficients[(k-1)*length(PMset@wellNames) + l,i], "\n", file=file, append=TRUE, sep=", ")
				} else {
					thisMat[i, (k-1)*noWells + l] <- 0
				}
			}
		}
	}
	
	#print(thisMat)
	
	pairs.breaks = seq(-100, 100, by=200/500)
	myPalette = colorpanel(n=500, low="blue",mid="black", high="yellow")
	#print(myPalette)
	heatmap.2(thisMat, Rowv=TRUE, Colv=FALSE, keysize=1, ylab="Biotypes", labRow=orgVec, col=myPalette, trace="none", density.info=c("none"), dendrogram=c("row"), margins=c(10,8), breaks=pairs.breaks, labCol="")
}
#	first = TRUE;
#	for ( i in 1:length(PMset@organisms)){
#		for(j in 1:PMset@replicates){
#			
#			thisWell <- modelWell(PMset, plate, replicateTag=paste(PMset@organisms[i],".",j,sep=""), well)
#			#print(paste(PMset@organisms[1],".",j," ",j," ",well,sep=""))
#			if(first){
#				plot.gcFitModel(thisWell, pch=i)
#				first = FALSE
#			} else {
#				plot.gcFitModel(thisWell, add=TRUE, pch=i)
#			}
#		}
#	}
#}


"getAreaVector" <- function(PMset="PMset", replicateTag="replicateTag"){
	aVec = vector(mode="numeric", length=length(PMset@plateTypes) * length(PMset@wellNames))
	z = 1
	for (x in PMset@organisms){
	for (y in 1:PMset@replicates){
		replicateTag = paste(x,"^",y,sep="")
	for ( i in PMset@plateTypes ){
		for (j in PMset@wellNames){
			#print(paste(i," ",splat[[1]][1]," ",splat[[1]][2],"\n"))
			thisWell = modelWell(PMset, i, replicateTag, j)
			val  <-  log(thisWell$parameters$integral, base=2);
			if(is.finite(val)){
				aVec[z] <- val
			}
			else {
				aVec[z] <- 0
			}
			#aVec[z] <- thisWell$parameters$integral
			z = z+1
		}
	}
	}
	}
	
	return(aVec)
}

"getSignalVector" <- function(PMset="PMset", plate = "character", replicateTag="character", lg=TRUE, zero=FALSE){
	aVec = vector(mode="numeric", length=length(PMset@wellNames))
	z = 1
	#A = getReplicate(PMset, plate, replicateTag)
	for (j in PMset@wellNames){
			#print(paste(i," ",splat[[1]][1]," ",splat[[1]][2],"\n"))
			
		thisWell = signalValue(PMset, plate, replicateTag, j, zero=zero)
		if(lg){
			if(thisWell < 0){
				thisWell = 0
			}
			val  <-  log(thisWell, base=2);
		} else {
			val <- thisWell
		}
		if(is.finite(val)){
			aVec[z] <- val
		}
		else {
			aVec[z] <- 0
		}
		z = z + 1
		#aVec[z] <- thisWell$parameters$integral
	}
	return(aVec)
}

"pmHeatmap" <- function(PMset="PMset", orgVec="vector", plateVec="vector", log=TRUE, zero=FALSE){
	noWells = length(PMset@wellNames) + 2
	thisMat = matrix(data = NA,nrow = length(orgVec)*PMset@replicates, ncol = (noWells *(length(plateVec))))
	rnames= c()
	for( i in 1:length(orgVec)){
		for( j in 1:PMset@replicates){
			rnames <- c(rnames, paste(orgVec[i], "^", j, sep=""))
		}
	}
	#print(rnames)
	rownames(thisMat) <- rnames
	for( i in 1:length(orgVec)){
		for( j in 1:PMset@replicates){
		#	print(j)
			for( k in 1:length(plateVec)){
				for(l in 1:(noWells - 2)){
					#print(paste(plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l]))
					#print(log(signalValue(PMset, plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l]), base=2))
					if(log){
						if(zero){
							val <- log(signalValue(PMset, plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l], zero=TRUE), base=2)
						} else{
							val <- log(signalValue(PMset, plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l]), base=2)
						}
					} else {
						if(zero){
							val <- signalValue(PMset, plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l], zero=TRUE)
						} else{
							val <- signalValue(PMset, plateVec[k], paste(orgVec[i], "^", j, sep=""), PMset@wellNames[l])
						}
					}
					if(!is.finite(val)){
						thisMat[(i-1)*PMset@replicates+j,(k-1)*noWells + l] <- 0
					} else {
						thisMat[(i-1)*PMset@replicates+j,(k-1)*noWells + l] <- val
					}
				}
				thisMat[(i-1)*PMset@replicates+j,(k-1)*noWells + l + 1] <- NA
				thisMat[(i-1)*PMset@replicates+j,(k-1)*noWells + l + 2] <- NA
			}
		}
	}
	#print(thisMat)
	#print(max(thisMat))
	pairs.breaks = seq(0,max(thisMat, na.rm=TRUE), by=max(thisMat,na.rm=TRUE)/500)
	#print(pairs.breaks)
	myPalette = colorpanel(n=500, low="blue",mid="black", high="yellow")
	#print(myPalette)
	heatmap.2(thisMat, Rowv=TRUE, Colv=FALSE, keysize=1, ylab="Biotypes", col=myPalette, trace="none", density.info=c("none"), dendrogram=c("row"), margins=c(10,8), breaks=pairs.breaks, labCol="")
}
	

"getSignalMatrix" <- function(PMset="PMset"){
	aVec = vector(mode="numeric", length=length(PMset@plateTypes) * length(PMset@wellNames))
	z = 1
	for (x in PMset@organisms){
	for (y in 1:PMset@replicates){
		replicateTag = paste(x,"^",y,sep="")
		print(replicateTag)
	for ( i in PMset@plateTypes ){
		for (j in PMset@wellNames){
			#print(paste(i," ",splat[[1]][1]," ",splat[[1]][2],"\n"))
			
			thisWell = signalValue(PMset, i, replicateTag, j)
			val  <-  log(thisWell, base=2);
			if(is.finite(val)){
				aVec[z] <- val
			}
			else {
				aVec[z] <- 0
			}
			#aVec[z] <- thisWell$parameters$integral
			z = z+1
		}
	}
	}
	}
	
	return(aVec)
}

			

"getESet" <- function(PMset="PMset", plateVec="vector", orgVec="vector", signal="TRUE", log="TRUE"){
	eMat <- matrix(data=NA, nrow = (length(plateVec) * length(PMset@wellNames)), ncol=(length(orgVec) * PMset@replicates))
	colNames = vector(mode="character", length=(length(orgVec) * PMset@replicates))
	x = 1
	for ( i in orgVec ){
		for( j in 1:PMset@replicates){
			colNames[x] = paste(i,"^",j,sep="")
			x = x +1
		}
	}
	rowNames = vector(mode="character", length=(length(plateVec) * length(PMset@wellNames)))
	x=1
	for ( i in plateVec){
		for( j in PMset@wellNames){
			rowNames[x] = paste(i,".",j,sep="")
			x = x + 1
		}
	}
	#print(length(colNames))
	#print(length(rowNames))
	colnames(eMat) <- colNames
	rownames(eMat) <- rowNames
	
	for ( i in 1:length(colNames) ){
		for (j in 1:length(rowNames)){
			splat = strsplit(rowNames[j], "\\.")
			#print(paste(i," ",splat[[1]][1]," ",splat[[1]][2],"\n"))
			if(signal){
				val = signalValue(PMset, splat[[1]][1], colNames[i], splat[[1]][2])
				if(log){
					eMat[j,i] <- log(val)
				} else {
					#print(val)
					eMat[j,i] <- val
				}
			} else{
				thisWell = modelWell(PMset, splat[[1]][1], colNames[i], splat[[1]][2])
				val <- log(thisWell$parameters$integral, base=2)
				if(is.finite(val)){
					eMat[j,i] <- val
				}
				else {
					eMat[j,i] <- 0
				}
			}
			#eMat[j,i] <- log(thisWell$parameters$integral, base=2)
			#eMat[j,i] <- thisWell$parameters$integral
		}
	}
	#print(eMat)
	colnames(eMat) <- make.names(colnames(eMat))
	return( new("ExpressionSet",exprs=eMat))
}
