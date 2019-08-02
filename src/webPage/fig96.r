#####   R-code for making Figure 9.6. It reads data from soilResidual.dat"

	dat   <- matrix(scan("soilResidual.dat"),186,5,byrow=T)

	id    <- seq(1,186)
	
	x1    <- dat[,1]
	x1    <- (x1-mean(x1))
	
	x2    <- dat[,2]
	x2    <- (x2-mean(x2))
	
	x3    <- dat[,3]
	x3    <- (x3-mean(x3))
	
	x4    <- dat[,4]
	x4    <- (x4-mean(x4))
	
	x5    <- dat[,5]
	x5    <- (x5-mean(x5))

	postscript("fig96.ps",width=6,height=7.5,horizontal=F)

	par(mfrow=c(3,2),mar=c(4,4,1,1))

	plot(id,x1,type="o",lty=1,pch=1,xlab="n",ylim=c(-3,3),
	     ylab=expression(X[1]),mgp=c(2,1,0),cex=0.6)
	plot(id,x2,type="o",lty=1,pch=1,xlab="n",ylim=c(-3,3),
	     ylab=expression(X[2]),mgp=c(2,1,0),cex=0.6)
	plot(id,x3,type="o",lty=1,pch=1,xlab="n",ylim=c(-3,3),
	     ylab=expression(X[3]),mgp=c(2,1,0),cex=0.6)
	plot(id,x4,type="o",lty=1,pch=1,xlab="n",ylim=c(-3,3),
	     ylab=expression(X[4]),mgp=c(2,1,0),cex=0.6)
	plot(id,x5,type="o",lty=1,pch=1,xlab="n",ylim=c(-3,3),
	     ylab=expression(X[5]),mgp=c(2,1,0),cex=0.6)

	graphics.off()
	