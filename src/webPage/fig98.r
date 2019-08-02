#####   R-code for making Figure 9.8. It reads data from files
#####   "fig98crit1.dat" and "fig98crit2.dat."

	dat1   <- matrix(scan("fig98crit1.dat"),,2,byrow=T)

	id1    <- dat1[,1]
	yn1    <- dat1[,2]

	dat2   <- matrix(scan("fig98crit2.dat"),,2,byrow=T)

	id2    <- dat2[,1]
	yn2    <- dat2[,2]

	
	postscript("fig98.ps",width=7,height=3.7,horizontal=F)

	par(mfrow=c(1,2),mar=c(4,4,1,1))

	plot(id1,yn1,type="l",lty=1,xlab="n",
	     ylab=expression(C[n*1]),mgp=c(2,1,0),cex=0.7)
	lines(id1,rep(12.32881,length(yn1)),lty=2)
	title(xlab="(a)",cex=0.9)

	plot(id2,yn2,type="l",lty=1,xlab="n",
	     ylab=expression(C[n*1*5]),mgp=c(2,1,0),cex=0.7)
	lines(id2,rep(29,length(yn2)),lty=2)
        title(xlab="(b)",cex=0.9)
	
	graphics.off()
	