#####   R-code for making Figure 9.7. It reads data from "soil.dat"

	dat   <- matrix(scan("soil.dat"),189,6,byrow=T)

	id    <- dat[,1]
	
	x1    <- dat[,2]
	x2    <- dat[,3]
	x3    <- dat[,4]
	x4    <- dat[,5]
	x5    <- dat[,6]

	postscript("fig97.ps",width=6.5,height=7,horizontal=F)

	par(mfrow=c(3,2),mar=c(4,4,2,2))

        tx1 <- ar.yw(x1)
	d1  <- density(tx1$resid[4:189],width=0.5)
        plot(d1$x,d1$y,xlim=c(-3,3),
	     type="l",lty=1,xlab=expression(X[1]),ylab="density",
             mgp=c(2,1,0),cex=0.7)

	tx2 <- ar.yw(x2)
	d2  <- density(tx2$resid[3:189],width=0.5)
	plot(d2$x,d2$y,xlim=c(-3,3),
	     type="l",lty=1,xlab=expression(X[2]),ylab="density",
             mgp=c(2,1,0),cex=0.7)

	tx3 <- ar.yw(x3)
	d3  <- density(tx3$resid[2:189],width=0.5)
        plot(d3$x,d3$y,xlim=c(-3,3),
	     type="l",lty=1,xlab=expression(X[3]),ylab="density",
             mgp=c(2,1,0),cex=0.7)

	tx4 <- ar.yw(x4,order.max=1)
	d4  <- density(tx4$resid[2:189],width=0.5)
        plot(d4$x,d4$y,xlim=c(-3,3),
	     type="l",lty=1,xlab=expression(X[4]),ylab="density",
             mgp=c(2,1,0),cex=0.7)

	tx5 <- ar.yw(x5)
	d5  <- density(tx5$resid[2:189],width=0.5)
        plot(d5$x,d5$y,xlim=c(-3,3),
	     type="l",lty=1,xlab=expression(X[5]),ylab="density",
             mgp=c(2,1,0),cex=0.7)


	graphics.off()

#	write.table(cbind(tx1$resid,tx2$resid,tx3$resid,tx4$resid,tx5$resid),
#	            "resid.dat",sep=" ")
	