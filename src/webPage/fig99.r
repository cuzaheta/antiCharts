##### R-code for making Figure 9.9

postscript("fig99.ps",width=7.5,height=7.5,horizontal=F)

par(mfrow=c(2,2),mar=c(5,4,1,1))

   a <- c(0,0.2,0.4,0.6,0.8,1)
arl1 <- c(200.0014,65.7688,24.2355,10.7574,5.7345,3.5602)
arl2 <- c(200.059,200.059,200.059,200.059,200.059,200.059)

plot(a,arl1,type="l",lty=1,xlab="a",xlim=c(0,1),ylab="Out-of-control ARL",
     ylim=c(0,201),mgp=c(2,1,0),cex=0.8,xaxt="n",yaxt="n")
lines(a,arl2,lty=2)
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4",
     "0.6","0.8","1.0"),cex=.6)
axis(2,at=c(0,50,100,150,200),labels=c("0","50","100","150","200"),cex=.6)
title(xlab="(a)",cex=0.8)


   a <- c(0,0.2,0.4,0.6,0.8,1)
arl1 <- c(7.9313,6.6393,5.5759,4.6759,4.0253,3.5602)
arl2 <- c(20.6668,32.935,57.0757,103.0173,168.6493,200.059)

plot(a,arl1,type="l",lty=1,xlab="a",xlim=c(0,1),ylab="Out-of-control ARL",
     ylim=c(0,201),mgp=c(2,1,0),cex=0.8,xaxt="n",yaxt="n")
lines(a,arl2,lty=2)
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4",
     "0.6","0.8","1.0"),cex=.6)
axis(2,at=c(0,50,100,150,200),labels=c("0","50","100","150","200"),cex=.6)
title(xlab="(b)",cex=0.8)


   a <- c(0,0.5,1,1.5,2,2.5,3)
arl1 <- c(200.0014,118.8366,96.0702,85.5468,79.9357,76.7794,76.1507)
arl2 <- c(200.0590,134.2838,100.8009,90.1457,84.9329,82.7324,81.8320)

plot(a,arl1,type="l",lty=1,xlab="a",xlim=c(0,3),ylab="Out-of-control ARL",
     ylim=c(0,201),mgp=c(2,1,0),cex=0.8,xaxt="n",yaxt="n")
lines(a,arl2,lty=2)
axis(1,at=c(0,0.5,1,1.5,2,2.5,3),labels=c("0","0.5","1.0","1.5","2.0",
     "2.5","3.0"),cex=.6)
axis(2,at=c(0,50,100,150,200),labels=c("0","50","100","150","200"),cex=.6)
title(xlab="(c)",cex=0.8)


   a <- c(-3,-2.5,-2,-1.5,-1,-0.5,0)
arl1 <- c(5.6078,6.4421,8.3083,13.1103,34.4693,184.7309,200.0014)
arl2 <- c(4.3864,4.9288,6.2143,9.4128,20.6668,77.3018,200.059)

plot(a,arl1,type="l",lty=1,xlab="a",xlim=c(-3,0),ylab="Out-of-control ARL",
     ylim=c(0,201),mgp=c(2,1,0),cex=0.8,xaxt="n",yaxt="n")
lines(a,arl2,lty=2)
axis(1,at=c(-3,-2.5,-2,-1.5,-1,-0.5,0),labels=c("-3.0","-2.5","-2.0","-1.5",
     "-1.0","-0.5","0"),cex=.6)
axis(2,at=c(0,50,100,150,200),labels=c("0","50","100","150","200"),cex=.6)
title(xlab="(d)",cex=0.8)

graphics.off()


