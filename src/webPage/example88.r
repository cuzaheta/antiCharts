##### R-code for Example 8.8. It writes the related data to the file
##### "example88.dat" and makes "fig810.ps"

set.seed(1000)

#### Generate 500 observations of the IC data from the
#### IC process distribution t(3)/sqrt(3)

y = rt(500,3)/sqrt(3)

#### Define the boundary q1, q2, ... q_p-1 where p=10

y1=sort(y)      # ordered y from the smallest to the largest

q1=(y1[50]+y1[51])/2
q2=(y1[100]+y1[101])/2
q3=(y1[150]+y1[151])/2
q4=(y1[200]+y1[201])/2
q5=(y1[250]+y1[251])/2
q6=(y1[300]+y1[301])/2
q7=(y1[350]+y1[351])/2
q8=(y1[400]+y1[401])/2
q9=(y1[450]+y1[451])/2

f0=rep(0.1,10)  # f0 is the IC distribution of the categorized data

#### Generate 100 phase II observations for online process monitoring,
#### among which 50 observations are from the t(3)/sqrt(3) distribution, 
#### and another 50 observations are from the t(3)/sqrt(3)+1 distribution.
#### So, there is a mean shift at tau=51.

x1=rt(50,3)/sqrt(3)
x2=rt(50,3)/sqrt(3)+1

x = c(x1,x2)
n = length(x)

kP=0.01
hP=11.377

gnl=matrix(0,n,10)    # Define g_nl

for(i in 1:n){
   if(x[i] <= q1){gnl[i,1]=1}
   if ((x[i] > q1) & (x[i] <= q2)){gnl[i,2]=1}
   if ((x[i] > q2) & (x[i] <= q3)){gnl[i,3]=1}
   if ((x[i] > q3) & (x[i] <= q4)){gnl[i,4]=1}
   if ((x[i] > q4) & (x[i] <= q5)){gnl[i,5]=1}
   if ((x[i] > q5) & (x[i] <= q6)){gnl[i,6]=1}
   if ((x[i] > q6) & (x[i] <= q7)){gnl[i,7]=1}
   if ((x[i] > q7) & (x[i] <= q8)){gnl[i,8]=1}
   if ((x[i] > q8) & (x[i] <= q9)){gnl[i,9]=1}
   if (x[i] > q9){gnl[i,10]=1}
}

# Define the charting statistic C_{n,P}

Bn = rep(0,n)
Unobs=matrix(0,n,10)
Unexp=matrix(0,n,10)
CnP = rep(0,n)

Bn[1]=sum((gnl[1,]-f0)^2/f0)

if(Bn[1] > kP){
   Unobs[1,]=gnl[1,]*(1-kP/Bn[1])
   Unexp[1,]=f0*(1-kP/Bn[1])
}
if(Bn[1] <= kP){
   Unobs[1,]=t(rep(0,10))
   Unexp[1,]=t(rep(0,10))
}

CnP[1]=sum((Unobs[1,]-Unexp[1,])^2/Unexp[1,])

for(i in 2:n){

   Bn[i]=sum(((Unobs[i-1,]-Unexp[i-1,])+(gnl[i,]-f0))^2/(Unexp[i-1,]+f0))

   if(Bn[i]>kP){
      Unobs[i,]=(Unobs[i-1,]+gnl[i,])*(1-kP/Bn[i])
      Unexp[i,]=(Unexp[i-1,]+f0)*(1-kP/Bn[i])
   }
   if(Bn[i]<=kP){
      Unobs[i,]=t(rep(0,10))
      Unexp[i,]=t(rep(0,10))
   }

   CnP[i]=sum((Unobs[i,]-Unexp[i,])^2/Unexp[i,])

}

write.table(round(cbind(x,Bn,CnP),digits=3), 
            file="example88.dat",col.names=T, row.names=F)

postscript("fig810.ps",width=7,height=3.5,horizontal=F)

par(mfrow=c(1,2),mar=c(5,4,1,1))

ii = seq(1,n)

plot(ii,x,type="p",pch=16,xlab="n",ylab=expression(X[n]),mgp=c(2,1,0),
     xlim=c(1,100), ylim=c(-3.1,3.5), cex=0.8)
title(xlab="(a)",cex=0.9)

plot(ii,CnP,type="o",lty=1,pch=16,xlab="n",ylab=expression(C[list(n,P)]),
     mgp=c(2,1,0),xlim=c(1,100), ylim=c(0,70),cex=0.8)
lines(ii,rep(hP,n),lty=2,cex=0.8)
title(xlab="(b)",cex=0.9)

graphics.off()
