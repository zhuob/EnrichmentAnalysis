data = read.csv("padog-data-analysis.csv")
index = which(data$meaca<0.05)

par(mfrow=c(1,2))
plot(-log(data$meaca[-index])/log(10), -log(data$GSEA[-index])/log(10),pch=19, 
	xlim=c(0, 6/log(10)), ylim=c(0,6/log(10)), 
	xlab="MEACA", ylab="GSEA", cex.lab=1.2)
points(-log(data$meaca[index])/log(10), -log(data$GSEA[index])/log(10), pch=19, col="blue")
abline(0,1)
#(h=-log(0.05))
#abline(v=-log(0.05))


plot(-log(data$meaca[-index])/log(10), -log(data$CAMERA_ModT[-index])/log(10),pch=19, 
	xlim=c(0, 6/log(10)), ylim=c(0,6/log(10)), 
	xlab="MEACA", ylab="CAMERA", cex.lab=1.2)
points(-log(data$meaca[index])/log(10), -log(data$CAMERA_ModT[index])/log(10), pch=19, col="blue")
abline(0,1)

