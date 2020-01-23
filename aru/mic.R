#Variable	Definition
#Euip_Type	Either Bioacoustic Recorder (SM3 or SM4 with built in Mics) or Microphone (separate SM2 mics)
#MODEL	SM3 versus SM4 versus Acoustic Microphones
#OWNER	Owner of equipment (can ignore)
#PurchaseYear	Year ARU or microphone was purchased in, in most cases, first used
#Unit_NAME	Unique name of microphone or ARU.  For SM3 and SM4 unit, this occurs twice, once for the right and once for the left channel.
#Date_Cal	Date that microphone sensitivity was measured
#dBDiff	Diference between the target value and the actual value.  This is the response variable
#Sensitivity	Sensitivity of microphoneo n a given date. The range of this number depends on model, Tone volume and calibrator ID
#TargetdB	Target dB for a given mic, tone and calibrator combination, dependson mic specs and ARU model
#Tone_Volume	Either 94 for calibrator or 70 for calibrator with adaptor
#Calibrator_ID	Number of calibrator. There are all supposed to be 94 dB but there is some variation
#MicDamage	NONE or BROKEN or unknown for blanks
#Test_Unit	For SM2 mics, can ignore for now.
#Observer	Person doing testing, should not be needed and not completely consistently filled in
#Channel	IMPORTANT for SM3 and SM4 units because each unit has two lines of data, one for left and one for right channel.

x <- read.csv("~/GoogleWork/abmi/aru/mic/BU-Microphone-Data-Take-1-30November2017xlsx.csv")

x$dyr <- as.integer(substr(as.character(x$Date_Cal),8,9)) + 2000 - x$PurchaseYear
x$dbd <- sqrt(abs(x$dBDiff))
ii <- !is.na(x$dyr) & !is.na(x$dbd) & x$MicDamage == "NONE"

m <- lm(dbd ~ dyr, x[ii,])
summary(m)

plot(dbd ~ jitter(dyr), x[ii,], pch=19, col="#00000020",
    ylab="sqrt(abs(dB difference))", xlab="Calibration year - Purchase year")
#lines(lowess(dbd ~ dyr), col=2, lwd=2)
abline(m, col=2, lwd=2, lty=2)
lines(lowess(x=x[ii,"dyr"], y=x[ii,"dbd"]), col=2, lwd=2)
legend("topright", lwd=2, lty=c(1,2), col=2, legend=c("lowess", "lm"), bty="n")

m1 <- lm(dbd ~ dyr * MODEL, x[ii,])
summary(m1)
AIC(m, m1)

op <- par(mfrow=c(1,3))
for (i in levels(x$MODEL)) {
    xx <- x[ii & x$MODEL==i,]
    mx <- lm(dbd ~ dyr, xx)
    plot(dbd ~ jitter(dyr), xx, pch=19, col="#00000020", main=i, ylim=c(0,8), xlim=c(0,6),
        ylab="sqrt(abs(dB difference))", xlab="Calibration year - Purchase year")
    abline(mx, col=2, lwd=2, lty=2)
    lines(lowess(x=xx[,"dyr"], y=xx[,"dbd"]), col=2, lwd=2)
    legend("topright", lwd=2, lty=c(1,2), col=2, legend=c("lowess", "lm"), bty="n")
}

par(op)
