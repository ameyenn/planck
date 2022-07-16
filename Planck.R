
# Investigation of Planck's equation using wavelength in micron
# I was interested from a Climate point of view
# reference: https://en.wikipedia.org/wiki/Planck%27s_law
# NOTE: check calculations at 
# https://www.spectralcalc.com/blackbody_calculator/blackbody.php

# there is a lot of confusing terminology involving angles etc.
# I have tried to keep it simple hopefully I have not introduce errors
# in calculations or theory. email address meyenna@gmail.com
# Date: 16/7/2022

# Define Constants

# Using SI units
      h= 6.62607015*10^(-34) #Joules seconds
      c=3*10^8               #meters seconds
      k=1.38*1E-23           #Joules seconds
      sB=5.6E-8              #W · m -2 · K -4
      weinConstant = 0.2898 * 1E4 #results in microns
      #NOTE: 1. wavelengths are input as micron(um) and converted to meters
      #      2. Spectral Radiance units are J·s-1·m-2·sr-1·m-1
      #      3. Total in W s m^2
      #      4. Equating Integrated Value to Total W.m^2, multiply by pi
      #         reference http://people.tamu.edu/~kevinkrisciunas/planck.pdf

PhotonEnergy<-function(WL){
  #input is in microns and convert to meters
  #output is Joules/sec, same as Watts
  h*c/(WL *1E-6)
}

weinLawMaxWL<-function(T){
  weinConstant/T #result in micons, eg T=288, maxWL = 19.0625
}

powerRadiance<- function(T, emisivity, A){
  #Total energy emitted per second, at all wavelengths, per Area
  #area = 1m^2 by default measure
  #Stephan Boltzman Law, units w m^2 for a BlackBody
  
  A*emisivity*sB*T^4
}

tempGivenPower<-function(P, e){
  #assume A = 1
  (P/(e*sB))^0.25
}

planck<-function(WL, T){
  WL <- WL*1E-6 #convert um to m, SI units
  #returned value is W.m^2
  ((2*h*c*c) / (WL^5))*(1/(exp((h * c) /(WL * k * T))-1))*1E-6
}
options(scipen = 999) #turn of S notation
options(scipen=-999)  #turn on S notation

SpectralRadiance<-planck(10, 288)
SpectralRadiance #check at https://ncc.nesdis.noaa.gov/data/planck.html

#Plot Planck Curve
wL= seq(0.1,40, by=0.1)
t=288
heading <- c("plank curve==> Temp (K) ", t)
data<- sapply(wL, planck, t)
sum(data)*pi #give approx to Stefan-Boltzman x 10^1 ie integral value
plot(data, type="l", 
     main = heading,
     xlim=c(0.01, 400), 
     ylim=c(min(data), max(data)),
     xlab="wavelength, divide by 10",
     ylab="Intensity ")
abline(v=weinLawMaxWL(288)*10)


#Calculate Area under Planck SprectralRadiance curve for Wavelength

  powerRadiance(288, 1, 1) #Expected value of Integration

#NOTE: area = pi*Integral of Planck(WL) over wavelength interval

#Integration method 1 - sum of rectangles under curve and x-axis
df<-as.numeric(as.matrix(data))
spRadSum<-0
n<- length(df)-1
for (i in 1:n) {
  spRadSum<- ((i+1)-i)*0.5*(df[i+1]+df[i]) + spRadSum
}
pi*spRadSum #gives 365.1156 w.m^2 approx.

#Integration method 2
#R has an integration function, and one using the Trapezium method
library(pracma)
n <- 400
x <- seq(0.1, 40, len = n)
y <- planck(x,288)
pi*trapz(x, y) #gives 365.1156 w.m^2 approx.

