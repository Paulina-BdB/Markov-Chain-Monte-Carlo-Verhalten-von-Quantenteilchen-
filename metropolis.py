# -*- coding: utf8 -*- #
from numpy import *
import matplotlib.pyplot as plt
from random import uniform

def action(x0, x1, x2):
    m0 = 1                          #Parameter m0 von V(x)
    u = 1.0                         #Parameter mü von V(x)
    l = 0                           #Parameter lambda von V(x)
    a = 1                        #Gitterweite
    v1 = 0.5*(u*x1)**2 + l*(x1**4)
    v0 = 0.5*(u*x0)**2 + l*(x0**4)
    s = (0.5/a)*m0*((x2-x1)**2 + (x1-x0)**2) + a*(v0+v1)
        
    return s

def metropolis(numSteps, n, d):
    xj = [0]*n
    xSequence = []
    akzeptanz =[]
    i=1
    for i in range(numSteps):
        xNeu = [0]*n
    
        for k in range(len(xj)):
            delta = uniform(-d, d)
            r = uniform(0.0, 1.0)            
            xNeu[k] = xj[k] + delta
            #print(delta)
            #print("--------------------")
            #print(xNeu[k])
            if k == 0:
                if action(xj[len(xj)-1], xNeu[k], xj[k+1])< action(xj[len(xj)-1], xj[k], xj[k+1]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[len(xj)-1], xj[k], xj[k+1])-action(xj[len(xj)-1], xNeu[k], xj[k+1]))>r:
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                else:
                    akzeptanz.append(0)
                
            elif k < (len(xj)-1):
                if action(xj[k-1], xNeu[k], xj[k+1])< action(xj[k-1], xj[k], xj[k+1]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[k-1], xj[k], xj[k+1])-action(xj[k-1], xNeu[k], xj[k+1]))>r:
                    xj[k]= xNeu[k]
                    akzeptanz.append(1)
                else:
                    akzeptanz.append(0)
                    
            else:
                if action(xj[k-1], xNeu[k], xj[0])< action(xj[k-1], xj[k], xj[0]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[k-1], xj[k], xj[0])-action(xj[k-1], xNeu[k], xj[0]))>r:
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                else:
                    akzeptanz.append(0)
                
        xSequence.append(list(xj))
    return xSequence, akzeptanz


def mittelwert(xj):
    n = len(xj)
    mw = 0
    mwSquared = 0
    
    for x in xj:
        mw = mw + x/float(n)
        mwSquared = mwSquared + (x**2)/float(n)
    #print(mw)
        
    return mw, mwSquared


def mittelwertMetropolis(numSteps, n, delta):
    mwSquared = []
    mw = []
    for x in metropolis(numSteps, n, delta)[0]:
        mw.append(mittelwert(x)[0])
        mwSquared.append(mittelwert(x)[1])
        
    return mw, mwSquared


#DATEI VON METROPOLS PFADEN
def dateiMittelwerte(numSteps, n, delta):
    mittelwerteAbgeschnitten = []
    mittelwerte = mittelwertMetropolis(numSteps, n, delta)[1]
    
    for i in range(200, numSteps):
        mittelwerteAbgeschnitten.append(mittelwerte[i])
        
    datei = savetxt("mw"+str(numSteps)+"_"+str(delta)+"_1.txt", mittelwerteAbgeschnitten)
    
    return datei
 
 
def dateiAkzeptanzRate(numSteps, n, delta):
    akzeptanz = metropolis(numSteps, n, delta)[1]
    akzeptanzAbgeschnitten = []
    for i in range(200, len(akzeptanz)):
        akzeptanzAbgeschnitten.append(akzeptanz[i])
    #print akzeptanzAbgeschnitten    
    akzeptanzRate = mittelwert(akzeptanzAbgeschnitten)
    print(akzeptanzRate[0])
    akzRate = savetxt("akzR"+str(numSteps)+"_"+str(delta)+"_1.txt", array(akzeptanzRate))
    
    #return akzRate


#DELTA OPTIMAL
def deltaOptimal(d):
    dateiAkzeptanzRate(100000, 640, d)
    akz = loadtxt("akzR100000_"+str(d)+"_0.03125.txt")[0]
    
    if akz <= 0.462 and akz >= 0.452:
        return d
    elif akz >0.46:
        d = d/4.0 + d
        return deltaOptimal(d)

    else:
        d = (2*d)/3.0
        return deltaOptimal(d)

    
    
#Mittelwert Plotten
def plottenMittelwerte(numSteps, n, d):
    
    xWerte = arange(numSteps)
    yWerte = mittelwertMetropolis(numSteps, n, d)[0]
    yWerteSquared = mittelwertMetropolis(numSteps, n, d)[1]
        
    plt.plot(xWerte, yWerte, label = 'x')
    plt.plot(xWerte, yWerteSquared, label = 'x**2')
    plt.xlabel('Pfad j')
    plt.ylabel('Mittelwert')
    plt.title('Mittelwert jedes Pfades')
    plt.legend()
    plt.show()


#print "Intervall Delta: {0}, Akzeptanz Rate: {1}".format([-g, g],akzeptanz[i])
        

#MITTELWERT ALLER PFADE UND FEHLER
def ausgabeMwFehler(numSteps, n):
    for delta in deltaInt:
        messpkt1, messpkt2 = mittelwertMetropolis(numSteps, n, delta)
        print("---------------------------------------------------------------")
        print(" delta aus [{0}, {1}]".format(-delta, delta))
        mwGesX, mwGesXSquared = mittelwertGesamt(messpkt1), mittelwertGesamt(messpkt2)
        print("Mittelwert gesamt X: {0},Mittelwert gesamt X**2: {1}".format(mwGesX, mwGesXSquared))
    
        #FEHLER SCHAETZEN
        fehlerX = sqrt(varianz(messpkt1, mwGesX)/(numSteps-1))
        fehlerXSquared = sqrt(varianz(messpkt2, mwGesXSquared, 30)/(numSteps-1))
        print("Standardabw. X: {0},Standardabw X**2: {1}".format(fehlerX, fehlerXSquared))
        
        
        
def varianz(messpkt, mw):
    n = len(messpkt)
    var=0    
    for i in range(len(messpkt)):
        var = var + ((messpkt[i] - mw)**2)/float(n)
        
    return var

#Erwartungswert von x^2
n = 20
a = 1
u = 1.0
r = 1 + ((a*u)**2)*0.5 - a*u*sqrt(1 + ((a*u)**2)*0.25)
erwartungsWert = (1+r**n)/(2*u*sqrt(1+((a*u)**2)*0.25)*(1-r**n))
print(erwartungsWert)

plottenMittelwerte(1000, 20, 1.92)
#print mittelwert(mittelwertMetropolis(1000, 20, 1.92)[1])

messpkt1, messpkt2 = mittelwertMetropolis(1000, 20, 1.92)
mwGesX, mwGesXSquared = mittelwert(messpkt1)[0], mittelwert(messpkt2)[0]

print(sqrt(varianz(messpkt1, mwGesX)/(1000-1)))
print(sqrt(varianz(messpkt2, mwGesXSquared)/(1000-1)))


dateiMittelwerte(100000,640, 0.42)
dateiAkzeptanzRate(10000, 20, 5)

#d = uniform(0.0625, 0.625)
#print deltaOptimal(d)


                            
