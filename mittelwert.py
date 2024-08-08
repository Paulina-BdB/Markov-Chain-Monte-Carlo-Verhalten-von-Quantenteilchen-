# -*- coding: utf8 -*- #
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from random import uniform


def action(x0, x1, x2):
    m0 = 1                          #Parameter m0 von V(x)
    u = 1.0                         #Parameter mü von V(x)
    l = 0                           #Parameter lambda von V(x)
    a = 0.01                         #Gitterweite
    v1 = 0.5*(u*x1)**2 + l*(x1**4)
    v0 = 0.5*(u*x0)**2 + l*(x0**4)
    s = (0.5/a)*m0*(x2-x1)**2 + a*v1+(0.5/a)*m0*(x1-x0)**2 + a*v0
        
    return s

n = 200
a = 0.01
u = 1.0
r = 1 + ((a*u)**2)*0.5 - a*u*sqrt(1 + ((a*u)**2)*0.25)
erwartungsWert = (1+r**n)/(2*u*sqrt(1+((a*u)**2)*0.25)*(1-r**n))
print(erwartungsWert)

def metropolis(xj, numSteps, n, d):
    xSequence = []
    akzeptanz =[]
    i=1
    for i in range(numSteps):
        xNeu = [0]*n
    
        for k in range(len(xj)):
            delta = uniform(-d, d)
            r = uniform(0.0, 1.0)            
            xNeu[k] = xj[k] + delta
            
            if k == 0:
                if action(xj[len(xj)-1], xNeu[k], xj[k+1])< action(xj[len(xj)-1], xj[k], xj[k+1]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[len(xj)-1], xj[k], xj[k+1])-action(xj[len(xj)-1], xNeu[k], xj[k+1]))>r:
                    xj[k]=xNeu[k]
                    akzeptanz.append(1)
                else:
                    akzeptanz.append(0)
                
            elif k < (len(xj)-1):
                if action(xj[k-1], xNeu[k], xj[k+1])< action(xj[k-1], xj[k], xj[k+1]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[k-1], xj[k], xj[k+1])-action(xj[k-1], xNeu[k], xj[k+1]))>r:
                    xj[k]=xNeu[k]
                    akzeptanz.append(1)
                else:
                    akzeptanz.append(0)
                    
            else:
                if action(xj[k-1], xNeu[k], xj[0])< action(xj[k-1], xj[k], xj[0]):
                    xj[k] = xNeu[k]
                    akzeptanz.append(1)
                elif exp(action(xj[k-1], xj[k], xj[0])-action(xj[k-1], xNeu[k], xj[0]))>r:
                    xj[k]=xNeu[k]
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
        mw = mw + x
        mwSquared = mwSquared + (x**2)
        
    return mw/float(n), mwSquared/float(n)


#Mittelwert jedes Pfades aus Metropolis
def mittelwertMetropolis(numSteps, n, delta):
    xStart = [0]*n
    metropXStart = metropolis(xStart, numSteps, n, delta)[0]           #liste mit Vektoren (messwerte der Pfade)
    
    mwSquared = []
    mw = []
    for x in metropXStart:
        mwPfad = mittelwert(x)[0]
        mwPfadSquared = mittelwert(x)[1]
        mw.append(mwPfad)
        mwSquared.append(mwPfadSquared)
    return mw, mwSquared


def autokorrelation(numSteps, n, delta):
    xStart = [1]*n
    mwXSquared = mittelwertMetropolis(numSteps, n, delta)[1]
    mwGesamt = mittelwertGesamt(mwXSquared)
    #pfade = metropolis(xStart, numSteps, n, delta)[0]
    #messpktPfad = mittelwertMetropolis(numSteps, n, delta)
    fehlermwXSquared = sqrt(varianz(mwXSquared, mwGesamt, 0)/(numSteps-1))


    gammaSchaetz = []
    gammaLog = []
    fehler = []
    fehlerLog = []
    autokorNenner = 0
    '''
    for i in range(len(mwXSquared)):
        faktor2 = mwXSquared[i] - mwGesamt
        autokorNenner = autokorNenner + faktor2*faktor2
    '''            
    for t in range(len(mwXSquared)-1):
        autokor = 0
        autokorZaehler = 0
        messpkte = [0]*len(mwXSquared)
        fehlerLogT = 0
        #print(messpkte)
        for i in range(len(mwXSquared)-t):
                faktor1 = mwXSquared[i+t] - float(mwGesamt)
                faktor2 = mwXSquared[i] - float(mwGesamt)
                messpkte[i] = faktor1*faktor2
                #messpkte.append(faktor1*faktor2)
                #print(messpkte)
                autokor = autokor + faktor1*faktor2/float(len(mwXSquared)-t)
                #autokor = autokorZaehler
                #Fehler von log(Gamma)        
                fehlerLogZaehler = faktor1
        print(messpkte)
        fehlerLogNenner = autokorZaehler
        fehlerLogT = fehlerLogT + fehlermwXSquared*fehlerLogZaehler/fehlerLogNenner
        fehlerLog.append(fehlerLogT)
        #Fehler von Gamma
        #print(messpkte)
        var = varianz(messpkte, autokor, 50)
        fehler.append(sqrt(var/(len(mwXSquared)-t-51)))
        #Gamma
        gammaSchaetz.append(autokor)
        #log(Gamma)
        gammaLog.append(log(autokor))
        
    #tauAuto = (len(mwXSquared)-50-50)/float((gammaSchaetz[len(mwXSquared)-50]-gammaSchaetz[50]))
    #print(tauAuto)

    return gammaSchaetz, gammaLog, fehler, fehlerLog

    
#LINEARE REGRESSION: find f(x) = a+bx
def bestFitParameters(sigma, x, y):
    s = 0
    sx = 0
    sy = 0
    sxx = 0
    sxy = 0
    #print(sigma)
    for i , j, k in zip(sigma, x, y):
        s = s + 1/float(i)**2
        sx = sx + j/(i**2)
        sxx = sxx + (j/i)**2
        sy = sy + k/(i**2)
        sxy = sxy + j*k/(i**2)
        
    delta = sxx - sx**2    
    a = (sxx-sxy)/delta
    b = (s*sxy-sx*sy)/delta
    #print(delta/(s*sxy-sx*sy))

        
    return a, b
    

def akzeptanzMW(n, d, numSteps):
    xStart = [0]*n
    akz = metropolis(xStart, numSteps, n, d)[1]
    akzeptanzRate = mittelwert(akz)[0]

    return akzeptanzRate

#TABELLE AKZEPTANZRATE
def tabelleAkzeptanzrate(intervallD):
    akzeptanz = []
    for i in intervallD:
        akzeptanz.append(akzeptanzMW(200, i, 5000))
                     
    for i in range(len(intervallD)):
        g = intervallD[i]
        print("Intervall Delta: {0}, Akzeptanz Rate: {1}".format([-g, g],akzeptanz[i]))


def mittelwertGesamt(messpkt):
    mwGesamt = 0
    for i in range(30, len(messpkt)):
        mwGesamt = mwGesamt + messpkt[i]

    return mwGesamt/float(len(messpkt))

def varianz(messpkt, mw, start):
    n = len(messpkt)
    var=0
    mwSquared=0
    
    for i in range(start, len(messpkt)):
        var = var + ((messpkt[i] - mw)**2)/float(n)
        
    return var

#PLOTS VON MITTELWERTE UND AUTOKORRELATION (GAMMA)
def plottenMittelwerte(numSteps, n, d):
    
    xWerte = arange(numSteps)
    yWerte = mittelwertMetropolis(numSteps, n, d)[0]
    yWerteSquared = mittelwertMetropolis(numSteps, n, d)[1]
        
    plt.plot(xWerte, yWerte, label = 'x')
    plt.plot(xWerte, yWerteSquared, label = 'x**2')
    plt.xlabel('Pfad j')
    plt.ylabel('mittelwert pro Pfad, delta='+str(d))
    plt.legend()
    plt.show()
    
def plotAutokorre(numSteps, n, delta):
    xAchse = range(numSteps-1)
    gammaSchaetz = autokorrelation(numSteps, n, delta)[0]
    gammaLog = autokorrelation(numSteps, n, delta)[1]
    gammaFehler = autokorrelation(numSteps, n, delta)[2]
    #upperlimits = array([1, 0] * 5)
    #lowerlimits = array([0, 1] * 5)
    fehler = gammaFehler
    #error = 0.1 + 0.2 * x
    #fig, ax0 = plt.subplots(nrows=2, sharex=True)    
    plt.errorbar(xAchse, gammaSchaetz, yerr = fehler, label = 'Gamma +- Fehler')
    #ax0.set_title('variable, symmetric error')
    
    #plt.plot(xAchse, gammaSchaetz,  label = 'Gamma Schaetzwert')
    #plt.plot(xAchse, gammaLog, label = 'log von Gamma')
    #plt.plot(xAchse, gammaFehler, label = 'Fehler Gamma')
    plt.xlabel('t')
    plt.legend()
    plt.show()


deltaInt = [100, 50, 10, 2, 1, 0.5, 0.25, 0.05, 0.025]

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
            

#ausgabeMwFehler(200, 200)
plottenMittelwerte(10000, 200, 2)
#plotAutokorre(1000, 200, 2)
#autokorrelation(400, 200, 2)
#print("-----------------------------------------------------------------------------")
#print(autokorrelation(500, 200, 2)[2])

#fehler = autokorrelation(200, 200, 2)[3]
#xWerte = range(200)
#yWerte = autokorrelation(200, 200, 2)[1]
#print(bestFitParameters(fehler, xWerte, yWerte))

#mwTauAuto = [plotAutokorre(200, 200, 2*sqrt(0.001)), plotAutokorre(1000, 200, 2*sqrt(0.001)), plotAutokorre(5000, 200, 2*sqrt(0.001)), plotAutokorre(3000, 200, 2*sqrt(0.001)), plotAutokorre( 2000, 200, 2*sqrt(0.001))]

#print(mittelwert(mwTauAuto))
#print((-13929.193051-159276.756138-112324.979615-190288.298584)/4)

#tabelleAkzeptanzrate(deltaInt)
