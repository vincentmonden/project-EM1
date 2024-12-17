import numpy as np
from matplotlib import pyplot as pl
import json
from scipy.spatial import distance
from scipy.special import hankel2
import gc
import sys

#####################################################################################################
#----------------------------------------------------------------------------------------------------
#INPUT
#----------------------------------------------------------------------------------------------------
#####################################################################################################

with open("init.json") as f:
    variables = json.load(f)

J0 = variables["J0"]/1000000 #/1000000 => change input to correct units
jsigma = variables["jsigma"]
tc = 8*jsigma               #units: coulomb
epsilon = 8.854*variables["epsr"]            #megaAmpere
mu = 1.2566/1000000 * variables["mur"]        #teravolt, microsecond

d = variables["dist"]
lambdamin = 2 * np.pi * jsigma / 3 / np.sqrt(mu * epsilon) 
deltax = lambdamin * variables["deltax"] #in meter?
deltay = lambdamin * variables["deltay"]
deltat = np.sqrt(epsilon*mu)/np.sqrt(1/deltax**2 + 1/deltay**2)*variables["deltat"]#Units: second?


Jpos = (int(np.floor(variables["jpos"][0]/deltax)), int(np.floor(variables["jpos"][1]/deltay)))

lenx = variables["xlen"] - variables["xlen"] % deltax #Grootte heel het spel[m]
leny = variables["ylen"] - variables["ylen"] % deltay

lengths = []
meetpos = []
for voltpos in variables["voltpos"]:
    meetpos.append((int(np.floor(voltpos[0]/deltax)), int(np.floor(voltpos[1]/deltay))))
    lengths.append(np.min([distance.euclidean(variables["jpos"], [2*lenx - voltpos[0], voltpos[1]]), distance.euclidean(variables["jpos"], [voltpos[0], 2*leny - voltpos[1]]), distance.euclidean(variables["jpos"], [voltpos[0], - voltpos[1]]), distance.euclidean(variables["jpos"], [- voltpos[0], voltpos[1]])]))

lent = np.min(lengths)*np.sqrt(epsilon*mu)- np.min(lengths)*np.sqrt(epsilon*mu) % deltat#[microsec]
#dichtste pad met bounce


caplist = variables["capa"]#[picofarad]
for cap in caplist:
    cap[0] = int(np.floor(cap[0]/deltax)) #change inputto correct units
    cap[1] = int(np.floor(cap[1]/deltax))

amxpoints = int(lenx/deltax) +1 #aantal xpunten
amypoints = int(leny/deltay) +1 #+1 om beide edges mee te rekenen
amtpoints = int(lent/deltat) -30
print(f"The simulation has {amtpoints} timesteps")

assert amxpoints*amypoints*amtpoints < 8000000000, f'These time- and space steps take at least more than 6GB of RAM, so watch out!! Only change this assert when confident in computing power!!'

elecz = np.zeros((2, amxpoints, amypoints), dtype=np.float64) #Hierin wordt alle info gestored: 3D-array voor de 3 veranderlijken
magnx = np.zeros((2, amxpoints, amypoints - 1), dtype=np.float64) #Enkel deze 3 componenten zijn != 0
magny = np.zeros((2, amxpoints - 1, amypoints), dtype=np.float64) #magnetisch veld -1 omdat de edges vastgepinned zijn op E -> H tussen de E's en dus 1 minder datapunt

Jsource = J0*np.exp(-1*(np.linspace(0, lent, amtpoints) - tc)**2/(2*jsigma**2)) #Alle berekeningen van de source vectorised voor de loop
#Jsource = np.zeros(amtpoints)
#Jsource[1] = 1000000000
#############################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################################################################################

#Deze waardes veranderen niet doorheen de simulatie -> Maar 1 keer berekenen

coef1 = (deltat/(deltax * mu))
coef2 = (deltat/(deltay * mu))
coef3 = np.full((amxpoints-2,amypoints-2),(deltat/(deltax * epsilon))) #niet ideaal om een array te maken voor alle coefficienten, maar ik zie niet drect een manier om dit te vereenvoudigen zonder for-loops te itroduceren
coef4 = np.full((amxpoints-2,amypoints-2),(deltat/(deltay * epsilon)))
coef5 = deltat / epsilon / deltax / deltay
plotvaluesE = np.zeros((len(meetpos), amtpoints))

for cap in caplist:
    coefhulp = epsilon - (cap[2]*d / deltax/deltay)
    coef3[cap[0]][cap[1]] = deltat/(deltax * coefhulp)
    coef4[cap[0]][cap[1]] = deltat/(deltay * coefhulp)
    print(deltax/deltat)
    print(coef4[cap[0]][cap[1]])
k=0
for tijdstip in range(1, amtpoints): 
    hulp=[]
    if tijdstip % 50==0:
        print(tijdstip)
    gc.collect()
    
    magny[1][:][:] = magny[0][:][:] + coef1 * (np.array([elecz[0][i+1,:] - elecz[0][i,:] for i in range(amxpoints - 1)]))#CHECK deze array nog eenz
    magnx[1][:][:] = magnx[0][:][:] - coef2 * (np.array([elecz[0][:,i+1] - elecz[0][:,i] for i in range(amypoints - 1)])).T 
    elecz[1][1:-1,1:-1] = elecz[0][1:-1,1:-1] + coef3 * np.array([magny[1][i+1,1:-1] - magny[1][i,1:-1] for i in range(amxpoints - 2)])
    elecz[1][1:-1,1:-1] -= coef4 * np.array([magnx[1][1:-1, i+1] - magnx[1][1:-1,i] for i in range(amypoints - 2)]).T
    elecz[1][Jpos] -= coef5 * Jsource[tijdstip]
    
    for i, pos in enumerate(meetpos):
        plotvaluesE[i,tijdstip] = elecz[1][pos]
    
    elecz[0][:][:] = elecz[1][:][:]
    magnx[0][:][:] = magnx[1][:][:]
    magny[0][:][:] = magny[1][:][:]

#############################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#DATA-ANALYSIS
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################################################################################
print(plotvaluesE)
timelist = deltat*np.linspace(1, amtpoints, amtpoints)
pl.rcParams['text.usetex'] = True


#energy = np.zeros(amtpoints)
#for tijdstip in range(len(energy)):
#    if tijdstip % 50==0:
#        print(tijdstip)
#    energy[tijdstip]+=mu/2*np.sum(np.array([magny[tijdstip][i+1,1:-1]/2 + magny[tijdstip][i,1:-1]/2 for i in range(amxpoints - 2)])**2 + np.array([magnx[tijdstip][1:-1, i+1]/2 + magnx[tijdstip][1:-1,i]/2 for i in range(amypoints - 2)]).T**2)
#    energy[tijdstip]+=epsilon/2*np.sum(elecz[tijdstip]**2)
#    energy[tijdstip]*=deltax*deltay
#print(energy[100])
#pl.plot(timelist[25:], energy[25:]*10e12, color='navy')
##pl.ylabel(r"E[J]")
#pl.xlabel(r"$t[\mu s]$")
#pl.yscale('log')
#pl.show()
#exit()

#print(magny[:,57,57][2])
#print(magny[:,57,56][2])
#pl.plot(magny[:,int(np.floor(amxpoints/4)),int(np.floor(amxpoints/4))])
#pl.plot(magny[:,int(np.floor(amxpoints/4)) - 1,int(np.floor(amxpoints/4))], linestyle='dotted')
what = variables['what']
if 'volt' in what:
    for i, pos in enumerate(meetpos):
        pl.plot(timelist,-d*plotvaluesE[i]*1000000000000, color='navy')
        pl.ylabel(r"$\hat{v}[V]$")
        pl.xlabel(r"$t[\mu s]$")
        #pl.yscale('log')
        #pl.plot(magny[:,57,57])
        #pl.xlim([0,30])
        pl.show()

if 'Z' in what:
    for i,pos in enumerate(meetpos):
        espec = np.fft.fft(plotvaluesE[i], n=100*amtpoints)*deltat
        vspec = - d * espec
        omegas = np.fft.fftfreq(100*amtpoints, deltat)*2*np.pi
        Jspec = np.fft.fft(Jsource, n=100*amtpoints)*deltat 
        omegas2 = np.linspace(0,3/jsigma,100)
        omega=100000
        print(len(omegas))
    
        pl.plot(omegas, np.real(- omegas * Jspec * mu *hankel2(0, omegas*np.sqrt(epsilon*mu)*distance.euclidean(pos, Jpos))/4))
        #pl.plot(omegas, np.real(espec))
        #print(np.fft.fft(plotvaluesE[i], n=100*amtpoints)[np.argmin(np.abs(omegas - omega))]*deltat)
        #print(- omega * mu *hankel2(0, omega*np.sqrt(epsilon*mu)*distance.euclidean(pos, Jpos))/4)
        Zspec = vspec/Jspec
        print(len(Jspec))
        #pl.xlim([0, 3/jsigma])
        pl.xlim([0,200000])
        #pl.ylim([2e-5,2e-1])
        #pl.ylim([np.min(np.real(Zspec))+0.5,np.max(np.real(Zspec))+100000])
        #pl.yscale('log')
        #pl.plot(omegas, np.real(omegas*np.sqrt(epsilon*mu)*d/4 * hankel2(0,omegas*np.sqrt(epsilon*mu)*distance.euclidean(pos, Jpos))))
        #pl.plot(omegas, np.real(Zspec)/np.sqrt(mu/epsilon))
        #pl.plot(omegas, np.imag(Zspec)/np.sqrt(mu/epsilon))
        pl.show()


# Amai wat een bevalling die arrays:
# Per tijdsstap doen we alleberekeningen ineens, omdat er 1 minder datapunt is voor het magnetsch veld doen we [:-1,blabla]
# Omdat er nog 1 minder datapunt voor het electrisch veld berekend kan worden (leapfrog heeft een datapunt links en rechts nodig) wordt er een matrix optelling
# gedaan van binnen uit. Op deze manier zijn de boundary counditions ook ineens vervuld!
# Check zeker zelf nog eens of er nergens x en y ofzo omgewisseld zijn. Het is allemaal nogal verwarrend en een fout is snel gemaakt!
#
#
# TO DO:
# -> source terms + check leapfrog uit verandering va stapgrootte)
# -> max t verbeteren (geen bouncy)
# -> surface transfer impedance
# -> 