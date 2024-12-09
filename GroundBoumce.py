import numpy as np
from matplotlib import pyplot as pl
import json
from scipy.spatial import distance

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

voltpos = (int(np.floor(variables["voltpos"][0]/deltax)), int(np.floor(variables["voltpos"][1]/deltay)))
Jpos = (int(np.floor(variables["jpos"][0]/deltax)), int(np.floor(variables["jpos"][1]/deltay)))

lenx = variables["xlen"] - variables["xlen"] % deltax #Grootte heel het spel[m]
leny = variables["ylen"] - variables["ylen"] % deltay
lengtht = np.min([distance.euclidean(variables["jpos"], [2*lenx - variables["voltpos"][0], variables["voltpos"][1]]), distance.euclidean(variables["jpos"], [variables["voltpos"][0], 2*leny - variables["voltpos"][1]]), distance.euclidean(variables["jpos"], [variables["voltpos"][0], - variables["voltpos"][1]]), distance.euclidean(variables["jpos"], [- variables["voltpos"][0], variables["voltpos"][1]])])
print(lengtht)
lent = lengtht*np.sqrt(epsilon*mu)- lengtht*np.sqrt(epsilon*mu) % deltat#[microsec]
#verbeter voor dichtste pad met bounce


caplist = variables["capa"]#[picofarad]
for cap in caplist:
    cap[0] = int(np.floor(cap[0]/deltax)) #change inputto correct units
    cap[1] = int(np.floor(cap[1]/deltax))

amxpoints = int(lenx/deltax) +1 #aantal xpunten
amypoints = int(leny/deltay) +1 #+1 om beide edges mee te rekenen
amtpoints = int(lent/deltat) -30
print(f"The simulation has {amtpoints} timesteps")

assert amxpoints*amypoints*amtpoints < 8000000000, f'These time- and space steps take at least more than 6GB of RAM, so watch out!! Only change this assert when confident in computing power!!'

elecz = np.zeros((amtpoints, amxpoints, amypoints)) #Hierin wordt alle info gestored: 3D-array voor de 3 veranderlijken
magnx = np.zeros((amtpoints, amxpoints, amypoints - 1)) #Enkel deze 3 componenten zijn != 0
magny = np.zeros((amtpoints, amxpoints - 1, amypoints)) #magnetisch veld -1 omdat de edges vastgepinned zijn op E -> H tussen de E's en dus 1 minder datapunt

Jsource = J0*np.exp(-1*(np.linspace(0, lent, amtpoints) - tc)**2/(2*jsigma**2)) #Alle berekeningen van de source vectorised voor de loop

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

for cap in caplist:
    coefhulp = epsilon - (cap[2]*d / deltax/deltay)
    coef3[cap[0]][cap[1]] = deltat/(deltax * coefhulp)
    coef4[cap[0]][cap[1]] = deltat/(deltay * coefhulp)
    print(deltax/deltat)
    print(coef4[cap[0]][cap[1]])

for tijdstip in range(1, amtpoints): 
    if tijdstip % 50==0:
        print(tijdstip)
        
    magny[tijdstip][:][:] = magny[tijdstip-1][:][:] + coef1 * (np.array([elecz[tijdstip-1][i+1,:] - elecz[tijdstip-1][i,:] for i in range(amxpoints - 1)]))#CHECK deze array nog eenz
    magnx[tijdstip][:][:] = magnx[tijdstip-1][:][:] - coef2 * (np.array([elecz[tijdstip-1][:,i+1] - elecz[tijdstip-1][:,i] for i in range(amypoints - 1)])).T 
    elecz[tijdstip][1:-1,1:-1] = elecz[tijdstip-1][1:-1,1:-1] + coef3 * np.array([magny[tijdstip][i+1,1:-1] - magny[tijdstip][i,1:-1] for i in range(amxpoints - 2)])
    elecz[tijdstip][1:-1,1:-1] -= coef4 * np.array([magnx[tijdstip][1:-1, i+1] - magnx[tijdstip][1:-1,i] for i in range(amypoints - 2)]).T
    elecz[tijdstip][Jpos] -= coef5 * Jsource[tijdstip]

#############################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#DATA-ANALYSIS
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################################################################################

pl.rcParams['text.usetex'] = True
#print(magny[:,57,57][2])
#print(magny[:,57,56][2])
#pl.plot(magny[:,int(np.floor(amxpoints/4)),int(np.floor(amxpoints/4))])
#pl.plot(magny[:,int(np.floor(amxpoints/4)) - 1,int(np.floor(amxpoints/4))], linestyle='dotted')
timelist = deltat*np.linspace(0, amtpoints, amtpoints)
pl.plot(timelist,d*elecz[:,voltpos[0],voltpos[1]]*1000000000000, color='navy')
pl.ylabel(r"$\hat{v}[V]$")
pl.xlabel(r"$t[\mu s]$")
#pl.plot(magny[:,57,57])
#pl.xlim([0,30])
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