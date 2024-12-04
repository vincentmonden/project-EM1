import numpy as np
from matplotlib import pyplot as pl
import json

#####################################################################################################
#----------------------------------------------------------------------------------------------------
#INPUT
#----------------------------------------------------------------------------------------------------
#####################################################################################################

with open("init.json") as f:
    variables = json.load(f)

J0 = variables["J0"]/1000000 #/1000000 => change input to correct units
jsigma = variables["jsigma"]
tc = 5*jsigma               #units: coulomb
epsilon = 8.854*variables["epsr"]            #megaAmpere
mu = 1.2566/1000000 * variables["mur"]        #teravolt, microsecond

lambdamin = 2 * np.pi * jsigma / 3 / np.sqrt(mu * epsilon) 
deltax = lambdamin * variables["deltax"] #in meter?
deltay = lambdamin * variables["deltay"]
deltat = np.sqrt(epsilon*mu)/np.sqrt(1/deltax**2 + 1/deltay**2)*variables["deltat"]#Units: second?

lenx = variables["xlen"] - variables["xlen"] % deltax #Grootte heel het spel[m]
leny = variables["ylen"] - variables["ylen"] % deltay
lent = min([lenx,leny])*np.sqrt(epsilon*mu)- min([lenx,leny])*np.sqrt(epsilon*mu) % deltat#[microsec]
#verbeter voor dichtste pad met bounce

Jpos = (int(np.floor(variables["jpos"][0]/deltax)), int(np.floor(variables["jpos"][1]/deltay)))
caplist = variables["capa"]#[picofarad]
for cap in caplist:
    cap[0] = int(np.floor(cap[0]/deltax)) #change inputto correct units
    cap[1] = int(np.floor(cap[1]/deltax))

amxpoints = int(lenx/deltax) +1 #aantal xpunten
amypoints = int(leny/deltay) +1 #+1 om beide edges mee te rekenen
amtpoints = int(lent/deltat) +1
#amtpoints = min([amxpoints, amypoints]) #Simulatie stopt op tijd: niet juist, gebruik phase velocity
print(f"The simulation has {amtpoints} timesteps")

assert amxpoints*amypoints*amtpoints < 8000000000, f'These time- and space steps take at least more than 6GB of RAM, so watch out!! Only change this assert when confident in computing power!!'

elecz = np.zeros((amtpoints, amxpoints, amypoints)) #Hierin wordt alle info gestored: 3D-array voor de 3 veranderlijken
magnx = np.zeros((amtpoints, amxpoints, amypoints - 1)) #Enkel deze 3 componenten zijn != 0
magny = np.zeros((amtpoints, amxpoints - 1, amypoints)) #magnetisch veld -1 omdat de edges vastgepinned zijn op E -> H tussen de E's en dus 1 minder datapunt

Jsource = J0*np.exp(-1*(np.linspace(0, lent, amtpoints) - tc)**2/(2*jsigma**2)) #Alle berekeningen van de source vectorised voor de loop

d = variables["dist"]

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
    coef3[cap[0]][cap[1]] = 1/(1/coef3[cap[0]][cap[1]] - deltax*cap[2]*d/deltat)
    coef4[cap[0]][cap[1]] = 1/(1/coef4[cap[0]][cap[1]] - deltax*cap[2]*d/deltat)
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


#print(magny[:,57,57][2])
#print(magny[:,57,56][2])
#pl.plot(magny[:,int(np.floor(amxpoints/4)),int(np.floor(amxpoints/4))])
#pl.plot(magny[:,int(np.floor(amxpoints/4)) - 1,int(np.floor(amxpoints/4))], linestyle='dotted')
pl.plot(elecz[:,127,127])
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
# -> calculate voltages
# -> surface transfer impedance