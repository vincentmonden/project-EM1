import numpy as np
from matplotlib import pyplot as pl

global epsilon, mu

J0 = 1
jsigma = 0.5
tc = 5*jsigma
epsilon = 1
mu = 1

lambdamin = 2 * np.pi * jsigma / 3 / np.sqrt(mu * epsilon) 
deltax = lambdamin/24 #in meter?
deltay = lambdamin/24
deltat = np.sqrt(epsilon*mu)/np.sqrt(1/deltax**2 + 1/deltay**2)#Units: second?

print(deltax, deltat)
 

lenx = 5 - 5 % deltax #Grootte heel het spel
leny = 5 - 5 % deltay
lent = 100 - 100 % deltat

Jpos = (int(np.floor(lenx/deltax/2)), int(np.floor(lenx/deltay/2)))

amxpoints = int(lenx/deltax) +1 #aantal xpunten
amypoints = int(leny/deltay) +1 #+1 om beide edges mee te rekenen
amtpoints = int(lent/deltat) +1

elecz = np.zeros((amtpoints, amxpoints, amypoints)) #Hierin wordt alle info gestored: 3D-array voor de 3 veranderlijken
magnx = np.zeros((amtpoints, amxpoints, amypoints - 1)) #Enkel deze 3 componenten zijn != 0
magny = np.zeros((amtpoints, amxpoints - 1, amypoints)) #magnetisch veld -1 omdat de edges vastgepinned zijn op E -> H tussen de E's en dus 1 minder datapunt

Jsource = J0*np.exp(-1*(np.linspace(0, lent, amtpoints) - tc)**2/(2*jsigma**2)) #Alle berekeningen van de source vectorised voor de loop
#pl.plot(Jsource)
#pl.show()
print(Jsource)
#Allemaal test prints:
#print(np.shape(elecz[0][1:-1,1:-1]))

#print(np.shape(np.tile(elecz[:][1][1] - elecz[:][2][2], (amxpoints,1))))
test = np.array([[1,2,3],[4,5,6],[7,8,9]])
#print(np.array([elecz[0].T[:][i] + elecz[0].T[:][i+1] for i in range(amxpoints - 1)]))
print(test[1,:])

#Hier weer voor echt

#Deze waardes veranderen niet doorheen de simulatie -> Maar 1 keer berekenen
coef1 = (deltat/(deltax * mu))
coef2 = (deltat/(deltay * mu))
coef3 = (deltat/(deltax * epsilon))
coef4 = (deltat/(deltay * epsilon))
print(amtpoints)
for tijdstip in range(1, 30):#amtpoints): 
    #if tijdstip < 50:
        #print(tijdstip)
        
    magny[tijdstip][:][:] = magny[tijdstip-1][:][:] + coef1 * (np.array([elecz[tijdstip-1][i+1,:] - elecz[tijdstip-1][i,:] for i in range(amypoints - 1)]))#CHECK deze array nog eenz
    magnx[tijdstip][:][:] = magnx[tijdstip-1][:][:] - coef2 * (np.array([elecz[tijdstip-1][:,i+1] - elecz[tijdstip-1][:,i] for i in range(amxpoints - 1)])).T 
    elecz[tijdstip][1:-1,1:-1] = elecz[tijdstip-1][1:-1,1:-1] + coef3 * np.array([magny[tijdstip][i+1,1:-1] - magny[tijdstip][i,1:-1] for i in range(amxpoints - 2)])
    elecz[tijdstip][1:-1,1:-1] -= coef4 * np.array([magnx[tijdstip][1:-1, i+1] - magnx[tijdstip][1:-1,i] for i in range(amypoints - 2)]).T
    elecz[tijdstip][Jpos] -= Jsource[tijdstip]
print(magny[:,57,57][2])
print(magny[:,57,56][2])
pl.plot(magnx[:,57,57])
pl.plot(-magnx[:,57,56], linestyle='dotted')
#pl.plot(magny[:,57,57])
pl.xlim([0,30])
pl.show()

# Amai wat een bevalling die arrays:
# Per tijdsstap doen we alleberekeningen ineens, omdat er 1 minder datapunt is voor het magnetsch veld doen we [:-1,blabla]
# Omdat er nog 1 minder datapunt voor het electrisch veld berekend kan worden (leapfrog heeft een datapunt links en rechts nodig) wordt er een matrix optelling
# gedaan van binnen uit. Op deze manier zijn de boundary counditions ook ineens vervuld!
# Check zeker zelf nog eens of er nergens x en y ofzo omgewisseld zijn. Het is allemaal nogal verwarrend en een fout is snel gemaakt!
#
#
# TO DO:
# -> check array optelling enzo
# -> Add source terms (+ check leapfrog uit verandering va stapgrootte)
# -> Add capacitors
# -> ...