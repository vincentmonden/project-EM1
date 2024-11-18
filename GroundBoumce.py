import numpy as np

global epsilon, mu

epsilon = 1
mu = 1

deltax = 0.1 #Units: meter??
deltay = 0.1
deltat = 0.1 #Units: second?

lenx = 100 #Grootte heel het spel
leny = 100
lent = 100

amxpoints = int(lenx/deltax) +1 #aantal xpunten
amypoints = int(leny/deltay) +1 #+1 om beide edges mee te rekenen
amtpoints = int(lent/deltat) +1

elecz = np.zeros((amtpoints, amxpoints, amypoints)) #Hierin wordt alle info gestored: 3D-array voor de 3 veranderlijken
magnx = np.zeros((amtpoints, amxpoints - 1, amypoints -1)) #Enkel deze 3 componenten zijn != 0
magny = np.zeros((amtpoints, amxpoints-1, amypoints-1))

#Allemaal test prints:
test = np.array(([1,2,3,4,5],[6,7,8,9,0]))
print(test[0][:])

elecz[0][1][1] = 7777
print(np.shape(elecz[0][1:-1,1:-1]))

#print(np.shape(np.tile(elecz[:][1][1] - elecz[:][2][2], (amxpoints,1))))
print(np.array([elecz[0].T[:][i] + elecz[0].T[:][i+1] for i in range(amxpoints - 1)]))
print(np.shape(np.array([elecz[0][:-1,i] + elecz[0][:-1,i+1] for i in range(amxpoints - 1)])))

#Hier weer voor echt

for tijdstip in range(2, 2*amtpoints): #maal 2 voor de halfjes in de leapfrog: real t is dus tijdstip/2
    #print(tijdstip)
    if tijdstip % 2 == 0: 
        magny[tijdstip][:][:] = magny[tijdstip - 2][:][:] + (deltat/(deltax * mu)) * (np.array([elecz[tijdstip - 1][:-1,i+1] - elecz[tijdstip-1][:-1,i] for i in range(amxpoints - 1)])) #CHECK deze array nog eenz
        magnx[tijdstip][:][:] = magnx[tijdstip - 2][:][:] - (deltat/(deltay * mu)) * (np.array([elecz[tijdstip-1][i+1,:-1] - elecz[tijdstip-1][i,:-1] for i in range(amypoints - 1)]))
    else:
        elecz[tijdstip][1:-1,1:-1] = elecz[tijdstip - 2][1:-1,1:-1] + (deltat/(deltax * epsilon)) * np.array([magny[tijdstip - 1][:-1,i+1] - magny[tijdstip-1][:-1,i] for i in range(amxpoints - 2)])
        elecz[tijdstip][1:-1,1:-1] -= (deltat/(deltay * epsilon)) * np.array([magny[tijdstip - 1][i+1,:-1] - magny[tijdstip-1][i,:-1] for i in range(amypoints - 2)])
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