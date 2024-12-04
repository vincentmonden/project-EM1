init.json form:

epsr: relative permittivity
mur: relative permeability

xlen: length plate along x-axis (in meter)
ylen: length plate along y-axis (in meter)

J0: amplitude source term (in ampere)
jsigma: width source (in microseconds)
jpos: position source ([xposition, yposition] in meter)

deltax, deltay: stepsize as a multiple of lambda_min (0.05 means stepsize=0.05*lambda_min)
deltat: stepsize as a multiple of courant limit (0.5 means stepsize=0.5*courant limit)

capa: input for capacitors, add lists inside this list to include more capacitors. 
    The form: [xposition(meter), yposition(meter), capacitance(pico farad)]

dist: distance between the plates (in meter)