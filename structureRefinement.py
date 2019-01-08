#!/usr/bin/python


from Bio.PDB import *
import os
import sys
import numpy as np
import matplotlib.pyplot as pl
import time
import formulas_values_mc as fv
import formulas_polynomial_mc as fp
from mpl_toolkits.mplot3d import Axes3D





accA = ["N7", "N1", "N3"]
donA = ["N6"]
A = accA + donA

accU = ["O4", "O2"]
donU = ["N3"]
U = accU + donU

accG = ["N7", "O6", "N3"]
donG = ["N1", "N2"]
G = accG + donG

accC = ["N3", "O2"]
donC = ["N4"]
C = accC + donC

allAt = accA+donA+accC+donC+accG+donG+accU+donU
allDon = donA+donC+donG+donU
allAcc = accA+accC+accG+accU



# findContacts() find contacts to refine; input - NeighborSearch instance, list of residues, distance
def findContacts(ns, resList, distance):


	# makeSingleContacts() looks for similat atom contact (pairs of contacts with the same atom) ant remove that contact which has largest distance between atoms; input - contact
	def makeSingleContact(contacts):

		newC = []
		for n, c in enumerate(contacts):
			for i in contacts[:n]+contacts[n+1:]:
				if c["acc"] == i["acc"]:

					if c["acc"] - c["don"] < i["acc"] - i["don"]:
						newC.append(i)
					elif c["acc"] - c["don"] > i["acc"] - i["don"]:
						newC.append(c)				
					else: pass

				elif c["don"] == i["don"]:

					if c["acc"] - c["don"] < i["acc"] - i["don"]:
						newC.append(i)
					elif c["acc"] - c["don"] > i["acc"] - i["don"]:
						newC.append(c)
					else: pass

				else: pass

		for i in newC:
			try:
				contacts.remove(i)
			except ValueError:
				#print "ValueError - line 67", i
				pass

		return contacts
			

	# find() looks for contacts between res atoms and others atoms	
	def find(res, resAcc, resDon, allAcceptors, allDonors, allAtoms, dist, output, aA, dA, aC, dC, aG, dG, aU, dU):	# funkction looks for contacts beetween particular atoms od nucleobases
		for a in res:
			closeAtoms = []
			if a.id in resAcc+resDon:
				closeAtoms = ns.search(a.coord, dist, level="A")
			else: pass

			for c in closeAtoms:
				if a.parent == c.parent or a.parent.id[1] == c.parent.id[1] + 1 or a.parent.id[1] == c.parent.id[1] - 1:
					pass
				else:

					if a.id in resAcc and c.id in allDonors:
						
						if c.parent.resname.split()[0] == "A" and c.id in dA and {"acc": a,"don": c} not in output:
							output.append({"acc": a,"don": c})
						elif c.parent.resname.split()[0] == "C" and c.id in dC and {"acc": a,"don": c} not in output:
							output.append({"acc": a,"don": c})
						elif c.parent.resname.split()[0] == "G" and c.id in dG and {"acc": a,"don": c} not in output:
							output.append({"acc": a,"don": c})
						elif c.parent.resname.split()[0] == "U" and c.id in dU and {"acc": a,"don": c} not in output:
							output.append({"acc": a,"don": c})
						else: pass 

					elif a.id in resDon and c.id in allAcceptors:

						if c.parent.resname.split()[0] == "A" and c.id in accA and {"acc": c,"don": a} not in output:
							output.append({"acc": c,"don": a})
						elif c.parent.resname.split()[0] == "C" and c.id in accC and {"acc": c,"don": a} not in output:
							output.append({"acc": c,"don": a})
						elif c.parent.resname.split()[0] == "G" and c.id in accG and {"acc": c,"don": a} not in output:
							output.append({"acc": c,"don": a})
						elif c.parent.resname.split()[0] == "U" and c.id in accU and {"acc": c,"don": a} not in output:
							output.append({"acc": c,"don": a})
						else: pass 
					else: pass

	pairs = []

	for r in resList:

		if r.resname.split()[0] == "A":
			find(r, accA, donA, allAcc, allDon, allAt, distance, pairs, accA, donA, accC, donC, accG, donG, accU, donU)	

		elif r.resname.split()[0] == "C":
			find(r, accC, donC, allAcc, allDon, allAt, distance, pairs, accA, donA, accC, donC, accG, donG, accU, donU)

		elif r.resname.split()[0] == "G":
			find(r, accG, donG, allAcc, allDon, allAt, distance, pairs, accA, donA, accC, donC, accG, donG, accU, donU)

		elif r.resname.split()[0] == "U":
			find(r, accU, donU, allAcc, allDon, allAt, distance, pairs, accA, donA, accC, donC, accG, donG, accU, donU)
	
	singlePairs = makeSingleContact(pairs)
	return singlePairs



# playClarna() classyfied contacts; inputs - contact name, RNA structure, ClaRNA thresh value
def playClaRNA(contacts, structure, thresh):
	
	CLARNA_PATH = "/home/mateusz/Studia/ClaRNA_play/clarna_run.py -ipdb "
	os.system(CLARNA_PATH + structure + " > " + structure.split(".")[0]+"_CLARNA -thresh " + str(thresh))
	
	clarna = open(structure.split(".")[0]+"_CLARNA", "r").readlines()[2:-1]
	for c in contacts:
		c["class"] = "nope"
		for l in clarna:
			l = l.split()
			if l[1] == str(c["acc"].parent.id[1]) and l[3] == str(c["don"].parent.id[1]) and l[0] == str(c["acc"].parent.parent.id.split()[0]) and l[2] == str(c["don"].parent.parent.id.split()[0]):
				clas = list(l[7])

				try:
					qind = l[7].index('?')
					clas[qind] = 'q'
				except ValueError:
					pass
				c["class"] = "".join(clas)
				
			elif l[1] == str(c["don"].parent.id[1]) and l[3] == str(c["acc"].parent.id[1]) and l[0] == str(c["don"].parent.parent.id.split()[0]) and l[2] == str(c["acc"].parent.parent.id.split()[0]):
				clas = list(l[7])
				try:
					qind = l[7].index('?')
					clas[qind] = 'q'
				except ValueError:
					pass
				c["class"] = "".join(clas)
			else: pass

	os.system("rm "+structure.split(".")[0]+"_CLARNA")
	
	return contacts



# createName() creates contact name; input - contact, angle or distance
def createName(contact, WHAT):

	name = [contact["acc"].parent.resname.split()[0], contact["acc"].id, contact["don"].id, contact["don"].parent.resname.split()[0], contact["class"], WHAT]
	return "_".join(name)



# callFormula() gives polinomial from ordered set of stored polinomials from value or function method; input - contact name, name of orderes method
def callFormula(name, method):

	if method == "-value":
		return fv.__dict__[name]
	elif method == "-poly":
		return fp.__dict__[name]
	else:
		raise KeyError



# makeFormula() creates polinomial equation 
def makeFormula(name, slow):
	
	polin = np.poly1d(slow["polinomial"])

	x = np.arange(slow["minx"], slow["maxx"], slow["delta"])

	def lookForMaxima(lstM, polinomial, minx, maxx):

		lr = 0
		if name.split("_")[-1] == "distance":
			lr = 0.1
		elif name.split("_")[-1] == "angle":
			lr = 1

		polinDe = polinomial.deriv()

		for i in list(np.roots(polinDe)):

			f = i-lr
			m = i
			l = i+lr

			if polinomial(f) < polinomial(m) and polinomial(l) < polinomial(m) and minx <= m <= maxx:
				maxima.append(i)
			else: pass
	
	maxima = []
	lookForMaxima(maxima, polin, slow["minx"], slow["maxx"])

	print "\t\t\t", maxima
	print

	return maxima


# comparison of maximas values from Monte Carlo 'new_mc' with calculated maximas from polynomial  
def compareMaximas(new_mc, old):

	mb = [float(i) for i in old]
	old = mb

	if len(new_mc) == 2 and len(mb) == 0:
		return []

	elif len(new_mc) == 1 and len(mb) == 1 or len(new_mc) == 1 and len(mb) == 0:
		return new_mc

	elif len(new_mc) == 2 and len(old) == 1:
		if abs(new_mc[0] - old[0]) < abs(new_mc[1] - old[0]):
			return [new_mc[0]]
		elif abs(new_mc[0] - old[0]) > abs(new_mc[1] - old[0]):
			return [new_mc[1]]
		else:
			raise TypeError

	elif len(new_mc) == 2 and len(old) == 2:
		new_mc.sort()
		return new_mc




# refineDistance() gets list of two atoms to refine and proper distance to be set between these atoms; input - atoms, name of contact, dictionary with polinomial
def refineDistanceValue(atoms, name, slow):		

	A, B = [atoms[0], atoms[1]]

	polin = np.poly1d(slow["polinomial"])
	maximas_err = slow["maximas"]	# it is list with maksimas and its errors from Monte Carlo simulation
	maximas_mc = [i[0] for i in maximas_err]

	def checkDifference(value, lst):
		final_dist = 0
		diff = 1000
		for m in lst:
			if abs(m - value) < diff:
				diff = abs(m - value)
				final_dist = m
		return final_dist


	def lookForMaxima(lstM, polinomial, minx, maxx):

		lr = 0
		if name.split("_")[-1] == "distance":
			lr = 0.1
		elif name.split("_")[-1] == "angle":
			lr = 1

		polinDe = polinomial.deriv()

		for i in list(np.roots(polinDe)):

			f = i-lr
			m = i
			l = i+lr

			if polinomial(f) < polinomial(m) and polinomial(l) < polinomial(m) and minx <= m <= maxx:
				lstM.append(i)
			else: pass

	maxima = []
	lookForMaxima(maxima, polin, slow["minx"], slow["maxx"])

	maxima = compareMaximas(maximas_mc, maxima)

	atomdist = A - B
	print "\tDISTANCE VALUE", atomdist

	
	dist = "None"
	if len(maxima) == 1:
		dist = maxima[0]
		print "\tONE MAXIMUM", dist
		return dist

	elif len(maxima) > 1:

		maxima.sort()

		if maxima[0] <= atomdist <= maxima[-1]:
			gap = slow["delta"]
			lv = atomdist - gap
			mv = atomdist
			rv = atomdist + gap

			if polin(lv) < polin(mv) < polin(rv):		# higher function values up on right
				while polin(mv) < polin(rv):
					mv += gap
					rv += gap
					if polin(mv) >= polin(rv):	# if polin(mv) < polin(rv):
						dist = mv
						print dist, maxima, atomdist
						break 
				print "\thigher values on the right"
				return checkDifference(dist, maxima)

			elif polin(lv) > polin(mv) > polin(rv):	# higher function values up on left
				while polin(lv) > polin(mv):
					lv -= gap
					mv -= gap
					if polin(lv) <= polin(mv):	# if polin(mv) > polin(lv):
						dist = mv
						print dist, maxima, atomdist
						break
				print "\thigher values on the left"
				return checkDifference(dist, maxima)

			elif polin(lv) > polin(mv) < polin(rv):	# higher function values up on both sides
				mx = [i for i in maxima]
				mx.append(mv)
				mx.sort()
				mvi = mx.index(mv)
				ls = polin(mx[mvi-1])
				rs = polin(mx[mvi+1])
				if ls > rs:
					dist = mx[mvi-1]
				elif ls < rs:
					dist = mx[mvi+1]
			
				print "\tvalue in the hole"
				return checkDifference(dist, maxima)

			elif polin(lv) < polin(mv) > polin(rv):	# actual atom distance is OK that's mean atom lies at the maximum
				dist = mv
				print "\tvalue at the maximim"
				return checkDifference(dist, maxima)

		elif atomdist < maxima[0]:
			dist = maxima[0]
			print "\tvalue on the left", dist
			return dist

		elif atomdist > maxima[-1]:
			dist = maxima[-1]
			print "\tvalue on the right", dist
			return dist
		else: pass
	print "\t\tfinal distance", dist



# refineDistanceAndSph() moves atoms in order to statistical data; input - atoms, statistical distance
def refineDistanceAndSph(atoms, dist):		# function gets list of two atoms to refine and proper distance to be set between these atoms

	# function checkes if x, y ,z point lies on a sphere 
	def sphereValues(x, y, z, mx, my, mz, r, xyz):
		if round((x - mx)**2 + (y - my)**2 + (z - mz)**2, 2) == round(r**2, 2):
			xyz.append((x,y,z))
		elif (x - mx)**2 + (y - my)**2 + (z - mz)**2 != r**2:
			pass


	A, B = [atoms[0], atoms[1]]

	xm, ym, zm = [A.coord[0] + (B.coord[0] - A.coord[0])*0.5, A.coord[1] + (B.coord[1] - A.coord[1])*0.5, A.coord[2] + (B.coord[2] - A.coord[2])*0.5]		# middle point coordinants between atoms A and B
	middle = np.array((xm, ym, zm), dtype=float)
	D = np.linalg.norm(middle - A.coord)			# absolute distance between atom A and middle point of atoms A and B

	proportion = (D - dist * 0.5) / D

	xA, yA, zA = [(xm - A.coord[0])*proportion, (ym - A.coord[1])*proportion, (zm - A.coord[2])*proportion]
	xB, yB, zB = [(B.coord[0] - xm)*proportion, (B.coord[1] - ym)*proportion, (B.coord[2] - zm)*proportion]

	newCoordA = [A.coord[0] + xA, A.coord[1] + yA, A.coord[2] + zA]
	newCoordB = [B.coord[0] - xB, B.coord[1] - yB, B.coord[2] - zB]
	
	A.coord = newCoordA
	B.coord = newCoordB

	return {"coord": A.coord, "sphMiddle": middle}



# createSphPoints() creates surface points of sphere; input - stastical distance, value of proportion, coordinates, middle coordinates 
def createSphPoints(distance, proportion, coord, middle):

	# function checkes if x, y ,z point lies on a sphere 
	def sphereValues(x, y, z, mx, my, mz, r, xyz):
		if round((x - mx)**2 + (y - my)**2 + (z - mz)**2, 2) == round(r**2, 2):
			xyz.append((x,y,z))
		elif (x - mx)**2 + (y - my)**2 + (z - mz)**2 != r**2:
			pass


	xyz = []
	step = 0.02		# step for range values

	xvalues = np.arange(coord[0] - distance * proportion, coord[0] + distance * proportion, step)
	yvalues = np.arange(coord[1] - distance * proportion, coord[1] + distance * proportion, step)
	zvalues = np.arange(coord[2] - distance * proportion, coord[2] + distance * proportion, step)


	for x in xvalues:
		for y in yvalues:
			for z in zvalues:
				sphereValues(x, y, z, middle[0], middle[1], middle[2], distance * 0.5, xyz)

	return {"sphPoints": xyz}



# getAngleValue() calculate actual contact angle; input - atoms
def getAngleValue(atoms):

	A, D = [atoms[0], atoms[1]]

	vacc = A.get_vector()
	vdon = D.get_vector()

	if D.parent.resname.split()[0] == "G":
		vddon = D.parent["C2"].get_vector()
		atdon = D.parent["C2"]
	elif D.parent.resname.split()[0] == "C" or D.parent.resname.split()[0] == "U":
		vddon = D.parent["C4"].get_vector()
		atdon = D.parent["C4"]
	elif D.parent.resname.split()[0] == "A":
		vddon = D.parent["C6"].get_vector()
		atdon = D.parent["C6"]
	else: pass

	angleR = calc_angle(vacc, vdon, vddon)
	angle = angleR * 180 / np.pi
	return angle



# refineAngleValue() returns statistical value of contact in order to actual contact angle; input - atoms, contact name, dictionary with polinomial
def refineAngleValue(atoms, name, slow):		# function takes list of two atoms to refine and gives proper angle value to be set between these atoms

	polin = np.poly1d(slow["polinomial"])
	maximas_err = slow["maximas"]	# it is list with maksimas and its errors from Monte Carlo simulation
	maximas_mc = [i[0] for i in maximas_err]

	def checkDifference(value, lst):
		final_angle = 0
		diff = 1000
		for m in lst:
			if abs(m - value) < diff:
				diff = abs(m - value)
				final_angle = m
		return final_angle


	def lookForMaxima(lstM, polinomial, minx, maxx):

		lr = 0
		if name.split("_")[-1] == "distance":
			lr = 0.1
		elif name.split("_")[-1] == "angle":
			lr = 1

		polinDe = polinomial.deriv()

		for i in list(np.roots(polinDe)):

			f = i-lr
			m = i
			l = i+lr

			if polinomial(f) < polinomial(m) and polinomial(l) < polinomial(m) and minx <= m <= maxx:
				lstM.append(i) # 				maxima.append(i)
			else: pass

	maxima = []
	lookForMaxima(maxima, polin, slow["minx"], slow["maxx"])

	maxima = compareMaximas(maximas_mc, maxima)

	angle = getAngleValue(atoms) 
	print "\tANGLE VALUE", angle

	ang = "nope"
	if len(maxima) == 1:
		ang = maxima[0]
		print "\tONE MAXIMUM", ang
		return ang

	elif len(maxima) > 1:
			
		maxima.sort()

		if maxima[0] <= angle <= maxima[-1]:
			gap = slow["delta"]
			lv = angle - gap
			mv = angle
			rv = angle + gap

			if polin(lv) < polin(mv) < polin(rv):		# higher function values up on right
				while polin(mv) < polin(rv):
					mv += gap
					rv += gap
					if polin(mv) >= polin(rv):
						ang = mv
						print ang, maxima
						break 
				print "\thigher values on the right"
				return checkDifference(ang, maxima)

			elif polin(lv) > polin(mv) > polin(rv):	# higher function values up on left
				while polin(lv) > polin(mv):
					lv -= gap
					mv -= gap
					if polin(mv) >= polin(lv):
						ang = mv
						print ang, maxima
						break
				print "\thigher values on the left"
				return checkDifference(ang, maxima)

			elif polin(lv) > polin(mv) < polin(rv):	# higher function values up on both sides
				mx = [i for i in maxima]
				mx.append(mv)
				mx.sort()
				mvi = mx.index(mv)
				ls = polin(mx[mvi-1])
				rs = polin(mx[mvi+1])
				if ls > rs:
					ang = mx[mvi-1]
				elif ls < rs:
					ang = mx[mvi+1]
			
				print "\tvalue in hole"
				return checkDifference(ang, maxima)

			elif polin(lv) < polin(mv) > polin(rv):	# actual atom distance is OK that's mean atom lies at the maximum
				dist = mv
				print "\tvalue at maximim"
				return checkDifference(ang, maxima)

		elif angle < maxima[0]:
			ang = maxima[0]
			print "\tvalue on the left"
			return ang

		elif angle > maxima[-1]:
			ang = maxima[-1]
			print "\tvalue on the right"
			return ang
		else: pass
	print "\t\tfinal distance", ang



# refineAngle() refines contact angle in order to statistical data; input - contact, spherical points, middle coordinates, statistical value of angle
def refineAngle(bond, sph, middle, refAngle):

	A = bond["acc"]
	D = bond["don"]

	### plotting part - start
	X=[]
	Y=[]
	Z=[]

	for i in sph:
		X.append(float(i[0]))
		Y.append(float(i[1]))
		Z.append(float(i[2]))

	fig = pl.figure(figsize=(8,8))
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(X, Y, Z, c='w')
	ax.scatter(middle[0], middle[1], middle[2], marker = ".", s=80, c="k")

	ax.scatter(A.coord[0], A.coord[1], A.coord[2], c='r', marker="o", s=80)
	ax.scatter(D.coord[0], D.coord[1], D.coord[2], c='b', marker="o", s=80)
	pl.title("input")
##	pl.show()
	pl.close()
	### plotting part - end

	nowAngle = getAngleValue([A, D])

	order = []
	for p in sph:
		point = np.array(p, dtype=float)
		tup = (np.linalg.norm(A.coord - point), sph.index(p))
#D		tup = (np.linalg.norm(D.coord - point), sph.index(p))
		order.append(tup)

	order.sort()
	for d in order:
		fig = pl.figure()
		ax = fig.add_subplot(111, projection='3d')

		ax.scatter(X, Y, Z, c='w')
		ax.scatter(middle[0], middle[1], middle[2], marker = ".")


		A.coord = np.array(sph[d[1]], dtype=float)
#D 	D.coord = np.array(sph[d[1]], dtype=float)

		xD, yD, zD = [A.coord[0] - middle[0], A.coord[1] - middle[1], A.coord[2] - middle[2]]
#D		xA, yA, zA = [D.coord[0] - middle[0], D.coord[1] - middle[1], D.coord[2] - middle[2]]

		D.coord = [A.coord[0] - xD * 2., A.coord[1] - yD * 2., A.coord[2] - zD * 2.]
#D		A.coord = [D.coord[0] + xA * 2., D.coord[1] + yA * 2., D.coord[2] + zA * 2.]

		ax.scatter(A.coord[0], A.coord[1], A.coord[2], c='r')
		ax.scatter(D.coord[0], D.coord[1], D.coord[2], c='b')
#		pl.show()
		pl.close()
#		print "now angle - ",nowAngle, "ref Angle - ", refAngle, "actual Angle -", getAngleValue([A, D])
		if round(refAngle, 0) == round(getAngleValue([A, D]), 0):
			print "now angle - ",nowAngle, "ref Angle - ", refAngle, "actual Angle -", getAngleValue([A, D])
			break


	fig = pl.figure(figsize=(8,8))
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(X, Y, Z, c='w')
	ax.scatter(middle[0], middle[1], middle[2], marker = ".", s=80, c="k")

	pl.title("output")
	ax.scatter(A.coord[0], A.coord[1], A.coord[2], c='r', marker="o", s=80)
	ax.scatter(D.coord[0], D.coord[1], D.coord[2], c='b', marker="o", s=80)
##	pl.show()
	pl.close()



# moveMainTool() main function that calls functions in proper order; input - name of the structure to refine, name of refined structure, used module with polinomials 
def moveMainTool(strName, strOut, methode):

	structure = PDBParser().get_structure(strName.split("/")[-1], strName)
	allRes = Selection.unfold_entities(structure, "R")
	allAtoms = Selection.unfold_entities(structure, "A")
	NS = NeighborSearch(allAtoms)

	contactsToRefine = findContacts(NS, allRes, 4.5)
		
	cR = playClaRNA(contactsToRefine, strName, 0.1)
	
	for i in contactsToRefine: print i, i["acc"].parent.id, i["don"].parent.id
	print

	for c in cR:
		print "\n\nNEW CONTACT MANIPULATION\n"
		print c
		print [c["acc"], c["don"]], c["acc"].parent.resname, c["don"].parent.resname, c["acc"].parent.id[1], c["don"].parent.id[1]
		print createName(c, "distance")

		try:
			m = callFormula(createName(c, "distance"), methode)
			properDistance = refineDistanceValue([c["acc"], c["don"]], createName(c, "distance"), m)
			coordAndMiddle = refineDistanceAndSph([c["acc"], c["don"]], properDistance)
			print "distance DONE"
		except KeyError:
			pass

		try:
			m = callFormula(createName(c, "angle"), methode)
			properAngle = refineAngleValue([c["acc"], c["don"]], createName(c, "angle"), m)
			actualAngle = getAngleValue([c["acc"], c["don"]])

			prop = abs(properAngle - actualAngle) / 100
	
			sph =  createSphPoints(properDistance, prop, coordAndMiddle["coord"], coordAndMiddle["sphMiddle"])
			refineAngle(c, sph["sphPoints"], coordAndMiddle["sphMiddle"], properAngle)
			print "angle DONE"

		except KeyError:
			pass

		print "FINISHED"
	

	io = PDBIO()
	io.set_structure(structure)
	io.save(strOut)


	


if __name__ == "__main__":
	
	err = """\n\t ERROR ERROR ERROR\n\nstructureRefinement.py requires path to the structure to refine and declaration of polinomials to use\n\n\t$./structureRefinement.py <path to RNA structure> -option\n\n\toptions: -value or -poly\n"""

	try:

		start = time.time()
		at = moveMainTool(sys.argv[1], sys.argv[1].split("/")[-1].split(".")[0]+"_ref_"+sys.argv[2][1:]+".pdb", sys.argv[2])
		print "\t\tTAKEN TIME WAS...", round(time.time() - start, 3), "s"
		print "\t\tTAKEN TIME WAS...", str(round((time.time() - start)/60, 2)).split(".")[0]+" min "+str(float(str(round((time.time() - start)/60, 2)).split(".")[1])*60)[0:2]+" s"

	except IndexError:
		print err





