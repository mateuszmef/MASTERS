#!/usr/bin/python
#! -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import scipy.optimize
import pylab as pl
import anglesDistancesList as adlp
from scipy.stats import norm
import time



# getContacts() creates list with names of lists stored in "angleDistanceList" module; input - file name of module
def getContacts(anglesDistances):		

	contacts = []
	for line in open(anglesDistances, "r"):
		if not line.split():
			pass
		else:
			contacts.append(line.split()[0])

	return contacts



# returns polynomial equation to selected points of maximum; input - contact name, value of distance delta, value of angle delta, float value 0 < value < 1 of percent of highest point,
#		condition if value method or polynomial method
def countRanges(name, deltaDist, deltaAngle, percent, condition):

	l = adlp.__dict__[name]
	l.sort()

	minv = min(l)
	maxv = max(l)

	if name.split("_")[-1] == "distance":
		delta = deltaDist
	elif name.split("_")[-1] == "angle":
		delta = deltaAngle
	else: pass	


	# "bins" creation
	x = []
	y = []

	binNumbers = int((maxv - minv) / delta) + 1

	types = []
	for n, i in enumerate(range(binNumbers)):
		dd = {} 
		x.append(minv+delta*i+delta*0.5)
		dd["range"] = str(minv+delta*i)+" - "+ str(minv+(i+1)*delta)
		dd["counted"] = 0
		for p in l:
			if p >= minv+delta*i and p < minv+(i+1)*delta:
				dd["counted"] += 1
 			else: pass

		y.append(dd["counted"])

		types.append(dd)

	def runValueMethod():
		# value method - restriction to only this points which are larger than value of some percent of maksimum value of y values 

		XY = []
		for i in range(len(x)):
			XY.append((x[i], y[i]))

		iksigrek = []
		for t in XY:
			if t[1] >= max(y) * percent:
				iksigrek.append(t)
			else: pass

		XX = []
		YY = []
		for i in iksigrek:
			XX.append(i[0])
			YY.append(i[1])

		if len(XX) <= 6:
			degree = 3
		elif len(XX) > 6:
			degree = 4

		fitSel = np.poly1d(np.polyfit(XX, YY, degree))

#		pl.plot(XX, YY, ".-", XX, fitSel(XX), "--")
#		pl.title("value")
#		pl.show()

		return {"polinomial": fitSel, "delta": delta, "maxx": max(XX), "minx": min(XX)}	# <============================


	def runPolynomialMethod():
		# function method - restricion to that x arguments which fits to fited polinomial

		dr = float(len(l)) / float(len(x))

		ddd = 0
		if len(l) <= 60:
			ddd = 4
		else:
			ddd = int(round(7./205 * dr + 8, 1))

		pddd = np.poly1d(np.polyfit(x, y, ddd))


		maks = 0
		for ix in x:			# loop checks for maximum value of x arguments
			if pddd(ix) > maks:
				maks = pddd(ix)
			else: pass

		def findPercentP(polinomial, ikses, ma):
			for i in ikses:
				if polinomial(i) >= (ma * percent):
					return ikses.index(i)
				else:
					pass
	
		ls = findPercentP(pddd, x, maks)
		x.reverse()
		rs = findPercentP(pddd, x, maks)
		x.reverse()
		
		if rs == 0:
			rs = 1

		xx = x[ls:(-rs)]
		yy = y[ls:(-rs)]

		if len(xx) <= 6:
			degree = 3
		elif len(xx) > 6:
			degree = 4

		wow = np.poly1d(np.polyfit(xx, yy, degree))

#		pl.plot(xx, yy, ".-", xx, wow(xx), "--")
#		pl.title("poly")
#		pl.show()

		return {"polinomial": wow, "delta": delta, "maxx": x[(-rs)], "minx": x[ls]}


	if condition == "-poly":
		return runPolynomialMethod()
	elif condition == "-value":
		return runValueMethod()
	elif condition == "-test":
		runPolynomialMethod()
		runValueMethod()
	else:
		raise ValueError



# countLenght() returns only that contacts which its list contain at least some NUMBER OF VALUES; input - contact name, NUMBER OF VALUES
def countLenght(cont, value):
	
	nope = 0
	l = []
	for i in cont: 
		if len(adlp.__dict__[i]) >= value:
			if i.split("_")[4] == "nope":
				nope += 1
			else: pass			
			l.append(i)

		else: pass

	return l



# lookForComplementaries() returns that contacts which are complementary; input - contact name
def lookForComplementaries(cont):


	def makeComp(i):
		if i[2] in donA and i[3] == "A" or i[2] in donU and i[3] == "U" or i[2] in donG and i[3] == "G" or i[2] in donC and i[3] == "C":
			return True
		else: pass


	accA = ["N7", "N1", "N3"]
	donA = ["N6"]

	accU = ["O4", "O2"]
	donU = ["N3"]

	accG = ["N7", "O6", "N3"]
	donG = ["N1", "N2"]

	accC = ["N3", "O2"]
	donC = ["N4"]

	res = ["A", "C", "G", "U"]

	outputContacts = []

	for i in cont:

		i = i.split("_")

		if i[0] == 'A' and i[1] in accA:
			if makeComp(i) == True:
				outputContacts.append("_".join(i))

		elif i[0] == "U" and i[1] in accU:
			if makeComp(i) == True:
				outputContacts.append("_".join(i))

		elif i[0] == "C" and i[1] in accC:
			if makeComp(i) == True:
				outputContacts.append("_".join(i))

		elif i[0] == "G" and i[1] in accG:
			if makeComp(i) == True:
				outputContacts.append("_".join(i))

		else:
			pass

	return outputContacts



# giveOnlyClassified() returns only classified contacts; input - contact name
def giveOnlyClassyfied(con):
	
	a = []
	for i in con:
		if i.split("_")[4] != "nope":
			a.append(i)
		else:
			pass

	return a



# exportPolinomial() saves key information about created polinomials; input - file with saved polinomials, contact name, dictionary with keys of polinomial
def exportPolinomial(filee, name, dictt):

	filee.writelines(name+" = {'polinomial': [")
	for n, i in enumerate(list(dictt["polinomial"])):
		if n == len(list(dictt["polinomial"])) -1:
			filee.writelines(str(i))
		else:
			filee.writelines(str(i)+", ")
	filee.writelines("], 'delta': "+str(dictt["delta"])+", 'maxx': "+str(dictt["maxx"])+", 'minx': "+str(dictt["minx"]))

	filee.writelines("}\n")
	filee.flush()



# moveMainTool() main function that call proper function in proper order; input - terminal arguments
def moveMainTool(system):

	# default values
	knife = 20
	percent = 0.2
	con = "-poly"
	dd = 0.025
	da = 1

	allowed = ["-k", "-p", "-poly", "-value", "-dd", "-da", "-h"]
	
	if "-k" in system:
		knife = int(system[system.index("-k") + 1])

	if "-p" in sys.argv:
		percent = float(system[system.index("-p") + 1])

	if "-poly" in system:
		con = "-poly"
	
	if "-value" in system:
		con = "-value"

	if "-dd" in system:
		dd = float(system[system.index("-dd") + 1])

	if "-da" in system:
		da = float(system[system.index("-da") + 1])

	if "-h" in system and system[1] == "-h":
			h =  """\nTo run this program user is able to use default values of 'knife', 'percent', delta angle, 
delta distance and used method which creates polinomials or user can type it in the command line. 

examples of COMMANDS:

	>>> ./makePoly.py -poly -k 25 -p 0.3 -dd 0.03 -da 2
	>>> ./makePoly.py -poly -k 25 -da 2 -p 0.3 -dd 0.03
	>>> ./makePoly.py -da 2 -dd 0.03 -value					
	>>> ./makePoly.py\t# dafault
	
There is no strict order in typed options. Only restriction is to type value right next to typed option for example "-p 0.3". 0.3 value refers to "-p" option. If option not in command program uses its default values

	# default values
	option\t\tvalue
	'-k'\tknife = 20\t\t- at least values on list with distances or angles
	'-p'\tpercent = 0.2\t\t- 0 < percent < 1 - value by which selection of point of maximas are able to reach
	'-poly' OR '-value'\tcon = "-poly"\t\t- '-poly' - polymials method, '-value' - values method
	'-dd'\tdd = 0.025\t\t- range of bin, delta for distance values
	'-da'\tda = 1\t\t- range of bin, delta for angle values

			"""
			print h
			return

	print
	print "-k\t", knife
	print "-p\t", percent
	print "-dd\t", dd
	print "-da\t", da
	print

	contacts = getContacts(adlp.__name__+".py")										# list with all contact names stored in adlp module 

	cuted = countLenght(contacts, knife)												# list with contact names that contained at least 'knife' values 

	complementaries = lookForComplementaries(cuted)									# list with contact names that are complementary

	noNopes = giveOnlyClassyfied(complementaries)									# list with contact names that are classified by ClaRNA

#	d = countRanges("A_N1_N3_U_WW_cis_angle", 0.025, 1, percen, con)
#	d = countRanges("A_N1_N3_U_WW_cis_distance", 0.025, 1, percen, con)

	if con == "-poly":
		f = open("TEST_formulas_polinomial.py", "w")
	elif con == "-value":
		f = open("TEST_formulas_absolute_values.py", "w")
	else:
		raise ValueError

	f.writelines("#\n#\tmaximum y / polinomial value cutoff = "+str(int(percent * 100))+" % \n#\tknife = "+str(knife)+"\n#\n")

	print "EXPORTING DATA TO "+f.name+" FILE..."
	for c in noNopes:
		d = countRanges(c, 0.025, 1, percent, con)
		exportPolinomial(f, c, d)
	print "\nPROCESS IS DONE\n"
	



if __name__ == "__main__":

	moveMainTool(sys.argv)


