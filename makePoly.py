#!/usr/bin/python
#! -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import scipy.optimize
import pylab as pl
import anglesDistancesList as adlp
from scipy.stats import norm, poisson
import time
import orcasValues as ov
from sklearn.neighbors import KernelDensity



# getContacts() creates list with names of lists stored in "angleDistanceList" module; input - file name of module
def getContacts(anglesDistances):		

	contacts = []
	for line in open(anglesDistances, "r"):
		if not line.split():
			pass
		else:
			contacts.append(line.split()[0])

	return contacts


def calculateMaxima(name_contact, function, min_val, max_val):

	# searching maximas of polinomial
#	print "MIEJSCA ZEROWE WIELOMIANU\t" #, np.roots(function.deriv())

	maximas = list()
	diff = 0.1 if name_contact.split("_")[-1] == "distance" else 1
	for i in np.roots(function.deriv()):

		f, m, l = i-diff, i, i+diff
		if function(f) < function(m) and function(m) > function(l) and min_val <= m <= max_val:
			maximas.append(i)
#	print "EKSTERM - MAKSIMA FUNKCJI\t", maximas

#	print "przed", maximas

	if len(maximas) == 1:
		if "+" in str(maximas[0]):
#			print maximas
			maximas = [str(maximas[0]).split("+")[0][1:]]
#			print maximas
		else:
			pass
	if len(maximas) == 2:
		if "+" in str(maximas[0]):
#				print maximas
				maximas = [str(maximas[0]).split("+")[0][1:]]
#				print maximas
				print
		else:
			pass
			
#	print "po", maximas
#	print
	return maximas



def monteCarloSimulation(values):
	
	simulated_values = list()
	for v in values:
		if v < 50:
			simulated_values.append( poisson.rvs(mu=v, loc=int(round(np.sqrt(v))))) 
		elif v >= 50:
			simulated_values.append( round(np.random.normal(loc=v, scale=np.sqrt(v))) )

	return simulated_values


def calculateRMS(list_with_values, option = "only_err"):
	
	try:
		if len(list_with_values) == 1:
			avg = list_with_values[0]
			rms = 0.0
		elif len(list_with_values) > 1:
			avg = np.average(list_with_values)	
			rms_list = [(i - avg)**2 for i in list_with_values]
			rms = np.sqrt(np.sum(rms_list) / len(rms_list))
			print avg, "+-", rms

		if option == "only_err":
			return rms
		elif option == "avg_and_err":
			return (avg, rms)

	except:
		return ValueError


def calculateError(in_maximas):

	print "\nCALCULATE ERROR - start"

	errors = list()

	maximas = [mx for mx in in_maximas if len(mx) >=1]
	
	option = 0
	numbers = []
	for vals in maximas:
		if len(vals) == 1:
			numbers.append(1)
		elif len(vals) == 2:
			numbers.append(2)
		elif len(vals) > 2: 
			raise ValueError


	if float(np.sum(numbers)) / len(numbers) == 1:
		print "\nOPTION == 1\n"
		values = [float(val[0]) for val in maximas]
			
		average = np.average(values)
		err = calculateRMS(values, option = "avg_and_err")
		errors.append(err)
	
	elif float(np.sum(numbers)) / len(numbers) == 2:
		print "\nOPTION == 2\n"
		first = [float(val[0]) for val in maximas]
		second = [float(val[1]) for val in maximas]

		average = np.average(first)
		err = calculateRMS(first, option = "avg_and_err")
		errors.append(err)

		average = np.average(second)
		err = calculateRMS(second, option = "avg_and_err")
		errors.append(err)
	
	elif 1.0 < float(np.sum(numbers)) / len(numbers) < 2.0:

		print "\n 1 < OPTION < 2\n"
		first = list()
		second = list()

		print "DLUGOSC MAXIMAS"
		print len(maximas)
		new_maximas = list()

		for vals in maximas:
			if len(vals) == 1:
				new_maximas.append(float(vals[0]))

			elif len(vals) == 2:
				first.append(float(vals[0]))
				second.append(float(vals[1]))

			else:
				raise ValueError

		for val in new_maximas:
			if abs(val - float(first[0])) < abs(val - float(second[0])):
				first.append(val)
			elif abs(val - float(first[0])) > abs(val - float(second[0])):
				second.append(val)
			else:
				raise ValueError
		
		average = np.average(first)
		err = calculateRMS(first, option = "avg_and_err")
		errors.append(err)		

		average = np.average(second)
		err = calculateRMS(second, option = "avg_and_err")
		errors.append(err)

	else:
		raise ValueError

	print errors
	return errors

			
		

	print "\nCALCULATE ERROR - stop\n"



# returns polynomial equation to selected points of maximum; input - contact name, value of distance delta, value of angle delta, float value 0 < value < 1 of percent of highest point,
#		condition if value method or polynomial method
def countRanges(name, deltaDist, deltaAngle, percent, condition):

	l = adlp.__dict__[name]
	l.sort()

	minv = min(l)
	maxv = max(l)

	dimension = float(len(l))

	if name.split("_")[-1] == "distance":
		delta = deltaDist
		suffix = "distance [Angstroms]"
	elif name.split("_")[-1] == "angle":
		delta = deltaAngle
		suffix = "angle [degree]"
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

		y.append(dd["counted"]) #/dimension)

		types.append(dd)

	print dimension
	print y
	print "!!!!!!!!!!!!!!"



	pl.plot(x, y, '.-')
	pl.title(" ".join(name.split('_')))
	pl.xlabel(suffix)
	pl.ylabel("number of values")
	pl.ylim(ymin=0)
	
	try:
		pl.axvline(ov.__dict__[name], color = 'r')
	except:
		pass

#	pl.show()
	pl.close()
		

	def runValueMethod():
		# value method - restriction to only this points which are larger than value of some percent of maksimum value of y values 

#		pl.plot(x, y, "-")

		print "\n"+"VALUE METHOD"*4+"\n"

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
		pl.plot(XX, YY, ".-", color = 'r')
		pl.plot(XX, fitSel(XX), "--", color = 'g')

		mc_maximas = list()

		# Monte Carlo part
		for i in range(100):
			new_yy = monteCarloSimulation(YY)
			new_poly = np.poly1d(np.polyfit(XX, new_yy, degree))
			mc_maximas.append(calculateMaxima(name, new_poly, min(XX) , max(XX)))

			pl.plot(XX, new_poly(XX), '-')

		mc_maximas.append(calculateMaxima(name, fitSel, min(XX) , max(XX)))
		to_write_errors = calculateError(mc_maximas)

		try:
			for error in to_write_errors:
				pl.axvline(float(error[0]), color = 'g')
#				pl.axvline(float(error[0])+float(error[1]), color = 'g')
#				pl.axvline(float(error[0])-float(error[1]), color = 'g')
			pass
		except ValueError:
			pass


		try:
			pass
#			pl.axvline(ov.__dict__[name], color = 'r')
		except:
			pass

		pl.title(" ".join(name.split('_')))
		pl.xlabel(suffix)
		pl.ylabel("number of values")
#		pl.show()
		pl.close()
		# PLOTTING PART >>>>>>>>>>> HISTOGRAM
		
		print {"polinomial": fitSel, "delta": delta, "maxx": max(XX), "minx": min(XX), "maximas": to_write_errors}	
		return {"polinomial": fitSel, "delta": delta, "maxx": max(XX), "minx": min(XX), "maximas": to_write_errors}	# <============================


	def runPolynomialMethod():
		# function method - restricion to that x arguments which fits to fited polinomial

		print "\n"+"POLY METHOD"*4+"\n"

#		pl.plot(x, y, "-")

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
		pl.plot(xx, yy, ".-", color = 'r')
		pl.plot(xx, wow(xx), "--", color = 'g')

		mc_maximas = list()
		# Monte Carlo part
		for i in range(100):
			new_yy = monteCarloSimulation(yy)
			new_poly = np.poly1d(np.polyfit(xx, new_yy, degree))
			mc_maximas.append(calculateMaxima(name, new_poly, min(xx) , max(xx)))

			pl.plot(xx, new_poly(xx), '-')

		mc_maximas.append(calculateMaxima(name, wow, min(xx) , max(xx)))
		to_write_errors = calculateError(mc_maximas)

		try:
			for error in to_write_errors:
				pl.axvline(float(error[0]), color = 'g')
#				pl.axvline(float(error[0])+float(error[1]), color = 'g')
#				pl.axvline(float(error[0])-float(error[1]), color = 'g')
			pass
		except ValueError:
			pass


		try:
			pass
#			pl.axvline(ov.__dict__[name], color = 'r')
		except:
			pass


		pl.title(" ".join(name.split('_')))
		pl.xlabel(suffix)
		pl.ylabel("number of values")
#		pl.show()
		pl.close()

		print {"polinomial": wow, "delta": delta, "maxx": x[(-rs)], "minx": x[ls], "maximas": to_write_errors}
		return {"polinomial": wow, "delta": delta, "maxx": x[(-rs)], "minx": x[ls], "maximas": to_write_errors}


	if condition == "-poly":
		return runPolynomialMethod()
	elif condition == "-value":
		return runValueMethod()
	elif condition == "-test":
		runValueMethod()
		runPolynomialMethod()
		return (x, y)
	else:
		raise ValueError



def kdeFit(contact_name): #, mc, les):

	np.random.seed(1)

	XXX = adlp.__dict__[contact_name]
	X = adlp.__dict__[contact_name]

	minv = min(X) # -3 if contact_name.split("_")[-1] == "angle" else min(X) - 0.075
	maxv = max(X) # -3 if contact_name.split("_")[-1] == "angle" else max(X) - 0.075

	bw = 0.3 if contact_name.split("_")[-1] == "distance" else 5
	arg = "distance [Angstroms]" if contact_name.split("_")[-1] == "distance" else "angle [degree]"

	X = np.array(X)[:, np.newaxis]

	X_plot = np.linspace(minv, maxv, 1000)[:, np.newaxis]


	fig, ax = pl.subplots()
	kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(X)
	
	log_dens = kde.score_samples(X_plot)

	X_list = [float(x[0]) for x in X_plot]

	poly = np.poly1d(np.polyfit(X_list, np.exp(log_dens), 8))
	maxi = calculateMaxima(contact_name, poly, minv, maxv)
	for m in maxi:
		#pl.axvline(float(m), color = 'r')
		pass

	pl.hist(XXX, 30, normed = True)	
	pl.plot(X_plot[:, 0], np.exp(log_dens), '-')

	if "distance" in arg:
		pl.axvline(float(maxi[-1]), color = 'r')
		print float(maxi[-1])

	if "angle" in arg:
		pl.axvline(float(maxi[1]), color = 'r')
		print float(maxi[1])

	pl.title(" ".join(contact_name.split('_')))
	pl.xlabel(arg)
	pl.ylabel("density of probability")
	pl.show()
	pl.close()


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

	def makeString(mxs):
		res = "["
		for tup in mxs:
			n = "("+str(tup[0])+", "+str(tup[1])+")"
			res = res + n
			if len(mxs) == 1:
				res = res + "]"			
			elif tup == mxs[-1]:
				res = res + "]"
			else:
				res = res + ", "
		return res


	filee.writelines(name+" = {'polinomial': [")
	for n, i in enumerate(list(dictt["polinomial"])):
		if n == len(list(dictt["polinomial"])) -1:
			filee.writelines(str(i))
		else:
			filee.writelines(str(i)+", ")
	filee.writelines("], 'delta': "+str(dictt["delta"])+", 'maxx': "+str(dictt["maxx"])+", 'minx': "+str(dictt["minx"])+", 'maximas': "+makeString(dictt["maximas"]))

	filee.writelines("}\n")
	filee.flush()



# moveMainTool() main function that call proper function in proper order; input - terminal arguments
def moveMainTool(system):

	# default values
	knife = 20
	percent = 0.2
#	con = "-poly"
#	con = "-value"
	con = "-test"
	dd = 0.02
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

	if "-test" in system:
		con = "-test"

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
	print len(contacts)
	print len(contacts)/2
	print

	cuted = countLenght(contacts, knife)												# list with contact names that contained at least 'knife' values 
	print len(cuted)
	print len(cuted)/2
	print

	complementaries = lookForComplementaries(cuted)									# list with contact names that are complementary
	print len(complementaries)
	print len(complementaries)/2
	print

	noNopes = giveOnlyClassyfied(complementaries)									# list with contact names that are classified by ClaRNA
	print len(noNopes)
	print len(noNopes)/2


	if con == "-poly":
		f = open("formulas_polynomial_mc.py", "w")
	elif con == "-value":
		f = open("formulas_values_mc.py", "w")
	elif con == "-test":
		f = open("TEST_formulas_values_mc.py", "w")
	else:
		raise ValueError

	f.writelines("#\n#\tmaximum y / polinomial value cutoff = "+str(int(percent * 100))+" % \n#\tknife = "+str(knife)+"\n#\n")

	print "EXPORTING DATA TO "+f.name+" FILE..."
	for c in noNopes:
		d = countRanges(c, dd, da, percent, con)
		exportPolinomial(f, c, d)

#		kdeFit(c)
#	print "\nPROCESS IS DONE\n"	



if __name__ == "__main__":

	moveMainTool(sys.argv)


