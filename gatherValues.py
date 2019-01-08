#!/usr/bin/python
#! -*- coding: utf-8 -*-

from Bio.PDB import *
import numpy as np
import os
import time


# getStrNames() creates list of structures in path "catalogPath"; input - path of catalog with structures
def getStrNames(catalogPath):			# metoda pobierajaca sciezke do katalogu w ktorym znajduja sie struktury RNA 
	
	print "function - getStrNames\n"

	return os.listdir(catalogPath)



# openStructure() creates BioPython structure, input - file name of structure
def openStructure(pdbFile):			
	
	print "function - openStructure\n"

	parser = PDBParser()
	structure = parser.get_structure("name_"+pdbFile, pdbFile)
	
	return structure



# parserPDB() parses structure into atoms and chains; input - structure
def parserPDB(structure):
	
	print "function - parserPDB\n"
	
	chainList = []
	atomList = []
	if len(structure) == 1:
		for chain in structure[0]:
			chainList.append(chain)
			for res in chain:
				for atom in res:
					atomList.append(atom)
		return (atomList, chainList)
		
	elif len(structure) > 1:

		raise KeyError
#			print "\n\tTHESE IS NO ONLY ONE MODEL OF THE STRUCTURE\n\tPROGRAM REQUIRES SINGLE MODEL STRUCTURES\n\tline of the program - 51\n"
	



# findHBond() finds "hydrogen bonds" between atoms of A, C, G, U, T nucleobases; input - structure, NeighbourSearch instance, radius of the sphere
def findHBond(structure, NSearch, distance):

	print "\nfunction - findHBond\n"
					
	def findingContactsNoHydrogens(atom, accAtoms, donAtoms, NS, dist, dictList):

		d = ["A", "C", "G", "U", "T"]

		if atom.get_name() in accAtoms + donAtoms and atom.parent.resname.split()[0] in d:
			close = NS.search(atom.coord, dist, level="A")
			for closeAtom in close:
				if atom.parent != closeAtom.parent and abs(atom.parent.id[1] - closeAtom.parent.id[1]) > 1 and \
					closeAtom.name in accAtoms + donAtoms and closeAtom.parent.resname.split()[0] in d:
					if atom.name in accAtoms and closeAtom.name in donAtoms and {"acc": atom, "don": closeAtom, "resacc": atom.parent.id[1], "resdon": closeAtom.parent.id[1]} not in dictList:
						dictList.append({"acc": atom, "don": closeAtom, "resacc": atom.parent.id[1], "resdon": closeAtom.parent.id[1]})

					elif atom.name in donAtoms and closeAtom.name in accAtoms and {"acc": closeAtom, "don": atom, "resacc": closeAtom.parent.id[1], "resdon": atom.parent.id[1]} not in dictList:
						dictList.append({"acc": closeAtom, "don": atom, "resacc": closeAtom.parent.id[1], "resdon": atom.parent.id[1]})

					else:
						pass
				else:
					pass				
		else:
			pass


	allAcc = ["N7", "N1", "N3", "O4", "O2", "O6"]
	allDonNoHydrogens = ["N6", "N3", "N1", "N2", "N4"]

	contactList = []			# list of dictionaries of found contacts

	for res in structure:

		if res.resname.split()[0] == "A":
			for a in res:
				findingContactsNoHydrogens(a, allAcc, allDonNoHydrogens, NSearch, distance, contactList)

		elif res.resname.split()[0] == "U":
			for a in res:
				findingContactsNoHydrogens(a, allAcc, allDonNoHydrogens, NSearch, distance, contactList)
		
		elif res.resname.split()[0] == "G":
			for a in res:
				findingContactsNoHydrogens(a, allAcc, allDonNoHydrogens, NSearch, distance, contactList)

		elif res.resname.split()[0] == "C":
			for a in res:
				findingContactsNoHydrogens(a, allAcc, allDonNoHydrogens, NSearch, distance, contactList)

		else:
			pass

	return contactList



# countAngleAndDistance() counts angle and distance of inputed contact; input - single contact, file with results
def countAngleAndDistance(bond, resultFile):

	vacc = bond["acc"].get_vector()
	vdon = bond["don"].get_vector()
	if bond["don"].parent.resname.split()[0] == "G":
		vddon = bond["don"].parent["C2"].get_vector()
		atdon = bond["don"].parent["C2"]
	elif bond["don"].parent.resname.split()[0] == "C" or bond["don"].parent.resname.split()[0] == "U":
		vddon = bond["don"].parent["C4"].get_vector()
		atdon = bond["don"].parent["C4"]
	elif bond["don"].parent.resname.split()[0] == "A":
		vddon = bond["don"].parent["C6"].get_vector()
		atdon = bond["don"].parent["C6"]
	else: pass

	dist = abs(bond["acc"] - bond["don"])

	angleR = calc_angle(vacc, vddon, vdon)
	angleD = angleR * 180 / np.pi
	angle = str(angleD)

	resultFile.writelines([str(bond["acc"].parent.parent.get_full_id()[0].split("/")[-1])+'\t'+str(bond["acc"].parent.parent.id)+'\t'+str(bond["don"].parent.parent.id)+'\t'+str(bond["resacc"])+'\t'+str(bond["acc"].parent.resname.split()[0])+'\t'+str(bond["acc"].id)+'\t'+str(bond["don"].id)+'\t'+str(atdon.id)+'\t'+str(bond["don"].parent.resname.split()[0])+'\t'+str(bond["resdon"])+'\t'+str(angleD)+'\t'+str(dist)+"\t"+bond["class"]+'\n'])

	resultFile.flush()
	
	return (angleD, dist)



# makeClaRNA() runs ClaRNA tool which classiefies contacts of the structure; function works independently and parallel to this program; result is a step of processing; input - structurs name, 
#		path to catalog with structures, path to catalog with semi-finished products, thresh for ClaRNA tool, path to ClaRNA tool
def makeClaRNA(name, pathPDBS, pathWorkFile, thresh, pathToClarna):

	print "\nfunction - makeClaRNA\n"

	command = pathToClarna+" -ipdb "+pathPDBS+name+" > "+pathWorkFile+" -thresh "+ str(thresh)
	os.system(command)

	classes = []

	for line in open(pathWorkFile, "r").readlines()[2:-1]:

		line = line.split()
		classes.append({"chainRes": [line[0]+"_"+line[1]+"_"+line[5], line[2]+"_"+line[3]+"_"+line[6]], "class": line[7]})
	
	plik = open(pathWorkFile, "a")
	plik.writelines([name])
	plik.close()

	os.system("cp "+pathWorkFile+" "+pathWorkFile+"_COPY")

	print "\tNumber of found classes by ClaRNA -", len(classes), "\n"

	return classes



# compareClarnaOnescontacts() compares output files of ClaRNA tool with list of found contacts previously
def compareClarnaOnesContacts(classes, onesContacts):

	print "compareClarnaOnesContacts\n"

	f = 0
	for ones in onesContacts:
		ones["class"] = "nope"
#		print ones
		one = [ones["acc"].parent.parent.id+"_"+str(ones["resacc"])+"_"+ones["acc"].parent.resname.split()[0], ones["don"].parent.parent.id+"_"+str(ones["resdon"])+"_"+ones["don"].parent.resname.split()[0]]

		for c in classes:
			if one[0] == c["chainRes"][0] and one[1] == c["chainRes"][1] or one[0] == c["chainRes"][1] and one[1] == c["chainRes"][0]:
				ones["class"] = c["class"]
#				print ones, c
#				print
				f += 1
			else:
				pass
		
	print "\tMy tool contanct which fits to ClaRNA contacts -", f, "\n"
	return onesContacts



# makeSimilar() adds values of measured distance and angle of contact to list of proper dictionary store in overallRes list; input - contact, file with results, list with tuples of suspected
#		characteristic contacts with non dupplicates, list of dictionaries with contact name and lists with angle and distance values
def makeSimilar(onCont, resultFile, contactT, overallRes):

	print "function - makeSimilar\n"

	for i in onCont:			# fragment, ktory tworzy liste z typami oddzialywan oraz liste z typami oddzialywan z wartosciami katow
#		print i
		try:
			angle = countAngleAndDistance(i, resultFile)
			distance = angle[1]
			angle = angle[0]

			if (i["acc"].parent.resname.split()[0], i["acc"].id, i["don"].id, i["don"].parent.resname.split()[0], i["class"]) in contactT:
				for di in overallRes:
					if di["type"] == (i["acc"].parent.resname.split()[0], i["acc"].id, i["don"].id, i["don"].parent.resname.split()[0], i["class"]):
						di["angles"].append(angle)
						di["distances"].append(distance)
				
			elif (i["acc"].parent.resname.split()[0], i["acc"].id, i["don"].id, i["don"].parent.resname.split()[0], i["class"]) not in contactT:
				contactT.append((i["acc"].parent.resname.split()[0], i["acc"].id, i["don"].id, i["don"].parent.resname.split()[0], i["class"]))
				overallRes.append({"type": (i["acc"].parent.resname.split()[0], i["acc"].id, i["don"].id, i["don"].parent.resname.split()[0], i["class"]), "distances": [distance], "angles": [angle]})
								
			else:
				pass
 			
		except TypeError:
			print "\n\tPROBLEM WITH makeSimilar() - line 204\n\t"
				


# makeLog() creates file with structures names which are explore; input - path to log file, structure name, condiction if process starts or ends, number of found contacts
def makeLog(filePath, string, where, con):

	if where == "in":
		f = open(filePath, "a")
		f.writelines([string+"\t"])
		f.close()
	elif where == "out":
		f = open(filePath, "a")
		f.writelines([string+"\t"+str(con)+"\n"])
		f.close()
	else: pass



# parseToFile() creates file with lists with values of found distances and angles 
def pasteToFile(directory, filePath):

	print "fanction - pasteToFile\n\n\tSAVING FOUND VALUES TO FILE...\n"

	g = open(filePath, "w")
	for i in directory:
		if "?" in i["type"][-1]:
			new = "A"
			for l in i["type"][-1]:
				if l == "?":
					new = new + "q"
				else:
					new = new + l
			brand = []
			for p in i["type"][:-1]:
				brand.append(p)
			brand.append(new[1:])
			i["type"] = brand
		else: pass

		g.writelines("_".join(i["type"])+"_angle"+" = ")
		g.writelines("[")
		for e in i["angles"]:
			l = len(i["angles"])
			if i["angles"].index(e) == l - 1:
				g.writelines(str(e))
			else:
				g.writelines(str(e)+", ")
		g.writelines("]"+"\n")

		g.writelines("_".join(i["type"])+"_distance"+" = ")
		g.writelines("[")
		for e in i["distances"]:
			l = len(i["distances"])
			if i["distances"].index(e) == l - 1:
				g.writelines(str(e))
			else:
				g.writelines(str(e)+", ")
		g.writelines("]"+"\n\n")

	g.close()



# moveTool() run all function in proper order
def moveTool():

	outs = "/home/mateusz/Studia/mgr/refinement/" # <<<--- to fill										# path to the directory which contains all output files
	g = "/home/mateusz/Studia/mgr/refinement/s/" # <<<--- to fill										# path to the directory which containes PDB files with RNA structures
	start = 1

	resultName = outs + "results"											# filename with results table
	resultFile = open(resultName, "w+")
	resultFile.writelines(["PDB name"+'\t'+"chain Acc"+'\t'+"chain Don"+'\t'+"no. Acc"+'\t'+"base Acc"+'\t'+"atom Acc"+'\t'+"atom Don"+'\t'+"atom Cx"+'\t'+"base Don"+'\t'+"no. Don"+'\t'+"angle on 'atom Nx' [degree]"+'\t'+"distance Acc - Don' [A]"+"\t"+"class"+'\n'])


	clarnaPath = "/home/mateusz/Studia/ClaRNA_play/clarna_run.py" # <<<--- to fill				# path to ClaRNA 
	clarnaWorkFile = outs + "ClaRNAwork"
	stop = len(os.listdir(g))

	log = outs + "toolLog"
	os.system("touch "+log)

	anglesDistancesFile = outs + "anglesDistancesLists.py"

	f = 0
	contactTypes = []
	overallRes = []
	bondClasses = []
	basePairs = []
	angleBasePairs = []

	distance = 4.5	# <<<--- to fill					# distance of potent contact

	for l in getStrNames(g):

		print "\n\t--------------------- NEW STRUCTURE ---------------------\n"
		print "\tWORKING\ton: "+l+", "+"\t"+str(start)+" / "+str(stop)+"\n"

		makeLog(log, l, "in", con = 0)									# makeLog() saves structure name into file when proces starts
		
		# looking for contacts section
		singleModelContacts = []											# list of found contacts between atoms of all structure
		try:
			structure = openStructure(g+"/"+l)							# openStructure() takes path to structure and returns structure as a representation of BioPython instance 
		except PDBExceptions.PDBConstructionException:
			pass

		try:
			atomsChainsTuple = parserPDB(structure)					# parserPDB() takes structure and returns tuple (all atoms list, all chains structures list)
			NS = NeighborSearch(atomsChainsTuple[0])					# NeighbourSearch() is a class which helps with searching potential contacts
			for chain in atomsChainsTuple[1]:							# loop iterates single chain list; looking for contacts of the chain
				for cont in findHBond(chain, NS, distance):			# loop iterates list from findHBond() ...
					if cont in singleModelContacts:						# and pass existing contact on singleModelContacts list ...
						pass	
					elif cont not in singleModelContacts:				# or add contact to the singleModelContacts list when contact isn't on it
						singleModelContacts.append(cont)
					else: pass
		except TypeError:
			pass
		
		singleModelContacts.sort()
		print "\tNumber of found contacts by this tool -", len(singleModelContacts), "\n"
	
		conta = len(singleModelContacts)

		# ClaRNA tool section
		classes = makeClaRNA(l, g, clarnaWorkFile, 0.4, clarnaPath)				# ClaRNA contact classifier tool; makeClaRNA() returns list with classified resiues contacts
		onesContacts = compareClarnaOnesContacts(classes, singleModelContacts)				# compareClarnaOnesContacts() compare classified contacts with singleModelContacts; returns 

		makeSimilar(onesContacts, resultFile, contactTypes, overallRes)						# makeSimilar() add values of the angle and distance to proper dictionary on ovarallRes list
		makeLog(log, l, "out", conta)																		# makeLog saves structure name into file when process ends

		start += 1
		o = time.time()
		print "\n\t^^^^^^^^^^^^^^^^^^^^^^ END ^^^^^^^^^^^^^^^^^^^^^^ \n"

	o = time.time()
	pasteToFile(overallRes, anglesDistancesFile)														# parseToFile() exports, saves values of the angles and distances into the file
	print "parseToFile time - ", round(time.time() - o, 2), "seconds \n"

	resultFile.close()




if __name__ == "__main__":

	ddd = time.time()
	moveTool()
	seconds = time.time() - ddd
	print "\t\t\t\tOverall time:"
	print "\t\t\t\t\t- in seconds\t", round(seconds, 2)
	print "\t\t\t\t\t- in munutes\t", round(seconds / 60, 2)


