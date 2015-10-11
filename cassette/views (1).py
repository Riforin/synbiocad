# Django imports 
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from .models import Results
from .forms import *

import os
from django.conf import settings
PROJECT_ROOT = settings.PROJECT_ROOT

# Imports for cassette generation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp
import copy
from intermine.webservice import Service
from pandas import *
from pandas import DataFrame, read_csv
import pandas as pd #this is how I usually import pandas
import matplotlib.pyplot as plt

# define global variables
HomologyLength = 100
PrimerMaxTm = 55
PrimerMaxLen = 60
OverhangMaxFrac = 1

def index(request):
	# if this is a POST request we need to process the form data
	if request.method == 'POST':
		# create a form instance and populate it with data from the request
		form = IndexForm(request.POST)
		# validate
		if form.is_valid():
			answer = request.POST['choices']

			if answer == '1':
				return HttpResponseRedirect('/cassette/charlocus/')
			elif answer == '2':
				return HttpResponseRedirect('/cassette/existinglocus/')
			elif answer == '3':
				return HttpResponseRedirect('/cassette/customlocus/')
	# If a GET (or any other method)we'll create a blank form) 
	else:
		form = IndexForm()

	return render(request, 'cassette/index.html', {'form': form})
# Pages possible from first branch
def charLocus(request):
	if request.method == 'POST':
		# create a form instance and populate it with data from the request
		form = CharLocusForm(request.POST)
		# validate
		if form.is_valid():
			answer = request.POST['choices']
			cutsite = request.POST['cutsite']
			request.session['cutsite'] = cutsite

			if answer == '1':
				return HttpResponseRedirect('/cassette/charlocus/prebuilt')
			elif answer == '2':
				return HttpResponseRedirect('/cassette/charlocus/standard')
			
	# If a GET (or any other method)we'll create a blank form) 
	else:
		form = CharLocusForm()

	return render(request, 'cassette/charloc.html', {'form': form})

def charPreBuilt(request):
	if request.method == 'POST':
		form = CharLocusPrebuiltForm(request.POST)
		if form.is_valid():
			name = request.POST['name']
			sequence = request.POST['sequence']
			cutsite = request.session.get('cutsite')

			Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editEmptyPrebuilt(name, sequence, cutsite)
			request.session['Lup'] = Lup
			request.session['Rup'] = Rup
			request.session['Ldown'] = Ldown
			request.session['Rdown'] = Rdown
			request.session['L'] = L
			request.session['R'] = R
			request.session['seqLen'] = seqLen
			request.session['donorSeq'] = donorSeq
			return HttpResponseRedirect('/cassette/results')
	else:
		form = CharLocusPrebuiltForm()

	return render(request, 'cassette/charlocprebuilt.html', {'form': form})

def results(request):
	Lup = request.session.get('Lup')
	Rup = request.session.get('Rup')
	Ldown = request.session.get('Ldown')
	Rdown = request.session.get('Rdown')
	L = request.session.get('L')
	R = request.session.get('R')
	seqLen = request.session.get('seqLen')
	donorSeq = request.session.get('donorSeq')

	return render(request, 'cassette/results.html', {'Lup': Lup, 'Rup': Rup, 'Ldown': Ldown, 'Rdown': Rdown, 'L': L, 'R': R, 'seqLen': seqLen, 'donorSeq': donorSeq})


def charStdORF(request):
	template_name = 'cassette/customloc/stdcustomorf.html'


# Pages possible from second branch
def existingLocus(request):
	template_name = 'cassette/existloc.html'

def existingDelete(request):
	template_name = 'cassette/existloc/delete.html'

def existingPremade(request):
	template_name = 'cassette/existloc/premade.html'

def existingStd(request):
	template_name = 'cassette/existloc/std.html'

def existingCustom(request):
	template_name = 'cassette/existloc/custom.html'

def existingNearby(request):
	template_name = 'cassette/existloc/nearby.html'

# Pages possiblefrom third branch
def customLocus(request):
	template_name = 'cassette/customloc.html'

def customNonVary(request):
	template_name = 'cassette/customloc/novary.html'

def customVary(request):
	template_name = 'cassette/customloc/vary.html'

#---------------  Helper Functions  ---------------------#
def editEmptyPrebuilt(name, sequence, cutname):
	labels=['208a', '1014a', '1114a', '607c', '308a', '1021b', '720a']
	ChrLetters=["Scer02", "Scer10", "Scer11", "Scer06", "Scer03", "Scer10", "Scer07"]
	ExpValues=[1.0, 1.2, 1.1, 1.1, 1.4, 1.0, 1.0]
	cutSeqs=["GTCCGCTAAACAAAAGATCT", "TTATGTGCGTATTGCTTTCA", "CTTGTGAAACAAATAATTGG", "CTATTTTTGCTTTCTGCACA", "TAGGATACAAGTGTCTGTGC", "CCTCTGTGTGGTGGTAATTG", "CAACAATTGTTACAATAGTA"]

	cutArray={'name' : Series(labels, index=labels),
				'exp. lev.' : Series(ExpValues, index=labels),
				'chrom. loc.' : Series(ChrLetters, index=labels),
				'sequence' : Series(cutSeqs, index=labels)
			 }

	cutFrame= DataFrame(cutArray)
	location=cutFrame.loc[cutname,'chrom. loc.']+".fasta"
	cutSequence=cutFrame.loc[cutname,'sequence']

	ChromosomeSeq=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\\" + location), "fasta").seq
	
	if ChromosomeSeq.find(cutSequence)==-1:
		ChromosomeSeq=ChromosomeSeq.reverse_complement()

	StartIndex=ChromosomeSeq.find(cutSequence)
	EndIndex=StartIndex+34
	
	UpSeq=ChromosomeSeq[StartIndex-HomologyLength:StartIndex]
	DownSeq=ChromosomeSeq[EndIndex:EndIndex+HomologyLength]
		
	UpHomRec = SeqRecord(UpSeq, id=cutname)
	DownHomRec = SeqRecord(DownSeq, id=cutname)

	orfRecord = SeqRecord(sequence, id=name)

	fragments = [UpHomRec, orfRecord, DownHomRec]

	return stitch(fragments)

def fetchGene(GeneName):
	
	service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
	template = service.get_template('Gene_GenomicDNA')

	rows = template.rows(
		E = {"op": "LOOKUP", "value": GeneName, "extra_value": "S. cerevisiae"}
	)
	
	# this service seems to return multiple similar genes but we want the first one only, so count
	# and it returns information about the gene you want
	count=0
	for row in rows:
		
		count=count+1
		if count==1:
			descr= row["description"]
			GeneSeq=Seq(row["sequence.residues"])
			GeneSysName=row["secondaryIdentifier"]
			#print(" ")
			#print("I think you want...... "+row["secondaryIdentifier"])
			#print(row["description"])
			#print(" ")
			#print(row["sequence.residues"])
			#print(" ")
			#print("Good choice! I have a feeling you're going to get lucky with this one.")
			#print(" ")
			#print("Give me a second to put some of my ducks in a circle...")
	   

			
	#let's create a record for the oldGene
	GeneRecord = SeqRecord(GeneSeq, id=GeneSysName)
	
	#now let's add some more information to make it useful
	GeneRecord.name=GeneName
	GeneRecord.features=GeneSysName
	GeneRecord.description=descr

	return GeneRecord 

def fetchNeighbor(NeighborRecord, direction, distance):


	# let's load the appropriate chromosome file. The record of the gene we looked up
	# contains in the "features" the systematic name, wherein the second letter
	# corresponds to chromosome number, e.g., 1=A etc
	if NeighborRecord.features[1]=="A":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer01.fasta"), "fasta")
	if NeighborRecord.features[1]=="B":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer02.fasta"), "fasta")
	if NeighborRecord.features[1]=="C":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer03.fasta"), "fasta")
	if NeighborRecord.features[1]=="D":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer04.fasta"), "fasta")
	if NeighborRecord.features[1]=="E":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer05.fasta"), "fasta")
	if NeighborRecord.features[1]=="F":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer06.fasta"), "fasta")
	if NeighborRecord.features[1]=="G":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer07.fasta"), "fasta")
	if NeighborRecord.features[1]=="H":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer08.fasta"), "fasta")
	if NeighborRecord.features[1]=="I":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer09.fasta"), "fasta")
	if NeighborRecord.features[1]=="J":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer10.fasta"), "fasta")
	if NeighborRecord.features[1]=="K":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer11.fasta"), "fasta")
	if NeighborRecord.features[1]=="L":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer12.fasta"), "fasta")
	if NeighborRecord.features[1]=="M":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer13.fasta"), "fasta")
	if NeighborRecord.features[1]=="N":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer14.fasta"), "fasta")
	if NeighborRecord.features[1]=="O":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer15.fasta"), "fasta")
	if NeighborRecord.features[1]=="P":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer16.fasta"), "fasta") 

	
	
	# let's explicitely name the sequences from the seq record
	NeighborSeq=NeighborRecord.seq
	ChromosomeSeq=ChromosomeRec.seq
	
	# flip the sequence to orient with respect to the old gene
	if ChromosomeSeq.find(NeighborSeq)==-1:
		ChromosomeSeq=ChromosomeSeq.reverse_complement()

	StartIndex=ChromosomeSeq.find(NeighborSeq)
	EndIndex=StartIndex+len(NeighborSeq)
	
	if direction=="upstream":
		DesiredSeq=ChromosomeSeq[StartIndex-distance:StartIndex]
	if direction=="downstream":
		DesiredSeq=ChromosomeSeq[EndIndex:EndIndex+distance]

	
	
	
	NeighborRec = SeqRecord(DesiredSeq, id=NeighborRecord.name)
	
	return NeighborRec

	#print(NeighborRec)

def getPrimer(currRecord):
	

	mp = 0
	length = 0
	primer = Seq("")

	seq=currRecord.seq
	
	while mp <= PrimerMaxTm and length <= PrimerMaxLen:
		primer = primer + seq[length]
		mp = MeltingTemp.Tm_staluc(primer)
		length += 1

	return primer           
		
def overhangPrimer(currRecord,prevSeq):
	#let's get the template-binding primer first
	primer=getPrimer(currRecord)
	
	
	#OK let's work on the overhang
	maxOhLen=PrimerMaxLen-len(primer)    
	maxFrac=1
	
	#let's decide on a max overhang length
	if round(len(primer)*(OverhangMaxFrac+1)) < 60:
			 maxOhLen=round(len(primer)*OverhangMaxFrac)
	
	#the index must be an integer!!!
	maxOhLen=int(maxOhLen)
	ohprimer=prevSeq.seq[-maxOhLen:]+primer #we add the .seq so that it returns a string
	
	return ohprimer      

def stitch(fragments):
	#this function takes seq records and prints primers

	#let's make an empty sequence file
	Nfrags=len(fragments)
	donor=Seq("")
	index=[]
	print("")
	for i in range (0, Nfrags):
		donor=donor+fragments[i]
	return fragments[1].reverse_complement(), 1, 2, 3,  4, 5, 6, 7
	# Dummy assignment setup to allow for compilation
	Lup = ""
	Rup = ""
	Ldown = ""
	Rdown = ""
	L = ""
	R = ""

	for i in range (0, Nfrags):
		if i==0:
			Lup = "Lup: "+ fragments[i].id + " " + getPrimer(donor)
			Rup = "Rup: "+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement())
		elif i==Nfrags-1:
			Ldown = "Ldown: "+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1])
			Rdown = "Rdown: "+ fragments[i].id + " " + getPrimer(donor.reverse_complement())
		else:
			L = "L: "+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1])
			R = "R: "+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement())

	sequenceLength = len(donor.seq)
	donorSequence = donor.seq

	return str(Lup), str(Rup), str(Ldown), str(Rdown), str(L), str(R), "Sequence Length: " + str(sequenceLength), "Sequence: " + str(donorSequence)

# Modified functions for input
def standardCassette(PromoterName,TerminatorName, orfName, orfSeq):
	
	#first, the promoter
	print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
	print("First, which PROMOTER do you want to use, e.g., TDH3")
	
	PromoterGeneRec=fetchGene(PromoterName)
	PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
	PromoterRec.id=PromoterRec.id+"ps"
	
	
	#second, the terminator
	print("Which TERMINATOR do you want to use, e.g., ADH1")
	TerminatorGeneRec=fetchGene(TerminatorName)
	TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
	TerminatorRec.id=TerminatorRec.id+"ts"
	
	
	#and last, the gene
	print("What is the name of your gene, e.g., KlGapDH")
	
	print("What's the sequence")
	
	orfRecord=SeqRecord(Seq(orfSeq), id=orfName)
	
	insertRec=[PromoterRec,orfRecord,TerminatorRec]
	return PromoterRec, orfRecord, TerminatorRe

def buildCassette(PromoterName, TerminatorName, orfName, orfSeq):
	
	#first, the promoter
	print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
	print("First, which PROMOTER do you want to use, e.g., TDH3")
	
	PromoterGeneRec=fetchGene(PromoterName)
	PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
	PromoterRec.id=PromoterRec.id+"ps"
	
	
	#second, the terminator
	print("Which TERMINATOR do you want to use, e.g., ADH1")
	TerminatorGeneRec=fetchGene(TerminatorName)
	TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
	TerminatorRec.id=TerminatorRec.id+"ts"
	
	#and last, the gene
	print("What is the name of your gene, e.g., KlGapDH")
	
	print("What's the sequence")
	
	orfRecord=SeqRecord(Seq(orfSeq), id=orfName)
	
	insertRec=[PromoterRec,orfRecord,TerminatorRec]
	return PromoterRec, orfRecord, TerminatorR

# This must be changed to be dynamic 
def variableCassette(N, toVary=0, variants=0):
	print("")
	print("Let's start building.")
	print("")
	
	# Store both name and sequence in a SeqRecord
	# Append them to a list
	# Return list as fragments to be stitched

	records = []
	for n in range(N):
		name = input("What is the name of sequence " + str(n+1) +":")
		sequence = input("What is the sequence of this fragment:")
		print("")
		Rec = SeqRecord(Seq(sequence), id = str(n+1))
		Rec.name = name
		records.append(Rec)

	variantRecords = []
	variantRecords.append(records)
	# This only happens if there are variants.
	if variants > 0:
		print("Time to make those variants you wanted.")
		for n in range(variants-1):
			name = input("What is the name of variant " + str(n+1) + ":")
			sequence = input("What is the sequence of this variant:")
			Rec = SeqRecord(Seq(sequence), id = str(n+1))
			Rec.name = name
			# Make a copy of the original, switch the fragments and add it to the list. 
			# Deep-copy ensures there are no pointer issues
			tempVariant = copy.deepcopy(records)
			tempVariant[toVary - 1] = Rec
			variantRecords.append(copy.deepcopy(tempVariant))
			print("")

	# Returns a list of lists of the SeqRecords of the fragments
	return variantRecords
