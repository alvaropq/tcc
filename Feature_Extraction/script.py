#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import collections
import argparse
from Bio import SeqIO
from itertools import product
import math,string,re,sys,getopt
import scipy.stats
import statistics as st
from scipy.fftpack import fft, ifft
import random
import pandas as pd
from igraph import *
import os
import operator
import scipy.stats


######################################    KMERS   #####################################################

def perms():
    global listTotal, caracteres
    caracteres = ["A", "C", "G", "T"]
    listTotal = []
    for k in range(1, ksize+1):
        permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
        # print (permsList)
        file = open(foutput, 'a')
        for perm in permsList:
            # print (perm)
            listTotal.append(perm)
            file.write("%s," % (str(perm)))
    file.write("class")
    file.write("\n")
    return


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return        
    

def chunksTwo(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def fileRecord():
    file = open(foutput, 'a')
    for x in probabilities:
        # print (x)
        file.write("%s," % (str(x[1])))
    file.write(labelDataset)
    file.write("\n")
    print ("Recorded Sequence!!!")
    return
    

def findKmers():
    global probabilities
    perms()
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        probabilities = []
        for k in range(1, ksize+1):
            kmer = {}
            totalWindows = (len(seq) - k) + 1 # (L - k + 1)
            permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
            for key in permsList:
                kmer[key] = 0
            kmer = collections.OrderedDict(sorted(kmer.items()))
            for subseq in chunksTwo(seq, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print (key)
                # print (value)
                probabilities.append([str(key), value/totalWindows])
        fileRecord()
    return

        
######################################################################################################   

########################################          KGAP               #################################


def kgap_record(name_seq):
    file = open(foutput, 'a')
    file.write('%s,' % name_seq)
    for x in number_kmers:
        # print (x)
        file.write('%s,' % (str(x[1])))
    file.write(labelDataset)
    file.write('\n')
    print('Recorded Sequence: %s' % name_seq)
    return


def chunks_kgap(seq, gap, before, after):
    seqlen = len(seq)
    chunks = []
    for i in range(0, seqlen):
        j = i + after + gap + before
        if j < seqlen:
            chunks.append(str(seq[i:j]))
            # print(str(seq[i:j]))
    return chunks


def header_kgap(permsList):
    file = open(foutput, 'a')
    file.write('%s,' % 'nameseq')
    for perm in permsList:
        file.write('%s,' % str(perm))
    file.write('label')
    file.write('\n')
    return


def compare(chunks, label):
    table = {}
    # total = len(chunks)
    for i in label:
        count = 0
        for j in chunks:
            aux = j[:before]
            aux2 = j[before + k:]
            aux3 = i[:before]
            aux4 = i[before + k:]
            if aux == aux3 and aux2 == aux4:
                count = count + 1

        # prob = count/total
        prob = count
        dic = {i: prob}
        table.update(dic)
    
    return table


def kgap():
    global number_kmers
    i = 0
    label = []
    gap = ''
    while i < k:
        gap = gap + '_'
        i = i + 1

    permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=before)]
    permsList2 = [''.join(str(i) for i in x) for x in product(caracteres, repeat=after)]

    for i in permsList:
        for j in permsList2:
            char = i + gap + j
            label.append(char)
    header_kgap(label)

    for seq_record in SeqIO.parse(finput, 'fasta'):
        number_kmers = []
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        # total = len(seq)
        chunks = chunks_kgap(seq, k, before, after)
        #print(chunks)
        table = compare(chunks, label)
        for key, value in table.items():
            number_kmers.append([str(key), value])
        kgap_record(name_seq)
    return

###############################################################################################################



####################################  PseKNC  #################################################################

def pseknc_record(listy,name_seq):
    i = 0
    z = 1
    listz = ''
    # Generating output file format
    while i < len(listy):
        if i < 1:
            listz = listz + str(listy[i])
            i = i + 1
        else:
            listz = listz + "," + str(listy[i])
            i = i + 1
    # Generating output file format
    output = name_seq +',' +listz +',' +labelDataset+'\n'
    # Writing to output file
    OutFile = open(foutput, 'a')
    OutFile.write(output)
    OutFile.close()
    print ('Recorded Sequence: %s' % (name_seq))
    return

def header_pseknc(permsList):
    file = open(foutput,'w')
    file.write('%s,' % ('nameseq'))
    i=0
        #print (permsList)
    while i < len(permsList):
            # print (perm)
        file.write('O'+str(i)+',')
        i = i+1
    
    file.write('label')
    file.write('\n')
    return

def header_pseknc2(permsList):
    file = open(foutput,'w')
    file.write('%s,' % ('nameseq'))
    i=0
        #print (permsList)
    while i < len(permsList):
            # print (perm)
        file.write('T'+str(i)+',')
        i = i+1
    
    file.write('label')
    file.write('\n')
    return


#  generate_permutations
# ----------------------------------------
# inputs: chars = int
# outputs: all possible oligonucleotides of length k

def generate_permutations(chars):
	allowed_chars = ['A','T','G','C']
	status = []
	for tmp in range(chars) :
		status.append(0)
	last_char = len(allowed_chars)
	rows = []
	for x in range(last_char ** chars) :
		rows.append("")
		for y in range(chars - 1 , -1, -1) :
				key = status[y]
				rows[x] = allowed_chars[key] + rows[x]
		for pos in range(chars - 1, -1, -1) :
			if(status[pos] == last_char - 1) :
					status[pos] = 0
			else :
					status[pos] += 1
					break;
	
	return rows

#  _mean: Calculates mean value of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the mean value of listy

def _mean(listy):
	return float(sum(listy))/len(listy)

#  _std(listy): Calculates standard deviation of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the standard deviation of listy
 
def _std(listy,ddof=1):
	mean = _mean(listy)
	# for i in listy, (i - mean)^2 => list
	temp=[math.pow(i-mean,2) for i in listy]
	# temp is added together, divided by length of listy,
	# and the square root is taken
	res = math.sqrt(float(sum(temp))/len(listy))
	return res


#  sepSequence
# -----------------------------------------
# inputs: seq, string and k, int
# output: list of k-tuples

def sepSequence(seq, k):
	i = k-1
	seqq = []
	while i < len(seq):
		j = 0
		nuc = ''
		while j < k:
			nuc = seq[i-j] + nuc
			j = j + 1
		seqq.append(nuc)
		i += 1
	return seqq

#  simplePseKNC: Calculates the frequencies of all possible
#  oligonucleotides in the input sequence
# ----------------------------------------
# inputs: inputFile = a string, outputFile = a string, k = int, 
#         formatt = string
# output: nothing, writing to a file

def simplePseKNC(inputFile,outputFile,k,formatt,geneticMaterial):
    listy = ''
    # Need to generate all possible oligonucleotides
    olinucz = generate_permutations(k)
    header_pseknc2(olinucz)
    for seq_record in SeqIO.parse(inputFile,'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        listy = sepSequence(seq,k)
        z = 1
        OutFileName = outputFile
        OutFile = open(OutFileName, 'a')
        aux = 0
        listz=''
        for olinuc in olinucz:
            freq = listy.count(olinuc)
            freq = int(freq/len(listy)*1000)/1000.0
            if olinucz.index(olinuc) < 1:
                #print(olinucz.index(olinuc))
                #print('olinuc: ',olinuc)
                listz = listz + str(freq)
            else:
                listz = listz + "," + str(freq)
        OutFile.write(name_seq+',')
        OutFile.write(listz)
        OutFile.write(','+labelDataset+'\n')
        print ('Recorded Sequence: %s' % (name_seq))
    OutFile.close()
    return

#  getValues: Returns a line of values for one property
#             from physicochemical property files
# -----------------------------------------
# input: prop = string of one property and supInfo = string
# output: values = a string representing a list of property values

def getValues(prop,supInfo):
	values = ""
	name = re.search(prop,supInfo)
	if name:
		strr = prop + '\s*\,(.+)'
		b = re.search(strr, supInfo)
		if b:
			values = b.group(1)
	return values

#  getSpecificValue: Returns a property value for a specific di- or tri-
#                    nucleotide
# -----------------------------------------
# input: olinuc = string, prop = string, supInfo = string
# output: value = an int that is property value of an olinuc

def getSpecificValue(olinuc,olinucs,prop,values,supInfo):
	values = values.split(",")
	#valueS = [float(x) for x in values.split(",")]
	count = olinucs.index(olinuc)
	value = values[count]
	return float(value)


#  hn: Hn function 
# -----------------------------------------
# inputs: olinuc = string, prop = string, supInfo = string
# output: temp = int

def hn(olinuc,olinucs,prop,supInfo,values):
    #values = getValues(prop,supInfo).rstrip()
    h0 = float(getSpecificValue(olinuc,olinucs,prop,values,supInfo)) 
    valueS = [float(x) for x in values.split(",")]
    temp = float((h0 - _mean(valueS)) / _std(valueS))
    return temp

#  theta2: Theta(i,i+j) function, for Type I
# -----------------------------------------
# input: seq = string, props = string, i = int,
#	 Type = int, k = int, j = int, supInfo = string
# output: summ = int

def theta2(seq,olinucs,props,i,k,j,supInfo):
    summ = 0
    values = ''
    for prop in props:
        values = getValues(prop,supInfo).rstrip()
        hn1 = hn(seq[i-1],olinucs,prop,supInfo,values)
        hn2 = hn(seq[i+j-1],olinucs,prop,supInfo,values)
        subsqr = math.pow(float(hn1-hn2),2)
        summ = summ + subsqr
    return float(summ)/len(props)

#  J: J(i,i+j) function, for Type II 
# --------------------------------------
# inputs: seq = string, prop = string, i = int, 
#         k = int, j = int, supInfo = string
# output: product = int

def J(seq,olinucs,prop,i,k,j,supInfo):
    values = getValues(prop,supInfo)
    hn1 = hn(seq[i-1],olinucs,prop,supInfo,values)
    hn2 = hn(seq[i+j-1],olinucs,prop,supInfo,values)
    return float(hn1*hn2)

#  theta1: Theta(j) and Tau(LamGam) function
# -----------------------------------------
# input: seq = string, props = string, Type = int
#        k = int, j = int, supInfo = string
# output: final = int

def theta1(seq,olinucs,props,Type,k,j,supInfo):
    k = int(k)
    gamma = len(props)
    seqq = sepSequence(seq,k)
    i = 1
    a = 0
    if Type == 1:
        var = len(seq) - int(k) - int(j) + 1
        while i <= var:
            b = 0
            b = theta2(seqq,olinucs,props,i,k,j,supInfo)
            a = a + b
            i = i + 1
    else:
        ii = 0
        var = len(seq) - int(k) - int(j/gamma)
        while i <= var:
            b = 0
            if ii == gamma:
                ii = 0
                b = J(seqq,olinucs,props[ii],i,k,int(j/gamma),supInfo) 
                a = a + b
            else:
                b = J(seqq,olinucs,props[ii],i,k,int(j/gamma),supInfo)
                a = a + b
            ii = ii + 1
            i = i + 1
    final = float(a)/var
    return final

#  pseKNCHelper: Creates list of adjusted frequency values for 4^k 
#  oligonucleotides and the lambda terms
# -----------------------------------------
# input: seq = string, props = string, Type = int, k = int, j = int, w = int
# output: freqs = list of ints

def pseKNCHelper(seq,olinucs,props,Type,k,j,w,geneticMaterial,supInfo):
    gamma = len(props)
    freqs = []
    seqq = sepSequence(seq,k)
    olinucs = olinucs.split(",")
    for olinuc in olinucs:
        freq = seqq.count(olinuc)
        freq = float(freq/len(seqq))
        total = 0
        i = 1
        if Type == 2:
            j = j/gamma
        while i <= j:
            total = total + theta1(seq,olinucs,props,Type,k,i,supInfo)
            i = i + 1
        total = float(freq/(1 + (float(w)*total)))
        total = int(total*1000) / 1000.0
        freqs.append(total)
    #Computing Lambda terms...
    fourK = math.pow(4,k)
    mu = fourK + 1
    while (fourK+1) <= mu <= (fourK + j):
        top = float(w) * theta1(seq,olinucs,props,Type,k,int(mu-fourK),supInfo)
        bottomTheta = 0
        bottom = 0
        i = 1
        while 1 <= i <= j:
            bottomTheta += theta1(seq,olinucs,props,Type,k,i,supInfo)
            i += 1
        bottom = 1 + (float(w) * bottomTheta)
        term = float(top / bottom)
        term = int(term * 1000) / 1000.0
        freqs.append(term)
        mu += 1
    return freqs

#  pseKNC: Opens input and output files, calls functions to calculate
#  values and writes outputs to the output file in appropriate format
# -----------------------------------------
# input: inputFile = string, outputFile = string, propNames = string,
#        Type = int, k = int, j = int, w = int, formatt = string
# output: nothing, write to output file

def pseKNC(inputFile,outputFile,props,Type,k,j,w,formatt,geneticMaterial):
    j = float(j)
    listy = ''
    output = ''
    outputs = ''
    gamma = len(props)
    #Getting supporting info from files
    if k == 2:
        SupFileName = 'Supporting Information S1 ' + geneticMaterial + '.txt'
    else:
        SupFileName = 'Supporting Information S3 DNA.txt'
    SupFile = open(SupFileName,'r')
    supInfo = SupFile.read()
    o = re.search('Physicochemical properties\,(.+)\n',supInfo)
    olinucs = ''
    if o:		
        olinucs = o.group(1).rstrip()
    SupFile.close()
    # Calculating frequencies
    aux = 0
    for seq_record in SeqIO.parse(inputFile,'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name	
        if Type == 1:
            listy = pseKNCHelper(seq,olinucs,props,Type,k,j,w,geneticMaterial,supInfo)
        else:
            listy = pseKNCHelper(seq,olinucs,props,Type,k,j*gamma,w,geneticMaterial,supInfo)
        if aux == 0:
            header_pseknc(listy)
            aux = 1
        pseknc_record(listy,name_seq)
    return

def parameters(x,t,k,j,w,s,finput,foutput,seq):
    inputFile = finput
    outputFile = foutput

    if seq == 2:
        geneticMaterial = 'RNA'
    else:
        geneticMaterial = 'DNA' 

    # Getting list of properties (ex. Tilt, shift)
    properties = ''
    props = []
    propFileName = x
    propFile = open(propFileName, 'r')
    for line in propFile:
        e = re.search('(.+)', line)
        if e:
            properties = properties + e.group(1) + ","
    propFile.close()
    props = properties.split(",")
    props = props[0:len(props)-1]
    
    # Checking if type argument is 1 or 2
    if t == 1 or t == 2:
        Type = t
    else:
        print('The Type argument must be 1 or 2.')
        sys.exit()

    # Checking if k argument is 2 or 3
    if k == 2 or k == 3:
        kay = k
    else:
        print('The K argument must be 2 or 3.')
        sys.exit()

    # Checking if weight argument is between 0 and 1.0
    if w > 0 and w <= 1:
	    weight = w
    else:
        print('The weight factor argument must be between 0.1 and 1.0.')
        sys.exit()
    
    # Checking if the lambda argument is a whole number and smaller than L-k
    InFileName = inputFile
    InFile = open(InFileName,'r')
    listy = ''
    for line in InFile:
        g = re.search('\>',line)
        if g:
            label = line
        else:
            listy = listy + line.rstrip()
    InFile.close()
    ell = len(listy)
    if (float(j) % 1) == 0 and (1 <= float(j) < (ell - kay)):
        lam = j
    else:
        print('Lambda must be a whole number and smaller than the length of the query sequence minus the k-tuple number (lambda < L-k).')
        print('Length of query sequence: ',ell)
        print('k: ',kay)
        sys.exit()
    
    formatt = 'csv'

    if kay == 2:
        SupFileName = 'Supporting Information S1 ' + geneticMaterial + '.txt'
    else:
        SupFileName = 'Supporting Information S3 DNA.txt'
	
    #SupFile = open(SupFileName,'r')
    
    tt = 0

    for prop in props:
        with open(SupFileName,'r')as SupFile:
            for line in SupFile:
                t = re.search(prop,line)
                if t:
                    tt=1
            if not(tt==1):
                print('"' + prop + ' was not found in the supporting information file."')
                print('Please check that the k value and the query properties correspond.')
                sys.exit()
        SupFile.close()

    if s == 1:
        simplePseKNC(inputFile,outputFile,int(kay),formatt,geneticMaterial)
    else:
        if kay == 2 or kay == 3:
            kay = int(kay)
        else:
            print('The k-tuple argument must be 2 or 3 (dinucleotides or trinucleotides).')
            sys.exit()
        pseKNC(inputFile,outputFile,props,Type,kay,lam,weight,formatt,geneticMaterial)

    return

########################################################################################################

#################################   Fourier-Class   ####################################################


def header_fourier():
	dataset = open(foutput, 'a')
	dataset.write("nameseq,average,median,maximum,minimum,peak,"
                  + "none_levated_peak,sample_standard_deviation,population_standard_deviation,"
                  + "percentile15,percentile25,percentile50,percentile75,amplitude,"
                  + "variance,interquartile_range,semi_interquartile_range,"
                  + "coefficient_of_variation,skewness,kurtosis,label")
	dataset.write("\n")
	return

        
def file_record_fourier():
    dataset = open(foutput, 'a')
    dataset.write("%s," % (str(name_seq)))
    for metric in features:
        dataset.write("%s," % (metric))
        # dataset.write("{0:.4f},".format(metric))
    dataset.write(label_dataset)
    dataset.write("\n")
    print("Sequence Analyzed!!")
    return


def feature_extraction_fourier():
    global features
    features = []
    average = sum(spectrum)/len(spectrum)
    features.append(average)
    ###################################
    median = np.median(spectrum)
    features.append(median)
	###################################
    maximum = np.max(spectrum)
    features.append(maximum)
    ###################################
    minimum = np.min(spectrum)
    features.append(minimum)
    ###################################
    peak = (len(spectrum)/3)/(average)
    features.append(peak)
    ###################################
    peak_two = (len(spectrumTwo)/3)/(np.mean(spectrumTwo))
    features.append(peak_two)
    ###################################
    standard_deviation = np.std(spectrum) # standard deviation
    features.append(standard_deviation)
    ###################################
    standard_deviation_pop = st.stdev(spectrum) # population sample standard deviation 
    features.append(standard_deviation_pop)
    ###################################
    percentile15 = np.percentile(spectrum, 15)
    features.append(percentile15)
    ###################################
    percentile25 = np.percentile(spectrum, 25)
    features.append(percentile25)
    ###################################
    percentile50 = np.percentile(spectrum, 50)
    features.append(percentile50)
    ###################################
    percentile75 = np.percentile(spectrum, 75)
    features.append(percentile75)
    ###################################
    amplitude = maximum - minimum
    features.append(amplitude)
    ###################################
    # mode = statistics.mode(spectrum)
    ###################################
    variance = st.variance(spectrum)
    features.append(variance)
    ###################################
    interquartile_range = np.percentile(spectrum, 75) - np.percentile(spectrum, 25)
    features.append(interquartile_range)
    ###################################
    semi_interquartile_range = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25))/2 
    features.append(semi_interquartile_range)
    ###################################
    coefficient_of_variation = standard_deviation/average
    features.append(coefficient_of_variation)
    ###################################
    skewness = (3 * (average - median))/standard_deviation
    features.append(skewness)   
    ###################################
    kurtosis = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25)) / (2 * (np.percentile(spectrum, 90) - np.percentile(spectrum, 10))) 
    features.append(kurtosis)
    ###################################
    return


def binary_fourier():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        A = []
        C = []
        T = []
        G = []
        for nucle in seq:
            if nucle == "A":
                A.append(1)
            else:
                A.append(0)
            if nucle == "C":
                C.append(1)
            else:
                C.append(0)
            if nucle == "T" or nucle =="U":
                T.append(1)
            else:
                T.append(0)
            if nucle == "G":
                G.append(1)
            else:
                G.append(0)
        FA = fft(A)
        FC = fft(C)
        FT = fft(T)
        FG = fft(G)
        for i in range(len(seq)):
        	specTotal = (abs(FA[i])**2) + (abs(FC[i])**2) + (abs(FT[i])**2) + (abs(FG[i])**2)
        	specTwo = (abs(FA[i])) + (abs(FC[i])) + (abs(FT[i])) + (abs(FG[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("A -- %s" % (A))
        # print("C -- %s" % (C))
        # print("T -- %s" % (T))
        # print("G -- %s" % (G))
        # print("\n")
        # print("A -- %s" % (abs(FA)))
        # print("C -- %s" % (abs(FC)))
        # print("T -- %s" % (abs(FT)))
        # print("G -- %s" % (abs(FG)))
        # print("\n")
        # print(spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def zcurve_fourier():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        ###################################
        ###################################
        R = 0 # x[n] = (An + Gn) − (Cn + Tn) ≡ Rn − Yn
        Y = 0
        M = 0 # y[n] = (An + Cn) − (Gn + Tn) ≡ Mn − Kn
        K = 0
        W = 0 # z[n] = (An + Tn) − (Cn + Gn) ≡ Wn − Sn
        S = 0
        ###################################
        ###################################
        x = []
        y = []
        z = []
        for nucle in seq:
            if nucle == "A" or nucle == "G":
            	R += 1
            	x.append((R)-(Y))
            else:
            	Y += 1
            	x.append((R)-(Y))
            if nucle == "A" or nucle == "C":
                M += 1
                y.append((M)-(K))
            else:
                K += 1
                y.append((M)-(K))
            if nucle == "A" or nucle == "T" or nucle == "U":
                W += 1
                z.append((W)-(S))
            else:
                S += 1
                z.append((W)-(S))
        FX = fft(x)
        FY = fft(y)
        FZ = fft(z)
        for i in range(len(seq)):
        	specTotal = (abs(FX[i])**2) + (abs(FY[i])**2) + (abs(FZ[i])**2)
        	specTwo = (abs(FX[i])) + (abs(FY[i])) + (abs(FZ[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print("\n")
        # print("X -- %s" % (x))
        # print("Y -- %s" % (y))
        # print("Z -- %s" % (z))
        # print("\n")
        # print("X -- %s" % (abs(FX)))
        # print("Y -- %s" % (abs(FY)))
        # print("Z -- %s" % (abs(FZ)))
        # print("\n")
        # print(spectrum)
        # print("\n")
        # print(spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def integer_fourier():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        integer = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
            	integer.append(0)
            elif nucle == "C":
            	integer.append(1)
            elif nucle == "A":
            	integer.append(2)
            else:
            	integer.append(3)
        FI = fft(integer)
        for i in range(len(seq)):
        	specTotal = (abs(FI[i])**2)
        	specTwo = (abs(FI[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("I -- %s" % (integer))
        # print("\n")
        # print("I -- %s" % (abs(FI)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def real_fourier():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        real = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
            	real.append(1.5)
            elif nucle == "C":
            	real.append(0.5)
            elif nucle == "A":
            	real.append(-1.5)
            else:
            	real.append(-0.5)
        FR = fft(real)
        for i in range(len(seq)):
        	specTotal = (abs(FR[i])**2)
        	specTwo = (abs(FR[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (real))
        # print("\n")
        # print("R -- %s" % (abs(FR)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def eiip_fourier():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        eiip = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
            	eiip.append(0.1335)
            elif nucle == "C":
            	eiip.append(0.1340)
            elif nucle == "A":
            	eiip.append(0.1260)
            else:
            	eiip.append(0.0806)
        Feiip = fft(eiip)
        for i in range(len(seq)):
        	specTotal = (abs(Feiip[i])**2)
        	specTwo = (abs(Feiip[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (eiip))
        # print("\n")
        # print("R -- %s" % (abs(Feiip)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def complex_number():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        complexNumber = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
            	complexNumber.append(1-1j)
            elif nucle == "C":
            	complexNumber.append(-1+1j)
            elif nucle == "A":
            	complexNumber.append(1+1j)
            else:
            	complexNumber.append(-1-1j)
        FcomplexNumber = fft(complexNumber)
        for i in range(len(seq)):
        	specTotal = (abs(FcomplexNumber[i])**2)
        	specTwo = (abs(FcomplexNumber[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (complexNumber))
        # print("\n")
        # print("R -- %s" % (abs(FcomplexNumber)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def atomic_number():
    header_fourier()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        atomicNumber = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
            	atomicNumber.append(66)
            elif nucle == "C":
            	atomicNumber.append(58)
            elif nucle == "A":
            	atomicNumber.append(70)
            else:
            	atomicNumber.append(78)
        FatomicNumber = fft(atomicNumber)
        for i in range(len(seq)):
        	specTotal = (abs(FatomicNumber[i])**2)
        	specTwo = (abs(FatomicNumber[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        feature_extraction_fourier()
        file_record_fourier()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (atomicNumber))
        # print("\n")
        # print("R -- %s" % (abs(FcomplexNumber)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return

#####################################################################################################################

#####################################  EntropyClass #################################################################

def header_entropy(ksize):
    file = open(foutput, 'a')
    file.write("nameseq,")
    for i in range(1, ksize+1):
        file.write("k" + str(i) + ",")
    file.write("label")
    file.write("\n")
    return


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return        
    

def chunks_two(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def file_record(information_entropy):
    file = open(foutput, 'a')
    file.write("%s," % (name_seq))
    for data in information_entropy:
        file.write("%s," % (str(data)))
    file.write(label_dataset)
    file.write("\n")
    print ("Recorded Sequence!!!")
    return
    

def entropy_equation():
    header_entropy(ksize)
    global name_seq, information_entropy
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        information_entropy = []
        for k in range(1, ksize+1):
            probabilities = []
            kmer = {}
            total_windows = (len(seq) - k) + 1 # (L - k + 1)
            for subseq in chunks_two(seq, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print(key)
                # print(value)
                probabilities.append(value/total_windows)
            if e == "Shannon" or e == "shannon":
                entropy_equation = [(p * math.log(p, 2)) for p in probabilities]
                entropy = -(sum(entropy_equation))
                information_entropy.append(entropy)
            else:
                q = 2
                entropy_equation = [(p ** q) for p in probabilities]
                entropy =  (1/(q - 1)) * (1 - sum(entropy_equation))
                information_entropy.append(entropy)
        file_record(information_entropy)
    return
#########################################################################################################

################################### Complex Network  ####################################################

def header_complex(file):
    """
    AB = avg_betweenness
    AEB = avg_edge_betweenness
    AD = avg_degree
    ASSD = assortativity_degree
    MAXD = maxdegree
    MIND = mindegree
    AS = avg_strength
    APL = average_path_length
    ASPL =average_shortest_paths
    TALU = transitivity_avglocal_undirected
    TU = transitivity_undirected
    AEC = avg_eigenvector_centrality
    NE = number_of_edges
    MOT3 = motifs_randesu_no_3 
    MOT4 = motifs_randesu_no_4
    """
    file = open(file, 'a')
    file.write("nameseq,")
    for i in range(1, threshold):
	    file.write("AB.{0},AD.{0},ASSD.{0},MAXD.{0},MIND.{0},AS.{0},APL.{0},TALU.{0},TU.{0},NE.{0},MOT3.{0},MOT4.{0},".format(i))
    file.write("label")
    file.write("\n")
    

def feature_extraction():
    metrics.append(mean(thresholdCN.betweenness(directed=False, weights=None, nobigint=True)))
    metrics.append(mean(thresholdCN.degree()))
    metrics.append(thresholdCN.assortativity_degree(directed=False)) # Returns the assortativity
    metrics.append(max(thresholdCN.degree()))
    metrics.append(min(thresholdCN.degree()))
    metrics.append(np.std(thresholdCN.degree())) # Returns the strength (weighted degree)
    metrics.append(thresholdCN.average_path_length(directed=False, unconn=False)) # Average path length
    metrics.append(thresholdCN.transitivity_avglocal_undirected()) # local transitivity (clustering coefficient) 
    metrics.append(thresholdCN.transitivity_undirected()) # global transitivity (clustering coefficient) 
    metrics.append(mean(cn.ecount())) # Counts the number of edges
    metrics.append(thresholdCN.motifs_randesu_no(size=3))
    metrics.append(thresholdCN.motifs_randesu_no(size=4))
    return


def patterns(seq, win):
    """
    Generate k-mers: subsequences of length k 
    contained in a biological sequence.
    """
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def file_record_complex(foutput):
    """Generates CSV with the extracted features"""
    file = open(foutput, 'a')
    file.write("%s," % (name_seq))
    metrics_preprocessing = np.nan_to_num(metrics)
    for x in metrics_preprocessing:
        # print (x)
        file.write("%s," % (x))
    file.write(label_dataset)
    file.write("\n")
    print ("Recorded Sequence!!!")
    return
    
    
def complex_network(finput,foutput,label,ksize,threshold):
    """Generates complex network"""
    global name_seq, cn, metrics, thresholdCN
    header_complex(foutput) # Header
    for seq_record in SeqIO.parse(finput, "fasta"): # Reading Sequences
        perm = []
        metrics = []
        cn = Graph()
        # print(summary(cn))
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        kmer = []
        for codon in patterns(seq, ksize): # Generates every chosen k pattern
            kmer.append(str(codon))
        # print(kmer)
        for i in range(len(kmer)-1): # Position -1 -- Build the Network
            cn.add_vertices(kmer[i]) if kmer[i] not in perm else kmer[i]
            cn.add_vertices(kmer[i+1]) if kmer[i+1] not in perm else kmer[i+1]
            cn.add_edges([(kmer[i],kmer[i+1])])
            perm.append(kmer[i]) if kmer[i] not in perm else kmer[i]
            perm.append(kmer[i+1]) if kmer[i+1] not in perm else kmer[i+1]
        # print(summary(cn))
        for t in range(1, threshold):
            """
            Extract features using a threshold scheme.
            Eliminates each edge of size < t and extract features again.
            """
            matrix_adj = pd.DataFrame(cn.get_adjacency())
            thresholdCN = np.where(matrix_adj < t, 0, matrix_adj)
            thresholdCN = Graph.Adjacency(thresholdCN.astype(int).tolist(), mode=ADJ_UNDIRECTED)
            # print(t)
            if thresholdCN.ecount() < 1:
                for i in range(t, threshold):
                    for i in range(1,13):
                        metrics.append(0)
                break
            else:
                feature_extraction()
        file_record_complex(foutput)
    # return cn.write_adjacency("cn")
    return
#########################################################################################################

#################################### Ficket Score #######################################################

def look_up_position_prob(value, base, position_para, position_prob, position_weight):

	"""look up positional probability by base and value"""

	if float(value) < 0:
		return None
	for idx, val in enumerate(position_para):
		if float(value) >= val:
			return float(position_prob[base][idx]) * float(position_weight[base])


def look_up_content_prob(value, base, content_para, content_prob, content_weight):

	"""look up content probability by base and value"""

	if float(value) < 0:
		return None
	for idx, val in enumerate(content_para):
		if float(value) >= val:
			return float(content_prob[base][idx]) * float(content_weight[base])


def fickett_value_orf(seq, type_seq):

	"""calculate Fickett value. Input is DNA sequence"""

	position_prob = {
		'A': [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
		'C': [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
		'G': [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
		'T': [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]}

	position_weight = {'A': 0.26, 'C': 0.18, 'G': 0.31, 'T': 0.33}
	position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]

	content_prob = {
		'A': [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
		'C': [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
		'G': [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
		'T': [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]}

	content_weight = {'A': 0.11, 'C': 0.12, 'G': 0.15, 'T': 0.14}
	content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]

	if len(seq) < 2:
		return 0
	fickett_score = 0
	seq = seq.upper()
	total_base = len(seq)
	A_content = float(seq.count('A')) / total_base
	C_content = float(seq.count('C')) / total_base
	G_content = float(seq.count('G')) / total_base
	if type_seq == 1:
		T_content = float(seq.count('T')) / total_base
	else:
		T_content = float(seq.count('U')) / total_base

	phase_0 = [seq[i] for i in range(0, len(seq)) if i % 3 == 0]
	phase_1 = [seq[i] for i in range(0, len(seq)) if i % 3 == 1]
	phase_2 = [seq[i] for i in range(0, len(seq)) if i % 3 == 2]
	
	A_position = max(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) / (min(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) + 1.0)

	C_position = max(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) / (min(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) + 1.0)

	G_position = max(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) / (min(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) + 1.0)

	if type_seq == 1:
		T_position = max(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) / (min(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) + 1.0)
	else:
		T_position = max(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) / (min(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) + 1.0)

	fickett_score += look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)
	
	fickett_score += look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)
	
	return fickett_score


def fickett_value_full_sequence(seq, type_seq):

	"""calculate Fickett from full sequence - CPC2"""

	position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
	content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]

	position_prob = {
		'A': [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
		'C': [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
		'G': [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
		'T': [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24]}
	
	position_weight = {'A': 0.062, 'C': 0.093, 'G': 0.205, 'T': 0.154}
	content_weight = {'A': 0.084, 'C': 0.076, 'G': 0.081, 'T': 0.055}

	content_prob = {
		'A': [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
		'C': [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
		'G': [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
		'T': [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]}

	if len(seq) < 2:
		return 0

	fickett_score = 0
	seq = seq.upper()
	total_base = len(seq)

	phase_0 = seq[::3]
	phase_1 = seq[1::3]
	phase_2 = seq[2::3]

	phase_0_A = phase_0.count('A')
	phase_1_A = phase_1.count('A')
	phase_2_A = phase_2.count('A')
	phase_0_C = phase_0.count('C')
	phase_1_C = phase_1.count('C')
	phase_2_C = phase_2.count('C')
	phase_0_G = phase_0.count('G')
	phase_1_G = phase_1.count('G')
	phase_2_G = phase_2.count('G')
	if type_seq == 1:
		phase_0_T = phase_0.count('T')
		phase_1_T = phase_1.count('T')
		phase_2_T = phase_2.count('T')
	else:
		phase_0_T = phase_0.count('U')
		phase_1_T = phase_1.count('U')
		phase_2_T = phase_2.count('U')

	A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
	C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
	G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
	T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
	A_position = max([phase_0_A, phase_1_A, phase_2_A]) / (min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
	C_position = max([phase_0_C, phase_1_C, phase_2_C]) / (min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
	G_position = max([phase_0_G, phase_1_G, phase_2_G]) / (min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
	T_position = max([phase_0_T, phase_1_T, phase_2_T]) / (min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)

	fickett_score += look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)

	fickett_score += look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)

	return fickett_score


def header_fickett_score(foutput):
	file = open(foutput, 'a')
	file.write('%s,' % 'nameseq')
	file.write('%s,' % 'fickett_score-ORF')
	file.write('%s,' % 'fickett_score-full-sequence')
	file.write('label')
	file.write('\n')
	return


def calculate_sequences(finput, foutput, label_dataset, type_seq):
	header_fickett_score(foutput)
	for seq_record in SeqIO.parse(finput, 'fasta'):
		seq = seq_record.seq
		seq = seq.upper()
		name_seq = seq_record.name
		measure_orf = fickett_value_orf(seq, type_seq)
		measure_full = fickett_value_full_sequence(seq, type_seq)
		file = open(foutput, 'a')
		file.write('%s,' % name_seq)
		file.write('%s,' % str(measure_orf))
		file.write('%s,' % str(measure_full))
		file.write('%s' % str(label_dataset))
		file.write('\n')
		print('Recorded Sequence: %s' % name_seq)
	return

#########################################################################################################

############################################ ORF ########################################################

def des_start_code(codes):
	if codes in ('ATG', 'AUG'):
		return True
	return False


def des_end_code(codes):
	if codes in ('TAA', 'UAA', 'TAG', 'UAG', 'TGA', 'UGA'):
		return True
	return False


def read_by_three(string, offset):
	flag = True
	length = len(string)
	start = end = -1
	i = 0
	result = set()
	while i < length-2:
		codes = string[i:i+3]
		if des_start_code(codes) and flag:
			start = i
			flag = False
		if des_end_code(codes) and not flag:
			end = i + 2
			flag = True
		if (end > start) and (start != -1):
			result.add((start + offset, end + offset))
		i = i + 3
	return result


def get_gc(string):
	gc = ((string.count('G') + string.count('C')) / len(string)) * 100
	return gc


def get_info(string, pos):
	length = pos[1] - pos[0] + 1
	gc = get_gc(string[pos[0]:pos[1]+1])
	return str(pos[0]), str(pos[1]), str(length), str(gc)


def orf(seq):
	result_info = []
	strings = [seq, seq[1:], seq[2:]]
	for index, string in enumerate(strings):
		# print(index)
		# print(string)
		positions = read_by_three(string, index)
		positions = sorted(positions, key=operator.itemgetter(0))
		# print(positions)
		for pos in positions:
			result_info.append(get_info(seq, pos))
	# print(result_info)
	# print(len(result_info))
	return result_info


def run(finput, foutput, label_dataset):
	file = open(foutput, 'a')
	file.write('nameseq,maximum_ORF_length,minimum_ORF_length,std_ORF_length,average_ORF_length,cv_ORF_length,' + ''
			   + 'maximum_GC_content_ORF,minimum_GC_content_ORF,std_GC_content_ORF,' + ''
			   + 'average_GC_content_ORF,cv_GC_content_ORF,label')
	file.write('\n')
	for seq_record in SeqIO.parse(finput, 'fasta'):
		seq = seq_record.seq
		seq = seq.upper()
		name_seq = seq_record.name
		file.write('%s,' % name_seq)
		measures = orf(seq)
		if len(measures) > 0:
			length_orf = []
			gc_mea = []
			for values in measures:
				length_orf.append(int(values[2]))
				gc_mea.append(float(values[3]))
			# print(length_orf)
			# print(gc_mea)
			file.write('%s,' % max(length_orf))
			file.write('%s,' % min(length_orf))
			file.write('%s,' % np.std(length_orf))
			file.write('%s,' % np.mean(length_orf))
			file.write('%s,' % scipy.stats.variation(length_orf))
			file.write('%s,' % max(gc_mea))
			file.write('%s,' % min(gc_mea))
			file.write('%s,' % np.std(gc_mea))
			file.write('%s,' % np.mean(gc_mea))
			file.write('%s,' % scipy.stats.variation(gc_mea))
			file.write('%s' % label_dataset)
			file.write('\n')
		else:
			file.write('0,0,0,0,0,0,0,0,0,0,')
			file.write('%s' % label_dataset)
			file.write('\n')
		print('Recorded Sequence: %s' % name_seq)
	return


#########################################################################################################


if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("########################## Script Feature Extraction: #############################")
    print("##########   Arguments: python3.5 -i input -o output -l label    ##################")
    print("##########               Author: Alvaro Pedroso Queiroz                ############")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-seq', '--seq', help='type of sequence, 1 = DNA, 2 = RNA and 3 = Protein')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    labelDataset = str(args.label)
    
    print("Select the feature extraction techniques: E.g., 1,2,3,4...\n")

    print("1 -  K-mers, k =3")
    print("2 - MonoMonoKGAP, k =1")
    print("3 - MonoDiKGAP, k =1")
    print("4 - DiMonoKGAP, k =1")
    print("5 - Type 2 Pseudo k-tuple nucleotide composition, k=3")
    print("6 - Fourier Transform with Complex Numbers Representation")
    print("7 - Shannon Entropy, k=12")
    print("8 - Complex Networks, k=3")
    print("9 - Ficket Score")
    print("10 - ORF")
    print("\n")

    tec = input()

    if(tec == ''):
        array = ['1','2','3','4','5','6','7','8','9','10']
    else:
        array = tec.split(",")

    diretorio = './'+foutput
    os.makedirs(diretorio)

    ### Configurando K-mer 3 ###
    if '1' in array:
        foutput = diretorio+'/Kmers3.csv'
        ksize = 3
        stepw = 1
        findKmers()

    ##### Configurando MonoMonoKGAP ###
    if '2' in array:
        foutput = diretorio+'/MonoMonoKGAP.csv'
        k = 1
        before = 1
        after = 1
        seq = int(args.seq)
        if seq == 1:
            caracteres = ["A", "C", "G", "T"]
        elif seq == 2:
            caracteres = ["A", "C", "G", "U"]
        elif seq == 3:
            caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        else:
            print('Check parameter: -seq, options: 1 = DNA, 2 = RNA and 3 = Protein')
        kgap()

    ##### Configurando MonoDiKGAP ###
    if '3' in array:
        foutput = diretorio+'/MonoDiKGAP.csv'
        k = 1
        before = 1
        after = 2
        seq = int(args.seq)
        if seq == 1:
            caracteres = ["A", "C", "G", "T"]
        elif seq == 2:
            caracteres = ["A", "C", "G", "U"]
        elif seq == 3:
            caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        else:
            print('Check parameter: -seq, options: 1 = DNA, 2 = RNA and 3 = Protein')
        kgap()
    
    ##### Configurando DiMonoKGAP ###
    if '4' in array:
        foutput = diretorio+'/DiMonoKGAP.csv'
        k = 1
        before = 2
        after = 1
        seq = int(args.seq)
        if seq == 1:
            caracteres = ["A", "C", "G", "T"]
        elif seq == 2:
            caracteres = ["A", "C", "G", "U"]
        elif seq == 3:
            caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        else:
            print('Check parameter: -seq, options: 1 = DNA, 2 = RNA and 3 = Protein')
        kgap()
    
    #### Configurando Type 2 Pseudo k-tuple nucleotide composition ####
    if '5' in array:
        foutput = diretorio+'/PseKNC.csv'
        x = 'propNames.txt'
        t = 2
        k = 3
        j = 1
        w = 0.5
        s = 2

        stepw = 1
        if seq == 1:
            caracteres = ["A", "C", "G", "T"]
        else:
            caracteres = ["A", "C", "G", "U"]
        parameters(x,t,k,j,w,s,finput,foutput,seq)
        
    #### Fourier + Complex ####
    if '6' in array:
        foutput = diretorio+'/FourierComplex.csv'
        label_dataset = str(args.label)
        representation = 6
        
        if representation == 1:
            binary_fourier()
        elif representation == 2:
            zcurve_fourier()
        elif representation == 3:
            real_fourier()
        elif representation == 4:
            integer_fourier()
        elif representation == 5:
            eiip_fourier()
        elif representation == 6:
            complex_number()
        elif representation == 7:
            atomic_number()
        else:
            print("This package does not contain this representation")

    #### Shannon ####
    if '7' in array:
        foutput = diretorio+'/Shannon.csv'
        ksize = 12
        stepw = 1
        e = "Shannon"
        if e == "Shannon" or e == "shannon" or e == "Tsallis" or e == "tsallis":
            entropy_equation()
        
    #### Complex Network ###
    if '8' in array:
        foutput = diretorio+'/ComplexNetwork.csv'
        ksize = 3
        threshold = 10 + 1 # always +1
        complex_network(finput,foutput,label_dataset,ksize,threshold)
        
    ### Ficket Score ###
    if '9' in array:
        foutput = diretorio+'/FicketScore.csv'
        label_dataset = str(args.label)
        type_seq = int(args.seq)
        calculate_sequences(finput, foutput, label_dataset, type_seq)
        
    ### ORF ###
    if '10' in array:
        foutput = diretorio+'/ORF.csv'
        label_dataset = str(args.label)
        run(finput, foutput, label_dataset)

    ### Formando Base de Dados ###
    if '1' in array:
        df1 = pd.read_csv(diretorio+'/Kmers3.csv')
    if '2' in array:
        df2 = pd.read_csv(diretorio+'/MonoMonoKGAP.csv')
    if '3' in array:
        df3 = pd.read_csv(diretorio+'/MonoDiKGAP.csv')
    if '4' in array:
        df4 = pd.read_csv(diretorio+'/DiMonoKGAP.csv')
    if '5' in array:
        df5 = pd.read_csv(diretorio+'/PseKNC.csv')
    if '6' in array:
        df6 = pd.read_csv(diretorio+'/FourierComplex.csv')
    if '7' in array:
        df7 = pd.read_csv(diretorio+'/Shannon.csv')
    if '8' in array:
        df8 = pd.read_csv(diretorio+'/ComplexNetwork.csv')
    if '9' in array:
        df9 = pd.read_csv(diretorio+'/FicketScore.csv')
    if '10' in array:
        df10 = pd.read_csv(diretorio+'/ORF.csv') 
    
    if '1' in array:
        df1 = df1.drop(['class'],axis=1)
    if '2' in array:
        df2 = df2.drop(['nameseq','label'],axis=1)
    if '3' in array:
        df3 = df3.drop(['nameseq','label'],axis=1)
    if '4' in array:
        df4 = df4.drop(['nameseq','label'],axis=1)
    if '5' in array:
        df5 = df5.drop(['nameseq','label'],axis=1)
    if '6' in array:
        df6 = df6.drop(['nameseq','label'],axis=1)
    if '7' in array:
        df7 = df7.drop(['nameseq','label'],axis=1)
    if '8' in array:
        df8 = df8.drop(['nameseq','label'],axis=1)
    if '9' in array:
        df9 = df9.drop(['nameseq','label'],axis=1)
    if '10' in array:
        df10 = df10.drop(['nameseq','label'],axis=1)

    df = pd.DataFrame()

    if '1' in array:
        df = pd.concat([df,df1],axis=1)
    if '2' in array:
        df = pd.concat([df,df2],axis=1)
    if '3' in array:
        df = pd.concat([df,df3],axis=1)
    if '4' in array:
        df = pd.concat([df,df4],axis=1)
    if '5' in array:
        df = pd.concat([df,df5],axis=1)
    if '6' in array:
        df = pd.concat([df,df6],axis=1)
    if '7' in array:
        df = pd.concat([df,df7],axis=1)
    if '8' in array:
        df = pd.concat([df,df8],axis=1)
    if '9' in array:
        df = pd.concat([df,df9],axis=1)
    if '10' in array:
        df = pd.concat([df,df10],axis=1)

    df['class'] = labelDataset
    df.to_csv(diretorio+'/Base_de_dados.csv',index= False, header = True)
    
#############################################################################