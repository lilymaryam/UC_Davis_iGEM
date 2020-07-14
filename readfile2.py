import sys
import fileinput
import os
import statistics

"""
Steps:
1. Input a directory into the command line
2. Read the first file in the directory
3. Output some stuff from that file
4. Store that stuff
5. Move on to the next file
6. Once the program has seen all the files, print out the stuff
"""

folder = sys.argv[1]
directory = os.listdir(folder)
#later make a thing where it'll do all the subdirectories (one for each genome)
#in a parent directory (ex: Promoters)
genes = []
hist = {}
for file in directory:
#need to get the program to open the file; right now it's just looking at the file name
	f = open(f'{folder}/{file}', 'r')
	genecount = 0
	for line in f:
			if line.startswith(">"):
				genecount += 1
	genes.append(genecount)
	if genecount not in hist:
		hist[genecount] = 1
	else: hist[genecount] += 1
print(statistics.mean(genes))
for i in range(1+max(hist.keys())):
	if i in hist:
		print(i, '#' * hist[i])
	else: print(i)
		


"""
	data = []
	for line in file:
		if line[0] == '#': continue # skip over comments
		if line.startswith('#'): continue # same as above
		line = line.rstrip() # remove newline (return character), often useful
		data.append(str(line)) # store the data
print(data)	

	genecount = 0
	for i in data:
			if ">" in i:
				genecount += 1
	print('GC:', genecount)



data = []
for line in fileinput.input():
	#if line[0] == '#': continue # skip over comments
	if line.startswith('#'): continue # same as above
	line = line.rstrip() # remove newline (return character), often useful
	data.append(str(line)) # store the data


genecount = 0
for i in data:
		if ">" in i:
			genecount += 1
print('GC:', genecount)





"""
