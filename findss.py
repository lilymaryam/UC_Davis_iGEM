#!/usr/bin/env python3
import sys
import gzip
import statistics

#assert(len(sys.argv) == 2)
#f = sys.argv[1]

def find_posdis(line):
	pos5distance = int(line[6])-int(line[4])
	pos3distance = int(line[5])-int(line[7])
	return pos5distance, pos3distance

def find_negdis(line):
	neg5distance = int(line[5])-int(line[7])
	neg3distance = int(line[6])-int(line[4])
	return neg5distance, neg3distance

def return_distances(f):
	neg5distances = []
	neg3distances = []
	pos5distances = []
	pos3distances = []
	fl = None
	if f.endswith('.gz'):
		fl = gzip.open(f,'rt')
	else:
		fl = open(f)
	for line in fl.readlines():
		line = line.split()
		if line[3] == '+':
			pos5distance,pos3distance = find_posdis(line)
			pos5distances.append(pos5distance)
			pos3distances.append(pos3distance)
		if line[3] == '-':
			neg5distance,neg3distance = find_negdis(line)
			neg5distances.append(neg5distance)
			neg3distances.append(neg3distance)
	fl.close()
	return pos5distances, pos3distances,neg5distances, neg3distances
	
def create_histograms(p5,p3,n5,n3):
	p5.sort()
	p3.sort()
	n5.sort()
	n3.sort()
	hists = []
	distances = [p5,p3,n5,n3]
	for d in distances:
		hist = {}
		for i in d:
			if i not in hist:
				hist[i] = 1
			else:
				hist[i] += 1
		hists.append(hist)
	return hists
	
def find_histstats(p5,p3,n5,n3):
	avgs = []
	stddevs = []
	distances = [p5,p3,n5,n3]
	for d in distances:
		avg = statistics.mean(d)
		stddev = statistics.stdev(d)
		avgs.append(avg)
		stddevs.append(stddev)
	return avgs
		
	
	
'''	
if __name__ == '__main__':
	p5 = [0,0,0,1,12,22,11]
	p3 = [0,0,0,1,1,1,4,5]
	n5 = [0,0,10]
	n3 = [0,1,2,3]
	print(find_histstats(p5,p3,n5,n3))
	print(create_histograms(p5,p3,n5,n3))
'''
		
	
	
	
		
#p5,p3,n5,n3 = return_distances(f)
#histograms = create_histograms(p5,p3,n5,n3)
'''
for h in histograms:
	print('distance','count',sep=', ')
	for d in h:
		print(d,h[d],sep=', ')
	#print(h)
	#for i in range(len(h)):
'''		
#print(find_histstats(p5,p3,n5,n3))		
	

	
			


'''
pos5distances = []
pos3distances = []
neg5distances = []
neg3distances = []
histograms = []
distances = [pos5distances, pos3distances, neg5distances, neg3distances]
step = 20
#print('distance','number of promoters',sep=', ')
with gzip.open(sys.argv[1], 'rt') as f:
	for line in f.readlines():
		line = line.strip()
		line = line.split()
		strand = line[3]
		if strand == '+':
			pos5distance, pos3distance = find(pos_di)
			pos5distances.append(pos5distance)
			pos3distances.append(pos3distance)			
		if strand == '-':
			
			neg5distances.append(neg5distance)
			neg3distances.append(neg3distance)			
	pos5distances.sort()
	pos3distances.sort()
	neg5distances.sort()
	neg3distances.sort()
	for d in distances:
		hist = {}
		for dis in d:
			dis = int(dis/step)
			if dis not in hist:
				hist[dis] = 1
			else: hist[dis] += 1
		histograms.append(hist)
for i in range(len(histograms)):
	for j in range(i+1,len(histograms)):
		d = distance(histograms[i],histograms[j])
		#write a function called distance, kl distance
		#throw out 0s, use t test from stats package 
		#print i j d
'''		
	
