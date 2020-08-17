#!/usr/bin/python3
import os
import sys
import motiflib
import readmeme
import argparse
import math

# setup
parser = argparse.ArgumentParser(
	description='Will test MEME parameters on fasta files of promoters.')
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='Path to jaspar directory')
parser.add_argument('--memepath', required=True, type=str,
	metavar='<str>', help='path to meme software')
parser.add_argument('--maxpromoterlength', required=False, type=int, default=400,
	metavar='<int>', help='maximum length of generated promoters [%(default)i]')
parser.add_argument('--promoterlengthstep', required=False, type=int, default=100,
	metavar='<int>', help='distance between promoter lengths being tested [%(default)i]')
parser.add_argument('--maxnumseq', required=False, type=int, default=20,
	metavar='<int>', help='maximum number of generated promoters [%(default)i]')
parser.add_argument('--numseqstep', required=False, type=int, default=5,
	metavar='<int>', help='distance between number of promoter seqs [%(default)i]')
parser.add_argument('--nummotifs', required=False, type=int, default=1,
	metavar='<int>', help='number of motifs for meme to find [%(default)i]')
parser.add_argument('--maxmarkov', required=False, type=int, default=1,
	metavar='<int>', help='highest markov order [%(default)i]')
parser.add_argument('--motiffrequency', required=False, type=float, default=0.9,
	metavar='<float>', help='background probability of A [%(default).3f]')
parser.add_argument('--PA', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of A [%(default).3f]')
parser.add_argument('--PC', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of C [%(default).3f]')
parser.add_argument('--PG', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of G [%(default).3f]')
parser.add_argument('--PT', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of T [%(default).3f]')
#parser.add_argument('--rfloat', required=True, type=float,
#	metavar='<float>', help='required floating point argument')
# optional arguments with default parameters
#parser.add_argument('--dstr', required=False, type=str, default='hello',
#	metavar='<str>', help='optional string argument [%(default)s]')
#parser.add_argument('--dint', required=False, type=int, default=1,
#	metavar='<int>', help='optional integer argument [%(default)i]')
#parser.add_argument('--dfloat', required=False, type=float, default=3.14,
#	metavar='<float>', help='optional floating point argument [%(default)f]')
# switches
parser.add_argument('--bothstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--negstrands', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()


#change to use math.isclose
assert(math.isclose(arg.PA + arg.PC + arg.PG + arg.PT, 1.0))
assert(arg.maxmarkov <= 4)
assert(arg.motiffrequency <= 1)



model = ['zoops', 'oops', 'anr']
background = {'A':arg.PA,'C':arg.PC,'G':arg.PG,'T':arg.PT}
freq = arg.motiffrequency

def convert_argtovar(maxpromlength, promoterstep, maxnumseq,numseqstep, maxmarkov):
	promoter = []
	num_seq = []
	markov_order = []
	for i in range(0,maxpromlength+1,promoterstep):
		if i != 0:
			promoter.append(i)
	for i in range(0,maxnumseq+1, numseqstep):
		if i != 0:
			num_seq.append(i)
	for i in range(maxmarkov+1):
		markov_order.append(i)	
	return promoter, num_seq, markov_order

promoter,numseq,markov_order = convert_argtovar(arg.maxpromoterlength, \
arg.promoterlengthstep,arg.maxnumseq,arg.numseqstep,arg.maxmarkov)

def generate_promoter(jasparfile, p, n, freq, background):
	tmpfile = f'/tmp/testmotif{os.getpid()}.fa' 
	cmd = f'python3 noahpalooza.py --jasparfile {arg.jasparfile} \
	--numseq {n} --seqlen {p} --freq {freq} --PA {background["A"]} --PC {background["C"]}\
	--PG {background["G"]} --PT {background["T"]} --bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile
	
def run_meme(promoterfile,m,o,nummotifs):
	meme = f'{arg.memepath} {promoterfile} -dna -markov_order {o} -mod {m} -nmotifs \
	{nummotifs} -revcomp'
	os.system(meme)
	meme_info = readmeme.read_memetxt('meme_out/meme.txt')
	motifs = readmeme.memepwm('meme_out/meme.txt')
	return motifs
	
#performance hands back one number 
def performance(motif,motifs): #(p,n,m,o,)
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motiflib.global_motcompare(motif,memepwm,background)
		scores.append(score)
	return scores
	
	
	
 




#is this ok?


#right thing to do, use /tmp 
#ask operating system for a temporary file , use temp file library practice on a separate script
#turn testmotif.fa into a temp file  
#ask for processid to add to test

#python tmp file import 
#ask in command line what they want to name it 



print('jasparfile','sequence','promoter length','number of sequences','jaspar \
position','jaspar strand','jaspar width','jaspar information','jaspar \
information percentage','markov model','meme model','meme motifid','motif e \
value','meme width','meme position','meme strand','p value','motif distance \
score','motif score percentage','positional distance','false negative','false\
 positive','fail',sep=', ')
#loops through files(should they be listed as strings or indexed thru as ints?)
#refernce arg directly 
#delete first loop
motif = motiflib.read_JASPAR(arg.jasparfile)
for p in promoter:
	for n in numseq:
		promoter_file = generate_promoter(arg.jasparfile,p,n,freq,background)
		for m in model:
			for o in markov_order:
				motifs = run_meme(promoter_file,m,o,arg.nummotifs)
				perf = performance(motif,motifs)
				print(p,n,m,o,perf)
				
				#print(promoter_file)
				#print(p,n,m,o, performance(promoter_file,m,o))
				
'''
	for p in promoter:
		for n in numseq:
			
			#put everything below here in 3ish functions
			print(p,n,performance(p,n,motif))
			#note that motifapalooza must be found in the same directory 
			
			#uses program noahpalooza (name can be changed)#
			#able to generate different files depending on cmd args
			#is this the best way to do it ?
			if arg.bothstrands:
			#contruct comannd append optional argument to end 
			#don't need to test this, do bothstrands all the time 
			
			#use this to test noahpalooza
				cmd = f'python3 noahpalooza.py --jasparfile {filepath} \
			 	--numseq {n} --seqlen {p} --PA {arg.PA} --PC {arg.PC} --PG\
			 	 {arg.PG} --PT {arg.PT} --bothstrands > testmotif{f}.fa '
			if arg.negstrands:
				cmd = f'python3 noahpalooza.py --jasparfile {filepath} \
			 	--numseq {n} --seqlen {p} --PA {arg.PA} --PC {arg.PC} --PG \
			 	{arg.PG} --PT {arg.PT} --negstrands > testmotif{f}.fa '
			else:
			 	cmd = f'python3 noahpalooza.py --jasparfile {filepath} \
			 	--numseq {n} --seqlen {p} --PA {arg.PA} --PC {arg.PC} --PG \
			 	{arg.PG} --PT {arg.PT}  > testmotif{f}.fa '
			os.system(cmd)
			
			#iterates through the three models of meme (probably hardcoded)
			for m in model:
				#loops through markov models (hardcoded,don't need more than 2)
				
				for o in markovorder:
						#meme is run in this loop and data is parsed
						#is there a better way to code optional arguments
						#(what is the method used in memerunner?)
						#should i add more parameters that will always be there?
						#modify so no repeat
						if arg.bothstrands or arg.negstrands:
							meme = f'{arg.memepath} testmotif{f}.fa -dna \
							-markov_order {o} -mod {m} -nmotifs \
							{arg.nummotifs} -revcomp'
						else: 
							meme = f'{arg.memepath} testmotif{f}.fa -dna \
							-markov_order {o} -mod {m} -nmotifs \
							{arg.nummotifs}'
						os.system(meme)
						
						
						#write a function to do all this 
						#site info
						
						#read_memetxt reads meme.txt and outputs every memesite
						meme_info = readmeme.read_memetxt('meme_out/meme.txt')
						#read_testmotif accepts file generated from 
						#motifapalooza (fasta form) and reads the > information 
						
						j_info = readmeme.read_testmotif(f'testmotif{f}.fa')
						#read_JASPAR takes a jaspar file and converts info into
						#a pwm in the form of an array of dictionaries
						
						jasparpwm = motiflib.read_JASPAR(filepath)
						#memepwm takes the meme.txt file and outputs two arrays
						#one contains a list of the motifs in pwm form
						#the other contains info regarding the motifs 
						#(is this a good way to do it ?)
						
						mememotifs, motif_info = \
						readmeme.memepwm('meme_out/meme.txt')
						#score_motifbit reads in a pwm in array/dictionary form
						#and outputs a tuple containing information of pwm 
						#and percentage of total possible info for that pwm* 
						#*(is that relevant?)
						
						bits = motiflib.score_motifbit(jasparpwm) 
						#distance scores will collect scores between each motif 
						#and the jaspar file in order to be referenced later 
						#(necessary?)
						
						distance_scores = []
						#fn collects all sequences found in meme.txt so that 
						#skipped sequences can be analyzed for false negatives
						#(necessary?)
						fn = []
						#this loop is explicitly for creating distance scores
						#(only needs to be done once)
						# is this the best place/best method?
						
						for i in range(0,len(mememotifs)):
							memepwm = mememotifs[i]
							#can distance scores contain motif_info?
							pwmdistance = motiflib.global_motcompare(jasparpwm\
							, memepwm,background)
							distance_scores.append(pwmdistance)
						#this loop is reading through all sites pulled from 
						#meme.txt
						
						
						#data will be printed as it is read
						for i in range(0, len(meme_info)):
							#better not to make variables?
							seq = meme_info[i][2]
							#used to catch false negatives
							fn.append(seq)
							motif = meme_info[i][0]
							m_strand = meme_info[i][1]
							m_pos = meme_info[i][3]
							p_val = meme_info[i][4]
							#this loop iterates through j info to find seq\
							#info related to location of jaspar motifs
							#must be done after sequences is identified
							#better way to do it?
							for j in range(len(j_info)):
								if seq == j_info[j][0]:
									j_pos = j_info[j][1]
									j_strand = j_info[j][2]
									break
							#this loop iterates through motifinfo to match \
							#motif to information
							#must be done after motif is identified
							#combine distance score? or more computiing power
							#better way to do it ?
							for j in range(len(motif_info)):
								if motif == motif_info[j][0]:
									distance_score = distance_scores[j]
									meme_wid = motif_info[j][1]
									meme_eval = motif_info[j][2]
									break	
							#outputs information about failure, false positive\
							#and distance between start points of motifs
							#doesn't run very well need advice
							positioninfo = motiflib.pos_accuracy(m_pos, j_pos,\
							len(memepwm), len(jasparpwm))
							#data display	
							#prints each iteration of meme as a single row of \
							#data with all parameters listed
							print(f,seq,p,n,j_pos,j_strand,len(jasparpwm),\
							bits[0],bits[1],o,m,motif,meme_eval,meme_wid, \
							m_pos,m_strand,p_val,distance_score[0],\
							distance_score[1], positioninfo[4],\
							positioninfo[2],positioninfo[3],positioninfo[5]\
							,sep=', ')
						#is this needed/done well?
						#catches false negatives
						for i in range(len(j_info)):
							if j_info[i][0] not in fn:
								if j_info[i][1] != '':
									print(f, j_info[i][0], '', '',j_info[i][1]\
									,j_info[i][2],len(jasparpwm),bits[0],\
									bits[1],'','', '', '', '', '', '', '', '',\
									'', 99,1,0,1,sep=', ')
'''
		

				
						






