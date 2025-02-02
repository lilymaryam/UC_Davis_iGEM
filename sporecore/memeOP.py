#!/usr/bin/python3
import os
import sys
import motiflib
import argparse
import math
import statistics


#dnafile optional, pa pc pg pt 
# setup
parser = argparse.ArgumentParser(
	description='Will test MEME parameters on fasta files of promoters.')
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='Path to jaspar directory')
parser.add_argument('--memepath', required=True, type=str,
	metavar='<str>', help='path to meme software')
parser.add_argument('--dnafile', required=False, type=str,
	metavar='<str>', help='fasta file containing promoter DNA')
parser.add_argument('--markov_order', required=False, type=int, default=0,
	metavar='<int>', help='chosen markov order [%(default)i]')
parser.add_argument('--minpromoterlength', required=False, type=int, 
	default=100,metavar='<int>', help='maximum length of generated promoters\
	 [%(default)i]')
parser.add_argument('--maxpromoterlength', required=False, type=int, 
	default=400,metavar='<int>', help='maximum length of generated promoters\
	 [%(default)i]')
parser.add_argument('--promoterlengthstep', required=False, type=int,
	default=100, metavar='<int>', help='distance between promoter lengths\
	 being tested [%(default)i]')
parser.add_argument('--minnumseq', required=False, type=int, default=5,
	metavar='<int>', help='maximum number of generated promoters \
	[%(default)i]')
parser.add_argument('--maxnumseq', required=False, type=int, default=20,
	metavar='<int>', help='maximum number of generated promoters \
	[%(default)i]')
parser.add_argument('--numseqstep', required=False, type=int, default=5,
	metavar='<int>', help='distance between number of promoter seqs\
	 [%(default)i]')
parser.add_argument('--nummotifs', required=False, type=int, default=1,
	metavar='<int>', help='number of motifs for meme to find [%(default)i]')
parser.add_argument('--numiterations', required=False, type=int, default=1,
	metavar='<int>', help='iterations for condensed data [%(default)i]')
parser.add_argument('--motiffrequency', required=False, type=float, 
	default=0.9, metavar='<float>', help='background probability of A \
	[%(default).3f]')
parser.add_argument('--PA', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of A [%(default).3f]')
parser.add_argument('--PC', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of C [%(default).3f]')
parser.add_argument('--PG', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of G [%(default).3f]')
parser.add_argument('--PT', required=False, type=float, default=0.25,
	metavar='<float>', help='background probability of T [%(default).3f]')
parser.add_argument('--bothstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--condenseddata', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()


assert(math.isclose(arg.PA + arg.PC + arg.PG + arg.PT, 1.0))
assert(arg.motiffrequency <= 1.0)
assert(arg.minpromoterlength <= arg.maxpromoterlength)

def convert_argtovar(minprom, maxprom, promstep, minseq,maxseq,seqstep):
	'''Accepts command line parameter arguments and converts them to ranges of
	parameters to be tested in MEME'''
	promoter = []
	num_seq = []
	for i in range(minprom,maxprom+1,promstep):
		if i != 0:
			promoter.append(i)
	for i in range(minseq,maxseq+1,seqstep):
		if i != 0:
			num_seq.append(i)
	return promoter, num_seq
	

def generate_promoter_mo(dnafile, jasparfile, p, n, o, freq,i): 
	'''Creates a fasta file of promoter sequences with embedded binding motifs.
	Background sequence is generated using markov models of DNA sequence from a
	fasta file input'''
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	cmd = f'python3 mbed.py --jasparfile {jasparfile} \
	--dnafile {dnafile} --markov_order {o} --numseq {n} --seqlen {p} --freq \
	{freq} --bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile

def generate_promoter_bg(jasparfile, p, n, freq, background,i):
	'''Creates a fasta file of promoter sequences with embedded binding motifs.
	Background sequence is generated using a nucleotide distribution frequency
	input from the command line.'''
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	cmd = f'python3 mbed.py --jasparfile {jasparfile} \
	--numseq {n} --seqlen {p} --freq {freq} --PA {background["A"]} --PC \
	{background["C"]} --PG {background["G"]} --PT {background["T"]} \
	--bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile


if not arg.dnafile:
	assert(arg.markov_order == 0)
	background = {'A':arg.PA,'C':arg.PC,'G':arg.PG,'T':arg.PT}
	
	


#initial conditions and parameter range setup
model = ['zoops','oops','anr']
freq = arg.motiffrequency
iterations = arg.numiterations

promoter,numseq = convert_argtovar(arg.minpromoterlength,\
arg.maxpromoterlength,arg.promoterlengthstep,arg.minnumseq,arg.maxnumseq,\
arg.numseqstep)
o = arg.markov_order


if arg.condenseddata:
	print('promoter_file','motif','number_of_sites','jaspar_file_info',\
	'jaspar_info_percentage','e_val','motif_score','motif_score_percentage',\
	'failure_rate','success_rate','false_positive_rate','false_negative',\
	'promoter_size','number_of_sequences','model','markov_order','iterations'\
	,'stdev(avg_eval)','stdev(avg_score)','stdev(avg_pscore)',\
	'stdev(avg_frate)','stdev(avg_srate)','stdev(avg_fprate)','stdev(avg_fn)'\
	,sep=', ')
else:
	print('file','sequence','motif','number_of_sites','meme_position', 'meme_strand'\
	,'meme_width','j_position','j_strand','j_width','p_value','score',\
	'score_percentage','e_value','false_pos','false_neg','positional_distance'\
	,'fail','overlap','overlap percentage', 'promoter length',\
	'number_of_sequences','model', 'markov_order','iterations',sep=', ')
jpwm = motiflib.read_JASPAR(arg.jasparfile)
bits, p_bits = motiflib.score_motifbit(jpwm)
for p in promoter:
	for n in numseq:
		for m in model:
			if arg.condenseddata:
				final = []
			for r in range(iterations):
				if arg.dnafile:
					#print('YEEEEEEEHHAAAWWWW')
					promoter_file = generate_promoter_mo(arg.dnafile, \
					arg.jasparfile, p, n, o, freq,r+1)
					background = motiflib.calcbg_frompromfile(promoter_file)
					#print('mo background', background)
				else:
					#print('NAHHHHHHHHHHHHHHH')
					promoter_file = generate_promoter_bg\
					(arg.jasparfile,p,n,freq,background,r+1)
				#print('background',background,'markov_order',arg.markov_order)
				j_info, numjsites = motiflib.read_testmotif(promoter_file)
				motifs, meme_info, motif_info = motiflib.run_meme\
				(arg.memepath,promoter_file,m,o,arg.nummotifs,\
				arg.bothstrands)
				scores = motiflib.performance_bg(jpwm,motifs,background)
				results, fn = motiflib.get_memedata(promoter_file,\
				meme_info,j_info,scores,motif_info,\
				jpwm,p,n,m,o,iterations)
				false_neg = motiflib.find_falsenegs(fn, j_info,p,n,\
				m,o,promoter_file,jpwm,iterations)
				if not arg.condenseddata:
					for i in range(len(results)):
						print(*results[i],sep = ',')
					for i in range(len(false_neg)):
						print(*false_neg[i],sep = ',')
					print('end')
				else:
					final.append(motiflib.find_condensedstats(results,\
					fn,motif_info,j_info,promoter_file,bits,p_bits,p,n,m,o\
					,false_neg,numjsites,iterations))	
			if arg.condenseddata:
				output = motiflib.avg_condensedstats(final,motif_info,p\
				,n,m,o,r)
				for i in range(len(output)):
					print(*output[i],sep = ',')
	
				


	






