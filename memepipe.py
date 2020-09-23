#!/usr/bin/python3
import os
import sys
import motiflib
import argparse
import math
import statistics
import memepipelib

# setup
parser = argparse.ArgumentParser(
	description='Will test MEME parameters on fasta files of promoters.')
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='Path to jaspar directory')
parser.add_argument('--memepath', required=True, type=str,
	metavar='<str>', help='path to meme software')
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
parser.add_argument('--maxmarkov', required=False, type=int, default=1,
	metavar='<int>', help='highest markov order [%(default)i]')
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
parser.add_argument('--negstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--condenseddata', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()


#change to use math.isclose
assert(math.isclose(arg.PA + arg.PC + arg.PG + arg.PT, 1.0))
assert(arg.maxmarkov <= 4)
assert(arg.motiffrequency <= 1)



model = ['zoops','oops','anr']
background = {'A':arg.PA,'C':arg.PC,'G':arg.PG,'T':arg.PT}
freq = arg.motiffrequency
iterations = arg.numiterations


def convert_argtovar(minprom, maxprom, promstep, minseq,maxseq,seqstep,\
maxmarkov):
	promoter = []
	num_seq = []
	markov_order = []
	for i in range(minprom,maxprom+1,promstep):
		if i != 0:
			promoter.append(i)
	for i in range(minseq,maxseq+1,seqstep):
		if i != 0:
			num_seq.append(i)
	for i in range(maxmarkov+1):
		markov_order.append(i)	
	return promoter, num_seq, markov_order


promoter,numseq,markov_order = convert_argtovar(arg.minpromoterlength,\
arg.maxpromoterlength,arg.promoterlengthstep,arg.minnumseq,arg.maxnumseq,\
arg.numseqstep,arg.maxmarkov)

'''
def generate_promoter(jasparfile, p, n, freq, background):
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}.fa' 
	cmd = f'python3 motifapalooza.py --jasparfile {arg.jasparfile} \
	--numseq {n} --seqlen {p} --freq {freq} --PA {background["A"]} --PC \
	{background["C"]} --PG {background["G"]} --PT {background["T"]} \
	--bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile
	
def run_meme(promoterfile,m,o,nummotifs):
	meme = f'{arg.memepath} {promoterfile} -dna -markov_order {o} -mod {m}\
	-nmotifs {nummotifs}'
	if arg.bothstrands:
		meme=meme+' -revcomp'
	#if arg.condenseddata:
	#	for i in range()
	os.system(meme)
	meme_info = motiflib.read_memetxt('meme_out/meme.txt')
	motifs, motif_info = motiflib.memepwm('meme_out/meme.txt')
	return motifs, meme_info, motif_info
	
	
#performance hands back one number 
def performance(motif,motifs,background):
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motiflib.global_motcompare(motif,memepwm,background)
		scores.append(score)
	return scores

	
#memeinfo is a single line of the meme.txt file

def get_memedata(meme_info, j_info, distance_scores, motif_info,jpwm,p,n,m,o):
	result = []
	fp = 0
	fn = []
	fail_c = 0
	for i in range(len(meme_info)):
		seq = meme_info[i][3]
		motif = meme_info[i][0]
		nsites = meme_info[i][1]
		m_strand = meme_info[i][2]
		m_pos = meme_info[i][4]
		p_val = meme_info[i][5]
		for j in range(len(j_info)):
			if seq == j_info[j][0]:
				j_pos = j_info[j][1]
				j_strand = j_info[j][2]
				#numjsites = j_info[j][3]
				break
		for j in range(len(motif_info)):
			if motif == motif_info[j][0]:
				score = distance_scores[j][0] 
				p_score = distance_scores[j][1]
				meme_wid = motif_info[j][1]
				meme_eval = motif_info[j][3]
				break	
		fpos, posdis, fl = motiflib.pos_accuracy(m_pos, j_pos,meme_wid,\
		len(jpwm))
		fp += fpos
		fn.append(seq)
		fail_c += fl
		result.append((seq,motif,nsites,m_strand,j_strand,p_val,meme_eval,\
		score,p_score,posdis, fpos, 0,fl,p,n,m,o))
	return result,fp,fn,fail_c
'''

'''
def find_falsenegs(fn,j_info):
	false_negs = []
	for i in range(len(j_info)):
		if j_info[i][0] not in fn:
			if j_info[i][1] != '':
					false_negs.append(( j_info[i][0],'','','',j_info[i][2],'',\
					'','','',99,0,1,1,p,n,m,o))
	return false_negs

			

def present_info(promoter_file,results, bits,p_bits,scores, fp,fn,fail_c,\
	j_info,motif_info,p,n,m,o):
	present = []
	false_neg = find_falsenegs(fn, j_info)
	if arg.condenseddata:
		for i in range(len(motif_info)):
			mot = motif_info[i][0]
			numsites = int(motif_info[i][2])
			success = numsites-fail_c
			e_val = float(motif_info[i][3])
			dis = scores[i][0]
			dis_p = scores[i][1]
			present = ((promoter_file,mot, numsites,bits,p_bits, e_val, dis,\
			dis_p,fail_c/numsites,success/numsites,fp/numsites,len(false_neg),p,n,m,o))		
	else:
		 for i in range(len(results)):
		 	present.append((results[i]))
		 for i in range(len(false_neg)):
		 	present.append(false_neg[i])
	return present


def get_sequenceinfo(jasparfile,p,n,freq,background,m,o,nummotifs,jpwm):
	seqinfo = []
	promoter_file = generate_promoter(jasparfile,p,n,freq,background)
	j_info = motiflib.read_testmotif(promoter_file)
	motifs, meme_info, motif_info = run_meme(promoter_file,m,o,nummotifs)
	scores = performance(jpwm,motifs,background)
	results, fp, fn, fail_c = memepipe.get_memedata(meme_info,j_info,scores,motif_info,\
	jpwm,p,n,m,o)
	false_negs = memepipelib.find_falsenegs(fn,j_info)
	for i in range(len(results)):
		seqinfo.append((results[i]))
	for i in range(len(false_negs)):
		seqinfo.append((false_negs[i]))
	return seqinfo	


def get_motifinfo(jasparfile,p,n,freq,background,m,o,nummotifs,jpwm,\
	numiterations):
	result = []
	final = []
	for i in range(numiterations):
		promoter_file = generate_promoter(jasparfile,p,n,freq,background)
		j_info = motiflib.read_testmotif(promoter_file)
		motifs, meme_info, motif_info = run_meme(promoter_file,m,o,nummotifs)
		scores = performance(jpwm,motifs,background)
		results, fp, fn, fail_c = memepipe.get_memedata(meme_info,j_info,scores,\
		motif_info,jpwm,p,n,m,o)
		present = memepipelib.present_info(promoter_file,results,bits,p_bits,scores,\
		fp,fn,fail_c,j_info,motif_info,p,n,m,o)
		result.append(present)
	if len(result) > 1:
		avg_nsites = []
		avg_score = []
		avg_eval = []
		avg_score = []
		avg_pscore = []
		avg_frate = []
		avg_srate = []
		avg_fprate = []
		avg_fn = []
		for i in range(len(result)):
			avg_nsites.append(result[i][2])
			avg_eval.append(result[i][5])
			avg_score.append(result[i][6])
			avg_pscore.append(result[i][7])
			avg_frate.append(result[i][8])
			avg_srate.append(result[i][9])
			avg_fprate.append(result[i][10])
			avg_fn.append(result[i][11])
		final.append((result[0][0],result[0][1],statistics.mean(avg_nsites),\
		result[0][3],result[0][4], statistics.mean(avg_eval),\
		statistics.mean(avg_score), statistics.mean(avg_pscore),\
		statistics.mean(avg_frate), statistics.mean(avg_srate),\
		statistics.mean(avg_fprate),statistics.mean(avg_fn),p,n,m,o,\
		numiterations,statistics.stdev(avg_eval),statistics.stdev(avg_score)\
		,statistics.stdev(avg_pscore),statistics.stdev(avg_frate),\
		statistics.stdev(avg_srate),statistics.stdev(avg_fprate)\
		,statistics.stdev(avg_fn)))
		return final
	else:
		return result

'''		




#right thing to do, use /tmp  
#ask for processid to add to test

#python tmp file import 
#ask in command line what they want to name it 
#(seq,motif,nsites,m_strand,j_strand,overlap,overlap_p,p_val,meme_eval,\
#		score,p_score,posdis, fpos, 0,fl,p,n,m,o)

if arg.condenseddata:
	print('promoter_file','motif','number_of_sites','jaspar_file_info',\
	'jaspar_info_percentage','e_val','distance','distance_percentage',\
	'failure_rate','success_rate','false_positive_rate','false_negative',\
	'promoter_size','number_of_sequences','model','markov_order','iterations'\
	,'stdev(avg_eval)','stdev(avg_score)','stdev(avg_pscore)',\
	'stdev(avg_frate)','stdev(avg_srate)','stdev(avg_fprate)','stdev(avg_fn)'\
	,sep=', ')
else:
	#seq, motif,nsites,m_pos, m_strand,meme_wid,j_pos,j_strand,len(jpwm),p_val,score,p_score,meme_eval,fpos,0,posdis,fl,overlap,overlap_p,p,n,m,o
	print('sequence','motif','number_of_sites','memem position', 'meme strand','meme width','j position','j_strand','j width','p value','score','score percentage','e value','false pos','false neg','positional distance','fail','overlap','overlap percentage', 'promoter length'\
	'number of sequences', 'model', 'markov_order',sep=', ')
jpwm = motiflib.read_JASPAR(arg.jasparfile)
bits, p_bits = motiflib.score_motifbit(jpwm)
for p in promoter:
	for n in numseq:
		for m in model:
			for o in markov_order:
				for r in range(iterations):
					print('	ALL INFO !!!!!!!!!!!!!!!!')
					print('r',r)
					promoter_file = memepipelib.generate_promoter_bg(arg.jasparfile,p,n,freq,background,arg.numiterations,arg.condenseddata)
					print('file',promoter_file)
					print('')
					j_info = motiflib.read_testmotif(promoter_file)
					print('j info',j_info)
					print('')
					motifs, meme_info, motif_info = memepipelib.run_meme(arg.memepath,promoter_file,m,o,arg.nummotifs,background)
					print('motifs',motifs)
					print('')
					print('meme info',meme_info)
					print('')
					print('motif info',motif_info)
					print('')
					scores = memepipelib.performance_bg(jpwm,motifs,background)
					print(scores)
					print('')
					results, fn = memepipelib.get_memedata(promoter_file,meme_info,j_info,scores,motif_info,\
					jpwm,p,n,m,o)
					print('results',results)
					print('')
					print('fn',fn)
					print('')
					false_neg = memepipelib.find_falsenegs(fn, j_info,p,n,m,o,promoter_file)
					print('false neg', false_neg)
					#print(memepipelib.get_memedata())
					if not arg.condenseddata:
						for i in range(len(results)):
							print(results[i])
						for i in range(len(false_neg)):
							print(false_neg[i])
						#present = memepipelib.present_info(promoter_file, results,bits,p_bits,\
						#scores, fp,fn,fail_c,j_info,motif_info,p,n,m,o,arg.condenseddata)
						#present = memepipelib.get_sequenceinfo_bg(arg.memepath,arg.jasparfile,p,n,freq,\
						#background,m,o,arg.nummotifs,jpwm,arg.numiterations,arg.condenseddata)
					else:
						#present = memepipelib.present_info(promoter_file,results,bits,p_bits, scores, fp,fn\
						#,fail_c,j_info,motif_info,p,n,m,o,arg.condenseddata)
						#result = memepipelib.get_variables_bg(arg.memepath,arg.jasparfile,p,n,freq,background,m,o,arg.nummotifs,jpwm,arg.numiterations,arg.condenseddata,arg.bothstrands,bits,p_bits)
						#present = memepipelib.get_motifinfo_bg(arg.memepath,arg.jasparfile,p,n,freq,background\
						#,m,o,arg.nummotifs,jpwm,arg.numiterations,arg.bothstrands,bits,p_bits,arg.condenseddata)
						final = memepipelib.find_condensedstats(results,fn,motif_info,j_info,promoter_file,bits,p_bits,p,n,m,o)	
						for i in range(len(final)):
							print(final[i],sep=',')
					'''
					for i in range(len(present)):
						row = str(present[i])[1:-1]
						print(row)
					if not arg.condenseddata:
						print('end')
					'''
				
				
				
#	promoter_file,results, bits,p_bits,scores, fp,fn,fail_c,\
#	j_info,motif_info,p,n,m,o,condenseddata		
				


	






