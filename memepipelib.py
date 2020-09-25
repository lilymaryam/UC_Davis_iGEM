import motiflib
import os
import statistics

#need to fix this
def convert_argtovar(minprom, maxprom, promstep, minseq,maxseq,seqstep):
	promoter = []
	num_seq = []
	for i in range(minprom,maxprom+1,promstep):
		if i != 0:
			promoter.append(i)
	for i in range(minseq,maxseq+1,seqstep):
		if i != 0:
			num_seq.append(i)
	return promoter, num_seq
	

#fix generate promter to create files for condenseddata
def generate_promoter_mo(dnafile, jasparfile, p, n, o, freq,i): 
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	cmd = f'python3 markov_motifapalooza.py --jasparfile {jasparfile} \
	--dnafile {dnafile} --markov_order {o} --numseq {n} --seqlen {p} --freq \
	{freq} > {tmpfile} '
	os.system(cmd)
	return tmpfile

def generate_promoter_bg(jasparfile, p, n, freq, background,i,condenseddata):
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	#tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}.fa' 
	cmd = f'python3 motifapalooza.py --jasparfile {jasparfile} \
	--numseq {n} --seqlen {p} --freq {freq} --PA {background["A"]} --PC \
	{background["C"]} --PG {background["G"]} --PT {background["T"]} \
	--bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile
	
def calcbg_fromdatafile(dnafile):
	A = 0
	C = 0
	G = 0
	T = 0
	with open(dnafile) as df:
		for line in df.readlines():
			if not line.startswith('>'): 
				line = line.strip()
				for i in range(len(line)):
					if line[i] == 'A': A += 1 
					if line[i] == 'C': C += 1 
					if line[i] == 'G': G += 1 
					if line[i] == 'T': T += 1 
	total = A+C+G+T
	#print(total)
	A = A/total
	C = C/total
	G = G/total
	T = T/total 
	background = {'A':A, 'C':C, 'G': G, 'T':T}
	#print(background)
	return background	
	
def run_meme(memepath,promoterfile,m,o,nummotifs,bothstrands):
	meme = f'{memepath} {promoterfile} -dna -markov_order {o} -mod {m}\
	-nmotifs {nummotifs}'
	if bothstrands:
		meme=meme+' -revcomp'
	#add palindrome and other options 
	#if arg.condenseddata:
	#	for i in range()
	os.system(meme)
	meme_info = motiflib.read_memetxt('meme_out/meme.txt')
	motifs, motif_info = motiflib.memepwm('meme_out/meme.txt')
	return motifs, meme_info, motif_info
	
#can combine these now	
#performance hands back one number 
def performance_bg(motif,motifs,background):
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motiflib.global_motcompare(motif,memepwm,background)
		scores.append(score)
	return scores
	
def performance_mo(motif,motifs):
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motiflib.local_motcompare(motif,memepwm)
		scores.append(score)
	return scores

def get_memedata(promoter_file, meme_info, j_info, distance_scores, motif_info\
,jpwm,p,n,m,o,iteration):
	result = []
	fp = 0
	fn = []
	#print('scores',distance_scores)
	#print('motif_info',motif_info)
	#print('meme info',meme_info)
	for i in range(len(meme_info)):
		seq = meme_info[i][3]
		#print('seq',seq)
		motif = meme_info[i][0]
		#print('motif',motif)
		nsites = meme_info[i][1]
		#print('nsites',nsites)
		m_strand = meme_info[i][2]
		#print('m strand',m_strand)
		m_pos = meme_info[i][4]
		#print('mpos',m_pos)
		p_val = meme_info[i][5]
		#print('meme info',meme_info)
		for j in range(len(j_info)):
			if seq == j_info[j][0]:
				j_pos = j_info[j][1][0]
				#print('j pos',j_pos)
				if j_pos != '':
					j_strand = j_info[j][1][1]
				else:
					j_strand = ''
				#print('j strand', j_strand) 
				break
		for j in range(len(motif_info)):
			if motif == motif_info[j][0]:
				score = distance_scores[j][0] 
				#print('score',score)
				p_score = distance_scores[j][1]
				#print('p_score',p_score)
				meme_wid = motif_info[j][1]
				#print('meme wid',meme_wid)
				meme_eval = motif_info[j][3]
				#print('meme eval', meme_eval)
				break	
		fpos,posdis,fl,overlap,overlap_p = motiflib.pos_accuracy(m_pos, j_pos,\
		meme_wid,\
		len(jpwm))
		#print('fpos',fpos)
		#print('posdis',posdis)
		#print('fl',fl)
		#print('overlap',overlap)
		#print('overlap_p', overlap_p)
		#print('fpos',fpos,'posdis',posdis,'fl',fl,'overlap',overlap_p)
		#print('seq','motif','nsites','mpos','m strand','meme wid','jpos'\
		#,'j_strand', 'j wid', 'pval', 'score','p score','eval', 'f pos', \
		#'f neg', 'posdis', 'fl', 'overlap', 'overlap_p', 'p len', 'numseq',\
		# 'model','order')
		result.append((promoter_file, seq, motif,nsites,m_pos, m_strand,\
		meme_wid,j_pos,j_strand,len(jpwm),p_val,score,p_score,meme_eval,fpos\
		,0,posdis,fl,overlap,overlap_p,p,n,m,o,iteration))
		#fp +=fpos
		fn.append(seq)
		#result.append((seq,motif,nsites,m_strand,j_strand,overlap,overlap_p\
		#,p_val,meme_eval,\
		#score,p_score,posdis, fpos, 0,fl,p,n,m,o))
	#print(fn)
	#print('getmemedata','result',result)
	#print('get meme data','fn',fn)
	#print('seq','motif','nsites','mpos','m strand','meme wid','jpos',\
	#'j strand', 'j wid', 'pval', 'score','p score','eval', 'f pos', 'f neg',\
	#'posdis','fl', 'overlap', 'overlap_p', 'p len', 'numseq', 'model','order')
	#for i in range(len(result)):
	#	print('result',i,result[i])
	return result,fn


#need to fix false negs everything 
def find_falsenegs(fn,j_info,p,n,m,o,promoter_file,jpwm,iteration):
	#print('j_info',j_info)
	#print('fn',fn)
	false_negs = []
	for i in range(len(j_info)):
		#print('j_info[i][0]',j_info[i][0])
		#print(fn[i])
		if j_info[i][0] not in fn:
			if j_info[i][1][0] != '':
				#print('is it ""?',j_info[i][1][0])
				false_negs.append(( promoter_file,j_info[i][0],'','','','','',\
				j_info[i][1][0],j_info[i][1][1],len(jpwm),'','','','',0,1,99,1\
				,'','',p,n,m,o,iteration))
	return false_negs	


#print('seq','motif','nsites','mpos','m strand','meme wid','jpos','j strand', \
#'j wid', 'pval', 'score','p score','eval', 'f pos', 'f neg', 'posdis', 'fl', \
#'overlap', 'overlap_p', 'p len', 'numseq', 'model','order')
#make sure this is correct
def find_condensedstats(result,fn,motif_info,j_info,promoter_file,bits,p_bits,\
p,n,m,o,false_neg,numjsites,iteration):
	condensed_stats = []
	#false_neg = find_falsenegs(fn, j_info,p,n,m,o,promoter_file)
	#print('motif_info', motif_info)
	#print('result',result)
	#print('numjsites',numjsites)
	#print('len fn', len(false_neg))
	for i in range(len(motif_info)):
		fl_ct = 0
		fp_ct = 0
		#ovl_sum = []
		#ovlp_sum = []
		#bits,p_bits
		#for k in range(len(fn)):
		#	fn_ct = 0
			#print(fn[k][1])
		#	if fn[k][1] == motif_info[i][0]:
				#print(fn[k][1])
		#		fn_ct += 1
		#print(motif_info[i][0],'fn',fn_ct)
		for j in range(len(result)):
			if result[j][2] == motif_info[i][0]:
				#print(result[j][1], motif_info[i][0],result[j][16],\
				#result[j][13])
				nsites = int(result[j][3])
				#print('nsites',nsites)
				fl_ct += (result[j][17])
				#print('fl_ct',fl_ct)
				fp_ct += (result[j][14])
				#print('fp_ct',fp_ct)
				score = (result[j][11])
				score_p = (result[j][12])
				evalue = (result[j][13])
				#ovl_sum.append(result[j][17])
				#ovlp_sum.append(result[j][18])
				#print('eval', result[j][12])
				#false_n = len(fn)
			#print('no!!!!!')
		#print('THIS IS CONDENSED DATA')
		#print('nsites',nsites)
		#print(motif_info[i][0])
		#print('fn',false_n,false_n/nsites)
		#print('fl_ct',fl_ct)
		#print('fp_ct',fp_ct)
		#print('score',score)
		#print('score_p',score_p) 
		#print('evals', evalue)
		#print('ovl_sum',ovl_sum)
		#print('ovlp_sum',ovlp_sum)
		condensed_stats.append((promoter_file,motif_info[i][0],bits,p_bits,\
		evalue,score,score_p, nsites,fl_ct/nsites,(nsites-fl_ct)/nsites,\
		fp_ct/nsites,len(false_neg)/numjsites,p,n,m,o,iteration))
	#print('file','motif','info','p_info','evalue','score','p score','nsites'\
	#,'failrate','suces rate','false p rate')
	#for i in range(len(condensed_stats)):
	#	print('condensed_stats',i, condensed_stats[i])
	#print('file','motif','info','p_info','evalue','score','p score','nsites',\
	#'failrate','suces rate','false p rate')
	#print('consensed_stats',condensed_stats)
	return condensed_stats	


#make sure this is good, make uncondensed have file name
def present_info(promoter_file,results, bits,p_bits,scores,fn,\
	j_info,motif_info,p,n,m,o,condenseddata):
	present = []
	false_neg = find_falsenegs(fn, j_info,p,n,m,o,promoter_file)
	if condenseddata:
		present = (find_condensedstats(results,fn,motif_info,j_info,\
		promoter_file,bits,p_bits,p,n,m,o))
		
		#for i in range(len(motif_info)):	
		#	mot = motif_info[i][0]
		#	numsites = int(motif_info[i][2])
		#	print('numsites',numsites)
		#	#print('fail_c', fail_c)
		#	#success = numsites-fail_c
		#	#print('success',success)
		#	e_val = float(motif_info[i][3])
		#	dis = scores[i][0]
		#	dis_p = scores[i][1]
		#	present.append((promoter_file,mot,numsites,bits,p_bits,e_val,dis,\
		#	dis_p,fail_c/numsites,success/numsites,fp/numsites,len(false_neg),\
		#p,n,m,o))		
	else:
		 for i in range(len(results)):
		 	present.append((results[i]))
		 for i in range(len(false_neg)):
		 	present.append(false_neg[i])
	return present
	
'''
#def get_variables_mo()

def get_variables_bg(memepath,jasparfile,p,n,freq,background,m,o,nummotifs,\
jpwm,i,condenseddata):
	promoter_file = generate_promoter_bg(jasparfile,p,n,freq,background,i,\
	condenseddata)
	j_info = motiflib.read_testmotif(promoter_file)
	motifs, meme_info, motif_info = run_meme(memepath,promoter_file,m,o,\
	nummotifs,background)
	scores = performance_bg(jpwm,motifs,background)
	results, fn = get_memedata(promoter_file, meme_info,j_info,scores,\
	motif_info,jpwm,p,n,m,o)
	false_negs = find_falsenegs(fn,j_info,p,n,m,o,promoter_file)
	return results, false_negs, promoter_file
	

def get_sequenceinfo_bg(memepath,jasparfile,p,n,freq,background,m,o,nummotifs,\
jpwm,numiterations,condenseddata):
	seqinfo = []
	#for i in range(numiterations):
	results, false_negs, promoter_file = get_variables_bg(memepath,jasparfile,\
	p,n,freq,background,m,o,nummotifs,jpwm,numiterations,condenseddata)
	
	#promoter_file = generate_promoter_bg(jasparfile,p,n,freq,background,i,\
	condenseddata)
	#j_info = motiflib.read_testmotif(promoter_file)
	#motifs, meme_info, motif_info = run_meme(memepath,promoter_file,m,o,\
	nummotifs,background)
	#scores = performance_bg(jpwm,motifs,background)
	#results, fn = get_memedata(meme_info,j_info,scores,motif_info,\
	#jpwm,p,n,m,o)
	#false_negs = find_falsenegs(fn,j_info,p,n,m,o)

	for i in range(len(results)):
		seqinfo.append((promoter_file,results[i]))
	for i in range(len(false_negs)):
		seqinfo.append((false_negs[i]))
	print(seqinfo)
	return seqinfo	

#can these be condensed

def get_sequenceinfo_mo(memepath,dnafile,jasparfile,p,n,freq,m,o,nummotifs,\
jpwm,bothstrands):
	seqinfo = []
	promoter_file = generate_promoter_mo(dnafile,jasparfile,p,n,o,freq)
	j_info = motiflib.read_testmotif(promoter_file)
	motifs, meme_info, motif_info = run_meme(memepath,promoter_file,m,o,\
	nummotifs,bothstrands)
	scores = performance_mo(jpwm,motifs)
	results, fp, fn, fail_c = get_memedata(meme_info,j_info,scores,motif_info,\
	jpwm,p,n,m,o)
	false_negs = find_falsenegs(fn,j_info,p,n,m,o)
	seqinfo.append(promoter_file)
	for i in range(len(results)):
		seqinfo.append((results[i]))
	for i in range(len(false_negs)):
		seqinfo.append((false_negs[i]))
	return seqinfo	
	


def get_motifinfo_bg(memepath,jasparfile,p,n,freq,background,m,o,nummotifs,\
jpwm,
	numiterations,bothstrands,bits,p_bits,condenseddata):
	result = []
	final = []
	for i in range(numiterations):
		promoter_file = generate_promoter_bg(jasparfile,p,n,freq,background,\
		numiterations,condenseddata)
		j_info = motiflib.read_testmotif(promoter_file)
		motifs, meme_info, motif_info = run_meme(memepath,promoter_file,m,o,\
		nummotifs,bothstrands)
		scores = performance_bg(jpwm,motifs,background)
		results, fn = get_memedata(promoter_file,meme_info,j_info,scores,\
		motif_info,jpwm,p,n,m,o)
		present = present_info(promoter_file,results,bits,p_bits,scores,fn,\
		j_info,motif_info,p,n,m,o,condenseddata)
		#print('present',present)
		result.append(present)
		print('motif_info', i, motif_info)
	#print('result',result)
	if numiterations > 1:
		for i in range(len(motif_info)):
			print('motif_info',motif_info)
			mot = motif_info[i][0]
			#final_info = []
			avg_nsites = []
			avg_score = []
			avg_eval = []
			avg_score = []
			avg_pscore = []
			avg_frate = []
			avg_srate = []
			avg_fprate = []
			avg_fn = []
			for j in range(len(result)):
				for k in range(len(result[j])):
					#print('MOT', mot,'result mot',result[i][j][1])
					if  mot == result[j][k][1]:
						#print('MOT', mot,'result mot',result[j][k][1])
						avg_nsites.append(result[j][k][7])
						avg_eval.append(result[j][k][4])
						avg_score.append(result[j][k][5])
						avg_pscore.append(result[j][k][6])
						avg_frate.append(result[j][k][8])
						avg_srate.append(result[j][k][9])
						avg_fprate.append(result[j][k][10])
						#avg_fn.append(result[j][k][11])	
			#print('FIND INFO HERE!!!!!!!!!!!!!!!!!!!')
			#print(mot)
			#print('avg_nsites',avg_nsites)
			#print('avg_eval',avg_eval)
			#print('avg_score',avg_score)
			#print('avg_pscore', avg_pscore)
			#print('avg_frate',avg_frate)
			#print('avg_srate',avg_srate)
			#print('avg_fprate',avg_fprate)
			#print('avg_fn',avg_fn)
			#'file','motif','info','p_info','evalue','score','p score',\
			'nsites','failrate','suces rate','false p rate
			final.append((result[j][k][0],mot,statistics.mean(avg_nsites),\
			result[j][k][3],result[j][k][4], statistics.mean(avg_eval),\
			statistics.mean(avg_score), statistics.mean(avg_pscore),\
			statistics.mean(avg_frate), statistics.mean(avg_srate),\
			statistics.mean(avg_fprate),statistics.mean(avg_fn),p,n,m,o,\
			numiterations,statistics.stdev(avg_eval),\
			statistics.stdev(avg_score),statistics.stdev(avg_pscore),\
			statistics.stdev(avg_frate),statistics.stdev(avg_srate),\
			statistics.stdev(avg_fprate),statistics.stdev(avg_fn)))\
		#print('final',final)
		return final
	else:
		return result
		




def get_motifinfo_mo(memepath,dnafile,jasparfile,p,n,freq,m,o,nummotifs,jpwm,\
numiterations,bothstrands,condenseddata,bits,p_bits):
	result = []
	final = []
	for i in range(numiterations):
		promoter_file = generate_promoter_mo(dnafile,jasparfile,p,n,o,freq,\
		numiterations,condenseddata)
		j_info = motiflib.read_testmotif(promoter_file)
		motifs, meme_info, motif_info = run_meme(memepath,promoter_file,m,o,\
		nummotifs,bothstrands)
		scores = performance_mo(jpwm,motifs)
		results, fn = get_memedata(promoter_file,meme_info,j_info,scores,\
		motif_info,jpwm,p,n,m,o)
		present = present_info(promoter_file,results,bits,p_bits,scores,fn,\
		j_info,motif_info,p,n,m,o,condenseddata)
		result.append(present)
	print('result',result)
	if numiterations > 1:
		for i in range(len(motif_info)):
			mot = motif_info[i][0]
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
		statistics.mean(avg_score),statistics.mean(avg_pscore),\
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
#promoter_file,motif_info[i][0],bits,p_bits,evalue,score,score_p, nsites,\
#fl_ct/nsites,(nsites-fl_ct)/nsites,fp_ct/nsites,p,n,m,o
def avg_condensedstats(final,motif_info,p,n,m,o,r):
	#print(motif_info)
	output = []
	#print('final',final)
	#print(len(final))
	if len(final) > 1:
		for i in range(len(motif_info)):
			#print('motif_info',motif_info)
			mot = motif_info[i][0]
			#final_info = []
			avg_nsites = []
			avg_score = []
			avg_eval = []
			avg_score = []
			avg_pscore = []
			avg_frate = []
			avg_srate = []
			avg_fprate = []
			avg_fnrate = []
			for j in range(len(final)):
				for k in range(len(final[j])):
					#print('final[j]',final[j])
					#print('final[j][k]',final[j][k])
					#print('final[j][k][1]',final[j][k][1])
					#print('mot',mot)
					#print('MOT', mot,'result mot',result[i][j][1])
					if  mot == final[j][k][1]:
						#print('MOT', mot,'final mot',final[j][k][1])
						avg_nsites.append(final[j][k][7])
						avg_eval.append(final[j][k][4])
						avg_score.append(final[j][k][5])
						avg_pscore.append(final[j][k][6])
						avg_frate.append(final[j][k][8])
						avg_srate.append(final[j][k][9])
						avg_fprate.append(final[j][k][10])
						avg_fnrate.append(final[j][k][11])
						#avg_fn.append(result[j][k][11])	
			#print('FIND INFO HERE!!!!!!!!!!!!!!!!!!!')
			#print(mot)
			#print('avg_nsites',avg_nsites)
			#print('avg_eval',avg_eval)
			#print('avg_score',avg_score)
			#print('avg_pscore', avg_pscore)
			#print('avg_frate',avg_frate)
			#print('avg_srate',avg_srate)
			#print('avg_fprate',avg_fprate)
			#print('avg_fn',avg_fnrate)
			#'file','motif','info','p_info','evalue','score','p score',\
			#'nsites','failrate','suces rate','false p rate
			output.append((final[j][k][0],mot,statistics.mean(avg_nsites),\
			final[j][k][2],final[j][k][3], statistics.mean(avg_eval),\
			statistics.mean(avg_score), statistics.mean(avg_pscore),\
			statistics.mean(avg_frate), statistics.mean(avg_srate),\
			statistics.mean(avg_fprate),statistics.mean(avg_fnrate),p,n,m,o,\
			r+1,statistics.stdev(avg_eval),statistics.stdev(avg_score)\
			,statistics.stdev(avg_pscore),statistics.stdev(avg_frate),\
			statistics.stdev(avg_srate),statistics.stdev(avg_fprate)\
			,statistics.stdev(avg_fnrate)))
		#print('final',final)
		return output
	else:
		for i in range(len(final)):
			for j in range(len(final[i])):
				output.append(final[i][j])
		#print(output)
		return output
		
		

if __name__ == '__main__':
	jasparfile = 'BS667.1.jaspar'
	dnafile = 'anidsteripromoters.fa'
	condenseddata = True
	i = 5
	freq = .9
	p = 400
	n = 20
	background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
	fl = generate_promoter_bg(jasparfile, p, n, freq, background,i,\
	condenseddata)
	print(calcbg_fromdatafile(dnafile))
	
	#print(fl)
	