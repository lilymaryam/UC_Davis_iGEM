#!/usr/bin/python3
import os
import sys
import motiflib
import readmeme
import argparse

# setup
parser = argparse.ArgumentParser(
	description='Will generate differing fasta files to test parameters of MEME.')
# required arguments
parser.add_argument('--jaspardirectory', required=True, type=str,
	metavar='<str>', help='File to jaspar directory')
parser.add_argument('--memepath', required=True, type=str,
	metavar='<str>', help='path to meme software')
#parser.add_argument('--rint', required=True, type=int,
#	metavar='<int>', help='required integer argument')
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
#parser.add_argument('--switch', action='store_true',
#	help='on/off switch')
# finalization
arg = parser.parse_args()

#jaspardirectory  #~/Downloads/JASPAR2020_CORE_fungi_non-redundant_pfms_jaspar/
#memepath = #~/scratch/meme/bin/meme

#should this be command line?
filenum = ['0265']
promoter = [100, 200]#, #400, 600]
numseq = [5, 10]#, 15]
model = ['zoops']#, #'oops', 'anr']
markovorder = ['0', '1']

print('jasparfile','sequence', 'motifs found', 'markov','motif','strand','meme position',\
'meme width', 'e value', 'jaspar position', 'jaspar width','motif distance score',\
'motif score percentage','false negative','false positive','positional distance','fail',\
sep=', ')
for f in filenum:
	filepath = f'{arg.jaspardirectory}/MA{f}.1.jaspar'
	for p in promoter:
		for n in numseq:
			#note that motifapalooza must be found in the same directory 
			cmd = f'python3 motifapalooza.py --jasparfile {filepath} \
			 --numseq {n} --seqlen {p} > testmotif{f}.fa '
			os.system(cmd)
			for m in model:
				for o in markovorder:
						os.system(f'{arg.memepath} testmotif{f}.fa -dna -markov_order {o}\
						 -mod {m}')
						memestats = readmeme.readmemeout('meme_out/meme.xml')
						#data parsing
						jpositions = readmeme.read_testmotif(f'testmotif{f}.fa')
						jasparpwm = motiflib.read_JASPAR(filepath)
						memepwm = readmeme.memepwm('meme_out/meme.txt')
						pwmdistance = motiflib.global_motcompare(jasparpwm, memepwm,\
						{'A':.25,'C':.25,'G':.25,'T':.25})
						for i in range(0,len(memestats)):
							memestat = memestats[i]
							jpos = jpositions[i]
							mpos = memestat[4]
							positioninfo = motiflib.pos_accuracy(mpos, jpos, \
							len(memepwm), len(jasparpwm))	
						#data display !! always at the very end 	
							if memestat[1] == 0:
								if len(jpositions[i]) > 0:
									print(f, memestat[0], memestat[1], o, '', '', '',\
									 '', '', jpositions[i],len(jasparpwm),pwmdistance[0],\
									pwmdistance[1],positioninfo[2],positioninfo[3],\
									positioninfo[4], positioninfo[5], sep =', ')
								else:
									print(f, memestat[0], memestat[1], o, '', '', '', \
									'', '', jpositions[i],'',pwmdistance[0], \
									pwmdistance[1], positioninfo[2],positioninfo[3],\
									positioninfo[4],positioninfo[5], sep =', ')
							else:
								if len(jpositions[i]) > 0: 
									print(f, memestat[0], memestat[1], o, memestat[2],\
									memestat[3], memestat[4], memestat[6], memestat[8],\
									jpositions[i], len(jasparpwm),pwmdistance[0],\
									pwmdistance[1], positioninfo[2],positioninfo[3],\
									positioninfo[4],positioninfo[5], sep=', ' )
								else: 
									print(f, memestat[0], memestat[1], o, memestat[2],\
									memestat[3], memestat[4], memestat[6], memestat[8],\
									jpositions[i], '',pwmdistance[0], pwmdistance[1], \
									positioninfo[2], positioninfo[3],positioninfo[4],\
									positioninfo[5], sep=', ' )
							

						
						






