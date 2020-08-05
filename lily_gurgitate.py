import argparse
import motiflib
import statistics

parser = argparse.ArgumentParser(
    description='Do motif instances look like motifs?')
parser.add_argument('--file', required=True, type=str,
    metavar='<str>', help='motif file in JASPAR format')
parser.add_argument('--seqs', required=False, type=int, default=5,
    metavar='<int>', help='sequences to generate [%(default)i]')
parser.add_argument('--iterations', required=False, type=int, default=5,
	metavar='<int>', help='sequences to generate [%(default)i]')
arg = parser.parse_args()

n = arg.iterations
num_seq = arg.seqs
#creates the initial pwm
motif = motiflib.read_JASPAR(arg.file)
print('number of seqs', 'distance', sep=', ')
for i in range(1,num_seq+1): 
	distances =[]
	for k in range(n):
		sites = []
		for l in range(i):
			#a single motif from the initial pwm
			site = motiflib.generate_site(motif)
			#sites is a list of all motifs generated 
			sites.append(site)
		pwm = motiflib.new_pwm(sites)
		distance = motiflib.compare_motifs(motif, pwm)
		distances.append(distance)
	avg_distance = statistics.mean(distances)
	print(i, f'{avg_distance:.3f}',sep=', ')
	


# if you collected the sites back into a motif, would it look like the motif?
# how would you measure the similarity?
# you might want to run this many times to look at the variance
# print(sites)