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
parser.add_argument('--end', required=False, type=int, default=100,
                    metavar='<int>', help='sequences to generate [%(default)i]')
arg = parser.parse_args()

n = arg.iterations
end = arg.end
num_seq = arg.seqs
distances = []
motif = motiflib.read_JASPAR(arg.file)
print()
for i in range(end):

    for j in range(n):
        sites = []
        for k in range(arg.seqs):
            site = motiflib.generate_site(motif)
            sites.append(site)
        pwm = motiflib.new_pwm(sites)
        distance = motiflib.compare_motifs(motif, pwm)
        distances.append(distance)
        #print(distances)
    #print(distances)
    avg_distance = statistics.mean(distances)
    print(i, avg_distance)

# if you collected the sites back into a motif, would it look like the motif?
# how would you measure the similarity?
# you might want to run this many times to look at the variance
# print(sites)