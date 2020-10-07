
import os
import argparse

## Makes arguments to input into speed_weeder

parser = argparse.ArgumentParser(
	description='Runs weeder on fasta file and reports results.')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<str>', help='fasta file of promoters')
parser.add_argument('--freq', required=True, type=str,
	metavar='<str>', help='two-letter freq code for background frequency')
parser.add_argument('--condenseddata', action='store_true',
	help='on/off switch')
arg = parser.parse_args()

## Runs weeder on the fasta file you input, via your terminal

os.system(f'weeder2.0/weeder2 -f {arg.fasta} -O {arg.freq}')
wr = f'{arg.fasta}.w2'
print(wr)

## Makes the weeder result into something readable and comparable to meme (copied from weeder_reader)

with open(wr) as wr:   #wr- weeder result
    motifs = []
    while True:

        line = wr.readline()
        if line == '': break
        if line.startswith('Matrix'):

            print(f'Motif')
            lineA = wr.readline()
            lineC = wr.readline()
            lineG = wr.readline()
            lineT = wr.readline()
            m = []
            lA = lineA.split()
            lC = lineC.split()
            lG = lineG.split()
            lT = lineT.split()
            size = len(lA)- 1
            for i in range(size):
                m.append({})
            for i in range(size):
                m[i]["A"] = float(lA[i + 1])
                m[i]["C"] = float(lC[i + 1])
                m[i]["G"] = float(lG[i + 1])
                m[i]["T"] = float(lT[i + 1])
            motifs.append(m)
            for x in m:
                print(x, '\n')



if arg.condenseddata:
    print('motifs', 'number of sites', 'distance', 'e-value')  # or p-value?











