
import json
import os
import sys

# Note: some hard-coded paths here


# Part 1: Run BLASTP

params = 'W=4 T=20 E=1e-5 mformat=3'
genomes = []
for file in os.listdir('Proteins'): genomes.append(file)

"""
n = 1
t = 786313 # Asp only
for s1 in genomes:
	if not s1.startswith('Asp'): continue
	print(s1)
		
	for file in os.listdir(f'Proteins/{s1}'):
		basename = os.path.basename(file)
		cluster, x = os.path.splitext(basename)
		for s2 in genomes:
			if s1 == s2: continue
			if not s2.startswith('Asp'): continue
			print(s1, cluster, s2, n, t, n/t)
			n += 1
			
			out = f'BLASTOUT/{s1}-{cluster}-{s2}.blastp'
			if os.path.isfile(out):
				print(f'{out} exists, skipping')
				continue
			cmd = f'blastp BLASTDB/{s2}.fa Proteins/{s1}/{file} {params}'
			os.system(f'{cmd} > {out}')
		#	sys.exit(0) # testing

"""
genome = {}  # indexed by species, contains all cids
cluster = {} # indexed by cid, contains list of all members
parent = {}  # indexed by pid, contains name of parent
species = {} # indexed by pid, contains the name of the species

for sid in os.listdir('Proteins'):
	if not sid.startswith('Asp'): continue
	for file in os.listdir(f'Proteins/{sid}'):
		cid, ext = file.split('.')
		if sid not in genome: genome[sid] = []
		genome[sid].append(cid)
		with open(f'Proteins/{sid}/{file}') as fp:
			for line in fp.readlines():
				if line.startswith('>'):
					line = line.rstrip()
					n, c, pid = line.split(':')
					parent[pid] = cid
					species[pid] = sid
					if cid not in cluster: cluster[cid] = []
					cluster[cid].append(pid)
					

# Part 2: Create JSON file of cluster info

data = {
	'species':  genome,
	'clusters': cluster,
	'parent_cluster':  parent,
	'parent_species': species,
	
}

with open('clusters.json', 'w') as fp:
	fp.write(json.dumps(data, indent=4))


