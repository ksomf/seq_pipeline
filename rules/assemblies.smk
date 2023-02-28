import os

assembly2params = {
	'hg38': { 'species'         :'homo_sapiens'
	        , 'fasta'           :'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz'
	        , 'gff'             :'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz'
	        , 'encode_blacklist':'https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz'
	        , 'tss'             :'https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz'
	        , 'promoters'       :'https://www.encodeproject.org/files/ENCFF140XLU/@@download/ENCFF140XLU.bed.gz'
	        , 'enhancers'       :'https://www.encodeproject.org/files/ENCFF212UAV/@@download/ENCFF212UAV.bed.gz' }
}

wildcard_constraints:
	assembly='|'.join(assembly2params.keys()),

rule download_assembly:
	output: os.path.join(config['reference_dir'],'{assembly}.fasta'),
	params: url=lambda w: assembly2params[w.assembly]['fasta']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_gff:
	output: os.path.join(config['reference_dir'],'{assembly}.gff'),
	params: url=lambda w: assembly2params[w.assembly]['gff']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_blacklist:
	output: os.path.join(config['reference_dir'],'{assembly}_blacklist.unsorted.bed'),
	params: url=lambda w: assembly2params[w.assembly]['encode_blacklist']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_tss:
	output: os.path.join(config['reference_dir'],'{assembly}_tss.bed'),
	params: url=lambda w: assembly2params[w.assembly]['tss']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_promoters:
	output: os.path.join(config['reference_dir'],'{assembly}_promoters.bed'),
	params: url=lambda w: assembly2params[w.assembly]['promoters']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_enhancers:
	output: os.path.join(config['reference_dir'],'{assembly}_enhancers.bed'),
	params: url=lambda w: assembly2params[w.assembly]['enhancers']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule generate_gtf:
	input:  os.path.join(config['reference_dir'],'{assembly}.gff'),
	output: os.path.join(config['reference_dir'],'{assembly}.gtf'),
	conda: '../envs/assemblies.yml'
	shell:
		'''
		#https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gff_to_gtf.md
		#https://github.com/NBISweden/AGAT
		agat_convert_sp_gff2gtf.pl --gff {input} -o {output}
		'''

rule generate_chromosome_sizes:
	input: os.path.join(config['reference_dir'],'{assembly}.fasta'),
	output:
		sizes=os.path.join(config['reference_dir'],'{assembly}_chromosome_sizes.tsv'),
		fai=os.path.join(config['reference_dir'],'{assembly}.fasta.fai'),
	params: url=lambda w: assembly2params[w.assembly]['chromosome_sizes']
	conda: '../envs/assemblies.yml'
	shell: 'samtools faidx {input} && cut -f1,2 {input}.fai > {output.sizes}'

rule generate_sorted_blacklist:
	input:
		blacklist=os.path.join(config['reference_dir'],'{assembly}_blacklist.unsorted.bed'),
		sizes    =os.path.join(config['reference_dir'],'{assembly}_chromosome_sizes.tsv'),
	output: os.path.join(config['reference_dir'],'{assembly}_blacklist.bed'),
	conda: '../envs/assemblies.yml'
	shell: 'bedtools sort -faidx {input.sizes} -i {input.blacklist} > {output}'

rule generate_whitelist:
	input:
		blacklist=os.path.join(config['reference_dir'],'{assembly}_blacklist.bed'),
		sizes    =os.path.join(config['reference_dir'],'{assembly}_chromosome_sizes.tsv'),
	output: os.path.join(config['reference_dir'],'{assembly}_whitelist.bed'),
	conda: '../envs/assemblies.yml'
	shell: 'bedtools complement -i {input.blacklist} -g {input.sizes} > {output}'

