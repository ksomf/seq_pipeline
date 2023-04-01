import os

assembly2params = {
	'hg38': { 'species'        : 'homo_sapiens'
	        , 'ensembl'  : { 'fasta' : 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' #According to https://www.biostars.org/p/342482/ is fine
	                      #, 'gff'   : 'https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz'
	                       , 'gtf'   : 'https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz'                     }
	        , 'ncbi'     : { 'fasta' : 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
	                      #, 'gff'   : 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz'
	                       , 'gtf'   : 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz' }
	        , 'ucsc'     : { 'fasta' : 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.masked.gz'
	                       , 'gtf'   : 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz' }
	        , 'chip_blacklist' : 'http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz'
	        , 'tss'            : 'https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz'
	        , 'promoters'      : 'https://www.encodeproject.org/files/ENCFF140XLU/@@download/ENCFF140XLU.bed.gz'
	        , 'enhancers'      : 'https://www.encodeproject.org/files/ENCFF212UAV/@@download/ENCFF212UAV.bed.gz' }
}

wildcard_constraints:
	assembly='|'.join(assembly2params.keys()),
	database='|'.join(['ensembl','ncbi','ucsc']),

rule download_assembly:
	output: os.path.join(config['reference_dir'],'{database}','{assembly}.fasta'),
	params: url=lambda w: assembly2params[w.assembly][w.database]['fasta']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_gff:
	output: os.path.join(config['reference_dir'],'{database}','{assembly}.gff'),
	params: url=lambda w: assembly2params[w.assembly][w.database]['gff']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

rule download_blacklist:
	output: os.path.join(config['reference_dir'],'{assembly}_blacklist.unsorted.bed'),
	params: url=lambda w: assembly2params[w.assembly]['chip_blacklist']
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

rule download_gtf:
	output: os.path.join(config['reference_dir'],'{database}','{assembly}.gtf'),
	params: url=lambda w: assembly2params[w.assembly][w.database]['gtf']
	conda: '../envs/assemblies.yml'
	shell: 'curl --location {params.url} | zcat > {output}'

#rule generate_gtf:
#	input:  os.path.join(config['reference_dir'],'{database}','{assembly}.gff'),
#	output: os.path.join(config['reference_dir'],'{database}','{assembly}.gtf'),
#	log:    os.path.join(config['reference_dir'],'{database}','{assembly}_agat_gff2gtf.log')
#	conda: '../envs/assemblies.yml'
#	shell:
#		'''
#		#https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gff_to_gtf.md
#		#https://github.com/NBISweden/AGAT
#		agat_convert_sp_gff2gtf.pl --gff {input} -o {output} > {log}
#		'''

rule generate_faidx:
	input: '{path}/{fa_name}.fasta',
	output: temp(local('{path}/{fa_name}.fasta.fai')),
	conda: '../envs/assemblies.yml'
	shell: 'samtools faidx {input}'

rule generate_chromosome_sizes:
	input:  os.path.join(config['reference_dir'],'{database}','{assembly}.fasta.fai'),
	output: os.path.join(config['reference_dir'],'{database}','{assembly}_chromosome_sizes.tsv'),
	conda: '../envs/assemblies.yml'
	shell: 'cut -f1,2 {input} > {output}'

rule generate_matching_names:
	input:
		blacklist=os.path.join(config['reference_dir'],'{assembly}_blacklist.unsorted.bed'),
	output:      os.path.join(config['reference_dir'],'{database}','{assembly}_blacklist.unsorted.bed'),
	conda: '../envs/assemblies.yml'
	shell:
		'''
			if [[ {wildcards.database} == "ensembl" ]]; then
				sed "s/^chr//g" {input.blacklist} > {output}
			else
				cp {input.blacklist} {output}
			fi
		'''

rule generate_sorted_blacklist:
	input:
		blacklist=os.path.join(config['reference_dir'],'{database}','{assembly}_blacklist.unsorted.bed'),
		sizes    =os.path.join(config['reference_dir'],'{database}','{assembly}_chromosome_sizes.tsv'),
	output:      os.path.join(config['reference_dir'],'{database}','{assembly}_blacklist.bed'),
	conda: '../envs/assemblies.yml'
	shell: 'bedtools sort -faidx {input.sizes} -i {input.blacklist} > {output}'

rule generate_whitelist:
	input:
		blacklist=os.path.join(config['reference_dir'],'{database}','{assembly}_blacklist.bed'),
		sizes    =os.path.join(config['reference_dir'],'{database}','{assembly}_chromosome_sizes.tsv'),
	output:      os.path.join(config['reference_dir'],'{database}','{assembly}_whitelist.bed'),
	conda: '../envs/assemblies.yml'
	shell: 'bedtools complement -i {input.blacklist} -g {input.sizes} > {output}'
