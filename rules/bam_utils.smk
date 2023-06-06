align_flags2flag = { 'read_is_mapped'                       : 1
                   , 'read_mapped_as_part_of_pair'          : 2
                   , 'read_is_unmapped'                     : 4
                   , 'mate_is_unmapped'                     : 8
                   , 'read_reverse_strand'                  : 16
                   , 'mate_reverse_strand'                  : 32
                   , 'first_in_pair'                        : 64
                   , 'second_in_pair'                       : 128
                   , 'not_primary_alignment'                : 256
                   , 'read_fails_platform_or_vendor_checks' : 512
                   , 'read_is_pcr_or_optical_duplicate'     : 1024 }

def create_align_flag(flags):
	return sum([ align_flags2flag[f] for f in flags])

rule sort_name_bam:
	input:             '{path}/{bamfile}.bam'
	output: temp(local('{path}/{bamfile}.sorted_by_name.bam'))
	threads: 8
	conda: '../envs/samtools.yml'
	shell: 'samtools sort -@ {threads} -n -o {output} {input}'

rule index_bam:
	input:  '{path}/{bamfile}.bam'
	output: '{path}/{bamfile}.bam.bai'
	threads: 4
	conda: '../envs/samtools.yml'
	shell: 'samtools index -@ {threads} {input} {output}'

wildcard_constraints:
	align_type='|'.join(['aligned', 'aligned.transcriptome']),

def should_use_whitelist( w ):
	if not config['use_whitelist']:
		return 'no'
	elif w.align_type == 'aligned':
		return 'yes'
	elif w.align_type == 'aligned.transcriptome':
		return 'no'

rule filter_bam:
	input:
		sorted_file=multiext(os.path.join(config['align_dir'], '{sample_id}.{aligner}_{align_type}.unfiltered'), '.bam', '.bam.bai'),
		whitelist=lambda w: os.path.join(config['reference_dir'],config['database'],f'{config["assembly"]}_whitelist.bed') if should_use_whitelist(w) == "yes" else [],
	output:
		filtered_file=os.path.join(config['align_dir'], '{sample_id}.{aligner}_{align_type}.bam'),
	params:
		phred_quality_cuttoff=30,
		pass_flag=create_align_flag(['read_mapped_as_part_of_pair']),
		fail_flag=create_align_flag(['read_is_unmapped', 'mate_is_unmapped', 'not_primary_alignment', 'read_fails_platform_or_vendor_checks', 'read_is_pcr_or_optical_duplicate']),
		whitelist=lambda wildcards: should_use_whitelist( wildcards ),
	threads: 8
	conda: '../envs/samtools.yml'
	shell:
		'''
			if [[ "{params.whitelist}" == "yes" ]]; then
				samtools view -@ {threads} -q {params.phred_quality_cuttoff} -F {params.fail_flag} -f {params.pass_flag} -L {input.whitelist} -o {output} {input.sorted_file[0]}
			elif [[ "{params.whitelist}" == "no" ]]; then
				samtools view -@ {threads} -q {params.phred_quality_cuttoff} -F {params.fail_flag} -f {params.pass_flag}                      -o {output} {input.sorted_file[0]}
			fi
		'''

rule bam_reads:
	input:  os.path.join('{path}', '{bam}.bam'),
	output: os.path.join('{path}', '{bam}.bam2counts.txt'),
	params:
	conda: '../envs/samtools.yml'
	shell: 'samtools view -c {input} > {output}'

rule bam_info:
	input:  os.path.join('{path}', '{bam}.bam'),
	output: os.path.join('{path}', '{bam}.baminfo.yml'),
	conda: '../envs/samtools.yml'
	shell: 'samtools stats {input} | grep ^SN | cut -f 2- | sed "s/\t/ /g" > {output}'

rule bam_info2tsv:
	input:  os.path.join('{path}', '{bam}.baminfo.yml'),
	output: os.path.join('{path}', '{bam}.baminfo.tsv'),
	run:
		obj = yaml.load(open(input[0],'r').read(), Loader=yaml.Loader)
		df  = pd.DataFrame([obj])
		df.to_csv( output[0], sep='\t', index=False )

rule bam2bed:
	input:             os.path.join('{path}', '{bamfile}.bam'),
	output: temp(local(os.path.join('{path}', '{bamfile}.bam2bed.bed'))),
	conda: '../envs/samtools.yml'
	threads: 2
	shell: 'bedtools bamtobed -i {input} > {output}'

rule bed2weirdsort:
	input:             os.path.join('{path}', '{bedfile}.bed'),
	output: temp(local(os.path.join('{path}', '{bedfile}.weirdsort.bed'))),
	conda: '../envs/samtools.yml'
	shell: 'sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 {input} > {output}' #sort by chrom, end, start, strand
