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
	threads: 16
	conda: '../envs/samtools.yml'
	shell: 'samtools sort -@ {threads} -n -o {output} {input}'

rule index_bam:
	input:  '{path}/{bamfile}.bam'
	output: '{path}/{bamfile}.bam.bai'
	threads: 4
	conda: '../envs/samtools.yml'
	shell: 'samtools index -@ {threads} {input} {output}'

rule filter_bam:
	input:
		sorted_file=multiext(os.path.join(config['align_dir'], '{sample_id}.{aligner}_aligned.unfiltered'), '.bam', '.bam.bai'),
		whitelist=expand(os.path.join(config['reference_dir'],'{assembly}_whitelist.bed'), assembly=config['assembly'], allow_missing=True),
	output:
		filtered_file=os.path.join(config['align_dir'], '{sample_id}.{aligner}_aligned.bam'),
	params:
		phred_quality_cuttoff=30,
		pass_flag=create_align_flag(['read_mapped_as_part_of_pair']),
		fail_flag=create_align_flag(['read_is_unmapped', 'mate_is_unmapped', 'not_primary_alignment', 'read_fails_platform_or_vendor_checks', 'read_is_pcr_or_optical_duplicate'])
	conda: '../envs/samtools.yml'
	shell: 'samtools view -q {params.phred_quality_cuttoff} -F {params.fail_flag} -f {params.pass_flag} -L {input.whitelist} -o {output} {input.sorted_file[0]}'
