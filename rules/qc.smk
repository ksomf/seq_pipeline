rule fastqc_fastq:
	input: sample_readlist
	output: [ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastqc.zip').replace('.fq.gz','_fastqc.zip').replace(' ','_') ) for f in sample_readlist ],
	params:
		output_dir=lambda wildcards, output: os.path.dirname(output[0]),
		safe_input=lambda wildcards, input: [ '"'+s+'"' for s in input ]
	threads: len(sample_readlist)
	conda: '../envs/qc.yml'
	shell: 'fastqc --noextract --threads {threads} --outdir {params.output_dir} {params.safe_input}'

rule make_multiqc_list:
	input: multiqc_inputs
	output: temp(local('multiqc_report_files.txt'))
	run:
		open(output[0], 'w+').write('\n'.join(input))

rule multiqc_report:
	input:  rules.make_multiqc_list.output
	output:
		report='multiqc_report.html',
		data=directory('multiqc_data/'),
	conda: '../envs/qc.yml'
	shell:
		'''
			test -d {output.data}   && rm -r {output.data}
			test -f {output.report} && rm    {output.report}
			multiqc --file-list {input}
		'''
