rule fastqc_fastq:
	input: sample_readlist
	output:
		html           =[ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastq.html') ) for f in sample_readlist ]
		html_compressed=[ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastq.zip' ) ) for f in sample_readlist ]
	params:
		output_dir=lambda wildcards, output: os.path.dirname(output[0]),
		initial_output_prefix=lambda wildcards, input: os.path.basename(input[0]).replace('.gz','').replace('.fastq','_fastq'),
	conda: '../envs/qc.yml'
	shell: 'fastqc --outdir {params.output_dir} {input}'

rule make_multiqc_list:
	input: [ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastq.zip' ) ) for f in sample_readlist ]
	output: temp(local('multiqc_report_files.txt'))
	run:
		open(output[0], 'w+').write('\n'.join(input))

rule run_multiqc:
	input:  'multiqc_report_files.txt'
	output: 'multiqc_report.html'
	params:
		unique_folders=lambda wildcards, input: os.path.dirname(qc_input_files)
	conda: '../envs/qc.yml'
	shell: 'multiqc --file-list {input}'
