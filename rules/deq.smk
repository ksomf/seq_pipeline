rule deq_diffbind:
	input:
		condition1_ips     = lambda wildcards: [ sample_id2bam[sample_id] for sample_id in condition2sample_ids[wildcards.condition1] ],
		condition1_inputs  = lambda wildcards: [ sample_id2bam[sample_id] for sample_id in condition2input_ids [wildcards.condition1] ],
		condition2_ips     = lambda wildcards: [ sample_id2bam[sample_id] for sample_id in condition2sample_ids[wildcards.condition2] ],
		condition2_inputs  = lambda wildcards: [ sample_id2bam[sample_id] for sample_id in condition2input_ids [wildcards.condition2] ],
		peaks         = lambda wildcards: [ os.path.join(config["peakcalling_dir"], f'{sample_id}_normal_peaks.narrowPeak') for sample_id in chain( condition2sample_ids[wildcards.condition1], condition2sample_ids[wildcards.condition2] ) ]
	output:
		counts = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}.counts.txt')
		peaks  = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}.txt')
	conda '../envs/deq.yml'
	script: '../scripts/deq.R'
