rule get_sailor:
	output: 'utils/sailor-1.0.4'
	params:
		url='https://s3-us-west-1.amazonaws.com/sailor-1.0.4/sailor-1.0.4'
	conda: '../envs/sailor.yml'
	shell: '''
		curl --location {params.url} --output {output}

		chmod u+x {output}
	'''
