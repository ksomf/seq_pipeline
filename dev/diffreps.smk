rule diffreps_install:
	output: TODO
	conda: '../envs/diffreps.yml'
	shell: 'cpanm PROFSHEN/diffReps-1.55.3.tar.gz'
