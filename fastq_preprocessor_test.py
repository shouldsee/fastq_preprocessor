def test():
	from path import Path
	import sys
	import subprocess
	subprocess.check_output("rm -rf test_out; python2 fastq_preprocessor.py $PWD/test_data/199R/S2 --newDir test_out",shell=True)
	import filecmp
	assert filecmp,filecmp("test_out.ICEFLAG-1H22C-DMSO-R1_S1_R1_raw.fastq.expect","test_out/ICEFLAG-1H22C-DMSO-R1_S1_R1_raw.fastq")

	with Path('test_out.synobio').makedirs_p():
		subprocess.check_output('''
		set -e
		BIN={sys.executable}
		ID=../test_data/199R/S2
		$BIN ../synobio.preprocessor.py --autoNewDir 1 --moveRaw 1 $ID
		'''.format(**locals()),shell=True,executable="bash")