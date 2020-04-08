def test():
	import subprocess
	subprocess.check_output("rm -rf test_out; python2 fastq_preprocessor.py test_data --newDIR test_out",shell=True)
	import filecmp
	assert filecmp,filecmp("test_out.ICEFLAG-1H22C-DMSO-R1_S1_R1_raw.fastq.expect","test_out/ICEFLAG-1H22C-DMSO-R1_S1_R1_raw.fastq")
	
