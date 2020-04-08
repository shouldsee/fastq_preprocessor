def test():
	import subprocess
	subprocess.check_output("rm -rf test_out; python2 fastq_preprocessor.py test_data --newDIR test_out",shell=True)
	
