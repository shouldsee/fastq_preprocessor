#!/usr/bin/env python
from setuptools import setup

# from distutils.core import setup

config = dict(
	name='fastq_preprocessor',
	version = '0.0.1', ### change in __init__.py in sync
	# package_dir={"": "."},
    packages=['.'],
	# include_package_data=True,
	license='MIT',
	author='Feng Geng',
	author_email='shouldsee.gem@gmail.com',
	# long_description=open('README.md').read(),
	# python_requires = '>=3.6',
	classifiers = [
	'Programming Language :: Python :: 2.7',
	'Programming Language :: Python :: 3.5',
	'Programming Language :: Python :: 3.6',
	'Programming Language :: Python :: 3.7',
	],
	install_requires=[
		x.strip() for x in open("requirements.txt","r")
        	if x.strip() and not x.strip().startswith("#")
	],
    entry_points={
        "console_scripts": [
            "preprocessor.py=fastq_preprocessor:main_entry",
            # "fastq_preprocessor=fastq_preprocessor:main_entry",
            "fastq_preprocess=fastq_preprocessor:main_entry",
            ]},	

)


if __name__ == '__main__':
	# from distutils.core import setup
	import os,glob,sys
	assert sys.version_info >= (2,7),('Requires python>=2.7, found python==%s'%('.'.join([str(x) for x in sys.version_info[:3]])))
	setup(**config)
