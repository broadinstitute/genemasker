from setuptools import setup
from codecs import open
import os
from genemasker import __version__

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst')) as f:
    long_description = f.read()

setup(
    name='genemasker',
	description='This application produces masked group/gene variant inclusion files for rare variant analysis',
	long_description=long_description, 
    version=__version__.version,
    url='',
    author='Ryan Koesterer / Trang Nguyen',
	author_email='',
	install_requires=['numpy',
		'pandas',
		'scikit-learn',
		'scipy',
		'psutil',
		'dask'],
    entry_points={
       'console_scripts': [
			'genemasker = genemasker.__main__:main',
           ],
       },
    packages=['genemasker'], 
	include_dirs = [], 
	package_data={},
    classifiers = [
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
		],
	zip_safe=False
)
