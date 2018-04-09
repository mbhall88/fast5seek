#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['h5py', 'pysam']

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Michael Hall",
    author_email='mbhall88@gmail.com',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Get paths for fast5 files contained in BAM, SAM, or fastq.",
    entry_points={
        'console_scripts': [
            'bam2fast5=bam2fast5.cli:cli',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='bam2fast5',
    name='bam2fast5',
    packages=find_packages(include=['bam2fast5']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mbhall88/bam2fast5',
    version='0.1.0',
    zip_safe=False,
)
