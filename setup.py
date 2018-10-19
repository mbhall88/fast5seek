#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['ont_fast5_api', 'pysam']

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
        'Programming Language :: Python :: 3.7',
    ],
    description="Get paths for fast5 files contained in BAM, SAM, or fastq.",
    entry_points={
        'console_scripts': [
            'fast5seek=fast5seek.cli:cli',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='fast5seek',
    name='fast5seek',
    packages=find_packages(include=['fast5seek']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mbhall88/fast5seek',
    version='0.1.1',
    zip_safe=False,
)
