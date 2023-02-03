#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup
import os
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

INSTALL_REQUIRES = ['numpy', 'pandas', 'bioframe', 'pomegranate', 'pybbi', 'natsort']

setup(
      name='hmm_bigwigs',
      packages=['hmm_bigwigs'],
      entry_points={
          'console_scripts': ['bigwig_hmm.py = hmm_bigwigs.__main__:main']},
      install_requires=INSTALL_REQUIRES,
      description='A tool to perform HMM analysis on genomic data in bigWig files',
      long_description=long_description,
      long_description_content_type='text/markdown',
      project_urls={'Source':'https://github.com/gspracklin/hmm_bigwigs',
                    'Issues':'https://github.com/gspracklin/hmm_bigwigs/issues'},
      author='George Spracklin',
      author_email='george.spracklin@gmail.com',
      classifiers=[
        "Programming Language :: Python",
    ],
    zip_safe=False
)
