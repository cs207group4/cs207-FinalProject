"""setup file

require libraries like numpy...

Some setup code borrowed from https://github.com/activescott/python-package-example
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pychemkin',

    version='1.1.0',

    description='Chemical Kinetics Library',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/cs207group4/cs207-FinalProject',

    # Author details
    author='CS207 Group4',
    author_email = 'zzy8200@gmail.com',

    # Choose your license
    license='GNU General Public License v3 (GPLv3)',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords= 'chemical kinetics',

    packages= ['pychemkin'],

    install_requires=['numpy>=1.13.3','scipy>=1.0.0','pandas>=0.20.3','tables>=3.4.2','matplotlib>=2.0.2'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    setup_requires = ['pytest-runner','numpy'],
    tests_require = ['pytest','coverage','coveralls','pytest-cov'],
    test_suite='tests',
    include_package_data=True
)
