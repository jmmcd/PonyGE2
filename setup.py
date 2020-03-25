from setuptools import setup

setup(
    name='PonyGE2',
    version='0.0.1',
    author='Raphael Fischer',
    author_email='raphael.fischer[at]tu-dortmund.de',
    packages=['datasets', 'grammars', 'parameters', 'seeds', 'src'],
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        "": ["*.txt", "*.bnf", "*.pybnf"]
    },
    include_package_data=True,
    license='LICENSE',
    description='Allows spatio-temporal gap filling with machine learning.',
    long_description=open('README.md').read(),
)