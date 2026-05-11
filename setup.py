from setuptools import setup, find_packages

with open(file='README.md', mode='r') as fh:
    long_description = fh.read()

setup(
    name='Metprofiler',
    version='0.1',
    author='Fabio Dorazio',
    description='Performs Differential methylation Analysis for a GEO dataset',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(exclude=('docs', 'tests')),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[ 
        "pandas==3.0.2",
        "numpy==2.4.4",
        "seaborn==0.13.2",
        "scipy==1.17.1",
        "statsmodels==0.14.6"],
    entry_points={
        'console_scripts': ['Metprofiler=Metprofiler.bin.main:main'],
    }
)