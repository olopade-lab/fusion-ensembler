from setuptools import setup, find_packages
from polyfuse.version import VERSION

setup(
    name='polyfuse',
    version=VERSION,
    description='A meta-caller for fusion transcript detection',
    url='https://github.com/olopade-lab/polyfuse',
    author='Anna Woodard',
    author_email='annawoodard@uchicago.edu',
    license='MIT',
    # download_url='https://github.com/Parsl/parsl/archive/{}.tar.gz'.format(VERSION),
    include_package_data=True,
    packages=find_packages(),
    install_requires=['parsl', 'seaborn', 'joblib', 'sklearn', 'tables', 'pyfaidx', 'upsetplot', 'twine'],
    keywords=['Workflows', 'Scientific computing', 'fusion transcript detection', 'bioinformatics'],
)
