from setuptools import setup, find_packages
from polyfuse.version import VERSION

setup(
    name='polyfuse',
    version=VERSION,
    description='A meta-caller for fusion transcript detection',
    url='https://github.com/olopade-lab/polyfuse',
    author='The Olopade Lab',
    author_email='annawoodard@uchicago.edu',
    license='MIT',
    include_package_data=True,
    packages=find_packages(),
    install_requires=['parsl', 'seaborn', 'joblib', 'sklearn', 'tables', 'pyfaidx', 'upsetplot', 'twine', 'pyensembl'],
    keywords=['Workflows', 'Scientific computing', 'fusion transcript detection', 'bioinformatics'],
    entry_points = {
        'console_scripts': ['polyfuse=polyfuse.polyfuse:run'],
    }
)
# conda install -c bioconda survivor=1.0.6
