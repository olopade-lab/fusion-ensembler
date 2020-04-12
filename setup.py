from setuptools import setup, find_packages

setup(
    name='polyfuse',
    # version=VERSION,
    # description='',
    # long_description='',
    url='https://github.com/olopade-lab/polyfuse',
    author='Anna Woodard',
    author_email='annawoodard@uchicago.edu',
    license='MIT',
    # download_url='https://github.com/Parsl/parsl/archive/{}.tar.gz'.format(VERSION),
    include_package_data=True,
    packages=find_packages(),
    install_requires=['parsl', 'multiprocessing'],
    # scripts = ['parsl/executors/high_throughput/process_worker_pool.py',
    #            'parsl/executors/extreme_scale/mpi_worker_pool.py',
    #            'parsl/executors/low_latency/lowlatency_worker.py',
    #            'parsl/executors/workqueue/workqueue_worker.py',
    # ],

    keywords=['Workflows', 'Scientific computing'],
    # entry_points={'console_scripts':
    #   [
    #    'parsl-globus-auth=parsl.data_provider.globus:cli_run',
    #    'parsl-visualize=parsl.monitoring.visualization.app:cli_run',
    #   ]}
)
