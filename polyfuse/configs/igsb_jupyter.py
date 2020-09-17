from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
from parsl.executors import HighThroughputExecutor

config = Config(
    executors=[
        HighThroughputExecutor(
            worker_debug=True,
            max_workers=7,
            address=address_by_hostname(),
            provider=SlurmProvider(
                'daenerys',
                nodes_per_block=1,
                init_blocks=10,
                max_blocks=10,
                worker_init='source activate gf',
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages',
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages; docker pull olopadelab/polyfuse',
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages',
                walltime='48:00:00'
            ),
        )
    ]
)
