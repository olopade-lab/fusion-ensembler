from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname, address_by_route
from parsl.executors import HighThroughputExecutor
from parsl.monitoring.monitoring import MonitoringHub
from parsl.utils import get_all_checkpoints

config = Config(
    executors=[
        HighThroughputExecutor(
            worker_debug=True,
            max_workers=5,
            address=address_by_hostname(),
            provider=SlurmProvider(
                'daenerys',
                nodes_per_block=1,
                init_blocks=13,
                max_blocks=13,
                scheduler_options='#SBATCH --exclude=kg15-8,kg15-11', # docker: Error response from daemon: Get https://registry-1.docker.io/v2/: dial tcp: lookup registry-1.docker.io on [::1]:53: dial udp [::1]:53: connect: cannot assign requested address.
                worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages',
                walltime='48:00:00'
            ),
        )
    ],
    checkpoint_mode='task_exit',
    checkpoint_files=get_all_checkpoints(),
   # monitoring=MonitoringHub(
   #     hub_address=address_by_hostname(),
   #     hub_port=55055,
   #     monitoring_debug=False,
   #     resource_monitoring_interval=10,
   # ),
   # retries=2,
)
