from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
from parsl.executors import WorkQueueExecutor
from parsl.utils import get_all_checkpoints

config = Config(
        executors=[
            WorkQueueExecutor(port=50056,
                autolabel=True,
                autocategory=True,
                shared_fs=True,
                worker_options='--cores 32',
                provider=SlurmProvider(
                    'daenerys',
                    nodes_per_block=1,
                    mem_per_node=0,
                    init_blocks=5,
                    max_blocks=5,
                    scheduler_options='#SBATCH --exclude=kg15-11,kg15-23', # docker: Got permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Post http://%2Fvar%2Frun%2Fdocker.sock/v1.37/containers/create: dial unix /var/run/docker.sock: connect: permission denied.
                    worker_init='source activate gf; docker stop $(docker ps -aq)',
                    walltime='240:00:00'
                    )
                )
            ],
        checkpoint_mode='task_exit',
        checkpoint_files=get_all_checkpoints(),
)
