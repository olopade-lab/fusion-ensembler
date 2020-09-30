from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
from parsl.executors import WorkQueueExecutor
from parsl.utils import get_all_checkpoints

config = Config(
        executors=[
            # WorkQueueExecutor(port=50056,
            WorkQueueExecutor(port=50057,
                autolabel=False,
                autocategory=False,
                # autolabel=True,
                # autocategory=True,
                shared_fs=True,
                # worker_options='--cores 38 --memory 190000 -dall',
                worker_options='-dall',
                provider=SlurmProvider(
                    'daenerys',
                    nodes_per_block=1,
                    # mem_per_node=0,
                    init_blocks=10,
                    max_blocks=10,
                    # scheduler_options='#SBATCH --exclude=kg15-11,kg15-23,kg15-12,kg15-13,kg15-14', # docker: Got permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Post http://%2Fvar%2Frun%2Fdocker.sock/v1.37/containers/create: dial unix /var/run/docker.sock: connect: permission denied.
                    # kg15-14: Error response from daemon: Cannot kill container: 347c9541f2d0: Cannot kill container 347c9541f2d0c84048209d1210d26540d61e26db7936a291958ca441d1bee197: grpc: the client connection is closing: failed precondition
                    # kg15-12: docker is stuck
                    worker_init='source activate gf; sh /cephfs/users/annawoodard/polyfuse/polyfuse/mem.sh >> /cephfs/users/annawoodard/polyfuse/polyfuse/logs/mem_$HOSTNAME.log &',
                    walltime='640:00:00'
                    )
                )
            ],
        checkpoint_mode='task_exit',
        checkpoint_files=get_all_checkpoints(),
        retries=3
)
# #SBATCH --exclude=kg15-13
