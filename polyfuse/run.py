import argparse
import importlib
import os
import glob
import subprocess
from polyfuse import apps

import parsl
parsl.set_stream_logger()

parser = argparse.ArgumentParser()
parser.add_argument("genome_lib", help="")
parser.add_argument("sample_dirs", help="")
parser.add_argument("left_fq", help="")
parser.add_argument("right_fq", help="")
parser.add_argument("--config", default=None, help="Parsl config to parallelize with")
parser.add_argument("--outdir", default=None, help="")
parser.add_argument("--container_type", default='docker',
        help="Container type to use; either 'singularity' or 'docker'")
args = parser.parse_args()
base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-2])
if args.config is None:
    args.config = os.path.join(base_dir, 'polyfuse', 'configs', 'local.py')
if args.outdir is None:
    args.outdir = os.path.join(base_dir, 'data', 'processed')
args.genome_lib = os.path.abspath(args.genome_lib)
args.sample_dirs = os.path.abspath(args.sample_dirs)
args.outdir = os.path.abspath(args.outdir)

spec = importlib.util.spec_from_file_location('', args.config)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
parsl.load(module.config)

assembly = os.path.join(args.genome_lib, 'ref_genome.fa')
annotation = os.path.join(args.genome_lib, 'ref_annot.gtf')
star_index = os.path.join(args.genome_lib, 'ref_genome.fa.star.idx')

if args.container_type == 'singularity':
    # TODO automate singularity hub builds from dockerhub
    for local, remote in [
                ('polyfuse.sif', 'docker://olopadelab/polyfuse:latest'),
                ('starseqr.sif', 'docker://eagenomics/starseqr:0.6.7')
            ]:
        image_path = '{base_dir}/docker/{local}'.format(base_dir=base_dir, local=local)
        # FIXME may require too much memory on some machines
        if not os.path.isfile(image_path):
            subprocess.call(
                'singularity build {image_path} docker://{remote}'.format(
                    image_path=image_path, remote=remote
                )
            )

sample_dirs = glob.glob(args.sample_dirs)
for sample_dir in sample_dirs:
    sample = os.path.split(sample_dir)[-1]
    output = os.path.join(args.outdir, sample)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    left_fq = apps.merge_lanes(
        os.path.join(sample_dir, args.left_fq),
        base_dir,
        sample,
        'R1'
    )
    right_fq = apps.merge_lanes(
        os.path.join(sample_dir, args.right_fq),
        base_dir,
        sample,
        'R2'
    )

    apps.run_arriba(
        os.path.join(output, 'arriba'),
        assembly,
        annotation,
        left_fq,
        right_fq,
        star_index,
        container_type=args.container_type
    )

    apps.run_starfusion(
        os.path.join(output, 'starfusion'),
        left_fq,
        right_fq,
        args.genome_lib,
        container_type=args.container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME
    )

parsl.wait_for_current_tasks()

print('finished processing!')
