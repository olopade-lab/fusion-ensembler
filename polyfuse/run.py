import argparse
import importlib
import os
import glob
import subprocess
from polyfuse import apps

import parsl
parsl.set_stream_logger()

parser = argparse.ArgumentParser()
parser.add_argument("sample_dirs", help="")
parser.add_argument("left_fq", help="")
parser.add_argument("right_fq", help="")
parser.add_argument("--config", default=None, help="Parsl config to parallelize with")
parser.add_argument("--out_dir", default=None, help="")
parser.add_argument("--library_dir", default=None, help="")
parser.add_argument("--ctat_release", default="GRCh38_gencode_v33_CTAT_lib_Apr062020", help="")
parser.add_argument("--ensemble_release", default=99, help="")
parser.add_argument("--container_type", default='docker',
        help="Container type to use; either 'singularity' or 'docker'")
args = parser.parse_args()
base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-2])
if args.config is None:
    args.config = os.path.join(base_dir, 'polyfuse', 'configs', 'local.py')
if args.out_dir is None:
    args.out_dir = os.path.join(base_dir, 'data')
if args.library_dir is None:
    args.library_dir = os.path.join(base_dir, 'data', 'library')
args.library_dir = os.path.abspath(args.library_dir)
args.sample_dirs = os.path.abspath(args.sample_dirs)
args.out_dir = os.path.abspath(args.out_dir)

spec = importlib.util.spec_from_file_location('', args.config)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
parsl.load(module.config)


# TODO split processing for model training and application, ensure no overwrite
if args.container_type == 'singularity':
    # TODO automate singularity hub builds from dockerhub
    # TODO pin versions for all images
    # TODO check everywhere to ensure no clobbering
    for local, remote in [
                ('polyfuse.sif', 'olopadelab/polyfuse:latest'),
                # ('pizzly.sif', 'olopadelab/pizzly:latest'),
                ('starseqr.sif', 'eagenomics/starseqr:0.6.7'),
                ('fusioncatcher.sif', 'olopadelab/fusioncatcher:latest'),
                ('starfusion.sif', 'trinityctat/starfusion:1.8.0')
                ('mapsplice2.sif', 'hiroko/mapsplice2-hg19')
            ]:
        image_path = '{base_dir}/docker/{local}'.format(base_dir=base_dir, local=local)
        # FIXME may require too much memory on some machines
        if not os.path.isfile(image_path):
            print('downloading {}'.format(image_path))
            subprocess.call(
                'singularity build {image_path} docker://{remote}'.format(
                    image_path=image_path, remote=remote
                ),
                shell=True
            )

fusion_catcher_build_dir = apps.download_fusioncatcher_build(
    output=os.path.join(args.library_dir, 'fusioncatcher')
)

ctat_dir = apps.download_ctat(
    args.library_dir,
    args.ctat_release
)

starseqr_star_index = apps.build_star_index(
    ctat_dir,
    'ref_genome.fa.starseqr.star.idx',
    container_type=args.container_type,
)

annotation = apps.download_ensemble_annotation(args.library_dir, args.ensemble_release)
assembly = apps.download_ensemble_assembly(args.library_dir, args.ensemble_release)
kallisto_index = apps.kallisto_index(assembly, args.container_type)

ref_split_by_chromosome_dir = apps.build_bowtie_index(
    apps.split_ref_chromosomes(
        ctat_dir,
        container_type=args.container_type
    ),
    container_type=args.container_type
)

gemtools_genome_index = apps.build_gemtools_genome_index(
    ctat_dir,
    container_type=args.container_type
)

gemtools_transcriptome_index_and_keys = apps.build_gemtools_transcriptome_index_and_keys(
    ctat_dir,
    gemtools_genome_index,
    container_type=args.container_type
)

sample_dirs = glob.glob(args.sample_dirs)
for sample_dir in sample_dirs:
    sample = os.path.split(sample_dir)[-1]
    output = os.path.join(args.out_dir, 'processed', sample)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    if (len(glob.glob(os.path.join(sample_dir, args.left_fq))) == 0) or  \
           (len(glob.glob(os.path.join(sample_dir, args.right_fq))) == 0):
        print('No fastqs found; skipping {}'.format(sample_dir))
        continue
    else:
        print('Fastqs found for {}'.format(os.path.join(sample_dir, args.left_fq)))

    left_fq = apps.gzip(
        apps.merge_lanes(
            os.path.join(sample_dir, args.left_fq),
            args.out_dir,
            sample,
            'R1'
        ),
        args.out_dir
    )
    right_fq = apps.gzip(
        apps.merge_lanes(
            os.path.join(sample_dir, args.right_fq),
            args.out_dir,
            sample,
            'R2'
        ),
        args.out_dir
    )

    arriba = apps.run_arriba(
        os.path.join(output, 'arriba'),
        ctat_dir,
        left_fq,
        right_fq,
        container_type=args.container_type
    )
    apps.parse_arriba(os.path.join(output, 'arriba'), inputs=[arriba])

    starfusion = apps.run_starfusion(
        os.path.join(output, 'starfusion'),
        left_fq,
        right_fq,
        ctat_dir,
        container_type=args.container_type
    )
    apps.parse_starfusion(os.path.join(output, 'starfusion'), inputs=[starfusion])

    starseqr = apps.run_starseqr(
        os.path.join(output, 'starseqr'),
        left_fq,
        right_fq,
        ctat_dir,
        starseqr_star_index,
        container_type=args.container_type
    )
    apps.parse_starseqr(os.path.join(output, 'starseqr'), inputs=[starseqr])

    fusioncatcher = apps.run_fusioncatcher(
        os.path.join(output, 'fusioncatcher'),
        left_fq,
        right_fq,
        fusion_catcher_build_dir,
        container_type=args.container_type
    )
    apps.parse_fusioncatcher(os.path.join(output, 'fusioncatcher'), inputs=[fusioncatcher])

    quant = apps.kallisto_quant(
        kallisto_index,
        assembly,
        os.path.join(output, 'pizzly'),
        left_fq,
        right_fq,
        container_type=args.container_type
    )

    pizzly = apps.run_pizzly(
        quant,
        annotation,
        assembly,
        os.path.join(output, 'pizzly'),
        container_type=args.container_type
    )
    apps.parse_pizzly(os.path.join(output, 'pizzly'), inputs=[pizzly])

    mapsplice2 = apps.run_mapsplice2(
        os.path.join(output, 'mapsplice2'),
        ctat_dir,
        ref_split_by_chromosome_dir,
        left_fq,
        right_fq,
        container_type=args.container_type
    )
    apps.parse_mapsplice2(os.path.join(output, 'mapsplice2'), inputs=[mapsplice2])

    chimpipe =  apps.run_chimpipe(
        os.path.join(output, 'chimpipe'),
        ctat_dir,
        gemtools_genome_index,
        gemtools_transcriptome_index_and_keys,
        left_fq,
        right_fq,
        container_type=args.container_type
    )


parsl.wait_for_current_tasks()
truth = apps.concatenate_true_fusions(args.sample_dirs, args.out_dir)

print('finished processing!')
