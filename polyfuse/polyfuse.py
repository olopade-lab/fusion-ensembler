import argparse
import importlib
import os
import glob
import subprocess
from polyfuse import calling

import parsl
parsl.set_stream_logger()

def run():
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
    module.config.run_dir = os.path.join(args.out_dir, 'rundir')
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
                    ('starfusion.sif', 'trinityctat/starfusion:1.8.0'),
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

    fusion_catcher_build_dir = calling.download_fusioncatcher_build(
        output=os.path.join(args.library_dir, 'fusioncatcher')
    )

    ctat_dir = calling.download_ctat(
        args.library_dir,
        args.ctat_release
    )

    starseqr_star_index = calling.build_star_index(
        ctat_dir,
        'ref_genome.fa.starseqr.star.idx',
        container_type=args.container_type,
    )

    annotation = calling.download_ensemble_annotation(args.library_dir, args.ensemble_release)
    assembly = calling.download_ensemble_assembly(args.library_dir, args.ensemble_release)
    kallisto_index = calling.kallisto_index(assembly, args.container_type)

    ref_split_by_chromosome_dir = calling.build_bowtie_index(
        calling.split_ref_chromosomes(
            ctat_dir,
            container_type=args.container_type
        ),
        container_type=args.container_type
    )

    sample_dirs = glob.glob(args.sample_dirs)
    for sample_dir in sample_dirs:
        sample = os.path.split(sample_dir)[-1]
        output = os.path.join(args.out_dir, 'processed', sample)
        os.makedirs(os.path.dirname(output), exist_ok=True)
        if (len(glob.glob(os.path.join(sample_dir, args.left_fq))) == 0) or  \
               (len(glob.glob(os.path.join(sample_dir, args.right_fq))) == 0):
            print('No fastqs found; skipping {} and {}'.format(
                os.path.join(sample_dir, args.left_fq),
                os.path.join(sample_dir, args.right_fq)
                )
            )
            continue
        else:
            print('Fastqs found for {} and {}'.format(
                os.path.join(sample_dir, args.left_fq),
                os.path.join(sample_dir, args.right_fq)
                )
            )

        left_fq = calling.gzip(
            calling.merge_lanes(
                os.path.join(sample_dir, args.left_fq),
                args.out_dir,
                sample,
                'R1'
            ),
            args.out_dir
        )
        right_fq = calling.gzip(
            calling.merge_lanes(
                os.path.join(sample_dir, args.right_fq),
                args.out_dir,
                sample,
                'R2'
            ),
            args.out_dir
        )

        arriba = calling.run_arriba(
            os.path.join(output, 'arriba'),
            ctat_dir,
            left_fq,
            right_fq,
            container_type=args.container_type
        )
        calling.parse_arriba(os.path.join(output, 'arriba'), inputs=[arriba])

        starfusion = calling.run_starfusion(
            os.path.join(output, 'starfusion'),
            left_fq,
            right_fq,
            ctat_dir,
            container_type=args.container_type
        )
        calling.parse_starfusion(os.path.join(output, 'starfusion'), inputs=[starfusion])

        starseqr = calling.run_starseqr(
            os.path.join(output, 'starseqr'),
            left_fq,
            right_fq,
            ctat_dir,
            starseqr_star_index,
            container_type=args.container_type
        )
        calling.parse_starseqr(os.path.join(output, 'starseqr'), inputs=[starseqr])

        fusioncatcher = calling.run_fusioncatcher(
            os.path.join(output, 'fusioncatcher'),
            left_fq,
            right_fq,
            fusion_catcher_build_dir,
            container_type=args.container_type
        )
        calling.parse_fusioncatcher(os.path.join(output, 'fusioncatcher'), inputs=[fusioncatcher])

        quant = calling.kallisto_quant(
            kallisto_index,
            assembly,
            os.path.join(output, 'pizzly'),
            left_fq,
            right_fq,
            container_type=args.container_type
        )

        pizzly = calling.run_pizzly(
            quant,
            annotation,
            assembly,
            os.path.join(output, 'pizzly'),
            container_type=args.container_type
        )
        calling.parse_pizzly(os.path.join(output, 'pizzly'), inputs=[pizzly])

        mapsplice2 = calling.run_mapsplice2(
            os.path.join(output, 'mapsplice2'),
            ctat_dir,
            ref_split_by_chromosome_dir,
            left_fq,
            right_fq,
            container_type=args.container_type
        )
        calling.parse_mapsplice2(os.path.join(output, 'mapsplice2'), inputs=[mapsplice2])

    parsl.wait_for_current_tasks()

    print('finished processing!')
