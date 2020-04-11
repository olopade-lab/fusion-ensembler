import parsl
from parsl.app.app import bash_app, python_app


@python_app(cache=True)
def build_star_index(
        assembly,
        annotation,
        path,
        image='uhrigs/arriba:1.1.0',
        container_type='docker',
        ):
    import os
    import subprocess

    if container_type == 'docker':
        command = [
            'docker run',
            '--rm',
            '-v {assembly}:/assembly:ro',
            '-v {annotation}:/annotation:ro',
            '{image}'
        ]
    elif container_type == 'singularity':
        command = [
            'singularity exec',
            '-B {assembly}:/assembly',
            '-B {annotation}:/annotation'
            ' docker://{image}'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'STAR',
        '--runMode genomeGenerate',
        '--genomeDir {path}',
        '--genomeFastaFiles {assembly}',
        '--sjdbGTFfile {annotation}',
        '--runThreadN {threads}',
        '--sjdbOverhang 200'
    ]

    subprocess.check_output(
        ' '.join(command).format(
            assembly=assembly,
            annotation=annotation,
            image=image,
            path=path,
            threads=os.environ.get('PARSL_CORES', 4)
        ),
        shell=True
    )

    return path


@python_app(cache=True)
def merge_lanes(fastq, base_dir, sample, tag='R1'):
    import glob
    import subprocess

    if len(glob.glob(fastq)) == 1:
        return fastq
    else:
        out_dir = '{base_dir}/data/interim/{sample}'.format(
            base_dir=base_dir,
            sample=sample
        )
        merged_fastq = '{out_dir}/merged.{tag}.fastq{ext}'.format(
            out_dir=out_dir,
            tag=tag,
            ext='.gz' if fastq.endswith('.gz') else ''
        )
        subprocess.check_output(
            'mkdir -p {out_dir}; cat {fastq} > {merged_fastq}'.format(
                out_dir=out_dir,
                fastq=fastq,
                merged_fastq=merged_fastq
                ),
            shell=True
        )
        return merged_fastq

@bash_app(cache=True)
def run_arriba(
        output,
        assembly,
        annotation,
        left_fq,
        right_fq,
        star_index,
        container_type='docker',
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os


    command = ['echo $HOSTNAME; mkdir -p {output}; cd {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {assembly}:/assembly:ro',
            '-v {annotation}:/annotation:ro',
            '-v {left_fq}:{left_fq}:ro'
            '-v {right_fq}:{right_fq}:ro',
            '-v {star_index}:/star_index:ro',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {assembly}:/assembly',
            '-B {annotation}:/annotation',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {star_index}:/star_index',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        '/usr/local/src/arriba_v1.2.0/run_arriba.sh',
        '/star_index',
        '/annotation',
        '/assembly',
        '/usr/local/src/arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz',
        '{left_fq}',
        '{right_fq}',
        '{threads}'
    ]

    return ' '.join(command).format(
            image=image,
            output=output,
            assembly=assembly,
            annotation=annotation,
            left_fq=left_fq,
            right_fq=right_fq,
            star_index=star_index,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
            threads=os.environ.get('PARSL_CORES', 4)
        )

