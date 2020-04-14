import parsl
from parsl.app.app import bash_app, python_app


@python_app(cache=True)
def build_star_index(
        assembly,
        annotation,
        output,
        image,
        container_type='docker',
        ):
    import os
    import subprocess
    import multiprocessing

    command = ['mkdir -p {output};']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {output}:/output',
            '-v {assembly}:/assembly:ro',
            '-v {annotation}:/annotation:ro',
            '{image}'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {output}:/output',
            '-B {assembly}:/assembly',
            '-B {annotation}:/annotation ',
            '{image}'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'STAR',
        '--runMode genomeGenerate',
        '--genomeDir /output',
        '--genomeFastaFiles /assembly',
        '--sjdbGTFfile /annotation',
        '--runThreadN {threads}',
        '--sjdbOverhang 200'
    ]

    subprocess.check_output(
        ' '.join(command).format(
            assembly=assembly,
            annotation=annotation,
            image=image,
            output=output,
            threads=os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
        ),
        shell=True
    )

    return output


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
            ext='.gz' if glob.glob(fastq)[0].endswith('.gz') else ''
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
    import multiprocessing


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
        '{cores}'
    ]

    return ' '.join(command).format(
            output=output,
            assembly=assembly,
            annotation=annotation,
            left_fq=left_fq,
            right_fq=right_fq,
            star_index=star_index,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
            cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 4)
        )

# FIXME Only run STAR once
@bash_app(cache=True)
def run_starfusion(
        output,
        left_fq,
        right_fq,
        genome_lib,
        container_type='docker',
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing


    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {left_fq}:{left_fq}:ro'
            '-v {right_fq}:{right_fq}:ro',
            '-v {genome_lib}:/genome_lib:ro',
            '-v {output}:/output',
            'trinityctat/starfusion:1.8.0'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {genome_lib}:/genome_lib',
            '-B {output}:/output',
            '{base_dir}/docker/starfusion.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        ' /usr/local/src/STAR-Fusion/STAR-Fusion ',
        '--left_fq {left_fq}',
        '--right_fq {right_fq}',
        '--genome_lib_dir /genome_lib',
        '-O /output',
        '--FusionInspector validate',
        '--examine_coding_effect',
        '--denovo_reconstruct',
        '--CPU {cores}'
    ]

    return ' '.join(command).format(
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        left_fq=left_fq,
        right_fq=right_fq,
        genome_lib=genome_lib,
        output=output,
        cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 8)
    )

@bash_app(cache=True)
def run_starseqr(
        output,
        left_fq,
        right_fq,
        genome_lib,
        star_index,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {left_fq}:{left_fq}:ro'
            '-v {right_fq}:{right_fq}:ro',
            '-v {genome_lib}:/genome_lib:ro',
            '-v {star_index}:/star_index:ro',
            '-v {output}:/output ',
            'eagenomics/starseqr:0.6.7'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {genome_lib}:/genome_lib',
            '-B {star_index}:/star_index:ro',
            '-B {output}:/output',
            '{base_dir}/docker/starseqr.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    import os

    command += [
        'starseqr.py ',
        '-1 {left_fq}',
        '-2 {right_fq}',
        '-p /output/ss',
        '-i /star_index',
        '-g /genome_lib/ref_annot.gtf',
        '-r /genome_lib/ref_genome.fa',
        '-m 1',
        '-vv',
        '-t {cores}'
    ]

    return ' '.join(command).format(
        star_index=star_index,
        output=output,
        genome_lib=genome_lib,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        left_fq=left_fq,
        right_fq=right_fq,
        cores=os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
    )

@bash_app(cache=True)
def download_fusioncatcher_build(output):
    import os

    if os.path.isfile(os.path.join(output, 'download.success')):
        return "echo 're-using existing fusioncatcher build'"
    command = """
    mkdir -p {output}
    cd {output}
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.aa
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ab
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ac
    wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ad
    cat human_v98.tar.gz.* | tar xz
    touch download.success
    ln -s human_v98 current
    cd -
    """

    return command.format(output=output)

@bash_app(cache=True)
def run_fusioncatcher(
        output,
        left_fq,
        right_fq,
        build_dir,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {left_fq}:{left_fq}:ro'
            '-v {right_fq}:{right_fq}:ro',
            '-v {build_dir}:/build_dir:ro',
            '-v {output}:/output ',
            'olopadelab/fusioncatcher:latest'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {build_dir}:/build_dir',
            '-B {output}:/output',
            '{base_dir}/docker/fusioncatcher.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    import os

    command += [
        'env PATH=/opt/fusioncatcher/v1.20/bin:$PATH fusioncatcher.py ',
        '-i {left_fq},{right_fq}',
        '-o /output',
        '-d /build_dir',
        '-p {cores}'
    ]

    return ' '.join(command).format(
        output=output,
        build_dir=build_dir,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        left_fq=left_fq,
        right_fq=right_fq,
        cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 16)
    )
