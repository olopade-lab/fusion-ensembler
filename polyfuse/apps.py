import parsl
from parsl.app.app import bash_app, python_app

@python_app
def assemble_training_data(sample, callers, out_dir):
    import os
    import pandas as pd

    # TODO: filter before loading, to save memory?
    caller_data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'), 'caller_data')
    caller_data = caller_data[caller_data['sample'] == sample]
    if not set(caller_data.caller.unique()).issubset(set(callers)):
        return [], []

    called_fusions = caller_data.fusion.unique() # TODO: deal with normalization
    true_fusions = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'true_fusions')
    true_fusions = true_fusions[true_fusions['sample'] == sample]

    X_train = []
    Y_train = []
    for fusion in called_fusions:
        row = []
        for c in callers:
            data = caller_data.loc[(caller_data.fusion == fusion) & (caller_data.caller == c), 'sum_J_S']
            if len(data) > 0:
                row += [data.values[0]]
                #row += [1]
            else:
                row += [0]
        X_train += [row]
        Y_train += [1 if any(true_fusions.fusion.isin([fusion])) else 0]

    return X_train, Y_train


@python_app
def concatenate_true_fusions(sample_dirs, out_dir):
    import pandas as pd
    import glob
    import os

    true_fusions = pd.concat([pd.read_pickle(path) for path in glob.glob(os.path.join(sample_dirs, 'truth.pkl'))])

    output = '{out_dir}/true_fusions.hdf'.format(out_dir=out_dir)
    true_fusions.to_hdf(output, 'true_fusions', mode='w')

    return output

@python_app
def concatenate_caller_data(out_dir):
    import pandas as pd
    import glob
    import os

    caller_data = pd.concat(
        [
            pd.read_pickle(f)[['fusion', 'spanning_reads', 'junction_reads', 'sample', 'caller']]
            for f in
            glob.glob(os.path.join(out_dir, '*', '*', 'fusions.pkl'))
        ]
    )
    caller_data['sum_J_S'] = caller_data['junction_reads'] + caller_data['spanning_reads']
    output = '{out_dir}/caller_data.hdf'.format(out_dir=out_dir)
    caller_data.to_hdf(output, 'caller_data', mode='w')

    return output

@python_app(cache=True)
def build_star_index(
        assembly,
        annotation,
        output,
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
            'eagenomics/starseqr:0.6.7'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {output}:/output',
            '-B {assembly}:/assembly',
            '-B {annotation}:/annotation ',
            '{base_dir}/docker/starseqr.sif'
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
            output=output,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
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
        return glob.glob(fastq)[0]
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

# FIXME This should write to the interim data directory, not the input data directory.
@python_app(cache=True)
def gzip(fastq):
    import subprocess

    if fastq.endswith('.gz'):
        return fastq
    else:
        subprocess.check_output('gzip {fastq}'.format(fastq=fastq), shell=True)
        return fastq + '.gz'

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


    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {assembly}:/assembly:ro',
            '-v {annotation}:/annotation:ro',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {star_index}:/star_index:ro',
            '-v {output}:/output',
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
            '-B {output}:/output',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        '/bin/bash -c "',
        'cd /output;',
        '/usr/local/src/arriba_v1.2.0/run_arriba.sh',
        '/star_index',
        '/annotation',
        '/assembly',
        '/usr/local/src/arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz',
        '{left_fq}',
        '{right_fq}',
        '{cores}"'
    ]

    return ' '.join(command).format(
            output=output,
            assembly=assembly,
            annotation=annotation,
            left_fq=left_fq,
            right_fq=right_fq,
            star_index=star_index,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
            cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 8)
        )


@python_app
def parse_arriba(out_dir, inputs=[]):
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'fusions.tsv')
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    data = pd.read_csv(path, sep='\t')
    # data = data[data.confidence.str.contains('medium|high')]
    data.rename(columns={'#gene1': 'gene1'}, inplace=True)
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['junction_reads'] = data.split_reads1 + data.split_reads2
    data['spanning_reads'] = data.discordant_mates
    data['caller'] = caller
    data['sample'] = sample

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output


# FIXME Only run STAR once
# FIXME add back validation?
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
            '-v {left_fq}:{left_fq}:ro',
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
        '/usr/local/src/STAR-Fusion/STAR-Fusion',
        '--left_fq {left_fq}',
        '--right_fq {right_fq}',
        '--genome_lib_dir /genome_lib',
        '-O /output',
        # '--FusionInspector validate',
        # '--examine_coding_effect',
        # '--denovo_reconstruct',
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

@python_app
def parse_starfusion(out_dir, inputs=[]):
    import os
    import re
    import pandas as pd

    path = os.path.join(out_dir, 'star-fusion.fusion_predictions.tsv')
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    data = pd.read_csv(path, sep='\t')

    data['breakpoint1'] = data['LeftBreakpoint'].str.rstrip('\:\+|\:\-').str.lstrip('chr')
    data['breakpoint2'] = data['RightBreakpoint'].str.rstrip('\:\+|\:\-').str.lstrip('chr')
    pattern = re.compile(r'\^.*')
    data['gene1'] = data.LeftGene.str.replace(pattern, '')
    data['gene2'] = data.RightGene.str.replace(pattern, '')
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['junction_reads'] = data['JunctionReadCount']
    data['spanning_reads'] = data['SpanningFragCount']
    data['caller'] = caller
    data['sample'] = sample

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output

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
            '-v {left_fq}:{left_fq}:ro',
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
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]), # FIXME allow custom outdir
        left_fq=left_fq,
        right_fq=right_fq,
        cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 8)
    )


@python_app
def parse_starseqr(out_dir, inputs=[]):
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'ss_STAR-SEQR/ss_STAR-SEQR_candidates.txt')
    sample = out_dir.split('/')[-2]
    caller = out_dir.split('/')[-1]
    data = pd.read_csv(path, sep='\t')

    data = data[data['DISPOSITION'] == 'PASS']
    data['gene1'] = data['LEFT_SYMBOL']
    data['gene2'] = data['RIGHT_SYMBOL']
    data['chromosome1'] = data.BRKPT_LEFT.str.extract(pat='(^\d.*)\:\d.*\:.*')
    data['chromosome2'] = data.BRKPT_RIGHT.str.extract(pat='(^\d.*)\:\d.*\:.*')
    # Transformation below is because starseqr is 0-indexed while the other callers are 1-indexed
    data['breakpoint1'] = data.BRKPT_LEFT.str.extract(pat='^\d.*\:(\d.*)\:.*').astype(float).astype('Int64') + 1
    data.breakpoint1 = data.chromosome1 + ':' + data.breakpoint1.astype(str)
    data['breakpoint2'] = data.BRKPT_RIGHT.str.extract(pat='^\d.*\:(\d.*)\:.*').astype(float).astype('Int64') + 1
    data.breakpoint2 = data.chromosome2 + ':' + data.breakpoint2.astype(str)
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['junction_reads'] = data.NREAD_JXNLEFT + data.NREAD_JXNRIGHT
    data['spanning_reads'] = data.NREAD_SPANS
    data['caller'] = caller
    data['sample'] = sample

    output = os.path.join(out_dir, 'fusions.pkl')
    data.to_pickle(output)

    return output

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
            '-v {left_fq}:{left_fq}:ro',
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

@python_app
def parse_fusioncatcher(out_dir, inputs=[]):
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'final-list_candidate-fusion-genes.txt')
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]

    data = pd.read_csv(path, sep='\t')
    data.rename(columns={'Gene_1_symbol(5end_fusion_partner)': 'gene1', 'Gene_2_symbol(3end_fusion_partner)': 'gene2'},
            inplace=True)
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['spanning_reads'] = data['Spanning_pairs']
    data['junction_reads'] = data['Spanning_unique_reads']

    data['caller'] = caller
    data['sample'] = sample

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output


@bash_app(cache=True)
def kallisto_index(
        genome_lib,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os

    if os.path.isfile(os.path.join(genome_lib, 'kallisto_index.idx')):
        return 'echo kallisto indexing completed'

    command = ['echo $HOSTNAME; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {genome_lib}:/genome_lib',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {genome_lib}:/genome_lib',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'kallisto index -i /genome_lib/kallisto_index.idx -k 31 /genome_lib/ref_annot.cdna.fa'
    ]

    return ' '.join(command).format(
        genome_lib=genome_lib,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )


@bash_app(cache=True)
def kallisto_quant(
        index,
        genome_lib,
        output,
        left_fq,
        right_fq,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {genome_lib}:/genome_lib',
            '-v {output}:/output ',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {genome_lib}:/genome_lib',
            '-B {output}:/output',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'kallisto quant',
        '-i /genome_lib/kallisto_index.idx',
        '--fusion',
        '-o /output',
        '{left_fq}',
        '{right_fq}'
    ]

    return ' '.join(command).format(
        output=output,
        left_fq=left_fq,
        right_fq=right_fq,
        genome_lib=genome_lib,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )


@bash_app(cache=True)
def run_pizzly(
        quant,
        genome_lib,
        output,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {genome_lib}:/genome_lib',
            '-v {output}:/output',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {genome_lib}:/genome_lib',
            '-B {output}:/output',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'pizzly ',
        '-k 31 ',
        '--gtf /genome_lib/ref_annot.gtf',
        '--cache /genome_lib/kallisto_index.cache.txt',
        '--align-score 2',
        '--insert-size 400',
        '--fasta /genome_lib/ref_annot.cdna.fa',
        '-o /output',
        '/output/fusion.txt'
    ]

    return ' '.join(command).format(
        output=output,
        genome_lib=genome_lib,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )
