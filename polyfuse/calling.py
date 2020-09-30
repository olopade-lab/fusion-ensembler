# TODO remove docker pulls once all versions are pinned
# TODO documentation
import parsl
from parsl.app.app import bash_app, python_app


@python_app(cache=True)
def split_ref_chromosomes(
        ctat_dir,
        container_type='docker'
        ):
    # TODO cite Shirley MD, Ma Z, Pedersen BS, Wheelan SJ. (2015) Efficient
    # "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196
    # https://dx.doi.org/10.7287/peerj.preprints.970v1
    # FIXME Handle case where these files exist
    import os
    import subprocess
    import shutil

    output = os.path.join(ctat_dir, 'ref_genome.fa.chromosomes')
    os.makedirs(output, exist_ok=True)
    shutil.copy(os.path.join(ctat_dir, 'ref_genome.fa'), output)
    shutil.copy(os.path.join(ctat_dir, 'ref_genome.fa.fai'), output)

    command = 'cd {output}; faidx --split-files ref_genome.fa'

    subprocess.check_output(
        command.format(
            output=output,
        ),
        shell=True
    )
    os.unlink(os.path.join(output, 'ref_genome.fa'))
    os.unlink(os.path.join(output, 'ref_genome.fa.fai'))

    return output

@python_app(cache=True)
def build_bowtie_index(
        ref_split_by_chromosome_dir,
        container_type='docker'
        ):
    import os
    import subprocess
    import multiprocessing
    import glob

    if os.path.isfile(
        '{ref_split_by_chromosome_dir}/indexed.success'.format(
            ref_split_by_chromosome_dir=ref_split_by_chromosome_dir)
        ):
        return ref_split_by_chromosome_dir

    chromosome_refs = glob.glob(
        os.path.join(ref_split_by_chromosome_dir, '*.fa')
    )

    command = []
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'bowtie-build',
        '--threads {threads}',
        '{chromosome_refs}',
        '{ref_split_by_chromosome_dir}/ref_genome.fa;'
        'touch {ref_split_by_chromosome_dir}/indexed.success'
    ]

    subprocess.check_output(
            ' '.join(command).format(
            wrap_dir=os.path.dirname(os.path.abspath(__file__)),
            ref_split_by_chromosome_dir=ref_split_by_chromosome_dir,
            chromosome_refs=','.join(chromosome_refs),
            threads=os.environ.get('PARSL_CORES', int(multiprocessing.cpu_count() / 1.5))
        ),
        shell=True
    )

    return ref_split_by_chromosome_dir

# def assemble_command(
#         container_type,
#         container,
#         bind_dirs,
#         command,
#         prefix=None,
#         postfix=None):
#     if prefix is None:
#         prefix = []
#     if postfix is None:
#         postfix = []
#     command = '; '.join(prefix) + '; '
#     if container_type == 'docker':
#         command += 'docker run --rm'

#         command += [
#             'docker run',
#             '--rm',
#             '-v {ctat_dir}:{ctat_dir}:ro',
#             '-v {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
#             '-v {output}:{output}',
#             'hiroko/mapsplice2-hg19'
#         ]
#     elif container_type == 'singularity':
#         command += [
#             'singularity exec',
#             '-B {ctat_dir}:{ctat_dir}',
#             '-B {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
#             '-B {output}:{output}',
#             '{base_dir}/docker/mapsplice2.sif'
#         ]

@bash_app(cache=True)
def run_mapsplice2(
        output,
        ctat_dir,
        ref_split_by_chromosome_dir,
        left_fq,
        right_fq,
        container_type='docker',
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing

    # MapSplice does not support gzipped reads, and further does not accept stdin for the reads
    # That is why we cannot use a pipe to provide the unzipped reads
    command = [
        'echo $HOSTNAME;',
        'mkdir -p {output};',
        'gunzip -c {left_fq} > {output}/reads_1.fq;',
        'gunzip -c {right_fq} > {output}/reads_2.fq;',
    ]
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {ctat_dir}:{ctat_dir}:ro',
            '-v {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            '-v {output}:{output}',
            'hiroko/mapsplice2-hg19'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {ctat_dir}:{ctat_dir}',
            '-B {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            '-B {output}:{output}',
            '{base_dir}/docker/mapsplice2.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'python mapsplice.py',
        '-c {ref_split_by_chromosome_dir}',
        '-x {ref_split_by_chromosome_dir}/ref_genome.fa',
        '-1 {output}/reads_1.fq',
        '-2 {output}/reads_2.fq',
        '--fusion',
        '--gene-gtf {ctat_dir}/ref_annot.gtf',
        '--output {output}',
        '--bam',
        '--threads {threads};'
        'rm {output}/reads_1.fq {output}/reads_2.fq'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        output=output,
        ctat_dir=ctat_dir,
        ref_split_by_chromosome_dir=ref_split_by_chromosome_dir,
        left_fq=left_fq,
        right_fq=right_fq,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        threads=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 10)
        # TODO Prefer not to hardcoode this max, which is tweaked for IGSB config.
        # Need to switch to WorkQueue Parsl executor which will optimize resource packing.
    )


@python_app(cache=True)
def build_star_index(
        ctat_dir,
        index_name,
        container_type='docker'
        ):
    import os
    import subprocess
    import multiprocessing

    output = os.path.join(ctat_dir, index_name)

    command = ['mkdir -p {output};']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {output}:/output',
            '-v {ctat_dir}/ref_genome.fa:/assembly:ro',
            '-v {ctat_dir}/ref_annot.gtf:/annotation:ro',
            'eagenomics/starseqr:0.6.7'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {output}:/output',
            '-B {ctat_dir}/ref_genome.fa:/assembly',
            '-B {ctat_dir}/ref_annot.gtf:/annotation',
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
            wrap_dir=os.path.dirname(os.path.abspath(__file__)),
            ctat_dir=ctat_dir,
            output=output,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
            threads=os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
        ),
        shell=True
    )

    return output


@python_app(cache=True)
def merge_lanes(fastq, out_dir, sample, tag='R1'):
    import glob
    import subprocess
    import os

    if len(glob.glob(fastq)) == 1:
        return glob.glob(fastq)[0]
    else:
        merged_fastq = '{out_dir}/interim/{sample}/merged.{tag}.fastq{ext}'.format(
            out_dir=out_dir,
            tag=tag,
            ext='.gz' if glob.glob(fastq)[0].endswith('.gz') else '',
            sample=sample
        )
        subprocess.check_output(
            'mkdir -p {dirname}; cat {fastq} > {merged_fastq}'.format(
                dirname=os.path.dirname(merged_fastq),
                out_dir=out_dir,
                fastq=fastq,
                merged_fastq=merged_fastq
                ),
            shell=True
        )
        return merged_fastq

@python_app(cache=True)
def gzip(fastq, out_dir):
    import subprocess

    if fastq.endswith('.gz'):
        return fastq
    else:
        dirname = '{out_dir}/interim/{sample}'.format(
            out_dir=out_dir,
            sample=sample
        )
        subprocess.check_output('mkdir -p {dirname}; gzip {fastq}; mv {fastq}.gz {dirname}'.format(
            fastq=fastq,
            dirname=dirname), shell=True
        )
        return os.path.join(dirname, os.path.basename(fastq) + '.gz')

@bash_app(cache=True)
def run_arriba(
        output,
        ctat_dir,
        left_fq,
        right_fq,
        container_type='docker',
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing


    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {ctat_dir}/ref_genome.fa:/assembly:ro',
            '-v {ctat_dir}/ref_annot.gtf:/annotation:ro',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {ctat_dir}/ref_genome.fa.star.idx:/star_index:ro',
            '-v {output}:/output',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {ctat_dir}/ref_genome.fa:/assembly',
            '-B {ctat_dir}/ref_annot.gtf:/annotation',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {ctat_dir}/ref_genome.fa.star.idx:/star_index',
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
            wrap_dir=os.path.dirname(os.path.abspath(__file__)),
            output=output,
            ctat_dir=ctat_dir,
            left_fq=left_fq,
            right_fq=right_fq,
            base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
            cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 8)
        )


@python_app
def parse_pizzly(out_dir, inputs=[]):
    import os
    import pandas as pd
    import json

    path = os.path.join(out_dir, 'out.json')
    if not os.path.exists(path):
        return None
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    with open(path, 'r') as f:
        json_data = json.load(f)['genes']

    data = pd.DataFrame(columns=[])
    data['gene1'] = [gf['geneA']['name'] for gf in json_data]
    data['gene2'] = [gf['geneB']['name'] for gf in json_data]
    data['junction_reads'] = [gf['paircount'] for gf in json_data]
    data['spanning_reads'] = [gf['splitcount'] for gf in json_data]
    data['caller'] = caller
    data['sample'] = sample

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'pizzly_' + c for c in data.columns if c not in common_columns}, inplace=True)

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output

@python_app
def parse_arriba(out_dir, inputs=[]):
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'fusions.tsv')
    if not os.path.exists(path):
        return None
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    data = pd.read_csv(path, sep='\t')
    # data = data[data.confidence.str.contains('medium|high')]
    data.rename(columns={
        '#gene1': 'gene1',
        'discordant_mates': 'spanning_reads',
        },
        inplace=True
    )
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['junction_reads'] = data.split_reads1 + data.split_reads2
    data['caller'] = caller
    data['sample'] = sample

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'arriba_' + c for c in data.columns if c not in common_columns}, inplace=True)

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output

@python_app
def parse_mapsplice2(out_dir, inputs=[]):
    # http://www.netlab.uky.edu/p/bioinfo/MapSplice2FusionJunctionFormat
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'fusions_well_annotated.txt')
    # path = os.path.join(out_dir, 'fusions_candidates.txt')
    if not os.path.exists(path):
        return None
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    columns = [
        'chrom', 'doner_end', 'acceptor_start', 'id', 'coverage', 'strand', 'rgb', 'block_count', 'block_size',
        'block_distance', 'entropy', 'flank_case', 'flank_string', 'min_mismatch', 'max_mismatch', 'ave_mismatch',
        'max_min_suffix', 'max_min_prefix', 'min_anchor_difference', 'unique_read_count', 'multi_read_count',
        'paired_read_count', 'left_paired_read_count', 'right_paired_read_count', 'multiple_paired_read_count',
        'unique_paired_read_count', 'single_read_count', 'encompassing_read_pair_count', 'doner_start',
        'acceptor_end', 'doner_iosforms', 'acceptor_isoforms', 'obsolete1', 'obsolete2', 'obsolete3', 'obsolete4',
        'minimal_doner_isoform_length', 'maximal_doner_isoform_length', 'minimal_acceptor_isoform_length',
        'maximal_acceptor_isoform_length', 'paired_reads_entropy', 'mismatch_per_bp', 'anchor_score',
        'max_doner_fragment', 'max_acceptor_fragment', 'max_cur_fragment', 'min_cur_fragment', 'ave_cur_fragment',
        'doner_encompass_unique', 'doner_encompass_multiple', 'acceptor_encompass_unique',
        'acceptor_encompass_multiple', 'doner_match_to_normal', 'acceptor_match_to_normal', 'doner_seq',
        'acceptor_seq', 'match_gene_strand', 'annotated_type', 'fusion_type', 'gene_strand', 'annotated_gene_donor',
        'annotated_gene_acceptor'
    ]
    data = pd.read_csv(path, sep='\t', names=columns, index_col=False)
    data['gene1'] = [d.strip(',') for d in data['annotated_gene_donor']]
    data['gene2'] = [d.strip(',') for d in data['annotated_gene_acceptor']]
    data['caller'] = caller
    data['sample'] = sample
    data.rename(columns={
        'coverage': 'junction_reads', # coverage: number of reads aligned to the fusion junction
        'encompassing_read_pair_count': 'spanning_reads',  # encompassing_read_pair_count: number of reads pairs surrounding (but not crossing) the fusion
        },
        inplace=True
    )

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'mapsplice2_' + c for c in data.columns if c not in common_columns}, inplace=True)

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
        ctat_dir,
        container_type='docker',
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing


    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {ctat_dir}:/ctat_dir:ro',
            '-v {output}:/output',
            'trinityctat/starfusion:1.8.0'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {ctat_dir}:/ctat_dir',
            '-B {output}:/output',
            '{base_dir}/docker/starfusion.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        '/usr/local/src/STAR-Fusion/STAR-Fusion',
        '--left_fq {left_fq}',
        '--right_fq {right_fq}',
        '--genome_lib_dir /ctat_dir',
        '-O /output',
        # '--FusionInspector validate',
        # '--examine_coding_effect',
        # '--denovo_reconstruct',
        '--CPU {cores}'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        left_fq=left_fq,
        right_fq=right_fq,
        ctat_dir=ctat_dir,
        output=output,
        cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 8)
    )

@python_app
def parse_starfusion(out_dir, inputs=[]):
    import os
    import re
    import pandas as pd

    path = os.path.join(out_dir, 'star-fusion.fusion_predictions.tsv')
    if not os.path.exists(path):
        return None

    sample = path.split('/')[-3]
    caller = path.split('/')[-2]
    data = pd.read_csv(path, sep='\t')

    data['breakpoint1'] = data['LeftBreakpoint'].str.rstrip('\:\+|\:\-').str.lstrip('chr')
    data['breakpoint2'] = data['RightBreakpoint'].str.rstrip('\:\+|\:\-').str.lstrip('chr')
    pattern = re.compile(r'\^.*')
    data['gene1'] = data.LeftGene.str.replace(pattern, '')
    data['gene2'] = data.RightGene.str.replace(pattern, '')
    data['caller'] = caller
    data['sample'] = sample
    data.rename(columns={
        'JunctionReadCount': 'junction_reads',
        'SpanningFragCount': 'spanning_reads'
        },
        inplace=True
    )

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'starfusion_' + c for c in data.columns if c not in common_columns}, inplace=True)

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output

@bash_app(cache=True)
def run_starseqr(
        output,
        left_fq,
        right_fq,
        ctat_dir,
        star_index,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os
    import multiprocessing

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {ctat_dir}:/ctat_dir:ro',
            '-v {star_index}:/star_index:ro',
            '-v {output}:/output ',
            'eagenomics/starseqr:0.6.7'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {ctat_dir}:/ctat_dir',
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
        '-g /ctat_dir/ref_annot.gtf',
        '-r /ctat_dir/ref_genome.fa',
        '-m 1',
        '-vv',
        '-t {cores}'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        star_index=star_index,
        output=output,
        ctat_dir=ctat_dir,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]), # FIXME allow custom outdir
        left_fq=left_fq,
        right_fq=right_fq,
        # too small of number of cores will give: [22167 rows x 6 columns]]'. Reason: 'IOError('bad message length',)'
        cores=min(os.environ.get('PARSL_CORES', multiprocessing.cpu_count()), 16)
    )


@python_app
def parse_starseqr(out_dir, inputs=[]):
    import os
    import pandas as pd

    path = os.path.join(out_dir, 'ss_STAR-SEQR/ss_STAR-SEQR_candidates.txt')
    if not os.path.exists(path):
        return None
    sample = out_dir.split('/')[-2]
    caller = out_dir.split('/')[-1]
    data = pd.read_csv(path, sep='\t')

    # data = data[data['DISPOSITION'] == 'PASS']
    data['gene1'] = data['LEFT_SYMBOL']
    data['gene2'] = data['RIGHT_SYMBOL']
    data['chromosome1'] = data.BRKPT_LEFT.str.extract(pat='(^\d.*)\:\d.*\:.*')
    data['chromosome2'] = data.BRKPT_RIGHT.str.extract(pat='(^\d.*)\:\d.*\:.*')
    # Transformation below is because starseqr is 0-indexed while the other callers are 1-indexed
    data['breakpoint1'] = data.BRKPT_LEFT.str.extract(pat='^\d.*\:(\d.*)\:.*').astype(float).astype('Int64') + 1
    data.breakpoint1 = data.chromosome1 + ':' + data.breakpoint1.astype(str)
    data['breakpoint2'] = data.BRKPT_RIGHT.str.extract(pat='^\d.*\:(\d.*)\:.*').astype(float).astype('Int64') + 1
    data.breakpoint2 = data.chromosome2 + ':' + data.breakpoint2.astype(str)
    # data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['junction_reads'] = data.NREAD_JXNLEFT + data.NREAD_JXNRIGHT
    data.rename(columns={'NREAD_SPANS': 'spanning_reads'}, inplace=True)
    data['caller'] = caller
    data['sample'] = sample

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'starseqr_' + c for c in data.columns if c not in common_columns}, inplace=True)

    output = os.path.join(out_dir, 'fusions.pkl')
    data.to_pickle(output)

    return output

@python_app(cache=True)
def download_ensemble_annotation(library_dir, release):
    import os
    import subprocess

    annotation = 'Homo_sapiens.GRCh38.{release}.gtf.gz'.format(release=release)

    if os.path.isfile(os.path.join(library_dir, 'ensemble', 'download.annotation.success')):
        return os.path.join(library_dir, 'ensemble', annotation)

    command = """
    mkdir -p {library_dir}/ensemble
    cd {library_dir}/ensemble
    wget ftp://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/{annotation}
    touch download.annotation.success
    cd -
    """

    subprocess.call(command.format(library_dir=library_dir, annotation=annotation, release=release), shell=True)

    return os.path.join(library_dir, 'ensemble', annotation)

@python_app(cache=True)
def download_ensemble_assembly(library_dir, release):
    import os
    import subprocess

    assembly = 'Homo_sapiens.GRCh38.cdna.all.fa.gz'

    if os.path.isfile(os.path.join(library_dir, 'ensemble', 'download.assembly.success')):
        return os.path.join(library_dir, 'ensemble', assembly)

    command = """
    mkdir -p {library_dir}/ensemble
    cd {library_dir}/ensemble
    wget ftp://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/cdna/{assembly}
    touch download.assembly.success
    cd -
    """

    subprocess.call(command.format(library_dir=library_dir, assembly=assembly, release=release), shell=True)

    return os.path.join(library_dir, 'ensemble', assembly)

@python_app(cache=True)
def download_ctat(library_dir, ctat_release):
    import os
    import subprocess

    if os.path.isfile(os.path.join(library_dir, ctat_release + '.plug-n-play', 'download.success')):
        return os.path.join(library_dir, ctat_release + '.plug-n-play', 'ctat_genome_lib_build_dir')

    command = """
    cd {library_dir}
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/{release}.plug-n-play.tar.gz
    tar xz {release}.plug-n-play.tar.gz
    touch {release}.plug-n-play/download.success
    cd -
    """
    subprocess.call(command.format(library_dir=library_dir, release=ctat_release), shell=True)

    return os.path.join(library_dir, ctat_release + '.plug-n-play', 'ctat_genome_lib_build_dir')

@python_app(cache=True)
def download_fusioncatcher_build(output):
    import os
    import subprocess

    if os.path.isfile(os.path.join(output, 'download.success')):
        return output

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
    subprocess.call(command.format(output=output), shell=True)

    return output

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
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {build_dir}/current:/build_dir:ro',
            '-v {output}:/output ',
            'olopadelab/fusioncatcher:latest'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {build_dir}/current:/build_dir',
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
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
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
    if not os.path.exists(path):
        return None
    sample = path.split('/')[-3]
    caller = path.split('/')[-2]

    data = pd.read_csv(path, sep='\t')
    data.rename(columns={
        'Gene_1_symbol(5end_fusion_partner)': 'gene1',
        'Gene_2_symbol(3end_fusion_partner)': 'gene2',
        'Spanning_pairs': 'spanning_reads',
        'Spanning_unique_reads': 'junction_reads'
        },
        inplace=True
    )
    data['fusion'] = data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1)
    data['caller'] = caller
    data['sample'] = sample

    common_columns = ['sample', 'caller', 'gene1', 'gene2', 'junction_reads', 'spanning_reads']
    data.rename(columns={c: 'fusioncatcher_' + c for c in data.columns if c not in common_columns}, inplace=True)

    output = os.path.join(os.path.dirname(path), 'fusions.pkl')
    data.to_pickle(output)

    return output


@bash_app(cache=True)
def kallisto_index(
        fasta_path,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os

    # if os.path.isfile(os.path.join(genome_lib, 'kallisto_index.idx')):
    #     return 'echo kallisto indexing completed'

    command = ['echo $HOSTNAME; ']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {fasta_dir}:/fasta_dir',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {fasta_dir}:/fasta_dir',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'kallisto index -i /fasta_dir/kallisto_index.idx -k 31 /fasta_dir/{fasta}'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        fasta_dir=os.path.dirname(fasta_path),
        fasta=os.path.basename(fasta_path),
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )


@bash_app(cache=True)
def kallisto_quant(
        index,
        fasta_path,
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
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {fasta_dir}:/fasta_dir',
            '-v {output}:/output ',
            'olopadelab/polyfuse'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {fasta_dir}:/fasta_dir',
            '-B {output}:/output',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'kallisto quant',
        '-i /fasta_dir/kallisto_index.idx',
        '--fusion',
        '-o /output',
        '{left_fq}',
        '{right_fq}'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        output=output,
        left_fq=left_fq,
        right_fq=right_fq,
        fasta_dir=os.path.dirname(fasta_path),
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )


@bash_app(cache=True)
def run_pizzly(
        quant,
        gtf,
        fasta,
        output,
        container_type,
        stderr=parsl.AUTO_LOGNAME,
        stdout=parsl.AUTO_LOGNAME):
    import os

    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'bash {wrap_dir}/wrap_docker.sh',
            '-v {gtf}:{gtf}',
            '-v {fasta}:{fasta}',
            '-v {output}:/output',
            'olopadelab/polyfuse:latest'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {gtf}:{gtf}',
            '-B {fasta}:{fasta}',
            '-B {output}:/output',
            '{base_dir}/docker/polyfuse.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        '/bin/bash -c "',
        'cd /output;',
        'pizzly ',
        '-k 31 ',
        '--gtf {gtf}',
        # '--cache /genome_lib/kallisto_index.cache.txt',
        '--align-score 2',
        '--insert-size 400',
        '--fasta {fasta}',
        '-o out',
        '/output/fusion.txt"'
    ]

    return ' '.join(command).format(
        wrap_dir=os.path.dirname(os.path.abspath(__file__)),
        output=output,
        gtf=os.path.abspath(gtf),
        fasta=os.path.abspath(fasta),
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )
