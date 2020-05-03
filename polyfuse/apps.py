import parsl
from parsl.app.app import bash_app, python_app


@python_app
def assemble_data_per_sample(sample, callers, out_dir):
    import os
    import pandas as pd
    import numpy as np

    caller_data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'), 'data')
    caller_data = caller_data[caller_data['sample'] == sample]
    # if not set(caller_data.caller.unique()) == set(callers):
    #     return [], []

    true_fusions = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'data')
    true_fusions = true_fusions[true_fusions['sample'] == sample]
    fusions = list(set(np.concatenate((caller_data.fusion.unique(), true_fusions.fusion.unique()))))

    x = []
    y = []
    encoded_feature_info = [
        ('confidence', 'arriba_confidence', ('high', 'medium', 'low')),
        # ('reading_frame', 'arriba_reading_frame', ('out-of-frame', 'in-frame', '.')),
        ('LargeAnchorSupport', 'starfusion_large_anchor_support', ('YES_LDAS', 'NO_LDAS'))
    ]
    extra_features = [
        'FFPM', 'LeftBreakEntropy', 'RightBreakEntropy',
        'coverage1', 'coverage2'
    ]

    for fusion in fusions:
        row = []
        for c in callers:
            cut = (caller_data.fusion == fusion) & (caller_data.caller == c)
            data = caller_data.loc[cut, ['spanning_reads', 'junction_reads']]
            if len(data) > 0:
                row += [1] + data.values[0].tolist()
            else:
                row += [0, 0, 0]
                # cut = (caller_data.fusion == fusion)
                # data = caller_data.loc[cut, ['spanning_reads', 'junction_reads']]
                # if len(data) > 0:
                #     row += [0] + data.mean().values.tolist()
                # else: # this is a true fusion that no callers identified
                #     row += [0, 0, 0]
        for feature in [f for f, _, _ in encoded_feature_info] + extra_features:
            view = caller_data.loc[caller_data.fusion == fusion, feature]
            index = view.first_valid_index()
            row += [view.loc[index] if index is not None else 0]

        x += [row]
        y += [1 if any(true_fusions.fusion.isin([fusion])) else 0]

    columns = []
    for c in callers:
        columns += [c + '_called', c + '_spanning_reads', c + '_junction_reads']
    columns += [feature for feature, _, _ in encoded_feature_info]
    columns += extra_features
    x = pd.DataFrame(x, columns=columns)

    for feature, prefix, categories in encoded_feature_info:
        #  transformation below is required so that one-hot encoding will still
        #  add a column for categories, even if they do not appear in this sample
        x[feature] = x[feature].astype(pd.api.types.CategoricalDtype(categories))
        one_hot_encoded_columns = pd.get_dummies(x[feature], prefix=prefix, dummy_na=True)
        x = pd.concat([x, one_hot_encoded_columns], axis=1).drop([feature], axis=1)


    return x, y, fusions

def assemble_data(samples, callers, out_dir):
    import pandas as pd

    data = [assemble_data_per_sample(sample, callers, out_dir) for sample in samples]

    x = pd.concat([d.result()[0] for d in data])
    y = sum([d.result()[1] for d in data], [])

    return x, y

@python_app
def concatenate_true_fusions(sample_dirs, out_dir):
    import pandas as pd
    import glob
    import os

    true_fusions = pd.concat([pd.read_pickle(path) for path in glob.glob(os.path.join(sample_dirs, 'truth.pkl'))])
    # true_fusions['fusion'] = true_fusions[['left_gene_name', 'right_gene_name']].apply(lambda x: '--'.join(sorted([str(i) for i in x])), axis=1) # FIXME
    true_fusions['fusion'] = true_fusions[['left_gene_name', 'right_gene_name']].apply(lambda x: '--'.join([str(i) for i in x]), axis=1) # FIXME

    output = '{out_dir}/true_fusions.hdf'.format(out_dir=out_dir)
    true_fusions.to_hdf(output, 'data', mode='w')

    return output

def parse_caller_data(out_dir, callers):
    import os
    import glob

    parsers = dict((c, globals()['parse_{}'.format(c)]) for c in callers)
    futures = []
    for c in callers:
        futures += [parsers[c](d) for d in glob.glob(os.path.join(out_dir, '*', c))]

    return futures

@python_app
def predict_per_sample(data, sample, out_dir, classifier_label, features, transformation, callers):
    import os
    import pandas as pd
    import numpy as np
    import joblib
    from polyfuse import transformations

    x, _, fusions = data

    classifier = joblib.load(os.path.join(out_dir, 'models', '{}.joblib'.format(classifier_label)))
    probabilities = classifier.predict_proba(getattr(transformations, transformation)(x[features]))[:, 1]
    predictions = classifier.predict(getattr(transformations, transformation)(x[features]))

    return pd.DataFrame(data={
            'sample': [sample for i in range(len(fusions))],
            'caller': ['Polyfuse' + classifier_label for i in range(len(fusions))],
            'fusion': fusions,
            'probability': probabilities.tolist(),
            'prediction': predictions.tolist()
        }
    )


def predict(samples, out_dir, classifiers, callers):
    import pandas as pd
    import os

    futures = []
    for sample in samples:
        for features, label, transformation in classifiers:
            sample_data = assemble_data_per_sample(sample, callers, out_dir)
            futures += [predict_per_sample(
                sample_data,
                sample,
                out_dir,
                label,
                features,
                transformation,
                callers)
            ]
    model_data = pd.concat([f.result() for f in futures])
    path = os.path.join(out_dir, 'model_data.hdf')
    model_data.to_hdf(path, 'data', mode='w')

    return path

@python_app
def score_caller(out_dir, sample, caller):
    import pandas as pd
    import os
    import numpy as np

    truth = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'data')
    cut_truth = truth[truth['sample'] == sample]
    data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'))
    cut_data = data.loc[(data.caller == caller) & (data['sample'] == sample)]
    fusions = set(np.concatenate(
        (data[data['sample'] == sample].fusion.unique(), cut_truth.fusion.unique()))
    )

    y_true = [1 if any(cut_truth.fusion.isin([f])) else 0 for f in fusions]
    y_pred = [1 if any(cut_data.fusion.isin([f])) else 0 for f in fusions]
    y_prob = [cut_data.loc[cut_data.fusion == f, 'sum_J_S'].values[0] if any(cut_data.fusion.isin([f])) else 0 for f in fusions]

    return y_true, y_pred, y_prob

@python_app
def score_model(out_dir, sample, model):
    import pandas as pd
    import os
    import numpy as np

    data = pd.read_hdf(os.path.join(out_dir, 'model_data.hdf'))
    cut_data = data.loc[(data.caller == model) & (data['sample'] == sample)]
    truth = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'))
    cut_truth = truth[truth['sample'] == sample]
    fusions = set(np.concatenate((cut_data.fusion.unique(), cut_truth.fusion.unique())))

    y_true = [1 if any(cut_truth.fusion.isin([f])) else 0 for f in fusions]
    y_pred = [cut_data.loc[cut_data.fusion == f, 'prediction'].values[0] for f in fusions]
    y_prob = [cut_data.loc[cut_data.fusion == f, 'probability'].values[0] for f in fusions]

    return y_true, y_pred, y_prob

def make_summary(out_dir, samples):
    import os
    import pandas as pd
    import numpy as np
    from sklearn import metrics

    futures = []
    model_data = pd.read_hdf(os.path.join(out_dir, 'model_data.hdf'))
    for caller in model_data.caller.unique():
        for sample in samples:
            futures += [(score_model(out_dir, sample, caller), caller, sample)]

    caller_data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'))
    for caller in caller_data.caller.unique():
        for sample in samples:
            futures += [(score_caller(out_dir, sample, caller), caller, sample)]

    samples = []
    callers = []
    tps = []
    fps = []
    fns = []
    tns = [] # FIXME not sure TN is well-defined
    recalls = []
    precisions = []
    accuracies = []
    f1s = []
    mccs = []
    for f, caller, sample in futures:
        y_true, y_pred, _ = f.result()
        cm = metrics.confusion_matrix(y_true, y_pred)
        tn = cm[0][0]
        fp = cm[0][1]
        fn = cm[1][0]
        tp = cm[1][1]

        samples += [sample]
        callers += [caller]
        tps += [tp]
        fps += [fp]
        fns += [fn]
        tns += [tn]
        recalls += [metrics.recall_score(y_true, y_pred)]
        precisions += [metrics.precision_score(y_true, y_pred)]
        accuracies += [metrics.accuracy_score(y_true, y_pred)]
        f1s += [metrics.f1_score(y_true, y_pred)]
        mccs += [metrics.matthews_corrcoef(y_true, y_pred)]

    return pd.DataFrame(data={
        'sample': samples,
        'caller': callers,
        'tp': tps,
        'fp': fps,
        'fn': fns,
        'tn': tns,
        'recall': recalls,
        'precision': precisions,
        'accuracy': accuracies,
        'f1': f1s,
        'mcc': mccs
        }
    )


@python_app
def concatenate_caller_data(out_dir, inputs=[]):
    import pandas as pd
    import glob
    import os

    columns = [
        'spanning_reads',
        'junction_reads',
        'sample',
        'caller',
        'gene1',
        'gene2',
        'confidence',
        'reading_frame',
        'LargeAnchorSupport',
        'FFPM',
        'LeftBreakEntropy',
        'RightBreakEntropy',
        'coverage1',
        'coverage2',
        'fusion',
        'sum_J_S'
    ]
    caller_data = pd.concat(
        [
            pd.read_pickle(f) #[['spanning_reads', 'junction_reads', 'sample', 'caller', 'gene1', 'gene2']]
            for f in
            glob.glob(os.path.join(out_dir, '*', '*', 'fusions.pkl'))
        ],
        ignore_index=True,
        sort=False
    )
    # caller_data['fusion'] = caller_data[['gene1', 'gene2']].apply(lambda x: '--'.join(sorted(x)), axis=1) # FIXME
    caller_data['fusion'] = caller_data[['gene1', 'gene2']].apply(lambda x: '--'.join(x), axis=1) # FIXME
    caller_data['sum_J_S'] = caller_data['junction_reads'] + caller_data['spanning_reads']
    output = os.path.join(out_dir, 'caller_data.hdf')
    # caller_data = caller_data.fillna(0)
    caller_data[columns].to_hdf(output, 'data', mode='w')

    return output

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

    command = 'cd {output}; faidx --split-files ref_genome.fa'

    subprocess.check_output(
        command.format(
            output=output,
        ),
        shell=True
    )

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
        os.path.join(ref_split_by_chromosome_dir, 'chr*.fa')
    )

    command = []
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
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
            ref_split_by_chromosome_dir=ref_split_by_chromosome_dir,
            chromosome_refs=','.join(chromosome_refs),
            threads=os.environ.get('PARSL_CORES', int(multiprocessing.cpu_count() / 2))
        ),
        shell=True
    )

    return ref_split_by_chromosome_dir


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


    command = ['echo $HOSTNAME; mkdir -p {output}; ']
    if container_type == 'docker':
        command += [
            'docker run',
            '--rm',
            '-v {ctat_dir}:{ctat_dir}:ro',
            '-v {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            '-v {left_fq}:{left_fq}:ro',
            '-v {right_fq}:{right_fq}:ro',
            '-v {output}:{output}',
            'hiroko/mapsplice2-hg19'
        ]
    elif container_type == 'singularity':
        command += [
            'singularity exec',
            '-B {ctat_dir}:{ctat_dir}',
            '-B {ref_split_by_chromosome_dir}:{ref_split_by_chromosome_dir}',
            '-B {left_fq}:{left_fq}',
            '-B {right_fq}:{right_fq}',
            '-B {output}:{output}',
            '{base_dir}/docker/mapsplice2.sif'
        ]
    else:
        raise RuntimeError('Container type must be either docker or singularity')

    command += [
        'python mapsplice.py',
        '-c {ref_split_by_chromosome_dir}',
        '-x {ref_split_by_chromosome_dir}/ref_genome.fa',
        '-1 {left_fq}',
        '-2 {right_fq}',
        '--fusion',
        '--gene-gtf {ctat_dir}/ref_annot.gtf',
        '--output {output}',
        '--bam',
        '--threads {threads}'
    ]

    return ' '.join(command).format(
        output=output,
        ctat_dir=ctat_dir,
        ref_split_by_chromosome_dir=ref_split_by_chromosome_dir,
        left_fq=left_fq,
        right_fq=right_fq,
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
        threads=os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
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
            'docker run',
            '--rm',
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

    if len(glob.glob(fastq)) == 1:
        return glob.glob(fastq)[0]
    else:
        merged_fastq = '{out_dir}/interim/{sample}/merged.{tag}.fastq{ext}'.format(
            out_dir=out_dir,
            tag=tag,
            ext='.gz' if glob.glob(fastq)[0].endswith('.gz') else ''
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
            'docker run',
            '--rm',
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
        ctat_dir,
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
            'docker run',
            '--rm',
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
        star_index=star_index,
        output=output,
        ctat_dir=ctat_dir,
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
    data['spanning_reads'] = data.NREAD_SPANS
    data['caller'] = caller
    data['sample'] = sample

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

    if os.path.isfile(os.path.join(library_dir, ctat_release, 'download.success')):
        return os.path.join(library_dir, ctat_release, 'ctat_genome_lib_build_dir')

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
            'docker run',
            '--rm',
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
            'docker run',
            '--rm',
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
            'docker run',
            '--rm',
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
            'docker run',
            '--rm',
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
        output=output,
        gtf=os.path.abspath(gtf),
        fasta=os.path.abspath(fasta),
        base_dir='/'.join(os.path.abspath(__file__).split('/')[:-2]),
    )
