# TODO df['sample'] -> df.sample_id
# TODO documentation
# TODO split up assemble_data_per_sample to get rid of 'assemble_truth' keyword arg
import parsl
from parsl.app.app import bash_app, python_app

from polyfuse.calling import parse_pizzly, parse_starseqr, parse_fusioncatcher, parse_mapsplice2, parse_starfusion, parse_arriba

@python_app
def assemble_data_per_sample(sample, callers, out_dir, encoded_features=None, extra_features=None, assemble_truth=True):
    import os
    import pandas as pd
    import numpy as np

    all_caller_data = pd.read_hdf(os.path.join(out_dir, 'all_caller_data.hdf'), 'data')
    # caller_data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'), 'data')
    sample_data = all_caller_data[all_caller_data['sample'] == sample]
    if not set(callers).issubset(set(sample_data.caller.unique())):
        for missing_caller in set(callers) - set(sample_data.caller.unique()):
            # if the pickle file is there, assume job finished successfully but no fusions were found
            if not os.path.isfile(os.path.join(out_dir, sample, missing_caller, 'fusions.pkl')):
                raise RuntimeError('problem processing sample {}: missing {}'.format(sample, missing_caller))

    if assemble_truth:
        true_fusions = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'data')
        true_fusions = true_fusions[true_fusions['sample'] == sample]
        fusions = list(set(np.concatenate((sample_data.fusion.unique(), true_fusions.fusion.unique()))))
    else:
        fusions = sample_data.fusion.unique().tolist()

    x = []
    y = []
    if encoded_features is None:
        encoded_features = [
            'arriba_confidence',
            'starfusion_LargeAnchorSupport',
            'fusioncatcher_Fusion_finding_method',
            'starseqr_FUSION_CLASS',
            'starseqr_SPLICE_TYPE',
            'mapsplice2_doner_match_to_normal',
            'starfusion_SpliceType',
            'starfusion_LeftBreakDinuc', # left break dinucleotide
            'starfusion_RightBreakDinuc',
            # https://arriba.readthedocs.io/en/latest/output-files/
            # arriba_site1/arriba_site2: These columns add information about the location of the
            # breakpoints. Possible values are: splice-site (at an exon boundary and oriented in way that the transcript
            # has likely been spliced), exon (inside an exon, but not at an exon boundary), intron, 5' UTR, 3' UTR, UTR
            # (overlapping with a 5' UTR as well as a 3' UTR), and intergenic.
            'arriba_site1',
            'arriba_site2',
            # https://arriba.readthedocs.io/en/latest/output-files/
            # arriba_type: Based on the orientation of the supporting reads and the coordinates of breakpoints, the type of
            # event can be inferred. Possible values are: translocation (between different chromosomes), duplication,
            # inversion, and deletion. If genes are fused head-to-head or tail-to-tail, this is indicated as 5'-5' or
            # 3'-3' respectively. Genes fused in such an orientation cannot yield a chimeric protein, since one of the
            # genes is transcribed from the wrong strand. This type of event is equivalent to the truncation of the genes.
            # Deletions with a size in the range of introns (<400kb) are flagged as read-through, because there is a high
            # chance that the fusion arises from read-through transcription rather than an underlying genomic deletion.
            # Intragenic duplications with both breakpoints at splice-sites are flagged as non-canonical-splicing, because
            # the supporting reads might originate from circular RNA, which are very abundant even in normal tissue, but
            # manifest as duplications in RNA-Seq data.
            'arriba_type'
        ]
    if extra_features is None:
        extra_features = [
            # https://github.com/STAR-Fusion/STAR-Fusion/wiki
            # The number of fusion-supporting reads depends on both the expression of the fusion transcript and the
            # number of reads sequenced. The deeper the sequenced data set, the greater the number of artifactual
            # fusions that will appear with minimal supporting evidence, and so taking into account the sequencing depth
            # is important to curtail overzealous prediction of fusion transcripts with ever-so-minimal supporting
            # evidence. We provide normalized measures of the fusion-supporting rna-seq fragments as FFPM (fusion
            # fragments per million total reads) measures.
            'starfusion_FFPM',
            # https://github.com/STAR-Fusion/STAR-Fusion/wiki
            # 'LeftBreakEntropy' and 'RightBreakEntropy' represent the Shannon entropy of the 15 exonic bases flanking
            # the breakpoint. The maximum entropy is 2, representing highest complexity. The lowest would be zero
            # (involving a 15 base mononucleotide run). Low entropy sites should generally be treated as less confident
            # breakpoints.
            'starfusion_LeftBreakEntropy',
            'starfusion_RightBreakEntropy',
            'arriba_coverage1',
            'arriba_coverage2',
            'starseqr_OVERHANG_BQ15', # Number of overhang fragmens with at least 15 base pairs
            'starseqr_TPM_FUSION', # Expression of the most abundant fusion transcript expressed in transcripts per million
            'starseqr_TPM_LEFT',
            'starseqr_TPM_RIGHT',
            'sum_J_S'
        ]

    for fusion in fusions:
        row = sample_data.loc[sample_data.fusion == fusion, ['spanning_reads', 'junction_reads']].mean().values.tolist()
        for c in callers:
            cut = (sample_data.fusion == fusion) & (sample_data.caller == c)
            data = sample_data.loc[cut, ['spanning_reads', 'junction_reads']]
            if len(data) > 0:
                row += [1] + data.values[0].tolist()
            else:
                row += [0, 0, 0]
                # Below was an attempt to impute missing data-- did not improve performance
                # cut = (sample_data.fusion == fusion)
                # data = sample_data.loc[cut, ['spanning_reads', 'junction_reads']]
                # if len(data) > 0:
                #     row += [0] + data.mean().values.tolist()
                # else: # this is a true fusion that no callers identified
                #     row += [0, 0, 0]
        for feature in encoded_features + extra_features:
            view = sample_data.loc[sample_data.fusion == fusion, feature]
            index = view.first_valid_index()
            row += [view.loc[index] if index is not None else 0]

        x += [row]
        if assemble_truth:
            y += [1 if any(true_fusions.fusion.isin([fusion])) else 0]

    columns = ['mean_spanning_reads', 'mean_junction_reads']
    for c in callers:
        columns += [c + '_called', c + '_spanning_reads', c + '_junction_reads']
    columns += encoded_features + extra_features
    x = pd.DataFrame(x, columns=columns)
    x = x.fillna(0) # fusions that only appear in truth set will have NaN means

    for feature in encoded_features:
        #  transformation below is required so that one-hot encoding will still
        #  add a column for categories, even if they do not appear in this sample
        categories = [c for c in all_caller_data[feature].unique() if str(c) != 'nan']
        x[feature] = x[feature].astype(pd.api.types.CategoricalDtype(categories))
        one_hot_encoded_columns = pd.get_dummies(x[feature], prefix=feature, dummy_na=True)
        x = pd.concat([x, one_hot_encoded_columns], axis=1).drop([feature], axis=1)

    sample_out_dir = os.path.join(out_dir, 'sample_data')
    os.makedirs(sample_out_dir, exist_ok=True)

    return x, y, fusions

def assemble_data(
        samples,
        callers,
        out_dir,
        encoded_features=None,
        extra_features=None,
        assemble_truth=True,
        tag=''):
    import os
    import pandas as pd

    data = [
        assemble_data_per_sample(
            sample,
            callers,
            out_dir,
            encoded_features,
            extra_features,
            assemble_truth=assemble_truth
        )
        for sample in samples
    ]

    x = pd.concat([d.result()[0] for d in data if d.result() is not None])
    y = sum([d.result()[1] for d in data if d.result() is not None], [])

    x = x.fillna(0)

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

def predict_consensus(samples, out_dir, callers, quorums):
    import pandas as pd
    import os
    from polyfuse.utils import get_consensus

    data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'))
    futures = []
    for sample in samples:
        cut_data = data.loc[data['sample'] == sample]

        callsets = []
        for caller in callers:
            callsets.append(set(cut_data[cut_data.caller == caller].fusion))
        for quorum in quorums:
            futures += [(get_consensus(callsets, quorum), sample, quorum)]

    samples = []
    callers = []
    fusions = []
    for f, sample, quorum in futures:
        consensus_fusions = f.result()
        samples += [sample for i in range(len(consensus_fusions))]
        callers += ['ConsensusQ{}'.format(quorum) for i in range(len(consensus_fusions))]
        fusions += list(consensus_fusions)

    consensus_data = pd.DataFrame(data={
            'sample': samples,
            'caller': callers,
            'fusion': fusions,
        }
    )

    path = os.path.join(out_dir, 'consensus_data.hdf')
    consensus_data.to_hdf(path, 'data', mode='w')

    return path

@python_app
def predict_per_sample(data, sample, out_dir, model_dir, classifier_label, features, transformation, callers, consensus=None):
    import os
    import pandas as pd
    import numpy as np
    import joblib
    from polyfuse import transformations

    if data is None:
        return None
    x, _, fusions = data

    classifier = joblib.load(os.path.join(model_dir, '{}.joblib'.format(classifier_label)))
    probabilities = classifier.predict_proba(getattr(transformations, transformation)(x[features]))[:, 1]
    predictions = classifier.predict(getattr(transformations, transformation)(x[features]))

    label = 'Polyfuse' + classifier_label
    if consensus is not None:
        label = consensus + label
        consensus_data = pd.read_hdf(os.path.join(out_dir, 'consensus_data.hdf'))
        cut_data = consensus_data.loc[(consensus_data['sample'] == sample) & (consensus_data.caller == consensus)]
        consensus_predictions = [1 if any(cut_data.fusion.isin([f])) else 0 for f in fusions]
        predictions = predictions | consensus_predictions


    return pd.DataFrame(data={
            'sample': [sample for i in range(len(fusions))],
            'caller': [label for i in range(len(fusions))],
            'fusion': fusions,
            'probability': probabilities.tolist(),
            'prediction': predictions.tolist()
        }
    )


def predict(samples, out_dir, model_dir, classifiers, callers, consensus=None, assemble_truth=True):
    import pandas as pd
    import os

    futures = []
    for sample in samples:
        for features, label, transformation in classifiers:
            sample_data = assemble_data_per_sample(sample, callers, out_dir, assemble_truth=assemble_truth)
            futures += [predict_per_sample(
                sample_data,
                sample,
                out_dir,
                model_dir,
                label,
                features,
                transformation,
                callers)
            ]
            if consensus is not None:
                for consensus_caller in consensus:
                    futures += [predict_per_sample(
                        sample_data,
                        sample,
                        out_dir,
                        model_dir,
                        label,
                        features,
                        transformation,
                        callers,
                        consensus=consensus_caller)
                    ]
    model_data = pd.concat([f.result() for f in futures if f.result() is not None])
    path = os.path.join(out_dir, 'model_data.hdf')
    model_data.to_hdf(path, 'data', mode='w')

    return path

@python_app
def score_consensus(out_dir, sample, caller):
    import pandas as pd
    import os
    import numpy as np

    truth = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'data')
    cut_truth = truth[truth['sample'] == sample]
    data = pd.read_hdf(os.path.join(out_dir, 'consensus_data.hdf'))
    cut_data = data.loc[(data.caller == caller) & (data['sample'] == sample)]
    if len(cut_data) == 0:
        return None
    fusions = set(np.concatenate(
        (data[data['sample'] == sample].fusion.unique(), cut_truth.fusion.unique()))
    )

    y_true = [1 if any(cut_truth.fusion.isin([f])) else 0 for f in fusions]
    y_pred = [1 if any(cut_data.fusion.isin([f])) else 0 for f in fusions]
    y_prob = [np.nan for f in fusions]

    return y_true, y_pred, y_prob

@python_app
def score_caller(out_dir, sample, caller):
    import pandas as pd
    import os
    import numpy as np

    truth = pd.read_hdf(os.path.join(out_dir, 'true_fusions.hdf'), 'data')
    cut_truth = truth[truth['sample'] == sample]
    data = pd.read_hdf(os.path.join(out_dir, 'caller_data.hdf'))
    cut_data = data.loc[(data.caller == caller) & (data['sample'] == sample)]
    if len(cut_data) == 0:
        return None
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
    if len(cut_data) == 0:
        return None
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

    consensus_data = pd.read_hdf(os.path.join(out_dir, 'consensus_data.hdf'))
    for caller in consensus_data.caller.unique():
        for sample in samples:
            futures += [(score_consensus(out_dir, sample, caller), caller, sample)]

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
        if f.result() is None:
            continue
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
        'arriba_confidence',
        'arriba_reading_frame',
        'starfusion_LargeAnchorSupport',
        'starfusion_FFPM',
        'starfusion_LeftBreakEntropy',
        'starfusion_RightBreakEntropy',
        'arriba_coverage1',
        'arriba_coverage2',
        'fusion',
        'sum_J_S'
    ]
    caller_data = pd.concat(
        [
            pd.read_pickle(f)
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

    output = os.path.join(out_dir, 'all_caller_data.hdf')
    caller_data.to_hdf(output, 'data', mode='w')

    output = os.path.join(out_dir, 'caller_data.hdf')
    caller_data[columns].to_hdf(output, 'data', mode='w')

    return output
