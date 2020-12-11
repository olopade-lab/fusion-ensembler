# took 4.5 hours to run
import os
import glob
import re
import itertools
import time
import joblib
import random
import itertools

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np
import seaborn as sns
sns.set_context('paper', font_scale=1.2)
import matplotlib.pyplot as plt
from matplotlib import container
from matplotlib.lines import Line2D

from sklearn.feature_selection import RFECV
from sklearn import svm, ensemble, datasets

import parsl
parsl.clear()

# from polyfuse.configs.local import config
from polyfuse.configs.igsb_jupyter import config
parsl.load(config)


from polyfuse import modeling, transformations

out_dir = '/cephfs/users/annawoodard/polyfuse/data/sim_50/processed'
training_fraction = 0.85
callers = ['starseqr', 'starfusion', 'arriba', 'fusioncatcher', 'pizzly', 'mapsplice2']

print('calculating recursive feature eliminatioin with cross-fold validation')

parsed_caller_data = modeling.parse_caller_data(out_dir, callers)

caller_data_path = modeling.concatenate_caller_data(out_dir, inputs=parsed_caller_data)
caller_data = pd.read_hdf(caller_data_path.result(), 'data')

samples = sorted(caller_data['sample'].unique())
training_samples = samples[:int(len(samples) * training_fraction)]
training_samples = samples[:30] # FIXME
print('will use {} samples'.format(len(training_samples)))

encoded_features = [
    'fusioncatcher_Fusion_finding_method',
    # 'fusioncatcher_Predicted_effect', # 32 features, probably not worth checking
    'starseqr_FUSION_CLASS',
    'starseqr_SPLICE_TYPE',
    'starseqr_ASSEMBLY_CROSS_JXN',
    # http://www.netlab.uky.edu/p/bioinfo/MapSplice2FusionJunctionFormat
    'mapsplice2_strand',
    'mapsplice2_flank_case', # non-zero for canonical and semi-canonical junctions
    'mapsplice2_flank_string', # the two basepairs after doner site combined the two basepairs before acceptor site
    'mapsplice2_doner_match_to_normal',
    'mapsplice2_acceptor_match_to_normal',
    # 'starfusion_SpliceType', # FIXME testing
    'starfusion_LargeAnchorSupport',
    'starfusion_LeftBreakDinuc',
    'starfusion_RightBreakDinuc',
    'arriba_strand1(gene/fusion)',
    'arriba_strand2(gene/fusion)',
    'arriba_site1',
    'arriba_site2',
    'arriba_type',
    'arriba_direction1',
    'arriba_direction2',
    'arriba_confidence',
    'arriba_reading_frame',
]

extra_features = [
    'spanning_reads',
    'junction_reads',
    'fusioncatcher_Counts_of_common_mapping_reads',
    'fusioncatcher_Longest_anchor_found',
    'starseqr_SPAN_CROSSHOM_SCORE',
    'starseqr_MINFRAG35',
    'mapsplice2_min_mismatch',
    'mapsplice2_max_mismatch',
    'mapsplice2_max_min_suffix',
    'mapsplice2_max_min_prefix',
    'mapsplice2_min_anchor_difference',
    'starseqr_DISTANCE',
    'starseqr_JXN_CROSSHOM_SCORE',
    'starseqr_OVERHANG_DIVERSITY',
    'starseqr_MINFRAG20',
    'starseqr_OVERHANG_MEANBQ',
    'starseqr_SPAN_MEANBQ',
    'starseqr_JXN_MEANBQ',
    'starseqr_OVERHANG_BQ15',
    'starseqr_SPAN_BQ15',
    'starseqr_JXN_BQ15',
    'starseqr_OVERHANG_MM',
    'starseqr_SPAN_MM',
    'starseqr_JXN_MM',
    'starseqr_OVERHANG_MEANLEN',
    'starseqr_SPAN_MEANLEN',
    'starseqr_JXN_MEANLEN',
    'starseqr_TPM_FUSION',
    'starseqr_TPM_LEFT',
    'starseqr_TPM_RIGHT',
    'mapsplice2_entropy',
    'mapsplice2_ave_mismatch',
    'mapsplice2_unique_read_count',
    'mapsplice2_paired_read_count',
    'mapsplice2_left_paired_read_count',
    'mapsplice2_right_paired_read_count',
    'mapsplice2_unique_paired_read_count',
    'mapsplice2_single_read_count',
    'mapsplice2_minimal_doner_isoform_length',
    'mapsplice2_maximal_acceptor_isoform_length',
    'mapsplice2_paired_reads_entropy',
    'mapsplice2_mismatch_per_bp',
    'mapsplice2_max_doner_fragment',
    'mapsplice2_max_cur_fragment',
    'mapsplice2_min_cur_fragment',
    'mapsplice2_ave_cur_fragment',
    'mapsplice2_doner_encompass_unique',
    'mapsplice2_doner_encompass_multiple',
    'mapsplice2_acceptor_encompass_unique',
    'mapsplice2_acceptor_encompass_multiple',
    'starfusion_FFPM',
    'starfusion_LeftBreakEntropy',
    'starfusion_RightBreakEntropy',
    'arriba_split_reads1',
    'arriba_split_reads2',
    'arriba_coverage1',
    'arriba_coverage2',
    # 'arriba_filters' # TODO regex to clean this up and pull out values
    'sum_J_S'
]

start = time.time()
x_train, y_train = modeling.extract_features(
    training_samples,
    callers,
    out_dir,
    encoded_features=encoded_features,
    extra_features=extra_features,
    tag='all_features'

)
print('extracted features in {:.1f}s'.format((time.time() - start)))
print('optimizing over {} features'.format(x_train.shape[1]))
print(x_train.columns)

start = time.time()
# classifier = ensemble.GradientBoostingClassifier(learning_rate=0.05, n_estimators=1250, subsample=1.)
classifier = ensemble.GradientBoostingClassifier()
classifier.fit(x_train, y_train)


rfecv = RFECV(classifier, cv=5, step=1, scoring='f1', n_jobs=-1)
rfecv = rfecv.fit(x_train, y_train)

rfecv_features = [feature for support, feature in zip(rfecv.get_support(), x_train.columns) if support]
print('optimal number of features:', rfecv.n_features_)
print('best features:', rfecv_features)

n_features = x_train.shape[1]
plt.figure(figsize=(8,8))
plt.barh(range(n_features), classifier.feature_importances_, align='center')
plt.yticks(np.arange(n_features), x_train.columns.values)
plt.xlabel('feature importance')
plt.ylabel('feature')
plt.savefig('../notebooks/plots/rfecv_feature_importances.pdf')

fig, ax = plt.subplots()
plt.xlabel("number of features selected")
plt.ylabel("f1")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig('../notebooks/plots/rfecv_feature_counts.pdf')
# fancier plotting option: https://www.scikit-yb.org/en/latest/api/model_selection/rfecv.html

print('calculated best features in {:.1f}s'.format(time.time() - start))
