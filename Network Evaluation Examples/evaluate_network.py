import sys

sys.path.append('..')

from network_evaluation_tools import data_import_tools as dit
from network_evaluation_tools import network_evaluation_functions as nef
from network_evaluation_tools import network_propagation as prop
import pandas as pd
import numpy as np

import pickle
import os
import sys

# get the output dirrectory
output_dir = sys.argv[3]
if output_dir[-1] != '/':
    output_dir += '/'

# creat the output dir
try:
    os.mkdir(output_dir)
except FileExistsError:
    print('Output dir already exists: ' + output_dir)

# Load network (We choose a smaller network here for the example's sake)
network = dit.load_network_file(sys.argv[1], verbose=True, delimiter='\t')
print(len(network.nodes))

# Load gene sets for analysis
genesets = dit.load_node_sets(sys.argv[2])

# Calculate geneset sub-sample rate
genesets_p = nef.calculate_p(network, genesets)

# Determine optimal alpha for network (can also be done automatically by next step)
alpha = prop.calculate_alpha(network)
print(alpha)

# Calculate network kernel for propagation
kernel = nef.construct_prop_kernel(network, alpha=alpha, verbose=True)

print(kernel.index)
print(genesets)

# Calculate the AUPRC values for each gene set
AUPRC_values = nef.small_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=30, cores=8, verbose=True)

# **Note about the above cell:** There are a several options for this particular step depending on the computational
# resources available and network size. If the network is sufficiently small (<250k edges), it is recommended to use
# the 'small_network_AUPRC_wrapper' function as it can be much faster, especially when run in parallel (at least 8G
# per core is recommended). If you would like to parallelize the AUPRC calculation with a larger network (between
# 250K and 2.5M edges), at least 16G per core is recommended, 32G per core if the network contains more than 2.5M
# edges. For larger networks, it is recommended to use the 'large_network_AUPRC_wrapper', which may be a slightly
# slower function, but more equipped to handle the larger memory footprint(required. To change the parllelization
# status of the function, change the 'cores' option to the number of threads you would like to utilize.)

# Construct null networks and calculate the AUPRC of the gene sets of the null networks
# We can use the AUPRC wrapper function for this
if os.path.exists('null_AUPRCs.pickle'):
    null_AUPRCs = pickle.load(open(output_dir + 'null_AUPRCs.pickle', 'rb'))
else:
    null_AUPRCs = []
    for i in range(10):
        shuffNet = nef.shuffle_network(network, max_tries_n=10, verbose=True)
        shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
        shuffNet_AUPRCs = nef.small_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=30, cores=8,
                                                          verbose=False)
        null_AUPRCs.append(shuffNet_AUPRCs)
        print('shuffNet', repr(i + 1), 'AUPRCs calculated')
    pickle.dump(null_AUPRCs, open(output_dir + 'null_AUPRCs.pickle', 'wb'))

# **Note about the above cell:** We use a small number to calculate the null AUPRC values, but a larger number of shuffled networks may give a better representation of the true null AUPRC value.  smaller number of networks here for this example, but larger numbers can be used, especially if the resulting distribution of null AUPRCs has a high variance relative to the actual AUPRC values, but we have found that the variance remains relatively small even with a small number of shuffled networks.

# Construct table of null AUPRCs
if os.path.exists('null_AUPRCs_table.pickle'):
    null_AUPRCs_table = pickle.load(open(output_dir + 'null_AUPRCs_table.pickle', 'rb'))
else:
    null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
    null_AUPRCs_table.columns = ['shuffNet' + repr(i + 1) for i in range(len(null_AUPRCs))]
    pickle.dump(null_AUPRCs_table, open(output_dir + 'null_AUPRCs_table.pickle', 'wb'))

# Calculate performance metric of gene sets
if os.path.exists('network_performance.pickle'):
    network_performance = pickle.load(open(output_dir + 'network_performance.pickle', 'rb'))
else:
    network_performance = nef.calculate_network_performance_score(AUPRC_values, null_AUPRCs_table, verbose=True)
    network_performance.name = 'Test Network'
    pickle.dump(network_performance, open(output_dir + 'network_performance.pickle', 'wb'))

# Calculate network performance gain over median null AUPRC
if os.path.exists('network_perf_gain.pickle'):
    network_performance = pickle.load(open(output_dir + 'network_perf_gain.pickle', 'rb'))
else:
    network_perf_gain = nef.calculate_network_performance_gain(AUPRC_values, null_AUPRCs_table, verbose=True)
    network_perf_gain.name = 'Test Network'
    pickle.dump(network_perf_gain, open(output_dir + 'network_perf_gain.pickle', 'wb'))

# Rank network on average performance across gene sets vs performance on same gene sets in previous network set
if os.path.exists('all_network_performance.pickle'):
    network_performance = pickle.load(open(output_dir + 'all_network_performance.pickle', 'rb'))
    all_network_performance_filt = pickle.load(open(output_dir + 'all_network_performance_filt.pickle', 'rb'))
    network_performance_rank_table = pickle.load(open(output_dir + 'network_performance_rank_table.pickle', 'rb'))
    network_performance_rankings = pickle.load(open(output_dir + 'network_perf_gain.pickle', 'rb'))
else:
    all_network_performance = pd.read_csv('~/Data/Network_Performance.csv', index_col=0)
    all_network_performance_filt = pd.concat(
        [network_performance, all_network_performance.ix[network_performance.index]], axis=1)
    network_performance_rank_table = all_network_performance_filt.rank(axis=1, ascending=False)
    network_performance_rankings = network_performance_rank_table['Test Network']
    pickle.dump(all_network_performance, open(output_dir + 'all_network_performance.pickle', 'wb'))
    pickle.dump(all_network_performance_filt, open(output_dir + 'all_network_performance_filt.pickle', 'wb'))
    pickle.dump(network_performance_rank_table, open(output_dir + 'network_performance_rank_table.pickle', 'wb'))
    pickle.dump(network_performance_rankings, open(output_dir + 'network_performance_rankings.pickle', 'wb'))

# Rank network on average performance gain across gene sets vs performance gain on same gene sets in previous network set
all_network_perf_gain = pd.read_csv('~/Data/Network_Performance_Gain.csv', index_col=0)
all_network_perf_gain_filt = pd.concat([network_perf_gain, all_network_perf_gain.ix[network_perf_gain.index]], axis=1)
network_perf_gain_rank_table = all_network_performance_filt.rank(axis=1, ascending=False)
network_perf_gain_rankings = network_perf_gain_rank_table['Test Network']

# Network Performance
network_performance_metric_ranks = pd.concat(
    [network_performance, network_performance_rankings, network_perf_gain, network_perf_gain_rankings], axis=1)
network_performance_metric_ranks.columns = ['Network Performance', 'Network Performance Rank',
                                            'Network Performance Gain', 'Network Performance Gain Rank']
network_performance_metric_ranks.sort_values(
    by=['Network Performance Rank', 'Network Performance', 'Network Performance Gain Rank', 'Network Performance Gain'],
    ascending=[True, False, True, False])

# Construct network summary table
network_summary = {}
network_summary['Nodes'] = int(len(network.nodes()))
network_summary['Edges'] = int(len(network.edges()))
network_summary['Avg Node Degree'] = np.mean(network.degree().values())
network_summary['Edge Density'] = 2 * network_summary['Edges'] / float(
    (network_summary['Nodes'] * (network_summary['Nodes'] - 1)))
network_summary['Avg Network Performance Rank'] = network_performance_rankings.mean()
network_summary['Avg Network Performance Rank, Rank'] = int(
    network_performance_rank_table.mean().rank().ix['Test Network'])
network_summary['Avg Network Performance Gain Rank'] = network_perf_gain_rankings.mean()
network_summary['Avg Network Performance Gain Rank, Rank'] = int(
    network_perf_gain_rank_table.mean().rank().ix['Test Network'])
for item in ['Nodes', 'Edges', 'Avg Node Degree', 'Edge Density', 'Avg Network Performance Rank',
             'Avg Network Performance Rank, Rank',
             'Avg Network Performance Gain Rank', 'Avg Network Performance Gain Rank, Rank']:
    print(item + ':\t' + repr(network_summary[item]))
