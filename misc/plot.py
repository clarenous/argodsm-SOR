import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


def all_array_sizes(im_list, jm_list, km_list):
    result = []
    for im in im_list:
        for jm in jm_list:
            for km in km_list:
                result.append((im, jm, km))
    return result


def get_nx_ny_list(node_count):
    if node_count == 1:
        return [(1, 1)]
    if node_count == 2:
        return [(1, 2), (2, 1)]
    if node_count == 4:
        return [(1, 4), (2, 2), (4, 1)]
    if node_count == 8:
        return [(1, 8), (2, 4), (4, 2), (8, 1)]
    if node_count == 16:
        return [(1, 16), (2, 8), (4, 4), (8, 2), (16, 1)]
    raise Exception('unsupported ndoe_count: %d' % node_count)


def combinitions_with_same_volume(im_list, jm_list, km_list):
    result_dict = {}
    sizes = all_array_sizes(im_list, jm_list, km_list)
    for size in sizes:
        volume = size[0]*size[1]*size[2]
        if volume not in result_dict:
            result_dict[volume] = []
        result_dict[volume].append((size[0], size[1], size[2]))
    return result_dict.items()


def parse_by_optimization_boost(df):
    df_non_optimized = df[df['optimized'] == 0]
    df_optimized = df[df['optimized'] == 1].copy()
    optimization_boost = np.array(
        df_non_optimized['time']) / np.array(df_optimized['time'])
    df_optimized.insert(
        df_optimized.shape[1], 'optimization_boost', optimization_boost)
    return df_optimized


def parse_by_cpu_time(df):
    cpu_time = np.array(df['time']) * np.array(df['node_count'])
    df_cpu_time = df.copy()
    df_cpu_time.insert(df.shape[1], 'cpu_time', cpu_time)
    return df_cpu_time


def select_by_im(df, im):
    return df[df['im'] == im]


def select_by_jm(df, jm):
    return df[df['jm'] == jm]


def select_by_km(df, km):
    return df[df['km'] == km]


def select_by_array_size(df, im, jm, km):
    return df[(df['im'] == im) & (df['jm'] == jm) & (df['km'] == km)]


def select_by_method(df, method):
    return df[df['method'] == method]


def select_by_node_count(df, node_count):
    return df[df['node_count'] == node_count]


def select_by_dsm_nx(df, dsm_nx):
    return df[df['dsm_nx'] == dsm_nx]


def select_by_dsm_ny(df, dsm_ny):
    return df[df['dsm_ny'] == dsm_ny]


def select_by_dsm_nx_ny(df, dsm_nx, dsm_ny):
    return df[(df['dsm_nx'] == dsm_nx) & (df['dsm_ny'] == dsm_ny)]


def fig_path(fig_name):
    return 'figs/'+fig_name


im_list = [128, 256, 512]
jm_list = [128, 256, 512]
km_list = [64, 128, 256]
node_count_list = [2, 8, 16]
fig_method_string_list = ['C', 'ArgoDSM V1', 'ArgoDSM V2', 'ArgoDSM V3', 'MPI']

# Create figure storage directory
os.makedirs('figs', exist_ok=True)

# Read from file
df = pd.read_csv('test-results.csv')
df = parse_by_cpu_time(df)
df = parse_by_optimization_boost(df)
print(df)

# Analyze the affect of the compiler '-O3' optimization flag
negative_boost = df[df['optimization_boost'] < 0.95]
neutral_boost = df[(df['optimization_boost'] >= 0.95) & (
    df['optimization_boost'] < 1.05)]
positive_boost = df[df['optimization_boost'] >= 1.05]
negative_counts = [
    len(select_by_method(negative_boost, 0)),
    len(select_by_method(negative_boost, 1)),
    len(select_by_method(negative_boost, 2)),
    len(select_by_method(negative_boost, 3)),
    len(select_by_method(negative_boost, 4)),
]
neutral_counts = [
    len(select_by_method(neutral_boost, 0)),
    len(select_by_method(neutral_boost, 1)),
    len(select_by_method(neutral_boost, 2)),
    len(select_by_method(neutral_boost, 3)),
    len(select_by_method(neutral_boost, 4)),
]
positive_counts = [
    len(select_by_method(positive_boost, 0)),
    len(select_by_method(positive_boost, 1)),
    len(select_by_method(positive_boost, 2)),
    len(select_by_method(positive_boost, 3)),
    len(select_by_method(positive_boost, 4)),
]
ax_optimization_bar = pd.DataFrame(
    np.array([negative_counts, neutral_counts, positive_counts]).T,
    index=fig_method_string_list,
    columns=['Negative Count', 'Neutral Count', 'Positive Count'],
).plot(kind='area', figsize=(12, 6), title='The effect of compiler \'-O3\' optimization flag')
for container in ax_optimization_bar.containers:
    ax_optimization_bar.bar_label(container)
for tick in ax_optimization_bar.get_xticklabels():
    tick.set_rotation(0)
ax_optimization_bar.set_ylabel('Number of different effects')
ax_optimization_bar.get_figure().savefig(fig_path('optimization_area.pdf'))
print('Negative Boost Count (Method 0 ~ 4):', negative_counts)
print('Neutral Boost Count (Method 0 ~ 4):', neutral_counts)
print('Positive Boost Count (Method 0 ~ 4):', positive_counts)

# Analyze the performance when array_size is fixed
array_size_best_index = []
array_size_best_data = []
for size in all_array_sizes(im_list, jm_list, km_list):
    im, jm, km = size[0], size[1], size[2]
    df_array_size = select_by_array_size(
        df, im, jm, km
    ).sort_values('time')
    array_size_best_index.append('(%d,\n %d,\n %d)' % (im, jm, km))
    array_size_best_data.append([
        select_by_method(df_array_size, 0).head(1).iloc[0].time,
        select_by_method(df_array_size, 1).head(1).iloc[0].time,
        select_by_method(df_array_size, 2).head(1).iloc[0].time,
        select_by_method(df_array_size, 3).head(1).iloc[0].time,
        select_by_method(df_array_size, 4).head(1).iloc[0].time,
    ])
ax_array_size_best_bar = pd.DataFrame(
    np.log10(np.array(array_size_best_data)),
    index=array_size_best_index,
    columns=fig_method_string_list
).plot(kind='bar', stacked=False, figsize=(18, 7), title='The SOR execution time (seconds) of each implementation in log10 scale')
for tick in ax_array_size_best_bar.get_xticklabels():
    tick.set_rotation(0)
ax_array_size_best_bar.set_xlabel('Array Shape (im, jm, km)')
ax_array_size_best_bar.set_ylabel(
    'SOR execution time (seconds) in log10 scale')
ax_array_size_best_bar.get_figure().savefig(fig_path('array_size_bar.pdf'))

# Analyze the performance for each volume
volumes_fig, volumes_axes = plt.subplots(
    nrows=2, ncols=3, sharex=False, figsize=(22, 12))
volumes_axes[1][0].set_xlabel('Array Shape (im, jm, km)')
volumes_axes[1][1].set_xlabel('Array Shape (im, jm, km)')
volumes_axes[1][2].set_xlabel('Array Shape (im, jm, km)')
volumes_axes = [
    volumes_axes[1][0], volumes_axes[0][0], volumes_axes[0][1],
    volumes_axes[0][2], volumes_axes[1][1], volumes_axes[1][2],
]
volumes_sub_fig_ordinals = ['d', 'a', 'b', 'c', 'e', 'f']
for item in combinitions_with_same_volume(im_list, jm_list, km_list):
    volume = int(math.log2(item[0]))
    if volume <= 20:
        continue
    volume_best_index = []
    volume_best_data = []
    for combinition in item[1]:
        im, jm, km = combinition[0], combinition[1], combinition[2]
        df_array_size = select_by_array_size(
            df, im, jm, km
        ).sort_values('time')
        volume_best_index.append('(%d,\n%d,\n%d)' % (im, jm, km))
        volume_best_data.append([
            select_by_method(df_array_size, 0).head(1).iloc[0].time,
            select_by_method(df_array_size, 1).head(1).iloc[0].time,
            select_by_method(df_array_size, 2).head(1).iloc[0].time,
            select_by_method(df_array_size, 3).head(1).iloc[0].time,
            select_by_method(df_array_size, 4).head(1).iloc[0].time,
        ])
    ax_volume_best_bar = pd.DataFrame(
        np.array(volume_best_data),
        index=volume_best_index,
        columns=fig_method_string_list
    ).plot(ax=volumes_axes[volume-21], kind='bar', stacked=False, title='(%s) Array Volume = 2 ^ %d' % (volumes_sub_fig_ordinals[volume-21], volume))
    for tick in ax_volume_best_bar.get_xticklabels():
        tick.set_rotation(0)
    print("Array Volume:", volume)
volumes_fig.suptitle(
    'The SOR execution time (secondes) of each Implementation under fixed array volumes'
)
volumes_fig.savefig(fig_path('volumes_bar.pdf'))

# Analyze the affect of node_count in MPI
mpi_boost_factor_index = []
mpi_speed_up_ratio_data = []
mpi_parallel_efficiency_data = []
for size in all_array_sizes(im_list, jm_list, km_list):
    im, jm, km = size[0], size[1], size[2]
    df_array_size = select_by_array_size(
        df, im, jm, km
    )
    df_mpi = select_by_method(df_array_size, 4).copy()
    df_c = select_by_method(df_array_size, 0)
    speed_up = np.array(df_c['time']) / np.array(df_mpi['time'])
    df_mpi.insert(df.shape[1], 'speed_up', speed_up)
    selected_data = [
        df_mpi.head(3).iloc[0],
        df_mpi.head(3).iloc[1],
        df_mpi.head(3).iloc[2]
    ]
    mpi_boost_factor_index.append('(%d,\n %d,\n %d)' % (im, jm, km))
    mpi_speed_up_ratio_data.append([
        selected_data[0]['speed_up'],
        selected_data[1]['speed_up'],
        selected_data[2]['speed_up']
    ])
    mpi_parallel_efficiency_data.append([
        selected_data[0]['speed_up']/selected_data[0]['node_count'],
        selected_data[1]['speed_up']/selected_data[1]['node_count'],
        selected_data[2]['speed_up']/selected_data[2]['node_count']
    ])
# plot speed_up_ratio
mpi_boost_factor_ax = pd.DataFrame(
    np.array(mpi_speed_up_ratio_data),
    index=mpi_boost_factor_index,
    columns=['Node Count = 2', 'Node Count = 8', 'Node Count = 16']
).plot(kind='bar', stacked=True, grid=False, figsize=(18, 7), title='MPI speedup ratio and parallel efficiency')
mpi_boost_factor_ax.set_ylabel('Speedup Ratio (Stacked Bar)')
# plot parappel_efficiency
mpi_boost_factor_line = pd.DataFrame(
    np.array(mpi_parallel_efficiency_data),
    index=mpi_boost_factor_index,
    columns=['Node Count = 2', 'Node Count = 8', 'Node Count = 16']
).plot(ax=mpi_boost_factor_ax.twinx(), kind='line', stacked=False, grid=False)
mpi_boost_factor_line.set_ylabel('Parallel Efficiency (Line)')
# figure
mpi_boost_factor_ax.set_xlabel('Array Shape (im, jm, km)')
mpi_boost_factor_ax.set_xticks(np.arange(len(mpi_boost_factor_index)))
mpi_boost_factor_ax.set_xticklabels(mpi_boost_factor_index)
for tick in mpi_boost_factor_ax.get_xticklabels():
    tick.set_rotation(0)
mpi_boost_factor_ax.get_figure().savefig(fig_path('mpi_boost_factor.pdf'))


def plot_argo_boost_factor_bar(method, ax):
    argo_boost_factor_index = []
    argo_speed_up_ratio_data = []
    argo_parallel_efficiency_data = []
    for size in all_array_sizes(im_list, jm_list, km_list):
        im, jm, km = size[0], size[1], size[2]
        df_array_size = select_by_array_size(
            df, im, jm, km
        )
        df_argo = select_by_method(df_array_size, method).copy()
        df_c = select_by_method(df_array_size, 0)
        speed_up = np.array(df_c['time']) / np.array(df_argo['time'])
        df_argo.insert(df.shape[1], 'speed_up', speed_up)
        df_argo = df_argo.sort_values('time')
        selected_data = [
            select_by_node_count(df_argo, 2).head(1).iloc[0],
            select_by_node_count(df_argo, 8).head(1).iloc[0],
            select_by_node_count(df_argo, 16).head(1).iloc[0],
        ]
        argo_boost_factor_index.append('(%d,\n%d,\n%d)' % (im, jm, km))
        argo_speed_up_ratio_data.append([
            selected_data[0]['speed_up'],
            selected_data[1]['speed_up'],
            selected_data[2]['speed_up']
        ])
        argo_parallel_efficiency_data.append([
            selected_data[0]['speed_up']/selected_data[0]['node_count'],
            selected_data[1]['speed_up']/selected_data[1]['node_count'],
            selected_data[2]['speed_up']/selected_data[2]['node_count']
        ])
    pd.DataFrame(
        np.array(argo_speed_up_ratio_data),
        index=argo_boost_factor_index,
        columns=['Node Count = 2', 'Node Count = 8', 'Node Count = 16']
    ).plot(
        ax=ax, kind='bar', stacked=True, grid=False, title='ArgoDSM Version %d' % method
    ).set_xticks([])
    pd.DataFrame(
        np.array(argo_parallel_efficiency_data),
        index=argo_boost_factor_index,
        columns=['Node Count = 2', 'Node Count = 8', 'Node Count = 16']
    ).plot(
        ax=ax.twinx(), kind='line', stacked=False, grid=False
    ).set_xticks([])
    ax.set_ylabel('Speedup Ratio (Stacked Bar)')
    # ax.twinx().set_ylabel('Parallel Efficiency (Line)')
    return argo_boost_factor_index


# Analyze the affect of node_count in ArgoDSM Versions
argo_boost_factor_fig, argo_boost_factor_axes = plt.subplots(
    nrows=3, ncols=1, sharex=True, figsize=(18, 15))
plot_argo_boost_factor_bar(1, argo_boost_factor_axes[0])
plot_argo_boost_factor_bar(2, argo_boost_factor_axes[1])
argo_boost_factor_xticks = plot_argo_boost_factor_bar(
    3, argo_boost_factor_axes[2]
)
argo_boost_factor_axes[2].set_xticks(np.arange(len(argo_boost_factor_xticks)))
argo_boost_factor_axes[2].set_xticklabels(argo_boost_factor_xticks)
argo_boost_factor_axes[2].set_xlabel('Array Shape (im, jm, km)')
for tick in argo_boost_factor_axes[2].get_xticklabels():
    tick.set_rotation(0)
argo_boost_factor_fig.suptitle(
    'ArgoDSM speedup ratio and parallel efficiency'
)
argo_boost_factor_fig.savefig(fig_path('argo_boost_factor.pdf'))


def plot_argo_nx_ny_bar(ax, im, jm, km):
    argo_nx_ny_index = []
    argo_nx_ny_data = []
    df_argo_nx_ny = select_by_array_size(df, im, jm, km)
    for node_count in node_count_list:
        for nx_ny in get_nx_ny_list(node_count):
            nx, ny = nx_ny[0], nx_ny[1]
            df_array_base = select_by_dsm_nx_ny(df_argo_nx_ny, nx, ny)
            argo_nx_ny_index.append('(%d, %d)' % (nx, ny))
            argo_nx_ny_data.append([
                select_by_method(df_array_base, 1).head(1).iloc[0].time,
                select_by_method(df_array_base, 2).head(1).iloc[0].time,
                select_by_method(df_array_base, 3).head(1).iloc[0].time,
            ])
    ax_argo_nx_ny_bar = pd.DataFrame(
        np.array(argo_nx_ny_data),
        index=argo_nx_ny_index,
        columns=fig_method_string_list[1:4]
    ).plot(ax=ax, kind='bar', stacked=False, title='Array Shape = (%d, %d, %d)' % (im, jm, km))
    for tick in ax_argo_nx_ny_bar.get_xticklabels():
        tick.set_rotation(0)


# Analyze the affect of dsm_nx and dsm_ny in ArgoDSM Versions
argo_nx_ny_fig, argo_nx_ny_axes = plt.subplots(
    nrows=2, ncols=2, sharex=True, figsize=(16, 8)
)
plot_argo_nx_ny_bar(argo_nx_ny_axes[0][0], 512, 512, 128)
plot_argo_nx_ny_bar(argo_nx_ny_axes[0][1], 512, 512, 256)
plot_argo_nx_ny_bar(argo_nx_ny_axes[1][0], 512, 256, 256)
plot_argo_nx_ny_bar(argo_nx_ny_axes[1][1], 256, 512, 256)
argo_nx_ny_axes[1][0].set_xlabel('(dsmNX, dsmNY)')
argo_nx_ny_axes[1][1].set_xlabel('(dsmNX, dsmNY)')
argo_nx_ny_fig.suptitle(
    'The SOR execution time (seconds) of ArgoDSM under different (dsmNX, dsmNY)'
)
argo_nx_ny_fig.savefig(fig_path('argo_nx_ny_bar.pdf'))
