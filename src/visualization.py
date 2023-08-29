import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import shutil  # remove non-empty directory
from sklearn.cluster import KMeans # kmeans regression
plt.rc('font', family="Times New Roman", size=24)


# convert xyz to phv
def xyz2phv(coord_xyz):
    point_s = np.linalg.norm(coord_xyz, axis=1)
    hor = np.arctan2(coord_xyz[:, 1], coord_xyz[:, 0])
    ver = np.arccos(coord_xyz[:, 2] / point_s)

    return np.asarray([point_s, hor, ver]).T


# generate horizontal and vertical map, plus k-means regression
def generate_hv_map(li_phv, rgb_phv, ref_values, label, results_folder):
    # figure 1: vertical-horizontal map
    plt.figure(figsize=(8, 6))
    fig, ax00 = plt.subplots()
    # labels = [label, 'point location']
    delta_hv = (li_phv - rgb_phv)[:, 1:] * 1000
    delta_hv_shift = delta_hv.copy()
    delta_hv_shift[:, 0] = np.multiply(delta_hv[:, 0], np.sin(li_phv[:, 2]))

    if label == 'BLK360' or label == 'RTC360':
        nCluster = 3
    else:
        nCluster = 1
    kmCluster = KMeans(n_clusters=nCluster, n_init='auto').fit(delta_hv_shift)
    categories = kmCluster.labels_
    colors = ['r', 'g', 'b']
    markers = ['o', 's', '^']
    # markers = ['s', '^', 'o']
    categories_name = ['T$_1$', 'T$_2$', 'T$_3$']
    for i in range(nCluster):
        regress = [index for index, category in enumerate(categories) if category == i]
        hor_reg = li_phv[regress, 1] * 180 / np.pi
        ver_reg = li_phv[regress, 2] * 180 / np.pi
        delta_hor_shift_reg = delta_hv_shift[regress, 0]
        delta_ver_shift_reg = delta_hv_shift[regress, 1]

        factor = 1 / 180 * np.pi * 1000
        Q = plt.quiver(hor_reg, ver_reg, np.multiply(delta_hor_shift_reg, factor),
                       np.multiply(delta_ver_shift_reg, factor), color=colors[0],
                       width=0.003, headwidth=3, headlength=4.5, scale=1, angles='xy', scale_units='xy')
        plt.quiverkey(Q, X=0.67, Y=0.955, U=1 * factor, coordinates='axes',
                       label=str(format(1, '.1f')) + ' mrad', labelpos='E', color='red')
        plt.scatter(hor_reg, ver_reg, label=categories_name[i],
                    c=colors[0], edgecolors=colors[0], s=15, marker=markers[0])

    plt.xlim(-200, 200)
    plt.ylim(155, -5)
    ax00.set_xticks([-202, -177, -90, 0, 90, 177, 202])
    ax00.set_xticklabels(('', '-180', '-90', '0', '90', '180', ''))
    ax00.set_yticks([-9, 3, 30, 60, 90, 120, 147, 157])
    ax00.set_yticklabels(('', '0', '', '60', '90', '', '150', ''))
    ax01 = plt.gca()
    if label == 'Faro_3296':
        plt.xlabel('Horizontal angle [$\degree$]')
        plt.ylabel('Vertical angle [$\degree$]')
    elif label == 'C10' or label == 'P50' or label == 'Faro_3287' or label == 'Faro_S120':
        ax01.axes.xaxis.set_ticklabels([])
        ax01.axes.yaxis.set_ticklabels([])
    elif label == 'BLK360' or label == 'RTC360':
        ax01.axes.xaxis.set_ticklabels([])
        plt.ylabel('Vertical angle [$\degree$]')
    else:
        ax01.axes.yaxis.set_ticklabels([])
        plt.xlabel('Horizontal angle [$\degree$]')
    plt.grid()
    plt.savefig(os.path.join(results_folder, f'ver2hor_{label}.pdf'), dpi=600, bbox_inches='tight')

    theta = np.sqrt(np.multiply(delta_hv_shift[:, 0], delta_hv_shift[:, 0]) +
                    np.multiply(delta_hv_shift[:, 1], delta_hv_shift[:, 1]))

    outlier_percent = np.sum((theta > max(ref_values)) |
                             (theta < -max(ref_values))) / theta.shape[0]
    outlier_percent = '%.1f%%' % (outlier_percent * 100)

    l1 = plt.axvline(-ref_values[0], color='g', linewidth=2, linestyle='-.', label='Beam divergence')
    l2 = plt.axvline(-ref_values[1], color='r', linewidth=2, linestyle='--', label='Image resolution')
    l3 = plt.axvline(-ref_values[2], color='purple', linewidth=2, linestyle=':', label='Scan resolution')

    # figure 2:
    plt.figure(figsize=(8, 6))
    fig, ax4 = plt.subplots()

    from matplotlib.patches import Ellipse, Circle
    plt.scatter(delta_hv_shift[:, 0], delta_hv_shift[:, 1], c='blue', s=10)
    cir1 = Circle(xy=(0.0, 0.0), radius=ref_values[0], alpha=1, color='g', fill=False, linestyle='--')
    cir2 = Circle(xy=(0.0, 0.0), radius=ref_values[1], alpha=1, color='red', fill=False, linestyle='-.')
    cir3 = Circle(xy=(0.0, 0.0), radius=ref_values[2], alpha=1, color='purple', fill=False, linestyle=':')
    ax4.add_patch(cir1)
    ax4.add_patch(cir2)
    ax4.add_patch(cir3)
    plt.axhline(0, color='gray', lw=1.0)
    plt.axvline(0, color='gray', lw=1.0)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    ax4.set_xticks([-2.1, -1.9, -1.0, 0, 1.0, 1.9, 2.1])
    ax4.set_xticklabels(('', '-2', '-1', '0', '1', '2', ''))
    ax4.set_yticks([-2.1, -1.9, -1.0, 0, 1.0, 1.9, 2.1])
    ax4.set_yticklabels(('', '-2', '-1', '0', '1', '2', ''))
    # plt.grid()

    if ref_values[1] >= np.max(ref_values):
        color_label = 'r'
    elif ref_values[0] >= np.max(ref_values):
        color_label = 'g'
    else:
        color_label = 'purple'
    if label == 'C10':
        plt.text(max(ref_values) - 0.25, -1.20, outlier_percent + ' out', color=color_label, va="center")
    else:
        plt.text(max(ref_values) - 0.25, -0.80, outlier_percent + ' out', color=color_label, va="center")

    ax4 = plt.gca()
    ax4.set_aspect(1)

    if label == 'Faro_3296':
        plt.xlabel('Hor. angle dev. [mrad]')
        plt.ylabel('Ver. angle dev. [mrad]')
        plt.legend(handles=[l1, l2, l3], loc='best', ncol=1, fontsize=18, framealpha=0.3)
    elif label == 'C10' or label == 'P50' or label == 'Faro_3287' or label == 'Faro_S120':
        ax4.axes.xaxis.set_ticklabels([])
        ax4.axes.yaxis.set_ticklabels([])
    elif label == 'BLK360' or label == 'RTC360':
        ax4.axes.xaxis.set_ticklabels([])
        plt.ylabel('Ver. angle dev. [mrad]')
    else:
        ax4.axes.yaxis.set_ticklabels([])
        plt.xlabel('Hor. angle dev. [mrad]')

    plt.savefig(os.path.join(results_folder, f'distribution_{label}.pdf'), dpi=600, bbox_inches='tight')

    return None


def parse_args():
    parser = argparse.ArgumentParser(
        description='automatic framework for generating visualization results')
    parser.add_argument('--input_dir', type=str, default='../data/results/estimation',
                        help='the directory to load data for visualization')
    parser.add_argument('--label', type=str, default='BLK360',
                        help='the instrument that collecting current data')
    parser.add_argument('--ref_value1', type=float, default=0.68,
                        help='the reference values (footprint size, mrad)')
    parser.add_argument('--ref_value2', type=float, default=0.40,
                        help='the reference values (pixel size, mrad)')
    parser.add_argument('--ref_value3', type=float, default=0.50,
                        help='the reference values (scan resolution, mrad)')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    # set a path and remove the old path for saving results for each individual scans
    work_dir = os.path.abspath(os.path.join(args.input_dir, '..', 'visualization'))
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    path_1 = os.path.join(args.input_dir, 'estimated_3Dcoordinates_XYZI.txt')
    path_2 = os.path.join(args.input_dir, 'estimated_3Dcoordinates_XYZRGB.txt')

    coord_1 = np.loadtxt(path_1, delimiter='	', skiprows=1, dtype=float)
    coord_2 = np.loadtxt(path_2, delimiter='	', skiprows=1, dtype=float)

    # 01. vertical / horizontal map
    generate_hv_map(xyz2phv(coord_1[:, 0:]), xyz2phv(coord_2[:, 0:]),
                    [args.ref_value1, args.ref_value2, args.ref_value3], args.label, work_dir)

    # coord_dev_vs_dist(coord_1[:, 0:], coord_2[:, 0:], [args.ref_value1, args.ref_value2], args.label, work_dir)

if __name__ == '__main__':
    main()
