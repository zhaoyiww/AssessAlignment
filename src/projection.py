import numpy as np
from scipy.interpolate import griddata  # for generating panoramic image
import pyransac3d as pyrsc  # for filtering outliers from detected planar targets
import os  # os.join.path
import random
import gc  # clear cache


# normalize the value to (0,255)
def scale_to_255(intensity, dtype=np.uint8):
    i_min = np.min(intensity)
    i_max = np.max(intensity)
    return ((intensity - i_min) / float(i_max - i_min) * 255).astype(dtype)


# convert xyz to phv
def xyz2shv(x, y, z):
    point_s = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    hor = np.arctan2(y, x)
    ver = np.arccos(z / point_s)
    return point_s, hor, ver


# downsample either randomly or w.r.t. the distance
def downsampling(raw_pc, downsampling_size, downsampling_type):
    """
    implement two ways for subsampling the raw large point cloud
    Input:
        raw_pc:
        sub_size:
        sub_type:

    Output:
        sub_index: the index of subsampled points in raw point cloud
    """
    if downsampling_type == 'random':
        sub_index = np.random.choice(raw_pc.shape[0], downsampling_size, replace=False)
    elif downsampling_type == 'dist_scope_ratio':
        # sub_index = []
        dist = np.linalg.norm(raw_pc[:, :3], axis=1)
        dist_diff = np.max(dist) - np.min(dist)
        index_close = np.asarray(np.where(dist <= dist_diff * 2 / 3))[0]
        index_far = np.asarray(np.where(dist > dist_diff * 2 / 3))[0]
        if index_far.shape[0] < int(downsampling_size / 2):
            index_close_idx = np.random.choice(index_close.shape[0], downsampling_size - index_far.shape[0], replace=False)
            sub_index = np.append(index_close[index_close_idx], index_far)
        else:
            index_close_idx = np.random.choice(index_close.shape[0], int(downsampling_size / 2), replace=False)
            index_far_idx = np.random.choice(index_far.shape[0], int(downsampling_size / 2), replace=False)
            sub_index = np.append(index_close[index_close_idx], index_far[index_far_idx])
    elif downsampling_type == 'dist_order_ascending':
        dist = np.linalg.norm(raw_pc[:, :3], axis=1)
        # suggest using meargesort or heapsort
        indx = np.argsort(dist, axis=-1, kind='mergesort', order=None)
        indx_fore = np.random.choice(indx[:int(len(indx) / 2)], int(downsampling_size / 2), replace=False)
        indx_back = np.random.choice(indx[int(len(indx) / 2):], int(downsampling_size / 2), replace=False)
        sub_index = np.append(indx_fore, indx_back)
    else:
        raise ValueError('The downsampling method is not defined')

    return sub_index


# self ransac plane fitting algorithm
def ransac_plane_fitting(coord_xyz, threshold=0.05, iterations=1000):
    """
    ransac plane fitting, a * x + b * y + c * z = d
    Input:
        coord_xyz: 3D point cloud, i.e., coordinates x, y, z
        threshold: the maximum distance of point to the plane
        iterations: the maximum iteration
    Output:
        plane_eq: plane equation, i.e., a, b, c, d
        plane_xyz: points that within the threshold
        inliers: point indexes that within the threshold
    """
    best_eq = []
    best_inliers = []
    for i in range(iterations):
        num_points = coord_xyz.shape[0]
        id_samples = random.sample(range(num_points), 3)
        vec_1 = coord_xyz[id_samples[1]] - coord_xyz[id_samples[0]]
        vec_2 = coord_xyz[id_samples[2]] - coord_xyz[id_samples[0]]

        normal = np.cross(vec_1, vec_2)
        # avoid the warning result
        np.seterr(invalid='ignore')
        normal = normal / np.linalg.norm(normal)
        d = np.dot(normal.T, coord_xyz[id_samples[1]].T)
        equation = np.hstack((normal, d))

        s = (np.dot(normal.T, coord_xyz.T) - d) / np.linalg.norm(normal)

        id_inliers = np.where(np.abs(s) <= threshold)[0]
        if len(id_inliers) > len(best_inliers):
            best_eq = equation
            best_inliers = id_inliers

    return best_eq, coord_xyz[best_inliers], best_inliers


def pc2img(point_cloud: str, scan_res: float):
    """
    this function goal: project 3D point cloud to 2D panoramic image

    Inputs:
        point_cloud:
        xx (dict): current
        xxx (str): path

    Returns:

    """

    # read raw point data from point cloud
    point_cloud = np.array(point_cloud)
    x, y, z, li = point_cloud[:, 0], point_cloud[:, 1], point_cloud[:, 2], point_cloud[:, 3]
    # r, g, b = point_cloud[:, 4], point_cloud[:, 5], point_cloud[:, 6]

    # calculate shv from xyz
    point_s, hor, ver = xyz2shv(x, y, z)

    # scan resolution, field-of-view (FoV) (hor, ver)
    res = scan_res
    # h = np.arange(-np.pi, np.pi, res, dtype=float)
    # v = np.arange(0, np.pi * 5 / 6, res, dtype=float)
    h = np.arange(min(hor), max(hor), res, dtype=float)
    v = np.arange(min(ver), max(ver), res, dtype=float)
    img_h, img_v = np.meshgrid(h, v)

    # image interpolation
    img = griddata((hor, ver), scale_to_255(li), (img_h, img_v), fill_value=0, method='linear', rescale=True)

    return img


def img2pc(point_cloud, bbox_result, scan_res, work_dir):
    """
    this function goal: reproject 2D panoramic image to 3D point cloud

    Inputs:
        point_cloud:
        xx (dict): current
        xxx (str): bbox_image_coord

    Returns:

    """

    x, y, z = point_cloud[:, 0], point_cloud[:, 1], point_cloud[:, 2]

    # calculate shv from xyz
    point_s, hor, ver = xyz2shv(x, y, z)

    # get the 9-dim raw point cloud, including H,V,X,Y,Z,I,R,G,B channels
    hor_ver = np.vstack((hor, ver)).T
    # pc_9_dim = np.concatenate((hor[:, None], ver[:, None], point_cloud), axis=1)

    res = scan_res
    # h = np.arange(-np.pi, np.pi, res, dtype=float)
    # v = np.arange(0, np.pi * 5 / 6, res, dtype=float)
    h = np.arange(min(hor_ver[:, 0]), max(hor_ver[:, 0]), res, dtype=float)
    v = np.arange(min(hor_ver[:, 1]), max(hor_ver[:, 1]), res, dtype=float)

    del point_s, hor, ver, x, y, z
    gc.collect()

    # create folders to save the results
    # os.makedirs(os.path.join(work_dir, 'extraction/', 'XYZIRGB'))
    os.makedirs(os.path.join(work_dir, 'extraction/', 'XYZI'))
    os.makedirs(os.path.join(work_dir, 'extraction/', 'XYZRGB'))
    addr = os.path.join(work_dir, 'extraction')

    coord_mean = []
    for i in range(bbox_result.shape[0]):
    # for i in np.arange(bbox_result.shape[0]):
        # get the range of each target in the sampled horizontal-vertical plane
        h_seg = h[int(bbox_result[i, 0]):int(bbox_result[i, 2])]
        v_seg = v[int(bbox_result[i, 1]):int(bbox_result[i, 3])]
        # wider edge for MAST project
        # h_min, h_max = min(h_seg) - 100 / 10000, max(h_seg) + 100 / 10000
        # v_min, v_max = min(v_seg) - 100 / 10000, max(v_seg) + 100 / 10000
        h_min, h_max = min(h_seg) - 5 / 10000, max(h_seg) + 5 / 10000
        v_min, v_max = min(v_seg) - 5 / 10000, max(v_seg) + 5 / 10000

        temp_indx = np.where((hor_ver[:, 0] <= h_max) & (hor_ver[:, 0] >= h_min)
                             & (hor_ver[:, 1] <= v_max) & (hor_ver[:, 1] >= v_min))

        power4 = point_cloud[temp_indx, :][0]

        # use ransac to filter the targets
        plane = pyrsc.Plane()
        # best_eq, best_inliers = plane.fit(power4[:, 2:5], 0.01)
        np.seterr(invalid='ignore')
        # _, best_inliers = plane.fit(power4[:, 2:5], 0.01)
        _, best_inliers = plane.fit(power4[:, :3], 0.01)
        # best_inliers = np.arange(power4)

        # np.savetxt(os.path.join(addr, 'XYZIRGB/', f'{i}.txt'), power4[best_inliers, 2:], fmt='%.5f')
        # save the approximate coordinates for each target center
        coord_mean.append(power4[best_inliers, :3].mean(axis=0))
        np.savetxt(os.path.join(addr, 'XYZI/', f'{str(i).zfill(3)}.txt'), power4[best_inliers, :4], fmt='%.5f')
        power4 = np.delete(power4, 3, axis=1)
        np.savetxt(os.path.join(addr, 'XYZRGB/', f'{str(i).zfill(3)}.txt'), power4[best_inliers, :], fmt='%.5f')

        # coord_mean = np.array(coord_mean)
    # coord_mean = np.hstack((np.array(range(coord_mean.shape[0]))[:, None], coord_mean))

    # with open(os.path.join(addr, f'Approximate_Coord_XYZIRGB.txt'), 'w') as f:
    #     f.write(coord_mean)
    coord_mean = np.hstack((np.array(range(np.array(coord_mean).shape[0]))[:, None], coord_mean))
    np.savetxt(os.path.join(addr, f'approximate_3Dcoordinates.txt'), coord_mean, fmt='%.5f')

    return bbox_result.shape[0], addr
