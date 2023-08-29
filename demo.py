import argparse
import h5py  # read mat file
import numpy as np
import os  # join path
import shutil  # remove non-empty directory
from pyntcloud import PyntCloud  # load point cloud
from src import projection
from src import detection
import cv2  # write projected image
from timeit import default_timer as time

start_time = time()


def parse_args():
    parser = argparse.ArgumentParser(
        description='automated target detection')
    parser.add_argument('--point_cloud', type=str,
                        default='./data/example.pts',
                        help='input raw point cloud \'mat\' or \'pts\'')
    parser.add_argument('--scan_res', type=float, default='0.00063',
                        help='input the scan resolution (i.g., 5mm/10m: 0.0005)')
    parser.add_argument('--work_dir', type=str, default='./data/demo/results/',
                        help='the directory to save the results')
    parser.add_argument('--ds_size', type=int, default=10000000,
                        help='the size of sub-sampled point cloud for projection')
    parser.add_argument('--ds_method', type=str, default='dist_scope_ratio',
                        help='downsampling methods, including random, dist_scope_ratio, or dist_order_ascending')
    parser.add_argument('--num_img', type=int, default=4,
                        help='the number of sub-images for detection, 1 or 4')
    parser.add_argument('--model_path', type=str,
                        default='./weights/yolox_x_C',
                        help='give the directory of the trained model and the config file')
    args = parser.parse_args()
    return args


def main():
    """
    The whole pipeline of automatic target detection, including loading, projection, detection, and extraction
    """

    args = parse_args()
    # set a path and remove the old path for saving results for each individual scans
    args.work_dir = os.path.abspath(os.path.join(args.point_cloud, "..", 'results'))
    if os.path.exists(args.work_dir):
        shutil.rmtree(args.work_dir)

    """ 00. Point cloud loading """

    print("Start running ...")
    # support loading different formats of input point cloud, and 'suffix' help judge the format of point cloud
    suffix = os.path.splitext(args.point_cloud)[-1]
    if suffix == '.mat':
        f = h5py.File(args.point_cloud, 'r')
        point_cloud = np.array(f['data'])  # 'dataloader', 'dataload', or 'data'
        point_cloud = point_cloud.T
    elif suffix == '.pts':
        point_cloud = PyntCloud.from_file(args.point_cloud, sep=" ", header=1,
                                          names=['x', 'y', 'z', 'i', 'r', 'g', 'b'])
        point_cloud = np.asarray(point_cloud.points)
    else:
        raise ValueError('The input point cloud format is not supported, please convert it into \'.pts\' or \'.mat\'')

    """ 01. projection """
    # subsample the raw point cloud for projection
    ds_size = args.ds_size
    ds_method = args.ds_method

    if point_cloud.shape[0] > ds_size:
        sub_index = projection.downsampling(point_cloud, downsampling_size=ds_size, downsampling_type=ds_method)
        sub_pc = point_cloud[sub_index, :]
    else:
        sub_pc = point_cloud

    time_spending = time() - start_time
    # project sub-sampled point cloud to a panoramic image
    img = projection.pc2img(sub_pc, args.scan_res)

    time_spending = time() - start_time - time_spending
    print(f'The running time for projection is: {time_spending} seconds')
    print('-----------------------------Projection is done----------------------------')

    time_spending = time() - start_time
    print(f'The running time before detection is: {time_spending} seconds')

    # generate a sub-folder for saving projection results
    os.makedirs(os.path.join(args.work_dir, 'projection'), exist_ok=True)
    img_folder = os.path.join(args.work_dir, 'projection')
    # split the whole panoramic image
    num_image = args.num_img
    if num_image == 1:
        cv2.imwrite(os.path.join(img_folder, 'projected_img.jpg'), img)
    else:
        hw_half = (img.shape / np.sqrt(num_image)).astype(int)
        h = range(0, img.shape[0], hw_half[0] - 1)
        w = range(0, img.shape[1], hw_half[1] - 1)
        for i in range(num_image):
            # overlap_ratio = 0.0
            curr_h = i // int(np.sqrt(num_image))
            curr_w = i % int(np.sqrt(num_image))
            sub_img = img[h[curr_h]: h[curr_h+1] + (1 - curr_h) * 200, w[curr_w]: w[curr_w+1] + (1 - curr_w) * 200]
            cv2.imwrite(os.path.join(img_folder, f'projected_img_{i}.jpg'), sub_img)

    """ 02. detection """
    # load training config file and trained model
    detector = "yolox_x_C"
    if detector == "frcnn_C":
        config = os.path.join(args.model_path, 'config.py')
        checkpoint = os.path.join(args.model_path, 'detector.pth')
        _, bbox_result = detection.detection_two_stages(config, checkpoint, args.work_dir)
    elif detector == "yolox_x_C":
        config = os.path.join(args.model_path, 'yolox_x.py')
        checkpoint = os.path.join(args.model_path, 'best_ckpt.pth')
        _, bbox_result = detection.detection_one_stage(config, checkpoint, args.work_dir)

    # get the bbox image coordinate results, and save it to '.txt' file

    bbox_result = np.array(bbox_result)
    # add multi-image situation
    if num_image != 1:
        hw_half = (img.shape / np.sqrt(num_image)).astype(int)
        for i in range(bbox_result.shape[0]):

            curr_h = int(bbox_result[i, 0]) // int(np.sqrt(num_image))
            curr_w = int(bbox_result[i, 0]) % int(np.sqrt(num_image))
            # derive the global image coordinates of all sub-images, height and width direction
            bbox_result[i, 1::2] += hw_half[1] * curr_w
            bbox_result[i, 2::2] += hw_half[0] * curr_h
    bbox_result = bbox_result[:, 1:]
    np.savetxt(os.path.join(args.work_dir, 'detection/', 'bbox_image_coordinates.txt'), bbox_result)

    time_spending = time() - start_time
    print(f'The running time after detection is: {time_spending} seconds')
    print('-----------------------------Detection is done-----------------------------')

    """ 03. Extraction """
    # reproject bbox image coordinates to raw point clouds
    num_detected, addr = projection.img2pc(point_cloud, bbox_result, args.scan_res, args.work_dir)
    print('{} targets are saved to {}'.format(num_detected, addr))

    time_spending = time() - start_time - time_spending
    print(f'The running time for extraction is: {time_spending} seconds')
    print('-------------Extraction is done, and target coordinates are saved--------------')

    """ 04. Estimation """
    time_spending = time() - start_time
    print(f'The running time for processing current scan is: {time_spending} seconds')
    print('----------------------------All done, great job!---------------------------')


if __name__ == '__main__':
    main()
 
