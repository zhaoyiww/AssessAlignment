# https://muyuuuu.github.io/2021/05/11/MMDetection-use/
import json
import os
from tqdm import tqdm
import numpy as np
import cv2
import sys
import time
from loguru import logger
import torch


def detection_two_stages(config: str, checkpoint: str, work_dir):
    sys.path.append('../mmdetection')
    from mmdet.apis import init_detector, inference_detector, show_result_pyplot

    # load config file and saved model
    config_file = config
    checkpoint_file = checkpoint
    device = 'cuda:0'
    model = init_detector(config_file, checkpoint_file, device=device)

    results = []
    path_files = os.listdir(os.path.join(work_dir, 'projection/'))
    for i, file in zip(range(len(path_files)), tqdm(path_files)):
        result = inference_detector(model=model, imgs=work_dir + '/projection/' + file)
        img = os.path.join(work_dir, 'projection/', file)
        # show_result_pyplot(model, img, result)  # can add score_thr=0.5    # visualization
        # save the detected result
        if not os.path.exists(os.path.join(work_dir, 'detection')):
            os.makedirs(os.path.join(work_dir, 'detection'))
        # os.makedirs(os.path.join(work_dir, 'detection'))
        model.show_result(img, result, out_file=os.path.join(work_dir, 'detection/', f'detected_result{i}.jpg'))

        # derive the necessary information from the result
        for cate_id, items in enumerate(result, 1):
            #
            for item in items:
                item = item.tolist()
                x, y, w, h, s = item[0], item[1], item[2], item[3], item[4]
                d = {'image_id': i,
                     'category_id': cate_id,
                     'bbox': [x, y, w, h],
                     'score': s}
                results.append(d)

    # save the detection result to a '.json' file
    with open(os.path.join(work_dir, 'detection/', 'detection_results.json'), 'w') as f:
        json.dump(results, f, indent=4)

    # only keep the image_id and bbox information
    bbox_result = []
    for i in np.arange(len(results)):
        # bbox_result.append(results[i]['bbox'])
        bbox_result.append(np.hstack((results[i]['image_id'], results[i]['bbox'])))

    return results, bbox_result


def get_image_list(path):
    IMAGE_EXT = [".jpg", ".jpeg", ".webp", ".bmp", ".png"]
    image_names = []
    for maindir, subdir, file_name_list in os.walk(path):
        for filename in file_name_list:
            apath = os.path.join(maindir, filename)
            ext = os.path.splitext(apath)[1]
            if ext in IMAGE_EXT:
                image_names.append(apath)
    return image_names


def image_demo(predictor, work_dir):
    results = []
    path_files = os.listdir(os.path.join(work_dir, 'projection/'))

    for i, file in zip(range(len(path_files)), tqdm(path_files)):
        outputs, img_info = predictor.inference(os.path.join(work_dir, "projection", file))
        img_ratio = img_info["ratio"]
        result_image = predictor.visual(outputs[0], img_info, predictor.confthre)

        # save the detected result
        os.makedirs(os.path.join(work_dir, 'detection'), exist_ok=True)
        cv2.imwrite(os.path.join(work_dir, "detection", f'detected_result{i}.jpg'), result_image)

        result = outputs[0].detach().cpu().numpy()
        # derive the necessary information from the result
        for cate_id, items in enumerate(result, 1):
            items = items.tolist()
            h0, w0, h1, w1, s = items[0] / img_ratio, items[1] / img_ratio, \
                                items[2] / img_ratio, items[3] / img_ratio, items[4] * items[5]
            d = {'image_id': i,
                 'category_id': cate_id,
                 'bbox': [h0, w0, h1, w1],
                 'score': s}
            results.append(d)

    # save the detection result to a '.json' file
    with open(os.path.join(work_dir, 'detection/', 'detection_results.json'), 'w') as f:
        json.dump(results, f, indent=4)

    # only keep the image_id and bbox information
    bbox_result = []
    for i in np.arange(len(results)):
        # bbox_result.append(results[i]['bbox'])
        bbox_result.append(np.hstack((results[i]['image_id'], results[i]['bbox'])))

    return results, bbox_result


# source: YOLOX
class Predictor(object):
    def __init__(
        self,
        model,
        exp,
        cls_names,
        trt_file=None,
        decoder=None,
        device="cpu",
        fp16=False,
        legacy=False,
    ):
        from yolox.data.data_augment import ValTransform
        self.model = model
        self.cls_names = cls_names
        self.decoder = decoder
        self.num_classes = exp.num_classes
        self.confthre = exp.test_conf
        self.nmsthre = exp.nmsthre
        self.test_size = exp.test_size
        self.device = device
        self.fp16 = fp16
        self.preproc = ValTransform(legacy=legacy)

    def inference(self, img):
        from yolox.utils import postprocess
        img_info = {"id": 0}
        if isinstance(img, str):
            img_info["file_name"] = os.path.basename(img)
            img = cv2.imread(img)
        else:
            img_info["file_name"] = None

        height, width = img.shape[:2]
        img_info["height"] = height
        img_info["width"] = width
        img_info["raw_img"] = img

        ratio = min(self.test_size[0] / img.shape[0], self.test_size[1] / img.shape[1])
        img_info["ratio"] = ratio

        img, _ = self.preproc(img, None, self.test_size)
        img = torch.from_numpy(img).unsqueeze(0)
        img = img.float()
        if self.device == "gpu":
            img = img.cuda()
            if self.fp16:
                img = img.half()  # to FP16

        with torch.no_grad():
            t0 = time.time()
            outputs = self.model(img)
            if self.decoder is not None:
                outputs = self.decoder(outputs, dtype=outputs.type())
            outputs = postprocess(
                outputs, self.num_classes, self.confthre,
                self.nmsthre, class_agnostic=True
            )
            logger.info("Infer time: {:.4f}s".format(time.time() - t0))
        return outputs, img_info

    def visual(self, output, img_info, cls_conf=0.35):
        from yolox.utils import vis
        ratio = img_info["ratio"]
        img = img_info["raw_img"]
        if output is None:
            return img
        output = output.cpu()

        bboxes = output[:, 0:4]

        # preprocessing: resize
        bboxes /= ratio

        cls = output[:, 6]
        scores = output[:, 4] * output[:, 5]

        vis_res = vis(img, bboxes, scores, cls, cls_conf, self.cls_names)
        return vis_res


def detection_one_stage(config: str, checkpoint: str, work_dir):
    sys.path.append('/scratch/zhawang/projects/AssessAlignment-2/YOLOX')
    from yolox.data.data_augment import ValTransform
    from yolox.utils import fuse_model, get_model_info, postprocess, vis
    # from mmdet.apis import init_detector, inference_detector, show_result_pyplot
    from yolox.exp import get_exp
    import torch

    # load config file and saved model
    device = 'cuda:0'
    exp = get_exp(config)
    model = exp.get_model()

    model.eval()
    model.to(device)

    ckpt = torch.load(checkpoint, map_location="cpu")
    model.load_state_dict(ckpt["model"])
    logger.info("loaded checkpoint done.")

    COCO_CLASSES = ('bota', 'checkerboard', 'leicaround', 'rotated_checker')
    predictor = Predictor(model, exp, COCO_CLASSES, None, None, 'gpu', False, False)
    results, bbox_result = image_demo(predictor, work_dir)

    return results, bbox_result


if __name__ == "__main__":
    config = "/scratch/zhawang/projects/AssessAlignment/weights/yolox_x_C/yolox_x.py"
    checkpoints = "/scratch/zhawang/projects/AssessAlignment/weights/yolox_x_C/best_ckpt.pth"
    work_dir = "/scratch/zhawang/projects/AssessAlignment/data/results"
    detection_one_stage(config, checkpoints, work_dir)
