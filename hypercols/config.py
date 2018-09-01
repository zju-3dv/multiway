"""Semantic segmentation config system.

This file specifies default config options for Fast R-CNN. You should not
change values in this file. Instead, you should write a config file (in yaml)
and use cfg_from_file(yaml_file) to load it and override the default options.

Most tools in $ROOT/tools take a --cfg option to specify an override file.
    - See tools/{train,test}_net.py for example code that uses cfg_from_file()
    - See experiments/cfgs/*.yml for example YAML config override files
"""

import os.path as osp
import numpy as np
from easydict import EasyDict as edict

__C = edict()
# Consumers can get config by:
#   from fast_rcnn_config import cfg
cfg = __C

# Minibatch size
__C.BATCH_SIZE = 4

# Use a prefetch thread in data_layer.
# Reduces data loading from 1.4 to 0.07 when loading 15 images with (321, 321)
__C.USE_PREFETCH = True

# Pair data loader layer settings.
__C.PAIR = edict()
__C.PAIR.RANDOM_NEGATIVE = True
__C.PAIR.SCALE = True
__C.PAIR.MAX_IMAGE_SIZE = 150
__C.PAIR.MIN_IMAGE_SIZE = 50
__C.PAIR.PCK_RADIUS = 0.05 * __C.PAIR.MAX_IMAGE_SIZE
__C.PAIR.MARGIN = 1 * __C.PAIR.PCK_RADIUS

# Flowweb Dataset
__C.FLOWWEB = edict()
__C.FLOWWEB.ROOT_PATH = '/home/qian/dataset/FlowWeb/pascal/data'
__C.FLOWWEB.CLASSMASK = ['aeroplane', 'bicycle', 'boat', 'bottle', 'bus', 'car',
                         'chair', 'table', 'motorbike', 'sofa', 'train',
                         'tvmonitor']
__C.FLOWWEB.TRAIN_PORTION = 0.8  # Portion of train set.

# Caltech-UCSD Birds-200-2011 Dataset
__C.CUB2011 = edict()
__C.CUB2011.ROOT_PATH = '/home/qian/dataset/CUB_200_2011'

# Caltech-UCSD Birds-200-2011 Dataset
__C.SINTEL = edict()
__C.SINTEL.ROOT_PATH = '/home/qian/dataset/Sintel'
__C.SINTEL.TRAIN_PORTION = 0.8  # Portion of train set.
__C.SINTEL.IMG_SUBSAMPLE = 10

# Flowweb Dataset
__C.VOC2011 = edict()
__C.VOC2011.KPT_ROOT_DIR = '/home/qian/dataset/VOC2011_Keypoints'
__C.VOC2011.ROOT_DIR = '/home/qian/dataset/VOC2011'
__C.VOC2011.TRAIN_PORTION = 0.8  # Portion of train set.

# PASCAL3D+ Dataset
__C.PASCAL3D = edict()
__C.PASCAL3D.ROOT_PATH = '/home/qian/dataset/PASCAL3D+_release1.1/'

# resize the bounding box to fit inside the 150 pixels
__C.PASCAL3D.SCALE = True
__C.PASCAL3D.MAX_IMAGE_SIZE = 150
__C.PASCAL3D.MIN_IMAGE_SIZE = 50
__C.PASCAL3D.PCK_RADIUS = 0.05 * __C.PASCAL3D.MAX_IMAGE_SIZE
__C.PASCAL3D.MARGIN = 1 * __C.PASCAL3D.PCK_RADIUS
__C.PASCAL3D.RANDOM_NEGATIVE = True
__C.PASCAL3D.SIMILAR_VIEW = True
__C.PASCAL3D.VIEW_THRESH = 30
__C.PASCAL3D.PHASE = 'train'

# KITTI_FLOW dataset
__C.KITTI_FLOW = edict()
__C.KITTI_FLOW.ROOT_PATH = '/home/qian/dataset/KITTI/flow'

# Catface dataset
__C.CATFACE = edict()
__C.CATFACE.ROOT_PATH = '/home/qian/dataset/Catface'

# Catface dataset
__C.CATFACEGT = edict()
__C.CATFACEGT.ROOT_PATH = '/home/qian/dataset/Catfacegt'

# Correspondence options
__C.CORRESPONDENCE = edict()

# max distance between query frame and reference frame
__C.CORRESPONDENCE.PAIR_SEARCH_DIST = 30

# minimum number of correspondences to be a valid correspondence pair for
# training
__C.CORRESPONDENCE.NUM_CORRES_THRESH = 200

# Ratio of negative correspondences in a batch
__C.CORRESPONDENCE.NEG_PAIR_RATIO = .5

# Limit the number of correspondence pairs per training.
__C.CORRESPONDENCE.MAX_NUM_CORRESPONDENCE = 1000

# For alignment to work, we set (we choose 32x so as to be able to evaluate
# the model for all different subsampling sizes):
# (1) input dimension equal to
# $n = 32 * k - 31$, e.g., 321 (for k = 11)
# Dimension after pooling w. subsampling:
# (16 * k - 15); (8 * k - 7); (4 * k - 3); (2 * k - 1); (k).
# For k = 11, these translate to
#           161;          81;          41;          21;  11
# height width
__C.IM_SHAPE = (321, 769)
# __C.IM_SHAPE = (321, 769)
# __C.IM_SHAPE = (353, 1217)

# Visualize correspondence using colored points
__C.CORRESPONDENCE.COLOR_CORRESPONDENCE = False

__C.CORRESPONDENCE.NUM_NEG_RND_POINTS = 50

__C.CORRESPONDENCE.DOWNSAMPLING_FACTOR = 8
__C.CORRESPONDENCE.KP_DOWNSAMPLING_FACTOR = 4

# Minimum pixel distance for random negative correspondence
# sqrt(8^2 + 8^2) = 11.31. Minimum diagonal distance when x8 subsampled
# __C.CORRESPONDENCE.MIN_PIXEL_DIST = 12
__C.CORRESPONDENCE.MIN_PIXEL_DIST = 20

__C.CORRESPONDENCE.PCK_RADII = [10]
# __C.CORRESPONDENCE.PCK_RADII = range(1, 101)
__C.CORRESPONDENCE.PCK_REFERENCE_THRESH = 3


# Hard negative mining
__C.CORRESPONDENCE.NEGATIVE_MINING = False
__C.CORRESPONDENCE.NEGATIVE_MINING_PIXEL_THRESH = 4

# (margin - dist)^2 as negative loss
__C.CORRESPONDENCE.MARGIN_DIST_SQ = True

# Use bilinear interpolation when evaluating
__C.CORRESPONDENCE.BILINEAR_INTERPOLATION = False

# Do post processing
__C.CORRESPONDENCE.PROCESS_HARD_NEGATIVE = True

# When the data layer only generate positive samples and hard negative
# IntraClassCorrespondenceDataLayer sets False
# VelodyneKITTIT sets True
__C.CORRESPONDENCE.NEGATIVE_PROCESSING_LAYER_POSITIVE_SAMPLES_ONLY = True

# Image Resize
__C.CORRESPONDENCE.SCALE = 1.


#
# Training options
#

__C.TRAIN = edict()



# Use horizontally-flipped images during training?
__C.TRAIN.USE_FLIPPED = True

# Iterations between snapshots
__C.TRAIN.SNAPSHOT_ITERS = 2000

# solver.prototxt specifies the snapshot path prefix, this adds an optional
# infix to yield the path: <prefix>[_<infix>]_iters_XYZ.caffemodel
__C.TRAIN.SNAPSHOT_INFIX = ''

# Kill the training process if the loss goes higher than the threshold
__C.TRAIN.LOSS_THRESHOLD = 1000


#
# MISC
#

# The mapping from image coordinates to feature map coordinates might cause
# some boxes that are distinct in image space to become identical in feature
# coordinates. If DEDUP_BOXES > 0, then DEDUP_BOXES is used as the scale factor
# for identifying duplicate boxes.
# 1/16 is correct for {Alex,Caffe}Net, VGG_CNN_M_1024, and VGG16
# __C.DEDUP_BOXES = 1./16.

# Pixel mean values (BGR order) as a (1, 1, 3) array
# These are the values originally used for training VGG16
# __C.PIXEL_MEANS = np.array([[[102.9801, 115.9465, 122.7717]]])
# mean value of used for Deeplab
# __C.PIXEL_MEANS = np.array([[[104.008, 116.669, 122.675]]])
# mean value of the KITTI images
__C.PIXEL_MEANS = np.array([[[94.95408014, 99.54709961, 94.48018654]]])

# For reproducibility
__C.RNG_SEED = 3

# A small number that's used many times
__C.EPS = 1e-14

# Root directory of project
__C.ROOT_DIR = osp.abspath(osp.join(osp.dirname(__file__), '..'))

# Place outputs under an experiments directory
__C.EXP_DIR = 'default'

__C.NUM_LABELS = 9 # number of labels

# ground truth segmentations are given only for the frames that are multiples
# of the following number. This is used in the correspondence.py
__C.SEGMENTATION_INTERVAL = 10 # must be an integer

# IMG scale. in the image transformation layer, scale the mean substracted raw
# image values by this value.
# -1.0: no scale
__C.IMG_SCALE = -1.0

# Use default interactive backend if True, Agg backend if False
__C.INTERACTIVE_PLOT = False

def get_output_dir(net, exp_dir=None):
    """Return the directory where experimental artifacts are placed.

    A canonical path is built using the name from an imdb and a network
    (if not None).
    """
    if exp_dir is None:
        path = osp.abspath(osp.join(__C.ROOT_DIR, 'output', __C.EXP_DIR))
    else:
        path = osp.abspath(osp.join(__C.ROOT_DIR, 'output', exp_dir))
    if net is None:
        return path
    else:
        return osp.join(path, net.name)

def _merge_a_into_b(a, b):
    """Merge config dictionary a into config dictionary b, clobbering the
    options in b whenever they are also specified in a.
    """
    if type(a) is not edict:
        return

    for k, v in a.iteritems():
        # a must specify keys that are in b
        if not b.has_key(k):
            raise KeyError('{} is not a valid config key'.format(k))

        # the types must match, too
        if type(b[k]) is not type(v):
            raise ValueError(('Type mismatch ({} vs. {}) '
                              'for config key: {}').format(type(b[k]),
                                                           type(v), k))

        # recursively merge dicts
        if type(v) is edict:
            try:
                _merge_a_into_b(a[k], b[k])
            except:
                print('Error under config key: {}'.format(k))
                raise
        else:
            b[k] = v

def cfg_from_file(filename):
    """Load a config file and merge it into the default options."""
    import yaml
    with open(filename, 'r') as f:
        yaml_cfg = edict(yaml.load(f))

    _merge_a_into_b(yaml_cfg, __C)

def cfg_from_list(cfg_list):
    """Set config keys via list (e.g., from command line)."""
    from ast import literal_eval
    assert len(cfg_list) % 2 == 0
    for k, v in zip(cfg_list[0::2], cfg_list[1::2]):
        key_list = k.split('.')
        d = __C
        for subkey in key_list[:-1]:
            assert subkey in d.keys()
            d = d[subkey]
        subkey = key_list[-1]
        assert subkey in d.keys()
        try:
            value = literal_eval(v)
        except:
            # handle the case when v is a string literal
            value = v
        assert type(value) == type(d[subkey]), \
            'type {} does not match original type {}'.format(
            type(value), type(d[subkey]))
        d[subkey] = value
