import glob
import numpy as np
import itertools
import tensorflow as tf
from config import cfg
from scipy.misc import imread, imresize
import scipy.io
import os
import cv2
from myalexnet_forward import network, run
from hypercol import compute_hypercols
import random
import os.path as osp

RESIZEFACTOR = 1.0
#PATH_NAME = 'dataset/WILLOW-ObjectClass/Car/'
testset = 'validation1'
PATH_NAME = '/home/qian/dataset/Catface/'+testset


def transform(img1, img2, coord1, coord2):
    # CROP and resize
    def crop_and_resize(im, coord):
        if im.ndim == 2:
            im = np.dstack([im] * 3)

        h, w, _ = im.shape

        max_size = float(cfg.PAIR.MAX_IMAGE_SIZE)
        min_size = float(cfg.PAIR.MIN_IMAGE_SIZE)

        if cfg.PAIR.SCALE:
            scale = False
            if max(w, h) > max_size:
                scale = max_size / max(w, h)
            elif min(w, h) < min_size:
                scale = min_size / min(w, h)
                if max(scale * w, scale * h) > max_size:
                    scale *= max_max(scale * w, scale * h)

            if scale > 0:
                im = cv2.resize(im, (int(scale * w), int(scale * h)),
                                interpolation=cv2.INTER_CUBIC)
                coord *= scale

        return im, coord

    img1, coord1 = crop_and_resize(img1, coord1)
    img2, coord2 = crop_and_resize(img2, coord2)

    return img1, img2, coord1, coord2

def ind2coord(ind, width):
    """ takse in num_batch x num_coords x k ind and genertes num_batch x
    num_coords x 2"""
    # k-th NN index: ind[:, k, :]
    x = np.floor(ind / width)
    y = ind % width
    xy_coords = np.concatenate((x[..., np.newaxis], y[..., np.newaxis]), axis=1)
    return xy_coords

if __name__ == '__main__':

    fnames = glob.glob(PATH_NAME+'/*.jpg')

    # load and resize images (will also save them afterwards)
    ims = [imread(f)[:,:,:3].astype('float32') for f in fnames]

    pairs = [c for c in itertools.combinations(range(len(fnames)), 2)]
    random.shuffle(pairs)

    layers = ['c1', 'c2', 'c3', 'c4', 'c5']
    tf.reset_default_graph()
    net = network()

    pck_all = np.zeros((100, 0))
    count = 0
    num_pts = []
    for pairidx in xrange(len(pairs)):
        fnm1 = fnames[pairs[pairidx][0]]
        fnm2 = fnames[pairs[pairidx][1]]
        img1 = imread(fnm1).astype(np.float64)
        img2 = imread(fnm2).astype(np.float64)

        coord1_fnm = os.path.splitext(fnm1)[0] + '.mat'
        coord1 = scipy.io.loadmat(coord1_fnm)

        coord2_fnm = os.path.splitext(fnm2)[0] + '.mat'
        coord2 = scipy.io.loadmat(coord2_fnm)

        mask1 = coord1['keypts_status'][0]
        mask2 = coord2['keypts_status'][0]
        mask = (mask1 & mask2).astype(bool)

        kp1, kp2 = coord1['keypts'].astype(np.float64), coord2['keypts'].astype(np.float64)
        kp1, kp2 = kp1[:, mask], kp2[:, mask]

        img1, img2, kp1, kp2 = transform(img1, img2, kp1, kp2)
	num_pts.append(len(kp1))
        # resize and convert from x,y to i,j
        kp1 = RESIZEFACTOR * kp1[[1, 0], :]
        kp2 = RESIZEFACTOR * kp2[[1, 0], :]
    
         # we can't send batches because each image has different dimensions
        out1 = run([img1], net, layers)
        hc1 = np.squeeze(compute_hypercols(out1, kp1.T, img1.shape[0]))
        hc1 = (np.transpose(hc1)).astype('float32')

        kp2_all = np.zeros([img2.shape[0]*img2.shape[1],2])
        for i in xrange(img2.shape[0]):
            kp2_all[i * img2.shape[1]:(i + 1) * img2.shape[1], 0] = i * np.ones(img2.shape[1]).T
            kp2_all[i * img2.shape[1]:(i + 1) * img2.shape[1], 1] = np.arange(img2.shape[1]).T
        out2 = run([img2], net, layers)
        hc2_all = np.squeeze(compute_hypercols(out2, kp2_all, img2.shape[0]))
        hc2_all = (np.transpose(hc2_all)).astype('float32')

        knn_ind = []
        for j in xrange(kp1.shape[1]):
            hc_j_dist = []
            for k in xrange(hc2_all.shape[1]):
                hc_j_dist.append(np.sum((hc2_all[:,k] - hc1[:,j])**2, axis=0))
            knn_ind.append(hc_j_dist.index(min((hc_j_dist))))

        kp2_knn = ind2coord(np.array(knn_ind), img2.shape[1]).T

        im_dist = np.sqrt(np.sum((kp2 - kp2_knn) ** 2, axis=0))

        pck_radii = range(1, 101)
        pcks = np.zeros((len(pck_radii), 1))
        for i, pck_radius in enumerate(pck_radii):
            pcks[i] = np.sum(im_dist < pck_radius).astype(np.float)/len(im_dist)

        pck_all = np.hstack((pck_all, pcks))

        count += 1
        pck_all = pck_all[:, :count]

        print(count, np.average(pck_all, axis=1, weights=num_pts))

    pck_avg = np.average(pck_all, axis=1, weights=num_pts)

    # Save the network parameters
    if not osp.exists(output_dir):
        os.makedirs(output_dir)

    np.savez('pck_hypercol_%s.npz'%(testset), pck_avg)







