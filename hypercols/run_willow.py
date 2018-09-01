import glob
import numpy as np
import tensorflow as tf
from scipy.misc import imread, imresize
import scipy.io
from myalexnet_forward import network, run
from hypercol import compute_hypercols

RESIZEFACTOR = 1.0
PATH_LIST = glob.glob('../dataset/WILLOW-ObjectClass/*')

for pname in PATH_LIST:

    fnames = glob.glob(pname+'/*.png')

    # load and resize images (will also save them afterwards)
    ims = [imresize(imread(f)[:,:,:3], RESIZEFACTOR).astype('float32') for f in fnames]
    
    # load all desired pixel positions
    kptsfiles = [f[0:-3]+'mat' for f in fnames]
    kpts = [scipy.io.loadmat(f)['pts_coord'] for f in kptsfiles]
    # resize and convert from x,y to i,j
    kpts = [RESIZEFACTOR*k[[1,0],:] for k in kpts]
    
    layers = ['c4', 'c5']
    
    tf.reset_default_graph()
    net = network()
    
    # we can't send batches because each image has different dimensions
    for im, k, f in zip(ims, kpts, fnames):
        out = run([im], net, layers)
        hc = np.squeeze(compute_hypercols(out, k.T, im.shape[0]))
        hc = (np.transpose(hc)).astype('float32')
        scipy.io.savemat(f + '.hypercols_kpts.mat',
                         {'desc': hc,
                          'frame': k[[1,0], :],
                          'nfeature': hc.shape[1],
                          'img': im.astype('uint8'),
                          'filename': f[len(pname)+1:]
                          })
