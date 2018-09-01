import itertools

import numpy as np
import scipy.interpolate


def compute_hypercols(out, pix, imdim):
    """ Compute hypercolumns.

    Args:
        out (list of ndarrays): list of activations of a neural network
                                ndarray are 4D: (batch, height, width, channels)
                                and correspond to one layer
        pix (n x 2 ndarray): list of pixels
        imdim (int): height  of input image

    Returns:
        hc (3D ndarray): format is (image, pix, hypercol dimensions)
    """
    # number of dimensions per layer
    hc_dims = [l.shape[-1] for l in out]
    # cummulative sum -> to keep track of where in vector we are adding stuff
    hc_cumdims = [0] + list(np.cumsum(hc_dims))
    # dimensions: image, pixel, hypercol
    hc = np.zeros([len(out[0]), len(pix), sum(hc_dims)])

    for l, im in itertools.product(range(len(out)),
                                   range(len(out[0]))):
        # i,j in layer coordinates
        ij = pix*(out[l].shape[1]-1)/float(imdim-1)
        grid = [range(x) for x in out[l].shape[1:3]]
        hc_layer = np.zeros((len(pix), hc_dims[l]))
        for i in range(hc_dims[l]):
            interpn = scipy.interpolate.RegularGridInterpolator(grid,
                                                                out[l][im, :, :, i],
                                                                bounds_error=False,
                                                                 fill_value=None)
            hc_layer[:, i] = interpn(ij)
            # TODO: transpose?
        hc[im, :, hc_cumdims[l]:hc_cumdims[l+1]] = np.array(hc_layer)  # copy    

    return hc
