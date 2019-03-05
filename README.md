# 

# Multi-Object Matching
This repository contains code for the following methods of multi-object matching:

1. **The proposed method**: [Multi-Image Semantic Matching by Mining Consistent Features](https://arxiv.org/abs/1711.07641), CVPR 2018 by Qianqian Wang, [Xiaowei Zhou](http://www.cad.zju.edu.cn/home/xzhou/) and [Kostas Daniilidis](http://www.cis.upenn.edu/~kostas/).
2. **Spectral method**: ["Solving the multi-way
matching problem by permutation synchronization"](http://people.cs.uchicago.edu/~risi/papers/PachauriKondorSinghNIPS2013.pdf), NIPS 2013.
3. **MatchLift**: ["Near-optimal joint object
matching via convex relaxation"](http://proceedings.mlr.press/v32/chend14.pdf), ICML 2014.
4. **MatchALS**: ["Multi-Image Matching via Fast Alternating Minimization"](http://www.cis.upenn.edu/~kostas/mypub.dir/xiaowei15iccv.pdf), CVPR 2015.

## Tutorial slides
To learn more about multi-object matching, please refer to [multiway-slides](https://www.dropbox.com/s/qsun5g4snw7jo5y/multiway-release.pdf?dl=0)

## Usage
1. Download [WILLOW-ObjectClass Dataset](http://www.di.ens.fr/willow/research/graphlearning/) at ```dataset/```
```
cd dataset/
wget http://www.di.ens.fr/willow/research/graphlearning/WILLOW-ObjectClass_dataset.zip
unzip WILLOW-ObjectClass_dataset.zip
# remove problematic image and annotation
rm -f WILLOW-ObjectClass/Face/image_0160.*
# there is an annotation error in Cars_030a.mat (coordinate swap between 6th and 7th keypoint), replace it with the correct one
mv Cars_030a.mat WILLOW-ObjectClass/Car/
```
2. Download [Alexnet Weights](http://www.cs.toronto.edu/~guerzhoy/tf_alexnet/bvlc_alexnet.npy) at ```hypercols/```, and then extract feature descriptor ***hypercolumn*** from AlexNet.
```
cd ../hypercols/
wget http://www.cs.toronto.edu/~guerzhoy/tf_alexnet/bvlc_alexnet.npy
python run_willow.py
```
3. run ```testWillow.m``` to test the code on WILLOW-ObjectClass dataset.

## Citation
If you find this code useful for your research, please cite the following paper:
```
@inproceedings{zhou2015multi,
  title={Multi-image matching via fast alternating minimization},
  author={Zhou, Xiaowei and Zhu, Menglong and Daniilidis, Kostas},
  booktitle={ICCV},
  year={2015}
}
@inproceedings{wang2018multi,
  title={Multi-Image Semantic Matching by Mining Consistent Features},
  author={Wang, Qianqian and Zhou, Xiaowei and Daniilidis, Kostas},
  booktitle={CVPR},
  year={2018}
}
```
If you have any questions, please contact qw246@cornell.edu




