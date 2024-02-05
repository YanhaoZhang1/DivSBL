# Sparse Bayesian Learning: A Diversified Scheme
This set of Matlab (Vision R2021b) functions contain the core code to reproduce the results of the above paper.

* DivSBL.m  ---> the main function for Diversified Block Sparse Bayesian Learning.

**For convenience, ```'DivSBL.m' is directly provided in each demo folders.``` **

    # Please note that the value of the `PRUNE_GAMMA` parameter in 'DivSBL.m' function needs to be manually adjusted according to the data.

    # We have already set appropriate sizes for the data style in the demo. If you need to change the data, you may need to adjust its values accordingly.


## Set path

-- The following 5 sub-files in DivSBL_public folder need to be manually read by the user into the MATLAB path, for use in comparative experiments in the demo folders. 

    * BCS_fast_rvm.m        ---> The function file for SBL(RVM) algorithm from [1],

    * CVX                   ---> Matlab software for disciplined convex programming from [2],

    * PCSBL.m               ---> The function file for PC-SBL algorithm from [3],

    * spgl1-2.1             ---> A solver for large-scale sparse reconstruction from [4],

    * StructOMP Released 2  ---> The function file for StructOMP algorithm from [5].




## Folder 1: demo_Synthetic signal 

    * DEMO_Synthetic_Signal.m   ---> generate a compressed sensing demo for synthetic signal using 7 different algorithms.
    * DivSBL.m                  ---> The main function for Diversified Block Sparse Bayesian Learning (DivSBL) algorithm.
    * BSBL.m                    ---> The function file for BSBL algorithm from [6].



## Folder 2: demo_Audio

    * DEMO_audio.m    ---> generate a compressed sensing demo for AudioSet from [7].
    (Some audio in WAV format for testing purposes are also provided here.）
  
  
  
## Folder 3: demo_Image

    * DEMO_image_public_parrots.m   ---> generate a compressed sensing demo for Parrot image.
    * DWTM.mat                      ---> A decrete wavelet transform matrix provided by [3].
    * CS_test_images                ---> Some classic test images.



##

**For bug reports, please contact me at email: yanhaozhang@buaa.edu.cn.**


Author: Yanhao Zhang.

Beihang University,  Jan, 27, 2024.

##

# References
[1] Ji, Shihao, Ya Xue, and Lawrence Carin. "Bayesian compressive sensing." IEEE Transactions on signal processing 56.6 (2008): 2346-2356.

[2] Grant, Michael, and Stephen Boyd. "CVX: Matlab software for disciplined convex programming, version 2.1." (2014).

[3] Fang, Jun, et al. "Pattern-coupled sparse Bayesian learning for recovery of block-sparse signals." IEEE Transactions on Signal Processing 63.2 (2014): 360-372.

[4] Van Den Berg, Ewout, and Michael P. Friedlander. "SPGL1: A solver for large-scale sparse reconstruction." (2007): 135.

[5] Huang, Junzhou, Tong Zhang, and Dimitris Metaxas. "Learning with structured sparsity." Proceedings of the 26th Annual International Conference on Machine Learning. 2009.

[6] Zhang, Zhilin, and Bhaskar D. Rao. "Extension of SBL algorithms for the recovery of block sparse signals with intra-block correlation." IEEE Transactions on Signal Processing 61.8 (2013): 2009-2015.

[7] Gemmeke, Jort F., et al. "Audio set: An ontology and human-labeled dataset for audio events." 2017 IEEE international conference on acoustics, speech and signal processing (ICASSP). IEEE, 2017.
