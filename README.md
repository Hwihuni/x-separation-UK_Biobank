
# x-separation application to UK biobank dataset
This repository contains the X-separation application for UK biobank datasets.
> H. Jeong, S. Oh, J. Lee, and H. Shin, _Application of susceptibility source separation (chi-separation) to UK Biobank protocol and clinical protocol using deep neural network_, ISMRM 2022,
> [[abstract]](https://archive.ismrm.org/2022/2364.html)



## UKB dataset

The whole data of UK biobank can be accessed by https://www.ukbiobank.ac.uk/

## Pretrained models

The checkpoints can be download form google drive link in /python/r2_r2star_mapping/checkpoints/download.txt

## Requirements
```
MATLAB: R2021b
PYTHON: use env.yml file in python/r2_r2star_mapping
FSL built with wsl
```

## ONNX
You can easily run the code in the Pipeline_onnx!

## In the case you do not have ONNX library
* change 'path' in Biobank_x_sep_step0_Dicomprocessing.m, Biobank_x_sep_step1_GREprocessing.m, Biobank_x_sep_step2_DLpreprocessing.m, and Biobank_x_sep_step3_reconstrction.m
* run Biobank_x_sep_step0_Dicomprocessing.m, Biobank_x_sep_step1_GREprocessing.m and Biobank_x_sep_step2_DLpreprocessing.m
* move Data_to_gpu_t2map.mat and Data_to_gpu_t2starmap.mat to python\r2_r2star_mapping\data
* run predict_r2.py and predict_r2star.py
* move inference file to 'inf_from_gpu' which is made in your directory
* run Biobank_x_sep_step3_reconstrction.m
