# ImSpiRE: <ins>Im</ins>age feauture-aided <ins>Sp</ins>at<ins>i</ins>al <ins>R</ins>esolution <ins>E</ins>nhancement method

![label1](https://img.shields.io/badge/version-v1.1.1-yellow)	![label2](https://img.shields.io/badge/license-MIT-green)



<!---**| [Overview](#overview) | [Installation](#installation) | [Quick Start](#quick-start) | [Parameter Details](#parameter-details) | [Run ImSpiRE in Python interface](#run-imspire-in-python-interface)  | [Citation](#citation) |** --->

**| [Overview](#overview) | [Installation](#installation) | [Quick Start](#quick-start) | [Parameter Details](#parameter-details) | [Run ImSpiRE in Python interface](#run-imspire-in-python-interface) |**



![Figure1](https://github.com/Yizhi-Zhang/MarkdownPicture/raw/master/ImSpiRE/Fig2_v1.1.2.png)

## Overview

ImSpiRE is an image feature-aided spatial resolution enhancement method for the ISC spatial transcriptome. Without any additional data except a parallel histology image, ImSpiRE constructs the subspot resolution transcriptional profiles while imputing gene expression of the unmeasured tissue regions. The basic idea of ImSpiRE is that the parallel histology images (usually H&E-stained images) provide valuable information (such as cell morphology, staining intensity and texture), and ImSpiRE leverages these image-derived features to enhance resolution. For spatial transcriptomes with high-pixel-density parallel histology images, ImSpiRE-enhanced profiles can achieve single-cell resolution.

ImSpiRE is a command line tool written in Python 3.8. It is available as a docker image, which eliminates the cumbersome installation and allows for portability across different computing systems. We also provide a Python library to meet the needs of other users.

## Installation

### Installation as a docker image

For convenience, we strongly recommend that users install ImSpiRE via the [Docker platform](https://www.docker.com/). Docker is available on a variety of [Linux platforms](https://docs.docker.com/desktop/install/linux-install/), [macOS](https://docs.docker.com/desktop/install/mac-install/) and [Windows 10](https://docs.docker.com/desktop/install/windows-install/) through Docker Desktop, and as a [static binary installation](https://docs.docker.com/engine/install/binaries/). You can find your preferred operating system [here](https://docs.docker.com/engine/install/).

After installing Docker, you can use the command below to pull the [Docker image of ImSpiRE](https://hub.docker.com/r/tongjizhanglab/imspire). It will take a few minutes to download the Docker image. 

   ```shell
   $ docker pull tongjizhanglab/imspire:1.1.2
   ```

You can then start a container based on this image and run ImSpiRE within it.

   ```shell
   $ docker run -it --name `whoami`_imspire -v ~:/home tongjizhanglab/imspire:1.1.2 /bin/bash
   ```

 You can also start a container with the IMAGE ID of ImSpiRE (`fbXXXXXXXX`).

   ```shell
   $ docker run -it --name `whoami`_imspire -v ~:/home fbXXXXXXX /bin/bash
   ```

### Installation as a Python library

ImSpiRE has been packaged and made available on [PyPI](https://pypi.org/project/imspire/). Before you begin, make sure you have pip3 on hand. Python's package installer is pip3. If you don't have pip3, you can get it [here](https://pip.pypa.io/en/stable/).

You may also need to prepare conda environment for ImSpiRE.

   ```shell
   $ conda create -n imspire python=3.8
   $ conda activate imspire 
   ```

ImSpiRE and other related packages can then be installed with a single command.

   ```shell
   $ pip3 install imspire 
   ```

Wait a minute and type the following command to check if ImSpiRE has been successfully installed.

   ```shell
   $ ImSpiRE -h
   ```

ImSpiRE includes an image feature extraction function. However, if users want to extract image features using CellProfiler, it must be installed separately. CellProfiler 4 should be pip installable in Python 3.8+ when a number of prerequisite packages are installed. More details can be found [here](https://github.com/CellProfiler/CellProfiler/wiki/Source-installation-(Linux)). 

   ```shell
   $ pip3 install cellprofiler==4.2.1
   ```

CellProfiler 4.2.1 can also be installed from source code.

   ```shell
   $ wget -O core-4.2.1.tar.gz https://github.com/CellProfiler/core/archive/refs/tags/v4.2.1.tar.gz
   $ wget -O CellProfiler-4.2.1.tar.gz https://github.com/CellProfiler/CellProfiler/archive/refs/tags/v4.2.1.tar.gz
   $ tar -xzvf CellProfiler-4.2.1.tar.gz
   $ tar -xzvf core-4.2.1.tar.gz
   ```

Due to problems with package dependencies, we recommend modifying the `setup.py` files of cellprofiler and cellprofiler-core to change the `pyzmq==18.0.1` to `pyzmq==18.1.1`, as they say too [here](https://github.com/CellProfiler/CellProfiler/issues/4462). You can now install cellprofiler-core and cellprofiler in order.

   ```shell
   $ cd core-4.2.1
   $ pip3 install .
   $ cd CellProfiler-4.2.1
   $ pip3 install .
   ```

You may need to use the command below to add the default installation path of pip3 to your system path.

   ```shell
   $ export PATH=~/.local/bin:$PATH
   ```

Then, type the command below to check whether CellProfiler has been installed successfully.

   ```shell
   $ cellprofiler -h
   ```

We note that you may encounter a failure to install wxPython 4.1.0 when installing CellProfiler. You can download the corresponding [WHL file](https://extras.wxpython.org/wxPython4/extras/) and then install it. As an example, for Ubuntu 18.0.4, you can install wxPython 4.1.0 with the following code.

   ```shell
   $ wget https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-18.04/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
   $ pip3 install wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
   ```

The CellProfiler piplines used in ImSpiRE can be available by the following commands.

   ```shell
   $ wget https://github.com/TongjiZhanglab/ImSpiRE/raw/master/CellProfiler_pipelines/cellprofiler_pipelines.zip
   $ unzip cellprofiler_piplines.zip
   ```

Furthermore, ImSpiRE employs two methods for extracting the foreground, one of which makes use of the Python package `backgroundremover`. Please keep in mind that when you run the program for the first time, it will check to see if you have the u2net models; if not, it will get them from u2net's google drive, as stated [here](https://github.com/xuebinqin/U-2-Net#usage-for-salient-object-detection). Download and place the pre-trained model `u2net.pth` in the directory `/.u2net`.

## Quick Start

### Step 1. Input preparation

ImSpiRE utilizes the count file in tab-delimited format or hierarchical-data format (HDF5 or H5) and the image file in TIFF format, as well as a file containing spot coordinates as input.

We provided a small [test dataset](https://github.com/TongjiZhanglab/ImSpiRE/tree/master/test/test_data) containing the raw count matrix, image and spot coordinates. A CellProfiler pipeline is also included in the test dataset for users if required. You can also type the command below to download.

   ```shell
   $ wget https://github.com/TongjiZhanglab/ImSpiRE/raw/master/test/test_data.zip
   $ unzip test_data.zip
   ```

### Step 2. Run ImSpiRE

Type the command below to run ImSpiRE.

   ```shell
   $ ImSpiRE -i test_data/ -c test_data/count_matrix.tsv -s test_data/image.tif -p ST -n test_output -O
   ```

Or use CellProfiler to extract image features.

   ```shell
   $ ImSpiRE -i test_data/ -c test_data/count_matrix.tsv -s test_data/image.tif -p ST -n test_output -m 2 --CellProfilerParam_Pipeline Cellprofiler_Pipeline_HE.cppipe -O
   ```

### Step 3. Output

The contents of the output directory in tree format will be displayed as described below, including the ImSpiRE-enhanced resolution profile stored in [Anndata h5ad file format](https://anndata.readthedocs.io/en/stable/index.html), the text file containing patch coordinates, the spot subimages and patch subimages if CellProfiler is used to extract image features, the image features, the matrices involved in OT and other supplementary results.

```shell
PATH/ProjectName
├── ProjectName_ResolutionEnhancementResult.h5ad  ## the ImSpiRE-enhanced resolution profile
├── ProjectName_PatchLocations.txt  ## the text file containing patch coordinates
├── ImageResults  ## the spot and patch subimages if CellProfiler is used to extract image features
│   ├── SpotImage
│   └── PatchImage
├── FeatureResults  ## the image features
│   ├── SpotFeature
│   └── PatchFeature
└── SupplementaryResults  ## the matrices involved in OT and other supplementary results
    ├── ot_C1_alpha_beta_epsilon.npy
    ├── ot_C2_alpha_beta_epsilon.npy
    ├── ot_M_alpha_beta_epsilon.npy
    ├── ot_T_alpha_beta_epsilon.npy
    └── ...
```

## Parameter Details

The parameter details of ImSpiRE are as follows:

```
usage: ImSpiRE [-h] <-i InputDir> <-c filtered_feature_bc_matrix.h5> <-s image.tif> [-n ProjectName] [-o OutputDir] [options]

ImSpiRE is an Image feature-aided Spatial Resolution Enhancenment method for in situ capturing (ISC) spatial transcriptome.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show version number of ImSpiRE and exit
  -i INPUT_DIR, --Input_Dir INPUT_DIR
                        This indicates the path to the directory for input datafiles. For a Visium
                        dataset, it is similar to the standard output format of Space Ranger,
                        which should contain a count file in the {Input_Dir} and a text file named
                        "tissue_positions_list.csv" that describes the spot locations in the
                        {Input_Dir}/spatial/. For a spatial transcriptomics (ST) dataset, it
                        should contain a count file and a text file named "pxl_pos.txt" that
                        describes the spot locations in the {Input_Dir}. Note that the file
                        "tissue_positions_list.csv" does not contain a header column, while the
                        file "pxl_pos.txt" includes a header column and contains two columns,
                        which represent the pixel coordinates of rows and columns of spots,
                        respectively.
  -c INPUT_COUNT_FILE, --Input_Count_File INPUT_COUNT_FILE
                        The input count file. It would typically be
                        "filtered_feature_bc_matrix.h5" and "count_matrix.tsv" for 10X Visium and
                        ST, respectively. Note that the count file of ST platform should be tab-
                        delimited.
  -s INPUT_SECTION_IMAGE_FILE, --Input_Section_Image_File INPUT_SECTION_IMAGE_FILE
                        The high-resolution tissue image, which should be a Tagged Image File
                        Format (TIF or TIFF) file.
  -t {H&E,IF}, --Input_Image_Type {H&E,IF}
                        The types of staining, including Haematoxylin & Eosin (H&E) staining and
                        Immunofluorescence (IF) staining. DEFAULT: "H&E".
  -o OUTPUT_DIR, --Output_Dir OUTPUT_DIR
                        The project directory, which is used to save all output files. DEFAULT:
                        "./".
  -n OUTPUT_NAME, --Output_Name OUTPUT_NAME
                        The project name, which is used to generate output file names as a prefix.
                        DEFAULT: "ImSpiRE".
  -p {Visium,ST}, --Platform {Visium,ST}
                        The platform that generates the dataset, Visium or ST. DEFAULT: "Visium".
  -m {1,2}, --Mode {1,2}
                        Two types of extracted image features. When this parameter is set to 1,
                        ImSpiRE will extract intensity and texture features of the image, which
                        are the objective features of the image itself. When this parameter is set
                        to 2, ImSpiRE will use CellProfiler to extract image features, which may
                        be more biologically significant. DEFAULT: 1.
  -O, --Overwriting     The switch of overwriting. If add the parameter "-O", ImSpiRE will
                        overwrite the exiting folders.
  --Verbose             The verbose flag. If add the parameter "--Verbose", ImSpiRE will verbose
                        the output.
  --Random_State RANDOM_STATE
                        Fix the seed for reproducibility. DEFAULT: 0.
  --Switch_Preprocess {ON,OFF}
                        The switch of basic filtering of spots and genes, ON or OFF. DEFAULT:
                        "ON".
  --Threshold_MinCounts THRESHOLD_MINCOUNTS
                        The minimum number of counts required for a spot to pass filtering, which
                        is enabled only if "--Switch_SpotFilter" is "ON". DEFAULT: 100.
  --Threshold_MaxCounts THRESHOLD_MAXCOUNTS
                        The maximum number of counts required for a spot to pass filtering, which
                        is enabled only if "--Switch_SpotFilter" is "ON". DEFAULT: 10000.
  --Threshold_MitoPercent THRESHOLD_MITOPERCENT
                        The maximum count percent of mitochondrial genes required for a spot to
                        pass filtering, which is enabled only if "--Switch_SpotFilter" is "ON".
                        DEFAULT: 20.
  --Threshold_MinSpot THRESHOLD_MINSPOT
                        The minimum number of spots expressed required for a gene to pass
                        filtering, which is enabled only if "--Switch_SpotFilter" is "ON".
                        DEFAULT: 10.
  --ImageParam_CropSize IMAGEPARAM_CROPSIZE
                        The pixel size of each patch subimage. For example, "--ImageParam_CropSize
                        100" means each patch subimage is 100*100 pixels. DEFAULT: 100.
  --ImageParam_PatchDist IMAGEPARAM_PATCHDIST
                        The pixel distance between adjacent patches. DEFAULT: 100.
  --ImageParam_TotalChannelNumber {3,4}
                        The total number of channels, which is needed only if "--Input_Image_Type"
                        is "IF". For example, "--ImageParam_TotalChannelNumber 3" means the TIFF
                        file contain three channels.
  --ImageParam_DAPIChannel {1,2,3,4}
                        The channel of the DAPI, which is needed only if "--Input_Image_Type" is
                        "IF". For example: "--ImageParam_DAPIChannel 1".
  --ImageParam_FiducialFrameChannel {1,2,3,4}
                        The channel of the fiducial frame, which is needed only if "--
                        Input_Image_Type" is "IF". For example: "--ImageParam_FiducialFrameChannel
                        3".
  --FeatureParam_ProcessNumber FEATUREPARAM_PROCESSNUMBER
                        The number of worker processes to create when extracting texture and
                        intensity features, which is used when "-m" is 1. DEFAULT: 8.
  --FeatureParam_FeatureType {0,1,2}
                        This determines which type of image features to use when "-m" is 1. 0 for
                        both texture and intensity features, 1 for texture features only and 2 for
                        intensity features only. DEFAULT: 0.
  --FeatureParam_ClipLimit FEATUREPARAM_CLIPLIMIT
                        The clipping limit, which is used when "-m" is 1. It is normalized between
                        0 and 1, with higher values representing more contrast. DEFAULT: 0.01.
  --FeatureParam_IterCount FEATUREPARAM_ITERCOUNT
                        Number of iterations Grabcut image segmentation algorithm should make
                        before returning the result. DEFAULT: 50.
  --CellProfilerParam_Pipeline CELLPROFILERPARAM_PIPELINE
                        The path to the CellProfiler pipline. It would be better to use different
                        piplines for H&E and IF samples. For H&E samples,
                        "Cellprofiler_Pipeline_HE.cppipe" is recommended. For IF samples,
                        "Cellprofiler_Pipeline_IF_C3/4.cppipe" is recommended based on the total
                        number of channels. In the docker image, the pipelines are stored in
                        "/root". You can also download them from
                        https://github.com/TongjiZhanglab/ImSpiRE to your own working directory.
                        DEFAULT: "Cellprofiler_Pipeline_HE.cppipe".
  --CellProfilerParam_KernelNumber CELLPROFILERPARAM_KERNELNUMBER
                        This option specifies the number of kernel to use to run CellProfiler.
                        DEFAULT: 8.
  --OptimalTransportParam_Alpha OPTIMALTRANSPORTPARAM_ALPHA
                        The weighting parameter between physical distance network and image
                        feature distance network, ranging from 0 to 1. For example, "--
                        OptimalTransportParam_Alpha 0.5" means ImSpiRE will equally consider the
                        weight of physical distance network and image feature distance network.
                        DEFAULT: 0.5.
  --OptimalTransportParam_Beta OPTIMALTRANSPORTPARAM_BETA
                        The constant interpolating parameter of Fused-gromov-Wasserstein transport
                        ranging from 0 to 1. DEFAULT: 0.5.
  --OptimalTransportParam_NumNeighbors OPTIMALTRANSPORTPARAM_NUMNEIGHBORS
                        The number of neighbors for nearest neighbors graph. DEFAULT: 5.
  --OptimalTransportParam_Epsilon OPTIMALTRANSPORTPARAM_EPSILON
                        The entropic regularization term with value greater than 0. DEFAULT:
                        0.001.

```

## Run ImSpiRE in Python interface

ImSpiRE can also be used step by step in the Python interface and easily integrated into custom scripts. [Here](https://nbviewer.org/github/TongjiZhanglab/ImSpiRE/blob/master/tutorial.ipynb) is a tutorial written using [Jupyter Notebook](https://jupyter.org/).

## Citation

