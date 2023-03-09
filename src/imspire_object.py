# -*- coding: utf-8 -*-
# @Author: huayuwei
# @Date:   2021-11-22 13:28:06
# @Last Modified by:   zhangyizhi
# @Last Modified time: 2023-03-09 17:26:41

###########
# imports #
###########

import numpy as np
import pandas as pd
import scanpy as sc
import tifffile as tif
import subprocess
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
from scipy import stats
from scipy.stats import pearsonr
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import os
import ot
import time
from utility import * 
import cv2
import scipy
import skimage
import skimage.morphology
import skimage.feature
import skimage.filters
from decimal import *
from skimage import io, img_as_float32, morphology, exposure
from skimage.feature import greycomatrix, greycoprops
from sklearn.preprocessing import minmax_scale
from skimage.color import separate_stains, hdx_from_rgb
from itertools import product
from multiprocessing import Pool
import warnings
warnings.filterwarnings("ignore")


class ImSpiRE_Parameters(object):
    """
    ImSpiRE_Parameters object is the container of parameters and default values.
    The values can be changed as needed.
    """
    def __init__(self):
        super(ImSpiRE_Parameters, self).__init__()
        #Basic Params#
        self.BasicParam_OutputDir="."
        self.BasicParam_OutputName="ImSpiRE"
        self.BasicParam_InputDir="./"
        self.BasicParam_InputCountFile="filtered_feature_bc_matrix.h5"  #OR count_matrix.txt
        self.BasicParam_InputImageFile="image.tif"
        self.BasicParam_InputImageType="H&E"    #["H&E","IF"]
        self.BasicParam_Verbose=False
        self.BasicParam_RandomState=0
        self.BasicParam_Overwriting=False
        self.BasicParam_PlatForm="Visium"   #["Visium","ST"]
        self.BasicParam_Mode=1  #[1,2]
        #Spot Preprocess#
        self.Switch_Preprocess=True
        self.Threshold_MinCounts=100
        self.Threshold_MaxCounts=10000
        self.Threshold_MitoPercent=20
        self.Threshold_MinSpot=10
        #Image segment#
        self.ImageParam_CropSize=100
        self.ImageParam_PatchDist=100
        self.ImageParam_TotalChannelNumber=None
        self.ImageParam_DAPIChannel=None
        self.ImageParam_FiducialFrameChannel=None
        ## Feature extraction#
        self.FeatureParam_ProcessNumber=8
        self.FeatureParam_FeatureType=0 #[0,1,2]
        self.FeatureParam_ClipLimit=0.01
        self.FeatureParam_IterCount=50
        #CellProfiler Params#
        self.CellProfilerParam_Pipeline="Cellprofiler_Pipeline_HE.cppipe"
        self.CellProfilerParam_KernelNumber=8
        #OT Solver Params#
        self.OptimalTransportParam_Alpha=0.5
        self.OptimalTransportParam_Beta=0.5
        self.OptimalTransportParam_NumNeighbors=5
        self.OptimalTransportParam_Epsilon=0.001
        self.OptimalTransportParam_NumIterMax=10


class ImSpiRE_Data(object):
    """
    ImSpiRE_Data object is the container to read 10x visium or ST files,
    and to do basic preprocessing if needed.

    Functions
    ---------
    read_10x_visium -- Read 10x-Genomics-formatted visum dataset.
    read_ST         -- Read ST dataset.
    preprocess      -- Basic filtering of spots and genes.
    """
    def __init__(self):
        super(ImSpiRE_Data, self).__init__()

    def read_10x_visium(self, path, count_file='filtered_feature_bc_matrix.h5'):
        """
        Read 10x-Genomics-formatted visum dataset.

        Parameters
        ----------
        path       -- Path to directory for visium datafiles.
        count_file -- Which file in the passed directory to use as the count file. 
                      Typically would be ‘filtered_feature_bc_matrix.h5’.

        Returns
        -------
        adata -- Annotated data matrix.
        """
        adata = sc.read_visium(path, count_file=count_file, load_images=False)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        tissue_positions_list_path=f"{path}/spatial/tissue_positions_list.csv"
        pos=pd.read_csv(tissue_positions_list_path,header=None)
        pos.columns=["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]
        pos.index=pos["barcode"].tolist()
        pos=pos.loc[adata.obs.index,]

        self.pos=pos 
        self.pos_not_in_tissue=pos.loc[pos["in_tissue"] == 0,:] 
        
        self.pos_in_tissue=pos.loc[pos["in_tissue"] == 1,:]
        self.pos_in_tissue_filter=pos.loc[pos["in_tissue"] == 1,:]

        
        adata.obs["pxl_row_in_fullres"]=pos["pxl_row_in_fullres"]
        adata.obs["pxl_col_in_fullres"]=pos["pxl_col_in_fullres"]

        adata.obs["array_row"]=pos["array_row"]
        adata.obs["array_col"]=pos["array_col"]

        adata.obs["barcode"]=pos["barcode"]

        self.adata = adata

    def read_ST(self, path, count_file):
        """
        Read ST dataset.

        Parameters
        ----------
        path       -- Path to directory for ST datafiles.
        count_file -- Which file in the passed directory to use as the count file. 
                      Note that it should be tab-delimited.

        Returns
        -------
        adata -- Annotated data matrix.
        """
        tissue_positions_list_path=f"{path}/pxl_pos.txt"
        pos=pd.read_csv(tissue_positions_list_path,sep="\t")
        pos.columns=["pxl_row_in_fullres","pxl_col_in_fullres"]
        
        exp_data=pd.read_csv(count_file,sep="\t",index_col=0)
        exp_data.index=[f"spot{i}" for i in range(exp_data.shape[0])]
        pos.index=exp_data.index.to_list()

        adata=sc.AnnData(exp_data)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        
        # pos.index=exp_data.index.to_list()
        self.pos_in_tissue=pos
        self.pos_in_tissue_filter=pos

        # pos=pos.loc[adata.obs.index,]
        
        adata.obs["pxl_row_in_fullres"]=pos["pxl_row_in_fullres"]
        adata.obs["pxl_col_in_fullres"]=pos["pxl_col_in_fullres"]

        self.adata = adata

    def preprocess(self, min_counts, max_counts, pct_counts_mt, min_cells):
        """
        Basic filtering of spots and genes.
        
        Parameters
        ----------
        min_counts      -- Minimum number of counts required for a spot to pass filtering.
        max_counts      -- Maximum number of counts required for a spot to pass filtering.
        pct_counts_mt   -- Maximum count percent of mitochondrial genes required for a spot to pass filtering.
        min_cells       -- Minimum number of cells expressed required for a gene to pass filtering.

        Returns
        -------
        pos_in_tissue_filter    -- Location matrices of filtered spots in tissue. 
        """
        sc.pp.filter_cells(self.adata, min_counts=min_counts)
        sc.pp.filter_cells(self.adata, max_counts=max_counts)
        self.adata = self.adata[self.adata.obs["pct_counts_mt"] < pct_counts_mt]
        sc.pp.filter_genes(self.adata, min_cells=min_cells)

        self.pos_in_tissue_filter=self.pos_in_tissue.loc[self.adata.obs.index,]


class ImSpiRE_HE_Image(object):
    """
    ImSpiRE_HE_Image object is the container for H&E image processing.

    Parameters
    ----------
    image_file_path -- Path to the high-resolution tissue image.
    platform        -- Visium or ST.
    output_dir      -- The project directory.
    output_name     -- The project name.
    iter_count      -- The Number of iterations Grabcut image segmentation algorithm should make 
                       before returning the result. Default is 50.

    Functions
    ---------
    segment_spot_image              -- Segment the spot subimages.
    generate_patch_locations        -- Determine the patch coordinates.
    generate_patch_locations_2      -- Another way to determine the patch coordinates.
    segment_patch_image             -- Segment the patch subimages.
    """
    def __init__(self, image_file_path, platform = "Visium", output_dir = "./", output_name = "ImSpiRE", iter_count = 50):
        super(ImSpiRE_HE_Image, self).__init__()
        self.image_file_path = image_file_path
        self.platform = platform
        self.output_dir = output_dir
        self.output_name = output_name
        self.iter_count = iter_count

    def segment_spot_image(self,
        pos_in_tissue_filter,
        output_path:str,
        crop_size:int = 100):
        """
        Segment the spot subimages. The image size is adjustable 
        by setting the parameter CROP_SIZE. Default pixel size is 100*100.

        Parameters
        ----------
        pos_in_tissue_filter    -- Location matrices of filtered spots in tissue.
        output_path             -- Path to save spot subimages.
        crop_size               -- Pixel size of each spot subimage. Default is 100.

        Outputs
        -------
        The output spot segmentation images will be saved in {output_path}.
        """
        im = read_image(self.image_file_path)
        # Check the exist of spot_image_output_path
        # spot_image_output_path = f"{output_path}/spot_image"
        spot_image_output_path = f"{output_path}"
        # mkdir(output_path,"spot_image")

        # if verbose:
        #     print("Generating spot images...")

        # Spot image segmentation
        for index, row in pos_in_tissue_filter.iterrows():
            # barcode = row["barcode"]
            # in_tissue = row["in_tissue"]
            # array_row = row["array_row"]
            # array_col = row["array_col"]
            pxl_row = row["pxl_row_in_fullres"]
            pxl_col = row["pxl_col_in_fullres"]

            spot_row_down = pxl_row - crop_size/2
            spot_row_up = pxl_row + crop_size/2
            spot_col_left = pxl_col - crop_size/2
            spot_col_right = pxl_col + crop_size/2
            
            spot_im = im[int(spot_row_down):int(spot_row_up+1), int(spot_col_left):int(spot_col_right+1),:]

            spot_path = f"{spot_image_output_path}/{index}_SpotImage_HE_{pxl_row}_{pxl_col}.tif"  
            tif.imwrite(spot_path, spot_im)

    def generate_patch_locations(self, 
        pos, 
        pos_in_tissue,
        dist:int = 100):
        """
        Determine the patch coordinates. The resolution-enhanced spots is defined as patches.
        The enhanced resolution can be specified by setting the parameter DIST.

        Parameters
        ----------
        pos             -- Location matrices of spots.
        pos_in_tissue   -- Location matrices of spots in tissue. 
        dist            -- This parameter specifies the pixel distance between patches.
        
        Returns
        -------
        patch_in_tissue -- A text file that contains patch coordinates.
        """ 
        im = read_image(self.image_file_path)

        scale, scale_row, scale_col = compressed_scale_HE(2000, im.shape)

        cimg = cv2.resize(im, (2000,2000), interpolation = cv2.INTER_AREA)
        params = [cv2.IMWRITE_JPEG_QUALITY, 80]
        cv2.imwrite(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cimg, params)
        
        cimg = cv2.imread(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cv2.IMREAD_UNCHANGED)
        if len(cimg.shape) == 3:
            cimg = cv2.cvtColor(cimg, cv2.COLOR_RGB2GRAY)
        cv2.normalize(cimg, cimg, 255, 0, cv2.NORM_MINMAX)
        if cimg.dtype == "uint16":
            cimg = cimg.astype("uint8")
        # remove_cimg = os.system("rm ./compressed_image.png")
        patch_locations = filter_patch_in_tissue(cimg, scale, scale_row, scale_col, pos, pos_in_tissue, dist)
        
        self.patch_in_tissue = patch_locations.loc[patch_locations["in_tissue"] == 1,:]

    def generate_patch_locations_2(self, 
        pos_in_tissue,
        dist:int = 100):

        im = read_image(self.image_file_path)

        scale, scale_row, scale_col = compressed_scale_HE(500, im.shape)

        cimg = cv2.resize(im, (500,500), interpolation = cv2.INTER_AREA)
        params = [cv2.IMWRITE_JPEG_QUALITY, 80]
        cv2.imwrite(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cimg, params)
        
        patch_locations = filter_patch_in_tissue_2(scale, scale_row, scale_col, pos_in_tissue, dist, self.platform, self.output_dir, self.output_name, self.iter_count)

        self.patch_in_tissue = patch_locations.loc[patch_locations["in_tissue"] == 1,:]

    def segment_patch_image(self,
        output_path,
        patch_in_tissue,
        crop_size:int = 100):
        
        im = read_image(self.image_file_path)

        # patch_image_output_path = f"{output_path}/patch_image"
        patch_image_output_path = f"{output_path}"
        # mkdir(output_path,"patch_image")

        # if verbose:
        #     print("Generating patch images...")

        # Patch image segmentation  
        for index, row in patch_in_tissue.iterrows():
            # patch_array_row = row["row"]
            # patch_array_col = row["col"]
            patch_pxl_row = row["pxl_row"]
            patch_pxl_col = row["pxl_col"]
            in_tissue = row["in_tissue"]

            patch_row_down = patch_pxl_row - crop_size/2
            patch_row_up = patch_pxl_row + crop_size/2
            patch_col_left = patch_pxl_col - crop_size/2
            patch_col_right = patch_pxl_col + crop_size/2

            patch_im = im[int(patch_row_down):int(patch_row_up+1), int(patch_col_left):int(patch_col_right+1),:]
            patch_path = f"{patch_image_output_path}/{index}_PatchImage_IF_{patch_pxl_row}_{patch_pxl_col}.tif"  
            tif.imwrite(patch_path, patch_im)


class ImSpiRE_IF_Image(object):
    """
    ImSpiRE_IF_Image object is the container for IF image processing.

    Parameters
    ----------
    image_file_path         -- Path to the high-resolution tissue image.
    DAPI_channel			-- The channel of the DAPI. Default is None.
    fiducial_frame_channel  -- The channel of the fiducial frame. Default is None.
    platform                -- Visium or ST.
    output_dir              -- The project directory.
    output_name             -- The project name.
    iter_count              -- The Number of iterations Grabcut image segmentation algorithm should make 
                               before returning the result. Default is 50.

    Functions
    ---------
    segment_spot_image              -- Segment the spot subimages.
    generate_patch_locations        -- Determine the patch coordinates.
    generate_patch_locations_2      -- Another way to determine the patch coordinates.
    segment_patch_image             -- Segment the patch subimages.
    """
    def __init__(self, image_file_path, DAPI_channel = None, fiducial_frame_channel = None, platform = "Visium", output_dir = "./", output_name = "ImSpiRE", iter_count = 50):
        super(ImSpiRE_IF_Image, self).__init__()
        self.image_file_path = image_file_path
        self.DAPI_channel = DAPI_channel
        self.fiducial_frame_channel = fiducial_frame_channel
        self.platform = platform
        self.output_dir = output_dir
        self.output_name = output_name
        self.iter_count = iter_count

    def segment_spot_image(self,
        pos_in_tissue_filter,
        output_path:str,
        crop_size:int = 100):

        im = read_image(self.image_file_path)

        if self.fiducial_frame_channel != None:
            im = modify_IF_image(im, self.DAPI_channel, self.fiducial_frame_channel)
        # Check the exist of spot_image_output_path
        # spot_image_output_path = f"{output_path}/spot_image"
        spot_image_output_path = f"{output_path}"
        # mkdir(output_path, "spot_image")

        # if verbose:
        #     print("Generating spot images...")

        # Spot image segmentation
        for index, row in pos_in_tissue_filter.iterrows():
            # barcode = row["barcode"]
            # in_tissue = row["in_tissue"]
            # array_row = row["array_row"]
            # array_col = row["array_col"]
            pxl_row = row["pxl_row_in_fullres"]
            pxl_col = row["pxl_col_in_fullres"]

            spot_row_down = pxl_row - crop_size/2
            spot_row_up = pxl_row + crop_size/2
            spot_col_left = pxl_col - crop_size/2
            spot_col_right = pxl_col + crop_size/2
                
            spot_im = im[:, int(spot_row_down):int(spot_row_up+1), int(spot_col_left):int(spot_col_right+1)]
            
            spot_path = f"{spot_image_output_path}/{index}_SpotImage_IF_{pxl_row}_{pxl_col}.tif"  
            tif.imwrite(spot_path, spot_im)

    def generate_patch_locations(self, 
        pos, 
        pos_in_tissue,
        dist:int = 100):

        im = read_image(self.image_file_path)
        
        scale, scale_row, scale_col = compressed_scale_IF(2000, im.shape)

        img = im.transpose(1,2,0)
        cimg = cv2.resize(img, (2000,2000), interpolation = cv2.INTER_AREA)
        params = [cv2.IMWRITE_JPEG_QUALITY, 80]
        cv2.imwrite(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cimg, params)
        
        cimg = cv2.imread(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cv2.IMREAD_UNCHANGED)
        if len(cimg.shape) == 3:
            cimg = cv2.cvtColor(cimg, cv2.COLOR_RGB2GRAY)
        cv2.normalize(cimg, cimg, 255, 0, cv2.NORM_MINMAX)
        if cimg.dtype == "uint16":
            cimg = cimg.astype("uint8")
        # remove_cimg = os.system("rm ./compressed_image.png")
        patch_locations = filter_patch_in_tissue(cimg, scale, scale_row, scale_col, pos, pos_in_tissue, dist)
        
        self.patch_in_tissue = patch_locations.loc[patch_locations["in_tissue"] == 1,:]

    def generate_patch_locations_2(self, 
        pos_in_tissue,
        dist:int = 100):

        im = read_image(self.image_file_path)
        
        scale, scale_row, scale_col = compressed_scale_IF(500, im.shape)

        img = im.transpose(1,2,0)
        cimg = cv2.resize(img, (500,500), interpolation = cv2.INTER_AREA)
        params = [cv2.IMWRITE_JPEG_QUALITY, 80]
        cv2.imwrite(f"{self.output_dir}/{self.output_name}/SupplementaryResults/compressed_image.png", cimg, params)
        
        patch_locations = filter_patch_in_tissue_2(scale, scale_row, scale_col, pos_in_tissue, dist, self.platform, self.output_dir, self.output_name, self.fiducial_frame_channel, self.iter_count)

        self.patch_in_tissue = patch_locations.loc[patch_locations["in_tissue"] == 1,:]

    def segment_patch_image(self,
        output_path,
        patch_in_tissue,
        crop_size:int = 100):

        im = read_image(self.image_file_path)
        
        if self.fiducial_frame_channel != None:
            im = modify_IF_image(im, self.DAPI_channel, self.fiducial_frame_channel)
        
        # patch_image_output_path = f"{output_path}/patch_image"
        patch_image_output_path = f"{output_path}"
        # mkdir(output_path, "patch_image")

        # if verbose:
        #     print("Generating patch images...")

        # Patch image segmentation  
        for index, row in patch_in_tissue.iterrows():
            # patch_array_row = row["row"]
            # patch_array_col = row["col"]
            patch_pxl_row = row["pxl_row"]
            patch_pxl_col = row["pxl_col"]
            # patch_filter=row["in_tissue"]

            patch_row_down = patch_pxl_row - crop_size/2
            patch_row_up = patch_pxl_row + crop_size/2
            patch_col_left = patch_pxl_col - crop_size/2
            patch_col_right = patch_pxl_col + crop_size/2

            patch_im = im[:, int(patch_row_down):int(patch_row_up+1), int(patch_col_left):int(patch_col_right+1)]

            patch_path = f"{patch_image_output_path}/{index}_PatchImage_IF_{patch_pxl_row}_{patch_pxl_col}.tif"  
            tif.imwrite(patch_path, patch_im)

            
class ImSpiRE_Image_Feature(object):
    """
    ImSpiRE_Image_Feature object is the container to extract texture and intensity features of subimages.

    Parameters
    ----------
    image_file_path -- Path to the high-resolution tissue image.
    clip_limit      -- Clipping limit, normalized between 0 and 1 (higher values give more contrast).
    n_processes     -- The number of worker processes to create.

    Functions
    ---------
    image_preprocess        -- Preprocess the high-resolution tissue image.
    run_extract_features    -- Extract texture and intensity features of subimages.

    """
    def __init__(self, image_file_path, clip_limit = 0.01, n_processes = 8):
        super(ImSpiRE_Image_Feature, self).__init__()
        self.image_file_path = image_file_path
        self.clip_limit = clip_limit
        self.n_processes = n_processes

    def image_preprocess(self):
        img = io.imread(self.image_file_path)

        img = separate_stains(img, hdx_from_rgb)
        img = minmax_scale(img.reshape(-1, 3)).reshape(img.shape)
        img = np.clip(img, 0, 1)
        img = exposure.equalize_adapthist(img, clip_limit=self.clip_limit)
        img = (255 * img).astype("uint8")

        self.processed_image=img
    
    def run_extract_features(self,
        processed_image,
        feature_set,
        img_meta,
        crop_size:int = 100):
        """
        Extract texture and intensity features of subimages.

        Parameters
        ----------
        feature_set     -- Features to be extracted.
        img_meta        -- Pixel coordinates of spots or patches in tissue. 
        crop_size       -- Pixel size of subimages. Default is 100.
        
        Returns
        -------
        texture_features    -- Texture features of subimages.
        intensity_features  -- Intensity features of subimages.
        """
        img=processed_image
        img_meta.columns = ["row","col"]

        # multiprocessing
        pool = Pool(self.n_processes)
        result = []

        for i in range(img_meta.shape[0]):
            ind = img_meta.index[i]
            row = img_meta.iloc[i]
            row, col= row[["row", "col"]].astype(int)
            spot_img = img[int(row - crop_size/2) : int(row + crop_size/2 + 1), int(col - crop_size/2) : int(col + crop_size/2 + 1),:]
            
            result.append(pool.apply_async(func = extract_image_features, args = (ind,spot_img,img.shape[2],feature_set)))

        pool.close()
        pool.join()

        # get values
        texture_features = pd.DataFrame()
        intensity_features = pd.DataFrame()    
        
        for i in range(img_meta.shape[0]):
            tmp = result[i].get()
            
            tmp_index = tmp['index']
            tmp_texture = tmp['texture_feature']
            tmp_intensity = tmp['intensity_feature']
            
            tmp_texture_df = pd.DataFrame(tmp_texture).T
            tmp_texture_df.index = [tmp_index]
            tmp_intensity_df = pd.DataFrame(tmp_intensity).T
            tmp_intensity_df.index = [tmp_index]
            
            texture_features = pd.concat([texture_features,tmp_texture_df])
            intensity_features = pd.concat([intensity_features,tmp_intensity_df])
        
        texture_features = texture_features.sort_index()
        intensity_features = intensity_features.sort_index()

        # Naming the features. f stands for channels, A stands for angles.
        # construct texture feature table
        int_low = 0.2
        int_high = 0.8
        int_step = 0.1
        q_bins = np.arange(int_low, int_high, int_step)

        channels = ["f" + str(i) for i in range(img.shape[2])]
        col_names = product(channels, feature_set, ["A1", "A2", "A3", "A4"])
        col_names = ["_".join(x) for x in col_names]

        texture_features.columns = col_names
        intensity_features.columns = ["_".join(x) for x in product(channels, ["{:.1f}".format(x) for x in q_bins])]

        
        self.texture_features = texture_features
        self.intensity_features = intensity_features
                

class ImSpiRE_HE_CellProfiler(object):
    """
    ImSpiRE_CellProfiler_HE object is the container for running CellProfiler to extract H&E image features.

    Parameters
    ----------
    image_file_path         -- Path to the high-resolution tissue images.
    cellprofiler_pipeline   -- Path to the CellProfiler pipeline.

    Functions
    ---------
    compute_image_features  -- Use CellProfiler to extract image features.
    filter_image_features   -- Filter out useless features from all image features.
    """
    def __init__(self, image_file_path, cellprofiler_pipeline):
        super(ImSpiRE_HE_CellProfiler, self).__init__()
        self.image_file_path = image_file_path
        self.cellprofiler_pipeline = cellprofiler_pipeline

    def compute_image_features(self, output_path, number_of_kernels):
        """
        Use CellProfiler to extract image features.

        Parameters
        ----------
        output_path         -- Path to save the image features file.
        number_of_kernels   -- This parameter specifies the number of kernels to use to run CellProfiler.

        Outputs
        -------
        A text file that contains image features will be saved in {output_path}.
        """
        n_kernel = number_of_kernels
        n_image = len(os.listdir(self.image_file_path))
        image_per_batch = int(Decimal(n_image/n_kernel).quantize(Decimal('0'), rounding = ROUND_HALF_UP))

        create_batch_data=os.system(f"cellprofiler -r -c -p {self.cellprofiler_pipeline} -o {output_path} -i {self.image_file_path}")
        initialize_batch_cmd=os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f 1 -l {image_per_batch} -o {output_path}/temp1 &" > {output_path}/batch_cmd.sh')
        n = n_kernel if (n_image/n_kernel)-int(n_image/n_kernel) >= 0.5 else n_kernel-1
        for i in range(1,n):
            n_first = int(i*image_per_batch+1)
            n_last = int((i+1)*image_per_batch) 
            n_last = n_last if n_last <= n_image else n_image
            create_batch_cmd = os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f {n_first} -l {n_last} -o {output_path}/temp{i+1} &" >> {output_path}/batch_cmd.sh')
        if n == n_kernel-1:
            finish_batch_cmd = os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f {(n_kernel-1)*image_per_batch+1} -l {n_image} -o {output_path}/temp{n_kernel} &" >> {output_path}/batch_cmd.sh')
        create_parallel = os.system(f'echo "wait;" >> {output_path}/batch_cmd.sh')
        run_batch = os.system(f"bash {output_path}/batch_cmd.sh")
        merge_files = os.system(f'cat {output_path}/temp*/Image_Features_Image.txt > {output_path}/Image_Features_Image.txt')
        remove_temp_files = os.system(f"rm -rf {output_path}/temp*")

    def filter_image_features(self,output_path):
        """
        Filter out useless features from all image features.

        Parameters
        ----------
        output_path     -- Path to the image features file.

        Returns
        -------
        image_features  -- Filtered features of subimages extracted by CellProfiler.
        
        Outputs
        -------
        A text file that contains filtered image features will be saved in {output_path}.
        """
        image_features = pd.read_csv(f"{output_path}/Image_Features_Image.txt", sep = "\t")
        image_features = image_features.drop_duplicates(keep = False)
        image_features_index = image_features["FileName_Image"].str.split("_", expand = True)[0].str.split(".", expand = True)[0].tolist()
        image_features.index = [i for i in image_features_index]
        image_features = image_features.sort_index()
        image_features.drop(columns = ["Channel_Image",
                                    "ExecutionTime_01Images",
                                    "ExecutionTime_02Metadata",
                                    "ExecutionTime_03NamesAndTypes",
                                    "ExecutionTime_04Groups",
                                    "ExecutionTime_06MeasureImageQuality",
                                    "ExecutionTime_05MeasureImageIntensity",
                                    "ExecutionTime_07MeasureTexture",
                                    "FileName_Image",
                                    "Frame_Image",
                                    "Group_Index",
                                    "Group_Number",
                                    "Height_Image",
                                    "ImageNumber",
                                    "ImageSet_ImageSet",
                                    "MD5Digest_Image",
                                    "Metadata_FileLocation",
                                    "Metadata_Frame",
                                    "Metadata_Series",
                                    "ModuleError_01Images",
                                    "ModuleError_02Metadata",
                                    "ModuleError_03NamesAndTypes",
                                    "ModuleError_04Groups",
                                    "ModuleError_06MeasureImageQuality",
                                    "ModuleError_05MeasureImageIntensity",
                                    "ModuleError_07MeasureTexture",
                                    "PathName_Image",
                                    "Scaling_Image",
                                    "Series_Image",
                                    "URL_Image",
                                    "Width_Image"], inplace = True)
        
        for col in image_features.columns:
            if len(set(image_features[col])) == 1:
                image_features.drop(columns = [col], inplace = True)

        self.image_features=image_features
        


class ImSpiRE_IF_CellProfiler(object):
    def __init__(self, image_file_path, cellprofiler_pipeline):
        super(ImSpiRE_IF_CellProfiler, self).__init__()
        self.image_file_path = image_file_path
        self.cellprofiler_pipeline = cellprofiler_pipeline

    def compute_image_features(self, output_path, number_of_kernels):
        n_kernel = number_of_kernels
        n_image = len(os.listdir(self.image_file_path))
        image_per_batch = int(Decimal(n_image/n_kernel).quantize(Decimal('0'), rounding = ROUND_HALF_UP))

        create_batch_data=os.system(f"cellprofiler -r -c -p {self.cellprofiler_pipeline} -o {output_path} -i {self.image_file_path}")
        initialize_batch_cmd=os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f 1 -l {image_per_batch} -o {output_path}/temp1 &" > {output_path}/batch_cmd.sh')
        n = n_kernel if (n_image/n_kernel)-int(n_image/n_kernel) >= 0.5 else n_kernel-1
        for i in range(1,n):
            n_first = int(i*image_per_batch+1)
            n_last = int((i+1)*image_per_batch) 
            n_last = n_last if n_last <= n_image else n_image
            create_batch_cmd = os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f {n_first} -l {n_last} -o {output_path}/temp{i+1} &" >> {output_path}/batch_cmd.sh')
        if n == n_kernel-1:
            finish_batch_cmd = os.system(f'echo "cellprofiler -c -r -p {output_path}/Batch_data.h5 -f {(n_kernel-1)*image_per_batch+1} -l {n_image} -o {output_path}/temp{n_kernel} &" >> {output_path}/batch_cmd.sh')
        create_parallel = os.system(f'echo "wait;" >> {output_path}/batch_cmd.sh')
        run_batch = os.system(f"bash {output_path}/batch_cmd.sh")
        merge_files = os.system(f'cat {output_path}/temp*/Image.txt > {output_path}/Image_Features_Image.txt')
        remove_temp_files = os.system(f"rm -rf {output_path}/temp*")

    def filter_image_features(self,output_path):
        image_features = pd.read_csv(f"{output_path}/Image_Features_Image.txt", sep = "\t")
        image_features = image_features.drop_duplicates(keep = False)
        image_features_index = image_features["FileName_DNA"].str.split("_", expand = True)[0].str.split(".", expand = True)[0].tolist()
        image_features.index = [i for i in image_features_index]
        image_features = image_features.sort_index()
        image_features.drop(columns = ["Channel_DNA",
                                    "ExecutionTime_01Images",
                                    "ExecutionTime_02Metadata",
                                    "ExecutionTime_03NamesAndTypes",
                                    "ExecutionTime_04Groups",
                                    "ExecutionTime_05ColorToGray",
                                    "ExecutionTime_06MeasureImageQuality",
                                    "ExecutionTime_07MeasureImageIntensity",
                                    "ExecutionTime_08IdentifyPrimaryObjects",
                                    "ExecutionTime_09MeasureImageAreaOccupied",
                                    "ExecutionTime_10MeasureTexture",
                                    "FileName_DNA",
                                    "Frame_DNA",
                                    "Group_Index",
                                    "Group_Number",
                                    "Height_DNA",
                                    "ImageNumber",
                                    "ImageSet_ImageSet",
                                    "MD5Digest_DNA",
                                    "Metadata_FileLocation",
                                    "Metadata_Frame",
                                    "Metadata_Series",
                                    "ModuleError_01Images",
                                    "ModuleError_02Metadata",
                                    "ModuleError_03NamesAndTypes",
                                    "ModuleError_04Groups",
                                    "ModuleError_05ColorToGray",
                                    "ModuleError_06MeasureImageQuality",
                                    "ModuleError_07MeasureImageIntensity",
                                    "ModuleError_08IdentifyPrimaryObjects",
                                    "ModuleError_09MeasureImageAreaOccupied",
                                    "ModuleError_10MeasureTexture",
                                    "PathName_DNA",
                                    "Scaling_DNA",
                                    "Series_DNA","URL_DNA","Width_DNA"], inplace = True)
        
        for col in image_features.columns:
            if len(set(image_features[col])) == 1:
                image_features.drop(columns = [col], inplace = True)

        self.image_features=image_features

   
class ImSpiRE_OT_Solver(object):
    """
    ImSpiRE_OT_Solver object is the container for solving OT problem.

    Parameters
    ----------
    spot_locations       -- Locations matrices of spots.
    patch_locations      -- Locations matrices of patches.
    spot_image_features  -- Image features matrices of spots.
    patch_image_features -- Image features matrices of patches.
    spot_gene_expression -- Gene expression features matrices of spots.

    Funtions
    --------
    setup_cost_matrices -- Compute cost matrices for OT.
    solve_OT            -- Solve OT problem and compute coupling matrices.
    """
    def __init__(self, spot_locations, patch_locations,
                       spot_image_features, patch_image_features,
                       spot_gene_expression,random_state=0):
        super(ImSpiRE_OT_Solver, self).__init__()
        self.spot_locations = spot_locations
        self.patch_locations = patch_locations
        self.spot_image_features = spot_image_features
        self.patch_image_features = patch_image_features
        self.spot_gene_expression = spot_gene_expression
        ot.utils.check_random_state(random_state)

    def setup_cost_matrices(self,alpha=0.5, num_neighbors=5):
        """
        Compute cost matrices for OT.

        Parameters
        ----------
        alpha         -- A constant interpolating between image features cost matrices and locations cost matirces.
        num_neighbors -- Number of neighbors for nearest neighbors graph.

        """
        self.alpha=alpha
        self.num_neighbors=num_neighbors
        locations_cost_matrices=compute_distance_matrices(self.spot_locations,self.patch_locations,'euclidean')
        image_features_cost_matrices=compute_distance_matrices(self.spot_image_features,self.patch_image_features,'correlation')
        M=(1-alpha)*locations_cost_matrices+alpha*image_features_cost_matrices
        M /= M.max()
        self.spot_gene_expression_pca=compute_pca(self.spot_gene_expression)
        C1=compute_graph_distance_matrices(self.spot_gene_expression_pca, num_neighbors=num_neighbors, metric='correlation')

        patch_locations_graph_distance_matrices=compute_graph_distance_matrices(self.patch_locations, num_neighbors=num_neighbors, metric='euclidean')
        patch_image_features_graph_distance_matrices=compute_graph_distance_matrices(self.patch_image_features, num_neighbors=num_neighbors, metric='correlation')
        C2=(1-alpha)*patch_locations_graph_distance_matrices+alpha*patch_image_features_graph_distance_matrices
        C2/=C2.max()

        self.M=M
        self.C1=C1
        self.C2=C2

    def solve_OT(self, beta=0.5, epsilon=0.001, **kwargs):
        """
        Solve OT problem and compute coupling matrices.

        Parameters
        ----------
        beta -- A trade-off parameter of Fused-gromov-Wasserstein transport.
        epsilon -- Entropic regularization term.
        """
        self.beta=beta
        self.epsilon=epsilon
        M=self.M
        C1=self.C1
        C2=self.C2
        source_p=ot.unif(M.shape[0])
        target_p=ot.unif(M.shape[1])
        T=entropic_fused_gromov_wasserstein(M, C1, C2, source_p, target_p, epsilon=epsilon, beta=beta, **kwargs)
        self.T=T


class ImSpiRE(object):
    """
    ImSpiRE object is the container to run ImSpiRE.

    Parameters
    ----------
    imspire_param   -- Parameters and default values used in ImSpiRE.

    Functions
    ---------
    run     -- Run ImSpiRE.
    """
    def __init__(self, imspire_param):
        super(ImSpiRE, self).__init__()
        self.imspire_param = imspire_param
    
    def run(self):
        """
        Run ImSpiRE.
        """
        #Create output folders#
        create_folder(self.imspire_param.BasicParam_OutputDir,
                      self.imspire_param.BasicParam_OutputName,
                      self.imspire_param.BasicParam_Overwriting)

        imdata=ImSpiRE_Data()
        if self.imspire_param.BasicParam_PlatForm=="Visium":
            imdata.read_10x_visium(self.imspire_param.BasicParam_InputDir, count_file=self.imspire_param.BasicParam_InputCountFile)
        else:
            imdata.read_ST(self.imspire_param.BasicParam_InputDir, count_file=self.imspire_param.BasicParam_InputCountFile)
        
        if self.imspire_param.Switch_Preprocess:
            imdata.preprocess(min_counts=self.imspire_param.Threshold_MinCounts, 
                              max_counts=self.imspire_param.Threshold_MaxCounts, 
                              pct_counts_mt=self.imspire_param.Threshold_MitoPercent, 
                              min_cells=self.imspire_param.Threshold_MinSpot)


        #ImageFeatureExtraction#
        if self.imspire_param.BasicParam_InputImageType=="H&E":
            im=ImSpiRE_HE_Image(self.imspire_param.BasicParam_InputImageFile,
                                self.imspire_param.BasicParam_PlatForm,
                                self.imspire_param.BasicParam_OutputDir,
                                self.imspire_param.BasicParam_OutputName,
                                self.imspire_param.FeatureParam_IterCount)
        else:
            im=ImSpiRE_IF_Image(self.imspire_param.BasicParam_InputImageFile,
                                self.imspire_param.ImageParam_DAPIChannel,
                                self.imspire_param.ImageParam_FiducialFrameChannel,
                                self.imspire_param.BasicParam_PlatForm,
                                self.imspire_param.BasicParam_OutputDir,
                                self.imspire_param.BasicParam_OutputName,
                                self.imspire_param.FeatureParam_IterCount)
        
        if self.imspire_param.BasicParam_Mode==1:
            im.generate_patch_locations_2(pos_in_tissue=imdata.pos_in_tissue_filter,
                                          dist=self.imspire_param.ImageParam_PatchDist)

            feature_set=["contrast","dissimilarity","homogeneity","ASM","energy","correlation"]

            if self.imspire_param.BasicParam_Verbose:
                print("Extracting features of spot images...")

            spot_ife = ImSpiRE_Image_Feature(self.imspire_param.BasicParam_InputImageFile,
                                             self.imspire_param.FeatureParam_ClipLimit,
                                             self.imspire_param.FeatureParam_ProcessNumber)
            spot_ife.image_preprocess()
            spot_ife.run_extract_features(processed_image=spot_ife.processed_image,
                                          feature_set=feature_set,
                                          img_meta=imdata.pos_in_tissue_filter[['pxl_row_in_fullres','pxl_col_in_fullres']],
                                          crop_size=self.imspire_param.ImageParam_CropSize)
            
            if self.imspire_param.BasicParam_Verbose:
                print("Extracting features of patch images...")

            patch_ife = ImSpiRE_Image_Feature(self.imspire_param.BasicParam_InputImageFile,
                                              self.imspire_param.FeatureParam_ClipLimit,
                                              self.imspire_param.FeatureParam_ProcessNumber)
            patch_ife.run_extract_features(processed_image=spot_ife.processed_image,
                                           feature_set=feature_set,
                                           img_meta=im.patch_in_tissue[['pxl_row','pxl_col']],
                                           crop_size=self.imspire_param.ImageParam_CropSize)
            
            spot_texture_features=spot_ife.texture_features.loc[imdata.pos_in_tissue_filter.index,]
            spot_intensity_features=spot_ife.intensity_features.loc[imdata.pos_in_tissue_filter.index,]

            patch_ife.texture_features.index=patch_ife.texture_features.index.astype('int')
            patch_ife.intensity_features.index=patch_ife.intensity_features.index.astype('int')
            patch_texture_features=patch_ife.texture_features.sort_index()
            patch_intensity_features=patch_ife.intensity_features.sort_index()
            patch_texture_features.index=list(range(patch_texture_features.shape[0]))
            patch_intensity_features.index=list(range(patch_intensity_features.shape[0]))

            spot_features=pd.concat([spot_texture_features,spot_intensity_features],axis=1)
            patch_features=pd.concat([patch_texture_features,patch_intensity_features],axis=1)

            spot_feature_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/FeatureResults/SpotFeature"
            patch_feature_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/FeatureResults/PatchFeature"
            
            spot_texture_features.to_csv(f"{spot_feature_output_path}/spot_texture_features.txt", sep = "\t")
            spot_intensity_features.to_csv(f"{spot_feature_output_path}/spot_intensity_features.txt", sep = "\t")
            spot_features.to_csv(f"{spot_feature_output_path}/spot_features.txt", sep = "\t")
    
            patch_texture_features.to_csv(f"{patch_feature_output_path}/patch_texture_features.txt", sep = "\t")
            patch_intensity_features.to_csv(f"{patch_feature_output_path}/patch_intensity_features.txt", sep = "\t")
            patch_features.to_csv(f"{patch_feature_output_path}/patch_features.txt", sep = "\t")
        else:
            if self.imspire_param.BasicParam_Verbose:
                print("Generating spot images...")

            spot_image_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/ImageResults/SpotImage"
            
            if self.imspire_param.BasicParam_Overwriting:
                if len(os.listdir(spot_image_output_path))!=0:
                    os.system(f"rm -rf {spot_image_output_path}/*")

            im.segment_spot_image(pos_in_tissue_filter=imdata.pos_in_tissue_filter,
                                  output_path=spot_image_output_path,
                                  crop_size=self.imspire_param.ImageParam_CropSize)
            
            if self.imspire_param.BasicParam_Verbose:
                print("Generating patch images...")

            patch_image_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/ImageResults/PatchImage"
            
            if self.imspire_param.BasicParam_Overwriting:
                if len(os.listdir(patch_image_output_path))!=0:
                    os.system(f"rm -rf {patch_image_output_path}/*")

            if self.imspire_param.BasicParam_PlatForm=="ST":
                im.generate_patch_locations_2(pos_in_tissue=imdata.pos_in_tissue_filter,
                                                 dist=self.imspire_param.ImageParam_PatchDist)
            else:
                im.generate_patch_locations(pos=imdata.pos, 
                                            pos_in_tissue=imdata.pos_in_tissue_filter,
                                            dist=self.imspire_param.ImageParam_PatchDist)
            im.segment_patch_image(patch_in_tissue=im.patch_in_tissue, 
                                   output_path=patch_image_output_path, 
                                   crop_size=self.imspire_param.ImageParam_CropSize)           
            #CellProfiler#
            if self.imspire_param.BasicParam_InputImageType=="H&E":
                spot_cp=ImSpiRE_HE_CellProfiler(image_file_path=f"{spot_image_output_path}", 
                                                cellprofiler_pipeline=self.imspire_param.CellProfilerParam_Pipeline)
                patch_cp=ImSpiRE_HE_CellProfiler(image_file_path=f"{patch_image_output_path}", 
                                              cellprofiler_pipeline=self.imspire_param.CellProfilerParam_Pipeline)
            else:
                spot_cp=ImSpiRE_IF_CellProfiler(image_file_path=f"{spot_image_output_path}", 
                                                cellprofiler_pipeline=self.imspire_param.CellProfilerParam_Pipeline)
                patch_cp=ImSpiRE_IF_CellProfiler(image_file_path=f"{patch_image_output_path}", 
                                                 cellprofiler_pipeline=self.imspire_param.CellProfilerParam_Pipeline)

            if self.imspire_param.BasicParam_Verbose:
                print("Extracting features of spot images...")

            spot_feature_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/FeatureResults/SpotFeature"
            
            if self.imspire_param.BasicParam_Overwriting:
                if len(os.listdir(spot_feature_output_path))!=0:
                    os.system(f"rm -rf {spot_feature_output_path}/*")

            spot_cp.compute_image_features(output_path=spot_feature_output_path,
                                           number_of_kernels=self.imspire_param.CellProfilerParam_KernelNumber)
            spot_cp.filter_image_features(output_path=spot_feature_output_path)

            if self.imspire_param.BasicParam_Verbose:
                print("Extracting features of patch images...")

            patch_feature_output_path=f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/FeatureResults/PatchFeature"
            
            if self.imspire_param.BasicParam_Overwriting:
                if len(os.listdir(patch_feature_output_path))!=0:
                    os.system(f"rm -rf {patch_feature_output_path}/*")

            patch_cp.compute_image_features(output_path=patch_feature_output_path,
                                            number_of_kernels=self.imspire_param.CellProfilerParam_KernelNumber)
            patch_cp.filter_image_features(output_path=patch_feature_output_path)

            spot_features=spot_cp.image_features.loc[imdata.pos_in_tissue_filter.index,]

            patch_cp.image_features.index=patch_cp.image_features.index.astype('int')
            patch_features=patch_cp.image_features.sort_index()
            patch_features.index=list(range(patch_features.shape[0]))

            spot_features.to_csv(f"{spot_feature_output_path}/Image_Features_Image_filter.txt", sep = "\t")
            patch_features.to_csv(f"{patch_feature_output_path}/Image_Features_Image_filter.txt", sep = "\t")

        im.patch_in_tissue.index=list(range(im.patch_in_tissue.shape[0]))
        im.patch_in_tissue.to_csv(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/{self.imspire_param.BasicParam_OutputName}_PatchLocations.txt", sep = "\t")

        
        #OT#
        spot_locations=np.array(imdata.pos_in_tissue_filter.loc[:,["pxl_row_in_fullres","pxl_col_in_fullres"]])
        patch_locations=np.array(im.patch_in_tissue.loc[:,["pxl_row","pxl_col"]])

        if self.imspire_param.BasicParam_Mode==1:
            spot_texture=pd.read_csv(f"{spot_feature_output_path}/spot_texture_features.txt",sep="\t",index_col=0)
            spot_intensity=pd.read_csv(f"{spot_feature_output_path}/spot_intensity_features.txt",sep="\t",index_col=0)

            patch_texture=pd.read_csv(f"{patch_feature_output_path}/patch_texture_features.txt",sep="\t",index_col=0)
            patch_intensity=pd.read_csv(f"{patch_feature_output_path}/patch_intensity_features.txt",sep="\t",index_col=0)

            if self.imspire_param.FeatureParam_FeatureType==0:
                spot_feature=pd.concat([spot_texture,spot_intensity],axis=1)
                patch_feature=pd.concat([patch_texture,patch_intensity],axis=1)
            elif self.imspire_param.FeatureParam_FeatureType==1:
                spot_feature=spot_texture.copy()
                patch_feature=patch_texture.copy()
            elif self.imspire_param.FeatureParam_FeatureType==2:
                spot_feature=spot_intensity.copy()
                patch_feature=patch_intensity.copy()
        else:
            spot_feature=pd.read_csv(f"{spot_feature_output_path}/Image_Features_Image_filter.txt",sep="\t",index_col=0)
            patch_feature=pd.read_csv(f"{patch_feature_output_path}/Image_Features_Image_filter.txt",sep="\t",index_col=0)
            
        
        spot_feature=spot_feature.dropna(axis=1)
        patch_feature=patch_feature.dropna(axis=1)

        commom_feature=list(set(spot_feature.columns).intersection(set(patch_feature.columns)))
        spot_feature = spot_feature.loc[:,commom_feature]
        patch_feature = patch_feature.loc[:,commom_feature]

        spot_feature=np.array(spot_feature.loc[imdata.pos_in_tissue_filter.index,])
        patch_feature=np.array(patch_feature.sort_index())
        
        exp_data=imdata.adata.to_df()
        exp_data=exp_data.loc[imdata.pos_in_tissue_filter.index,]

        if self.imspire_param.BasicParam_Verbose:
                print("Solving OT...")

        ot_solver=ImSpiRE_OT_Solver(spot_locations,patch_locations,
                                    spot_image_features=spot_feature,
                                    patch_image_features=patch_feature,
                                    spot_gene_expression=exp_data,
                                    random_state=self.imspire_param.BasicParam_RandomState)

        ot_solver.setup_cost_matrices(alpha=self.imspire_param.OptimalTransportParam_Alpha,
                                      num_neighbors=self.imspire_param.OptimalTransportParam_NumNeighbors)
        ot_solver.solve_OT(beta=self.imspire_param.OptimalTransportParam_Beta, 
                           epsilon=self.imspire_param.OptimalTransportParam_Epsilon,
                           numItermax=self.imspire_param.OptimalTransportParam_NumIterMax,
                           verbose=self.imspire_param.BasicParam_Verbose)
        
        np.save(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/SupplementaryResults/ot_M_alpha{self.imspire_param.OptimalTransportParam_Alpha}_beta{self.imspire_param.OptimalTransportParam_Beta}_epsilon{self.imspire_param.OptimalTransportParam_Epsilon}.npy",ot_solver.M)
        np.save(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/SupplementaryResults/ot_C1_alpha{self.imspire_param.OptimalTransportParam_Alpha}_beta{self.imspire_param.OptimalTransportParam_Beta}_epsilon{self.imspire_param.OptimalTransportParam_Epsilon}.npy",ot_solver.C1)
        np.save(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/SupplementaryResults/ot_C2_alpha{self.imspire_param.OptimalTransportParam_Alpha}_beta{self.imspire_param.OptimalTransportParam_Beta}_epsilon{self.imspire_param.OptimalTransportParam_Epsilon}.npy",ot_solver.C2)
        np.save(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/SupplementaryResults/ot_T_alpha{self.imspire_param.OptimalTransportParam_Alpha}_beta{self.imspire_param.OptimalTransportParam_Beta}_epsilon{self.imspire_param.OptimalTransportParam_Epsilon}.npy",ot_solver.T)
        
        if self.imspire_param.BasicParam_Verbose:
                print("Computing high resolution expression profiles...")

        exp_data_hr=compute_high_resolution_expression_profiles(exp_data,ot_solver.T)
        # exp_data_hr.to_csv(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/{self.imspire_param.BasicParam_OutputName}_ResolutionEnhancementResult.txt",sep="\t")
        adata_hr=sc.AnnData(exp_data_hr)
        adata_hr.write_h5ad(f"{self.imspire_param.BasicParam_OutputDir}/{self.imspire_param.BasicParam_OutputName}/{self.imspire_param.BasicParam_OutputName}_ResolutionEnhancementResult.h5ad")
        
        self.imdata=imdata
        self.ot_solver=ot_solver

    def output_supplementary_results(self):
        pass

