# -*- coding: utf-8 -*-
# @Author: huayuwei
# @Date:   2021-11-25 14:37:11
# @Last Modified by:   zhangyizhi
# @Last Modified time: 2023-01-04 17:01:03

###########
# imports #
###########

import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import os
import ot
import time
import tifffile as tif
import sys
import cv2
import scipy
import skimage
import skimage.morphology
import skimage.feature
import skimage.filters
from skimage import io, img_as_float32, morphology, exposure
from skimage.feature import greycomatrix, greycoprops
from sklearn.preprocessing import minmax_scale
from skimage.color import separate_stains, hdx_from_rgb
from itertools import product
from multiprocessing import Pool

#############
# functions #
#############

def read_image(image_file_path):
    image = tif.imread(image_file_path)
    return image

def mkdir(directory, name):
    if not os.path.exists(directory):
        raise ValueError(f'{directory} : No such directory!')
    elif not os.path.exists(f"{directory}/{name}"):
        os.makedirs(f"{directory}/{name}")
    else:
        print(f"The folder {directory}/{name} exists!")

def create_folder(output,name,overwriting=False):

    if not os.path.exists(output):
        raise ValueError(f'{output} : No such directory!')

    if os.path.exists(f"{output}/{name}"):
        if not overwriting:
            print(f"The project {output}/{name} exists! To avoid the overwriting, ImSpaRE exits. If overwriting is required, try the -O parameter. ")
            sys.exit(0)
        else:
            print(f"The project {output}/{name} exists! ImSpaRE overwrite it.")
    mkdir(output,name)
    mkdir(f"{output}/{name}","ImageResults")
    mkdir(f"{output}/{name}/ImageResults","SpotImage")
    mkdir(f"{output}/{name}/ImageResults","PatchImage")
    mkdir(f"{output}/{name}","FeatureResults")
    mkdir(f"{output}/{name}/FeatureResults","SpotFeature")
    mkdir(f"{output}/{name}/FeatureResults","PatchFeature")
    mkdir(f"{output}/{name}","SupplementaryResults")

def modify_IF_image(im, DAPI_channel, FiducialFrameChannel):
    """
    Modify the IF image to ensure the first channel is DAPI and to remove the fiducial frame channel 
    
    Parameters
    ----------
    im                      -- The original IF image.
    DAPI_channel            -- The channel of the DAPI.
    FiducialFrameChannel    -- The channel of the fiducial frame.

    Return
    ------
    The array of the modified IF image.
    """
    channels = list(set(range(im.shape[0]))-{DAPI_channel-1, FiducialFrameChannel-1})
    im_DAPI = im[int(DAPI_channel-1),:,:]
    if im.shape[0] == 3:
        channel2 = channels[0]
        im_channel2 = im[channel2,:,:]
        new_im = np.zeros(shape=(3, im.shape[1], im.shape[2]), dtype = 'uint16')
        new_im[0,:,:] = im_DAPI
        new_im[1,:,:] = im_channel2
    else:
        channel2 = channels[0]
        channel3 = channels[1]
        im_channel2 = im[channel2,:,:]
        im_channel3 = im[channel3,:,:]
        new_im = np.array([im_DAPI, im_channel2, im_channel3])

    return new_im

def filter_patch_in_tissue(cimg, scale, scale_row, scale_col, pos, pos_in_tissue, dist):
    """
    Generate patch locations and determine which patches are located on the tissue.

    Parameters
    ----------
    cimg            -- Compressed image (array).
    scale           -- The scale of compression.
    scale_row 		-- The compression scale of rows.
    scale_col		-- The compression scale of columns.
    pos             -- Location matrices of spots.
    pos_in_tissue   -- Location matrices of spots in tissue. 
    dist            -- This parameter specifies the interval between each patch.
    
    Returns
    ------
    patch_locations -- Location matrices of patches.
    """ 
    # Detect the tissue foreground and get the mask
    mask = tissue_detection(cimg, pos, pos_in_tissue, scale, scale_row, scale_col)

    # Generation patch locations
    sta = pos_in_tissue.describe()
        
    PXL_ROW_MIN = int(sta['pxl_row_in_fullres']['min'])
    PXL_ROW_MAX = int(sta['pxl_row_in_fullres']['max'])

    PXL_COL_MIN = int(sta['pxl_col_in_fullres']['min'])
    PXL_COL_MAX = int(sta['pxl_col_in_fullres']['max'])

    row_list = range(PXL_ROW_MIN, PXL_ROW_MAX, dist)
    col_list = range(PXL_COL_MIN, PXL_COL_MAX, dist)

    len_row = len(row_list)
    len_col = len(col_list)

    patch_locations = pd.DataFrame(
        {"row": np.repeat(range(len_row), len_col, axis=0),
        "col": list(range(len_col)) * len_row,
        "pxl_row": np.repeat(row_list, len_col, axis=0),
        "pxl_col": list(col_list) * len_row,
        "in_tissue": np.repeat([0], len_row * len_col)})

    for index, row in patch_locations.iterrows():
        patch_pxl_row = row['pxl_row']
        patch_pxl_col = row['pxl_col']
        
        srowpos = int(patch_pxl_row * scale_row)
        scolpos = int(patch_pxl_col * scale_col)

        height, width = mask.shape[:2]
        # Determine which patches are located on the tissue
        if 0 <= srowpos < height and 0 <= scolpos < width:
            tissue = int(mask[srowpos, scolpos] == 1)
        else:
            tissue = 0
            
        patch_locations['in_tissue'][index] = tissue
        
    return patch_locations

def filter_patch_in_tissue_2(scale, scale_row, scale_col, pos_in_tissue, dist, platform, output_dir, output_name, ff_channel=None, iter_count=50):
    """
    Generate patch locations and determine which patches are located on the tissue.
    """ 
    # Detect the tissue foreground and get the mask
    # os.system("backgroundremover -i compressed_image.png -a -ae 15 -o bgremover_output.png")
    
    ## mask1
    os.system(f"backgroundremover -i {output_dir}/{output_name}/SupplementaryResults/compressed_image.png -o {output_dir}/{output_name}/SupplementaryResults/bgremover_output.png")
    
    rembg = cv2.imread(f"{output_dir}/{output_name}/SupplementaryResults/bgremover_output.png", cv2.IMREAD_UNCHANGED)
    if len(rembg.shape) == 3:
        rembg = cv2.cvtColor(rembg, cv2.COLOR_RGB2GRAY)
    cv2.normalize(rembg, rembg, 255, 0, cv2.NORM_MINMAX)
    if rembg.dtype == "uint16":
        rembg = rembg.astype("uint8")

    # ret,th = cv2.threshold(cimg,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    ret, th = cv2.threshold(rembg, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_TRIANGLE)

    mask1 = np.where(th==0,0,1)

    ## mask2
    cimg = cv2.imread(f"{output_dir}/{output_name}/SupplementaryResults/compressed_image.png", cv2.IMREAD_UNCHANGED)
    if cimg.dtype == "uint16":
        cimg=cv2.convertScaleAbs(cimg, alpha=(255.0 / 65535.0))
    if cimg.shape[2]==4:
        cimg=np.delete(cimg,int(ff_channel-1),2)
    
    mask_tmp = np.zeros((cimg.shape[:2]),np.uint8)
    
    bgdModel=np.zeros((1,65),np.float64)
    fgdModel=np.zeros((1,65),np.float64)

    ROW_MIN=int(min(np.where(mask1==1)[0]))
    ROW_MAX=int(max(np.where(mask1==1)[0]))

    COL_MIN=int(min(np.where(mask1==1)[1]))
    COL_MAX=int(max(np.where(mask1==1)[1]))

    x=COL_MIN
    y=ROW_MIN
    w=COL_MAX-COL_MIN
    h=ROW_MAX-ROW_MIN

    rect=(x,y,w,h)
    mask_tmp, bgdModel, fgdModel = cv2.grabCut(cimg,mask_tmp,rect,bgdModel,fgdModel,iter_count,cv2.GC_INIT_WITH_RECT)
    mask2 = np.where((mask_tmp==2)|(mask_tmp==0),0,1).astype('uint8')

    mask=mask1*mask2
    
    cimg_new=cimg*mask[:,:,np.newaxis]
    cv2.imwrite(f"{output_dir}/{output_name}/SupplementaryResults/foreground.png", cimg_new)

    ## Generation patch locations
    if platform=="ST":
        PXL_ROW_MIN=int(min(np.where(mask==1)[0])/scale_row)
        PXL_ROW_MAX=int(max(np.where(mask==1)[0])/scale_row)
        
        PXL_COL_MIN=int(min(np.where(mask==1)[1])/scale_col)
        PXL_COL_MAX=int(max(np.where(mask==1)[1])/scale_col)
    else:
        sta = pos_in_tissue.describe()
            
        PXL_ROW_MIN = int(sta['pxl_row_in_fullres']['min'])
        PXL_ROW_MAX = int(sta['pxl_row_in_fullres']['max'])

        PXL_COL_MIN = int(sta['pxl_col_in_fullres']['min'])
        PXL_COL_MAX = int(sta['pxl_col_in_fullres']['max'])


    row_list = range(PXL_ROW_MIN, PXL_ROW_MAX, dist)
    col_list = range(PXL_COL_MIN, PXL_COL_MAX, dist)

    len_row = len(row_list)
    len_col = len(col_list)

    patch_locations = pd.DataFrame(
        {"row": np.repeat(range(len_row), len_col, axis=0),
        "col": list(range(len_col)) * len_row,
        "pxl_row": np.repeat(row_list, len_col, axis=0),
        "pxl_col": list(col_list) * len_row,
        "in_tissue": np.repeat([0], len_row * len_col)})

    for index, row in patch_locations.iterrows():
        patch_pxl_row = row['pxl_row']
        patch_pxl_col = row['pxl_col']
        
        srowpos = int(patch_pxl_row * scale_row)
        scolpos = int(patch_pxl_col * scale_col)

        height, width = mask.shape[:2]
        # Determine which patches are located on the tissue
        if 0 <= srowpos < height and 0 <= scolpos < width:
            tissue = int(mask[srowpos, scolpos] == 1)
        else:
            tissue = 0
            
        patch_locations['in_tissue'][index] = tissue
        
    return patch_locations

def extract_image_features(index, spot_img, img_shape2, feature_set):    
        features_dict={}
        
        ith_texture_f = []
        for c in range(img_shape2):
            glcm = greycomatrix(
                spot_img[:, :, c],
                distances=[1],
                # Angles are arranged in a counter clockwise manner, in radian.
                angles=[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
                levels=256,
                symmetric=True,
                normed=False,
            )
            glcm = glcm[1:, 1:]
            glcm = glcm / np.sum(glcm, axis=(0, 1))
            for feature_name in feature_set:
                ith_texture_f += greycoprops(glcm, feature_name)[0].tolist()
        # The first 6 features are intensity features, and the rest are Haralicks.

        # extract intensity features
        spot_mask = np.ones_like(spot_img[:, :, 0], dtype="bool")

        int_low = 0.2
        int_high = 0.8
        int_step = 0.1
        q_bins = np.arange(int_low, int_high, int_step)
        ith_int_f = []
        for c in range(img_shape2):
            for t in q_bins:
                ith_int_f.append(np.quantile(spot_img[:, :, c][spot_mask == True], t))
        
        features_dict["index"] = index
        features_dict["texture_feature"] = ith_texture_f
        features_dict["intensity_feature"] = ith_int_f
        
        return features_dict

def compressed_scale_HE(compressed_size, shape):
    origianl_size = max(shape)
    scale = 1.0 * compressed_size/origianl_size
    scale_row = 1.0 * compressed_size/shape[0]
    scale_col = 1.0 * compressed_size/shape[1]
    
    return min(1.0, scale), min(1.0, scale_row), min(1.0, scale_col)

def compressed_scale_IF(compressed_size, shape):
    origianl_size = max(shape)
    scale = 1.0 * compressed_size/origianl_size
    scale_row = 1.0 * compressed_size/shape[1]
    scale_col = 1.0 * compressed_size/shape[2]
    
    return min(1.0, scale), min(1.0, scale_row), min(1.0, scale_col)

def tissue_detection(cimg, pos, pos_in_tissue, scale, scale_row, scale_col):
    """
    Obtain the tissue foreground. 
    When the input is the IF stained image, we will caculate mask using two different ways
    to ensure that the tissue foreground is extracted as accurately as possible.

    Parameters
    ----------
    cimg                -- A numpy arrary of the compressed image.
    pos                 -- Location matrices of spots.
    pos_in_tissue       -- Location matrices of spots in tissue.
    scale               -- The float used to scale the spot coordinates to the (usu reduced-size) input image.
    scale_row           -- The float used to scale the row coordinates to the (usu reduced-size) input image.
    scale_col           -- The float used to scale the column coordinates to the (usu reduced-size) input image.

    Returns
    -------
    mask_new            -- A binary array where 1 = tissue and 0 = background.
    """
    spots = pos.drop(['barcode', 'in_tissue'], axis = 1)
    spots = [tuple(x) for x in spots.values]

    spot_diameter, average_spot_to_spot_distance = estimate_spot_diameter(spots)
    bounding_box = get_bounding_box(spots, scale, spot_diameter)
    # Caculate mask using the opencv grabcut algorithm
    mask = get_mask(cimg, bounding_box)

    # Caculate mask again based on spot locations
    mask2 = get_mask2(mask, average_spot_to_spot_distance, pos_in_tissue, scale_row, scale_col)
    mask_new = mask2*mask

    return mask_new

def estimate_spot_diameter(spots):
    """
    Function that takes a list of spots centre coordinates (row, col, y_pixel, x_pixel).
    The spot radius (in pixels) is estimated as a fixed proportion of the average spot to spot 
    distance between spots with even column coordinates.
    This estimate will only be valid for 5k arrays.
    If input spot coordinates is from a 1k array (both even and odd cordinates on same row) 
    the function will return an estimated spot radius of 1px (falling back to spot centers).

    Parameters
    ----------
    spots                           -- A list of spots centre coordinates (row, col, y_pixel, x_pixel).

    Returns
    -------
    spot_diameter                   -- The estimated spot radius in full-resolution pixels.
    average_spot_to_spot_distance   -- The distance between spots in full-resolution pixels.
    """
    # make a dictionary to sort spots by row and column coordinates
    positions_by_row_col = {}
    for spot in spots:
        row, col, rowpos, colpos = spot
        try:
            positions_by_row_col[row][col] = (colpos, rowpos)
        except KeyError:
            positions_by_row_col[row] = {col: (colpos, rowpos)}

    # calculate center to center distances of spots, for each spot calculate the distance to the
    # neighbour straight to the right i.e. the spot on the same row 2 columns later
    distances = []
    for row, cols in positions_by_row_col.items():
        for col in cols:
            try:
                x1, y1 = positions_by_row_col[row][col]
                x2, y2 = positions_by_row_col[row][col + 2]
                distances.append(((x2 - x1) ** 2 + (y2 - y1) ** 2) ** (1 / 2.0))
            except KeyError:
                pass  # End of row or missing spots

            # if both even and odd columns are present on the same row this is a old array
            # i.e fall back to no spot areas (estimated radius = 1 pixel)
            if col + 1 in positions_by_row_col[row]:
                return 1

    average_spot_to_spot_distance = np.mean(np.array(distances))

    spot_radius = average_spot_to_spot_distance * 0.65 * 0.5

    if spot_radius > 1:
        return 2.0 * spot_radius, average_spot_to_spot_distance
    else:
        return 1.0, average_spot_to_spot_distance

def get_bounding_box(spots, scale, padding):
    """
    Generates a rectangular bounding box based on a set of points WITHOUT the assumption that the box is axis-aligned.  
    The box is padded on each side by the amount padding specified in original (unscaled) image units.  
    The box returned is scaled.  
    Returns the (x,y) coordinates of four corners.  
    Makes heavy use of OpenCV convenience functions.

    Parameters
    ----------
    spots           -- A list of spots centre coordinates (row, col, y_pixel, x_pixel).
    scale           -- The float used to scale the spot coordinates to the (usu reduced-size) input image.
    padding         -- The padding to be used (generally spot diameter).

    Returns
    -------
    bounding_box    -- Box coordinates [ (x1,y1)...(x4,y4) ] as float.
    """
    x = np.array(
        [x for _, _, _, x in spots]
    )  # spots are row,col,y,x so we're taking x,y for openCV
    y = np.array([y for _, _, y, _ in spots])

    # compute the "rotated rectangle" surrounding the (padded) array - no assumption of axis alignment
    rrect = cv2.minAreaRect(np.column_stack((x, y)))
    # given an opencv rotated rectangle, add "padding" pixels to all sides and return a reconstructed RotatedRect
    (centerx, centery), (width, height), angle = rrect
    width += 2 * padding
    height += 2 * padding
    rrect = ((centerx, centery), (width, height), angle)

    # compute the corners of the rotated rectangle
    #
    # Given points 0, 1, 2, 3, these define four sides 0:1, 1:2, 2:3, and 3:0.
    # OpenCV code guarantees that 0:1 is opposite 2:3, 1:2 is opposite 3:0, and
    # 0:1 is adjacent to 1:2.  This implies that the lengths defined by 0:1 and 1:2
    # give the two side lengths of the rectangle.  Also, drawing the sides "in order"
    # will trace out a continuous contour.
    bbox = cv2.boxPoints(rrect)

    # scale the corners
    sbox = np.round(np.multiply(scale, bbox))

    return sbox

def get_mask(original, bounding_box):
    """
    The get_mask function takes an image (array) in grayscale and uses the opencv grabcut
    algorithm to detect tissue section(s) (foreground) on the glass slide surface (background).
    Markers for initialization of the grabcut algorithm are based on otsu tresholding of the
    grayscale and the gradient (local max-min) of the image.

    Parameters
    ----------
    original        -- The compressed image in grayscale.
    bounding_box    -- Box coordinates [ (x1,y1)...(x4,y4) ] as float.

    Returns
    -------
    mask            -- A numpy array that distinguishes tissue from background.
    """
    if len(original.shape) != 2:
        raise RuntimeError(
            "non-2D image (color?) passed to get_mask nD={}".format(len(original.shape))
        )

    gray = box(original, bounding_box)

    longest_side = max(original.shape[:2])
    # set sizes of objects and holes (assumes approx square input image)
    small_holes_size = int(longest_side / 2.0)
    small_objects_size = int(longest_side / 2.0)
    large_holes_size = int(longest_side * 50.0)

    # Calculate grayscale intensity Otsu treshold, often good proxy for where tissue is
    # (not as good when background and tissue has similar gray scale histograms)
    otsu_treshed = gray <= skimage.filters.threshold_otsu(gray)
    # Remove holes and tissue debris
    otsu_treshed = skimage.morphology.remove_small_objects(otsu_treshed, small_objects_size)
    otsu_treshed = skimage.morphology.remove_small_holes(otsu_treshed, small_holes_size)

    # Get the gradient (local max - local min) of the gray scale image,
    # high gradient usually indicates tissue (not as good for when larger
    # areas are out of focus or tissue is really smooth potentially
    # affected by low resolution images)
    gradient = skimage.filters.rank.gradient(gray, skimage.morphology.disk(5))

    # Binarize the gradient into two classes 1=FG and 0=BG using Otsu treshold
    inverted_grad = skimage.util.invert(gradient, signed_float=False)
    otsu_of_gradient = inverted_grad <= skimage.filters.threshold_otsu(inverted_grad)
    otsu_of_gradient = skimage.morphology.remove_small_objects(otsu_of_gradient, small_objects_size)
    otsu_of_gradient = skimage.morphology.remove_small_holes(otsu_of_gradient, small_holes_size)

    # Detect canny edges on the grayscale image (many edges usually indicate tissue)
    canny_edges = skimage.feature.canny(gray)
    closed_canny = skimage.morphology.closing(canny_edges)
    closed_canny = scipy.ndimage.distance_transform_edt(~closed_canny) <= longest_side * 0.01

    # Sum upp the two estimates of tissue placement
    # (gradient based and Outsu on grayscale intensity)
    BACKGROUND = 0
    DETECTED_BY_ONE_METHOD = 1
    DETECTED_BY_TWO = 2
    DETECTED_BY_ALL = 3
    otsu_sum = np.add(
        np.add(otsu_of_gradient.astype("uint8"), otsu_treshed.astype("uint8")),
        closed_canny.astype("uint8"),
    )

    # Start making markers for the grabcut
    markers_gc = np.zeros(gray.shape).astype("uint8")
    # to keep track of not yet classed vs obvious background
    classed = np.zeros(otsu_sum.shape).astype("uint8")

    ##### below is order dependent based on priority, pixels may be assign GC_BGD early and then
    ##### the same pixel may get GC_FGD later

    # If classed as background by both methods add a margin of 1% image longest side (in pixels) and
    # set to an obvious background pixels
    background = np.zeros(otsu_sum.shape).astype("uint8")
    background[otsu_sum == BACKGROUND] = 1
    background = scipy.ndimage.distance_transform_edt(background) >= longest_side * 0.01
    markers_gc[background == 1] = cv2.GC_BGD
    classed[background == 1] += 1

    # Take the two estimates (otsu_sum) fill all holes and set everything detected by at least one
    # method to be probable Background (instead of obvious background)
    # This is done so no holes will be classed as obvious background
    no_holes = np.zeros(otsu_sum.shape).astype("bool")
    no_holes[otsu_sum >= DETECTED_BY_ONE_METHOD] = True
    # remove_small_holes treats 0/1 mask different than false/true - use boolean
    no_holes = skimage.morphology.remove_small_holes(no_holes, large_holes_size)
    markers_gc[no_holes >= 1] = cv2.GC_PR_BGD
    classed[no_holes >= 1] += 1

    # If detected by at least one method set to be a possible foreground pixel
    markers_gc[otsu_sum >= DETECTED_BY_ONE_METHOD] = cv2.GC_PR_FGD
    classed[otsu_sum >= DETECTED_BY_ONE_METHOD] += 1

    # If detected by two methods add a margin of 5% (inward) image longest side (in pixels)
    # basically make the estimate smaller by some amount around the boundaries
    # set as an obvious foreground (object) pixel
    foreground = np.zeros(otsu_sum.shape).astype("uint8")
    foreground[otsu_sum == DETECTED_BY_TWO] = 1
    foreground = scipy.ndimage.distance_transform_edt(foreground) >= longest_side * 0.05
    markers_gc[foreground == 1] = cv2.GC_FGD
    classed[foreground == 1] += 1

    # If detected by all methods add a margin of 2.5% image longest side (in pixels)
    # same as above, but with smaller margin of safety due to greater certainty
    # set as an obvious foreground (object) pixel
    foreground = np.zeros(otsu_sum.shape).astype("uint8")
    foreground[otsu_sum == DETECTED_BY_ALL] = 1
    foreground = scipy.ndimage.distance_transform_edt(foreground) >= longest_side * 0.025
    markers_gc[foreground == 1] = cv2.GC_FGD
    classed[foreground == 1] += 1

    # Within tissue estimates (no_holes) but zero in outsu sum should be probable background
    # Essentially encourage interior holes when no method has indicated foreground.
    otsu_background = np.zeros(otsu_sum.shape).astype("uint8")
    otsu_background[otsu_sum == BACKGROUND] = 1
    probable_foreground_hole = np.add(otsu_background.astype("uint8"), no_holes.astype("uint8"))
    temp = np.ones(otsu_sum.shape).astype("uint8")
    temp[probable_foreground_hole == 2] = 0
    temp = scipy.ndimage.distance_transform_edt(temp) >= longest_side * 0.015
    temp = skimage.util.invert(temp, signed_float=False)
    probable_foreground_hole = temp.astype("uint8")
    markers_gc[temp == 1] = cv2.GC_PR_BGD

    # Set any unclassed pixels to be possible background
    markers_gc[classed == 0] = cv2.GC_PR_BGD

    # make the image compatible with the grabcut (color)
    im = cv2.cvtColor(gray.astype("uint8"), cv2.COLOR_GRAY2RGB)

    # run the opencv grabcut
    bgmodel = np.zeros((1, 65), dtype = "float64")
    fgmodel = np.zeros((1, 65), dtype = "float64")
    mask, bgmodel, fgmodel = cv2.grabCut(
        im, markers_gc, None, bgmodel, fgmodel, 6, cv2.GC_INIT_WITH_MASK
    )
    mask = np.where((mask == 2) | (mask == 0), 0, 1)

    # reembed the mask in the original image shape
    mask = unbox(
        np.zeros(original.shape[:2], dtype = np.uint8),
        bounding_box,
        mask,
        border_thickness=0,
        interp=cv2.INTER_NEAREST,
    )

    return mask

def box(original, bounding_box, interp = cv2.INTER_LINEAR):
    """
    Given a bounding box comprising coordinates of vertices (NOT x,y,w,h), extract a possibly rotated
    rectangle from the original image and return the subimage
    """
    trans, cols, rows = _get_bbox_transform(bounding_box)
    return cv2.warpAffine(
        original, trans, (cols, rows), borderMode=cv2.BORDER_REPLICATE, flags=interp
    )

def _get_bbox_transform(bounding_box):
    len1 = np.linalg.norm(bounding_box[1] - bounding_box[0])
    len2 = np.linalg.norm(bounding_box[2] - bounding_box[1])

    cols, rows = int(len1 + 0.5), int(len2 + 0.5)
    srcpts = np.float32(bounding_box[:3])
    dstpts = np.float32(np.array([[0, 0], [cols - 1, 0], [cols - 1, rows - 1]]))

    trans = cv2.getAffineTransform(src=srcpts, dst=dstpts)

    return trans, cols, rows

def unbox(original, bounding_box, boxed_image, border_thickness = None, interp = cv2.INTER_LINEAR):
    """
    Function for re-embedding a cropped bounding boxed image back in the original image (or shape).
    If border thickness is > 0, then we return a color version of what otherwise MUST be grayscale inputs.
    """
    if len(boxed_image.shape) == 3 and len(original.shape) == 2:
        original = cv2.cvtColor(original, cv2.COLOR_GRAY2RGB)
        output = np.zeros(original.shape, dtype="uint8")
        for i in range(3):
            output[:, :, i] = _unbox1(original[:, :, i], bounding_box, boxed_image[:, :, i], interp)
    else:
        output = _unbox1(original, bounding_box, boxed_image, interp)

    # add a blue border
    if border_thickness:
        if border_thickness < 0:
            raise ValueError("negative border thickness passed to unbox")

        if len(output.shape) == 2:
            output = cv2.cvtColor(output, cv2.COLOR_GRAY2RGB)

        cv2.drawContours(
            output, [np.intp(bounding_box)], 0, color=(29, 67, 122), thickness=2
        )

    return output

def _unbox1(original, bounding_box, boxed_image, interp = cv2.INTER_LINEAR):
    sentinel = np.iinfo(boxed_image.dtype).max
    img_max = np.max(boxed_image)
    if img_max == sentinel:
        cv2.normalize(boxed_image, boxed_image, sentinel - 1, 0, cv2.NORM_MINMAX)

    trans, _, _ = _get_bbox_transform(bounding_box)
    output = cv2.warpAffine(
        src=boxed_image,
        M=trans,
        flags=cv2.WARP_INVERSE_MAP | cv2.BORDER_CONSTANT | interp,
        dsize=(original.shape[1], original.shape[0]),
        borderValue=sentinel,
    )
    output[output == sentinel] = original[output == sentinel]
    return output

def get_mask2(mask, average_spot_to_spot_distance, pos_in_tissue, scale_row, scale_col):
    """
    For the IF stained image, the mask is calculated again based on spot positions, considering that using
    the opencv algorithm alone cannot always extract the foreground of IF images completely correctly.

    Parameters
    ----------
    mask                            -- The mask obtained by opencv grabcut algorithm.
    average_spot_to_spot_distance   -- The distance between spots in full-resolution pixels.
    pos_in_tissue                   -- A data frame of spots in tissue.
    scale_row                       -- The float used to scale the row coordinates to the (usu reduced-size) input image.
    scale_col                       -- The float used to scale the column coordinates to the (usu reduced-size) input image.
    Returns
    -------
    mask2                           -- A binary array where 1 = tissue, 0 = background.
    """
    # Convert mask to a 0-1 matrix
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if mask[i,j] > 0:
                mask[i,j] = 1
            else:
                mask[i,j] = 0
    
    pxl_dis = int(average_spot_to_spot_distance)

    # Caculate mask2
    mask2 = np.zeros(shape = mask.shape)
    for index, row in pos_in_tissue.iterrows():
        pxl_row = row['pxl_row_in_fullres']
        pxl_col = row['pxl_col_in_fullres']
        
        ssquare_row_down = (pxl_row - pxl_dis)*scale_row
        ssquare_row_up = (pxl_row + pxl_dis)*scale_row
        ssquare_col_left = (pxl_col - pxl_dis)*scale_col
        ssquare_col_right = (pxl_col + pxl_dis)*scale_col
        
        mask2[int(ssquare_row_down):int(ssquare_row_up), int(ssquare_col_left):int(ssquare_col_right)]  = 1
    
    return mask2
    


def compute_pca(X, n_components=30):
    """
    PCA from sklearn.

    Parameters
    ----------
    X            -- A matrix, obs as columns and vars as rows.
    n_components -- Number of principal components to compute.

    Returns
    -------
    Y -- PCA matrix
    """
    pca = PCA(n_components=n_components)
    Y=pca.fit_transform(X)
    return Y

def compute_distance_matrices(x1, x2=None, metric='euclidean'):
    """
    Compute distance between samples in x1 and x2 using function scipy.spatial.distance.cdist

    Parameters
    ----------
    x1 -- A matrix of size d.
    x2 -- Another matrix of size d.
    metric -- The distance metric used to compute.

    Returns
    -------
    M -- Normalized distance matrix computed with given metric.
    """
    if x2 is None:
        x2 = x1
    M = cdist(x1, x2, metric=metric)
    M /= M.max()
    return M

def compute_graph_distance_matrices(X, num_neighbors=5, metric='euclidean'):
    """
    Compute graph-based distance matrices for points in X.
    
    Parameters
    ----------
    X             -- The matrix of sample.
    num_neighbors -- Number of neighbors for nearest neighbors graph.
    metric        -- The distance metric used to calculate the k-Neighbors for each sample point.

    Returns
    -------
    shortest_paths_matrices -- The normalized shortest paths matrices.
    """
    # Shortest paths matrices
    kneighbors_graph_matrices = kneighbors_graph(X, num_neighbors, mode='connectivity', include_self=True,
                                                 metric=metric)
    shortest_paths_matrices = dijkstra(csgraph=csr_matrix(kneighbors_graph_matrices), directed=False,
                                      return_predecessors=False)
    shortest_paths_matrices_max = np.nanmax(shortest_paths_matrices[shortest_paths_matrices != np.inf])
    shortest_paths_matrices[shortest_paths_matrices > shortest_paths_matrices_max] = shortest_paths_matrices_max #set threshold for shortest paths

    # Set normalized shortest paths matrices
    shortest_paths_matrices = shortest_paths_matrices / shortest_paths_matrices.max()

    return shortest_paths_matrices

def entropic_fused_gromov_wasserstein(M, C1, C2, p, q, loss_fun='square_loss', epsilon=0.001, beta=0.5, **kwargs):
    '''
    Computes the entropic FGW transport between two graphs.

    Parameters
    ----------
    M        -- Metric cost matrix between features across domains.
    C1       -- Metric cost matrix representative of the structure in the source space.
    C2       -- Metric cost matrix representative of the structure in the target space.
    p        -- Distribution in the source space.
    q        -- Distribution in the target space.
    loss_fun -- Loss function used for the solver.
    beta     -- A trade-off parameter of Fused-gromov-Wasserstein transport.
    epsilon  -- Entropic regularization term.
    
    Returns
    -------
    T -- Optimal transportation coupling matrix.
    '''
    constC, hC1, hC2 = ot.gromov.init_matrix(C1, C2, p, q, loss_fun)

    G0 = p[:, None] * q[None, :]

    def f(G):
        return ot.gromov.gwloss(constC, hC1, hC2, G)

    def df(G):
        return ot.gromov.gwggrad(constC, hC1, hC2, G)

    T = ot.optim.gcg(p, q, (1-beta)*M, reg1=epsilon, reg2=beta, f=f, df=df, G0=G0, **kwargs)
    return T

def compute_high_resolution_expression_profiles(original_expression_profiles,T):
    '''
    Compute high resolution expression profiles.
    
    Parameters
    ----------
    original_expression_profiles -- Original gene expression profiles.
    T                            -- Optimal transportation coupling matrix.

    Returns
    -------
    high_resolution_expression_profiles -- High resolution expression profiles.
    '''
    gene_name=original_expression_profiles.columns
    original_expression_profiles=np.asarray(original_expression_profiles)
    high_resolution_expression_profiles=np.dot(original_expression_profiles.T,T)
    high_resolution_expression_profiles=pd.DataFrame(high_resolution_expression_profiles.T)
    high_resolution_expression_profiles.columns=gene_name
    # high_resolution_expression_profiles=high_resolution_expression_profiles.apply(lambda x: (x - np.min(x))*1e6 / np.sum(x))

    return high_resolution_expression_profiles



