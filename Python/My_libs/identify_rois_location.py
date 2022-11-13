"""
Objective
    To anatomically identify a region of interest (ROI) in a nifty brain image mask.

Description
    This script compares an input brain image mask, in the nifty format (.nii), with the AAL3v1 brain parcellation [1], to anatomically name the ROIs of the input image. The result will be an CSV file, with four columns:
        roi -> The index number of the ROI in the input brain image;
        Anatomical Description -> The name of the region the ROI is a part of;
        Hemisphere -> If this region is located at the right or left hemisphere or non-lateralized;
        Percentage -> how much of the region are contained in the ROIs
    
    In the input image, null values are considered empty areas (or no brain areas).

    Parameters Descriptions
        The parameters you'll have to change are in the Head section:
            baseatlas_path -> The location of the reference brain image file [string];
            baseatlasinfo_path -> The location of the .mat file, with the anatomical information of the reference image [string];
            targetaltas -> The input brain image mask [string];
            outputfolder_csv -> The output csv file [string];
    
    More
        The reference files are located on the page https://github.com/LexC/Brain_ROIS_ID/. And, as an example, I added a functional brain atlas (Shen_2mm.nii)[2] on the github page, so you can test out the script.


Limitations:
    For this specific script, the input image must be in the MNI space, with an isotropic voxel size of 1 or 2mm. But, by changing the reference Atlas and it's information file, this limitation can be overcome - just make sure that the information file are in the same format.

References:
    [1] Rolls, E.T., Huang, C.C., Lin, C.P., Feng, J. and Joliot, M., 2020. Automated anatomical labelling atlas 3. Neuroimage, 206, p.116189.
    [2] Shen, X., Tokoglu, F., Papademetris, X. and Constable, R.T., 2013. Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. Neuroimage, 82, pp.403-415.

Last Updated (yyy/mm/dd): 2021/11/02 
Contact for any doubts: carvalhoacastro@gmail.com
"""
#%%
""" Libraries ================================= """
import numpy as np
import pandas as pd
import nibabel as nib
import scipy.io as sio

#%% 
""" Head ================================= """

baseatlas_path = 'AAL3v1.nii'
baseatlasinfo_path = 'ROI_MNI_V7_List.mat'

targetaltas =  'Shen_2mm.nii'

outputfolder_csv = 'test.csv'

#%% 
""" Party ================================= """

# Loading Files
baseatlas = nib.load(baseatlas_path)
baseatlasimg = baseatlas.get_fdata()

taltas = nib.load(targetaltas)
tatlasimg = taltas.get_fdata()

# Checking files compatibility 
baseatlas_affine = baseatlas.affine
taltas_affine = taltas.affine

compt = list(baseatlas_affine == taltas_affine)
compt_aux = [] 
for i in compt:
    compt_aux += list(i)

if False in compt_aux:
    print('\nERROR: Files are not compatible.\nThe base and the target atlases must have the same affine tranformation matrix\n')
else:

    # Comparing the base Atlas with the target one
    voxelsid = [int(i) for i in np.unique(tatlasimg) if i != 0]

    rois = {}
    for roi_id in voxelsid:

        voxels_loc = np.where(tatlasimg == roi_id)

        corresp = []
        for vox,p in enumerate(voxels_loc[0]):
            x = voxels_loc[0][vox]
            y = voxels_loc[1][vox]
            z = voxels_loc[2][vox]

            corresp += [baseatlasimg[x,y,z]]

        corresp = np.unique(corresp,return_counts=True)
        
        rois[roi_id] = {
            'baseids': [int(i) for i in corresp[0]],
            'Percentage': corresp[1]/len(voxels_loc[0])
        }
    
    # Defining base atlas info
    baseatlasinfo = sio.loadmat(baseatlasinfo_path)['ROI']

    batlasinfodict = {}
    for i in range(baseatlasinfo.shape[1]):
        vals = baseatlasinfo[0,i] 
        keys = baseatlasinfo[0,i].dtype.descr

        batlasinfodict[vals[0][0][0]] = {}
        for j in range(1,5):
            batlasinfodict[vals[0][0][0]][keys[j][0]] = vals[j][0]

    batlasinfodict[0] = {
        keys[1][0]:'NB',
        keys[2][0]:'NoBrainArea',
        keys[3][0]:'No Brain Area',
        keys[4][0]:'N',
    }

    # Giving names to the locations
    for x in rois:
        rois[x]['Anatomical Description'] = []
        rois[x]['Hemisphere'] = []
        
        for y in rois[x]['baseids']:
            rois[x]['Anatomical Description'] += [batlasinfodict[y]['Anatomical_Description']]
            rois[x]['Hemisphere'] += [batlasinfodict[y]['Hemisphere']]

    # Saving results em CSV
    resultable = {
        'roi':[],
        'Anatomical Description': [],
        'Hemisphere': [],
        'Percentage': [],
    }

    hemps = {
        'L':'Left',
        'R':'Right',
        'N':'Non-Lateral'
    }

    for roi in rois:
        for i,x in enumerate(rois[roi]['baseids']):
            
            resultable['roi'] += [roi]
            for z in ['Anatomical Description','Percentage']:
                resultable[z] += [rois[roi][z][i]]
            resultable['Hemisphere'] += [hemps[rois[roi]['Hemisphere'][i]]]
            
    resultable = pd.DataFrame(resultable)
    resultable.to_csv(outputfolder_csv,index=False)
