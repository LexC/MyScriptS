""" Libraries ================================= """

from cgi import test
import os
#import h5py
import shutil

import numpy as np
import nibabel as nib
import networkx as nx

from scipy import stats


"""----------------------------------------------------------------------------
Files and Folders
----------------------------------------------------------------------------"""

def makefolder(path,isfile=False):
    """Creates a folder and it's subfolders if they do not exist

    Parameters
    ----------
    path : str
        The target folder's path
    isfile : bool, optional
        If Trues, it understands that the last entry of the path variable is a file and ignores it. By default False
    """

    if isfile:
        paths = path.split('/')[:-1]
    else:
        paths = path.split('/')

    path = ''
    for p in paths:
        path = f'{path}/{p}'
        if not os.path.isdir(path):
            os.makedirs(path)
            print(f'\nFolder created: {path}')

def nofileflag(filename,override=False):
    """Returns False if the file 'filename' exists and True if doesn't. It will always return True if overide is True.

    Parameters
    ----------
    filename : str
        file's path
    override : bool, optional
        by default False

    Returns
    -------
    bool
        ...
    """

    if override:
        flag = True
    else:
        flag = not os.path.isfile(filename)
    return flag

def copyfile(src,dst,override=True):    
    """Copies the file src to the file dst. 

    Parameters
    ----------
    src : str
        source file's path
    dst : str
        destination file's path
    override : bool, optional
        If true, and there is a dst file, it overrides; by default True
    """

    makefolder(dst,isfile = True)

    if nofileflag(src,override):    
        shutil.copyfile(src, dst)


# ---------- \ .txt file / ----------

def to_txt(filepath,input,param):
    """Writes or append a .txt file.

    Parameters
    ----------
    filepath : str
        file's path
    input : str
        The text that will be writen in the file
    param : str
        one of the 2 options bellow
            'w': write
            'a': append
    """

    if param in ['w','a']:
        with open(filepath, param) as file:
            file.write(input)
    else:
        print('The param variable must be "w" (for write) or "a" (for append)' )

def changestr_intxt(infile,oldstr,newstr,outfile=False):
    """Changes strings in a .txt file

    Parameters
    ----------
    infile : str
        the path for the input file
    oldstr : str
        the string to be changed
    newstr : str
        the new string
    outfile : str, optional
        the path for the output file
    """

    with open(infile, "r") as file:
        changedtxt = file.read().replace(oldstr, newstr)
    
    if type(outfile) == bool:
        outfile = infile
    
    to_txt(outfile,changedtxt,'w')
    

# ---------- \ HDF5 / ----------

def makeorupdate_attrs(h,name,value,override=True):
    """Creates or updates the attribute 'name' of the 'h' hdf5 file.

    Parameters
    ----------
    h : _type_
        h5py group or dataset pointer
    name : _type_
        _description_
    value : _type_
        _description_
    override : bool, optional
        _description_, by default True
    """
    if not h.attrs.__contains__(name):
        h.attrs.create(name,value)
    elif override:
        h.attrs.__delitem__(name)
        h.attrs.create(name,value)

def makeorupdate_dset(h,name,dset,override=True):
     
    if not h.__contains__(name):
        h.create_dataset(name,shape=dset.shape,dtype=dset.dtype, data=dset,compression="gzip")
    elif override:
        h.__delitem__(name)
        h.create_dataset(name,shape=dset.shape,dtype=dset.dtype, data=dset,compression="gzip")

def cc_attrs(inh5,outh5,name,override=True):
    '''
    check and copy Atribute
    '''
    attr = inh5.attrs[name]
    makeorupdate_attrs(outh5,name,attr,override)
    return attr


# ---------- \ Nifty / ----------

def nii2gz(file, delfile = False):
    
    print(f'FILE: {file}')
    niifile = nib.load(file)
    print('Compressing...')
    nib.save(niifile,f'{file}.gz')
    if delfile:
        print('Deleting the original file...')
        os.remove(file)
    print('DONE \n')

def nii2gz_folder(folder, delfiles = False):
    
    for root, _, files in os.walk(folder):
        for f in files:
            if f[-4:] == '.nii':
                
                file = root+'/'+f
                nii2gz(file, delfiles)
                
    print(' ------------> ALL DONE \n')

def atlas2rois(atlasdir,roisfolder,bin = True):

    atlas = nib.load(atlasdir)
    data = np.array(atlas.get_fdata())

    makefolder(roisfolder)
    
    for x in range(1,int(np.max(data))+1):
        
        roi = np.zeros(data.shape)
        
        if bin:
            roi[data==x] = 1

        roifile = nib.Nifti1Image(roi, atlas.affine, atlas.header)
        
        filename = f'{roisfolder}{str(x)}.nii.gz'

        nib.save(roifile,filename)

def AtlasNewROISValues(niidir, netvec, savename):
    """
    Change the ROI's IDvalue to the network's list values and save a .nii file
    with the new information.

    niidir = the location of a nifti file, with ROI'S enumerated from 1 to n.
    network = a list | len(network)=n
    savename = a string with the save name of the new nifti file
    """

    atlas = nib.load(niidir)
    img = atlas.get_fdata()

    base = np.zeros(img.shape)
    for i,x in enumerate(netvec):
        base[img == i+1] = x

    img = nib.Nifti1Image(img, atlas.affine, atlas.header)
    nib.save(img, savename)


"""----------------------------------------------------------------------------
Math
----------------------------------------------------------------------------"""

# ---------- \ Transformations / ----------
#asdfasdf

def ztransform(r,inv=False):
    """Fisher's z-Transformation.

    Args:
        r: int, float or numpy array

        inv: bool
            Makes the inverse Fisher's z Transformation

    Returns:
        float or numpy array : The output type depends on the input.
    """
    
    if not inv:
        z = np.log((1 + r) / (1 - r)) * (1/2)
    else:
        #z = (np.exp(2 * r) - 1)/(np.exp(2 * r) + 1)
        z = np.tanh(r)
    

    return z

def fisher_mean(array, axis=None,ignorenan = False):
    """Calculate the mean value(s), at the fisher's space, of one of it's dimensions

    Args:
        array: array_like
            A (n)-dimension array    

        axis: {int, tuple of int, None}, optional
            the dimension in which the calculation will be done.

        ignorenan: bool, opitional
            ignore NaNs values

    Returns:
        numpy.array
            A (n-1)-dimension array
    """

    array = ztransform(array)

    if ignorenan:
        array = np.nanmean(array, axis=axis)
    else:
        array = np.mean(array, axis=axis)

    array = ztransform(array,inv=True)

    return array

def NonlinearTransformation3D(vec, tmtx, inv=False):
    """ Calculates a nonlinear Transform of a 3d vector

    Args:
        vec: array like
            coordenates, len(vec) = 3

        tmtx: array like:
            Transformation matrix, .shape =(4, 4)

        inv: bool, optional
            if True, invertes the transformation matrix.

    Returns:
        list
            new coordenates, len(list) = 3.
    """

    vec = np.array(vec)
    agvec = np.array(np.append(vec, 1))

    if inv:
        tmtx = np.linalg.inv(np.array(tmtx))
        
    newvec = np.dot(tmtx,agvec)

    return newvec[:3]


"""----------------------------------------------------------------------------
Statstics
----------------------------------------------------------------------------"""

def t_test(a,b,p):
    """_summary_

    Args:
        a (_type_): _description_
        b (_type_): _description_
        p (_type_): _description_

    Returns:
        _type_: _description_

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
    """
    
    c = a+b
    c = stats.shapiro(c) 
    cp = c.pvalue #pvalue>p means it is parametric
    
    if cp>p:
        d = stats.ttest_ind(a,b)
    else:
        d = stats.mannwhitneyu(a,b)

    results = {
        'shapiro': {
            'parametric':cp>p,
            'values':c,
            'pthreshold':p
        },
        'test':d
    }
    
    return results


"""----------------------------------------------------------------------------
Graph Theory
----------------------------------------------------------------------------"""

def nparray2nxgraph(cmat,wethed=True,threshold=np.nan,thressign=np.nan,directional=False,delselfloop=True):

    cmat = np.array(cmat)

    if delselfloop:
        for i in range(cmat.shape[0]):
            cmat[i,i] = 0
    
    if not directional:
        for i in range(cmat.shape[0]):
            for j in range(cmat.shape[1]):
                if i<j:
                    cmat[i,j]=0

    flag = False
    if not wethed:
        if np.isnan(threshold):
            flag = True
        else:
            if thressign == 1 or thressign == -1:
                cmat = cmat*thressign
                aux = np.zeros(cmat.shape)
                aux[cmat>=threshold] = 1
                cmat = aux
            else:
                flag = True

    if flag:
        print('ERROR: You must define a threshold and the thressign (1 or -1)')
    else:
        return nx.convert_matrix.from_numpy_array(cmat)

"""----------------------------------------------------------------------------
MRI
----------------------------------------------------------------------------"""

def altas_seed(seed_mni,atlasdir):

   
    atlas = nib.load(atlasdir)
    atlasimg = atlas.get_fdata()

    x,y,z = np.round(NonlinearTransformation3D(seed_mni,atlas.affine))
    roinumb = int(atlasimg[x,y,z])
    seed_xyz = [x,y,z]

    return roinumb, seed_xyz

def MRISpaceTransf(vec,opt,inv=False):
    """
    
    https://brainmap.org/icbm2tal/
    """
    tmtx = {
    'icbm_spm2tal':np.array([
        [0.9464, 0.0034, -0.0026, -1.068],
        [-0.0083, 0.9479, -0.058, -1.0239],
        [0.0053, 0.0617,  0.901,  3.1883],
        [0.0, 0.0,  0.0,  1.0]]),    
    
    'icbm_fsl2tal':np.array([
        [0.9254, 0.0024, -0.0118, -1.0207],
        [-0.0048, 0.9316, -0.0871, -1.7667],
        [0.0152, 0.0883,  0.8924,  4.0926],
        [0.0, 0.0,  0.0,  1.0]])
    }

    return NonlinearTransformation3D(vec, tmtx[opt],inv)

 




