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
    """_summary_

    Parameters
    ----------
    h : _type_
        _description_
    name : _type_
        _description_
    dset : _type_
        _description_
    override : bool, optional
        _description_, by default True
    """
    if not h.__contains__(name):
        h.create_dataset(name,shape=dset.shape,dtype=dset.dtype, data=dset,compression="gzip")
    elif override:
        h.__delitem__(name)
        h.create_dataset(name,shape=dset.shape,dtype=dset.dtype, data=dset,compression="gzip")

def cc_attrs(inh5,outh5,name,override=True):
    """Check and copy Atribute
    

    Parameters
    ----------
    inh5 : _type_
        _description_
    outh5 : _type_
        _description_
    name : _type_
        _description_
    override : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    """
    attr = inh5.attrs[name]
    makeorupdate_attrs(outh5,name,attr,override)
    return attr


# ---------- \ Nifty / ----------

def nii2gz(file, delfile = False):
    """_summary_

    Parameters
    ----------
    file : _type_
        _description_
    delfile : bool, optional
        _description_, by default False
    """
    print(f'FILE: {file}')
    niifile = nib.load(file)
    print('Compressing...')
    nib.save(niifile,f'{file}.gz')
    if delfile:
        print('Deleting the original file...')
        os.remove(file)
    print('DONE \n')

def nii2gz_folder(folder, delfiles = False):
    """_summary_

    Parameters
    ----------
    folder : _type_
        _description_
    delfiles : bool, optional
        _description_, by default False
    """
    for root, _, files in os.walk(folder):
        for f in files:
            if f[-4:] == '.nii':
                
                file = root+'/'+f
                nii2gz(file, delfiles)
                
    print(' ------------> ALL DONE \n')

def atlas2rois(atlasdir,roisfolder,bin = True):
    """_summary_

    Parameters
    ----------
    atlasdir : _type_
        _description_
    roisfolder : _type_
        _description_
    bin : bool, optional
        _description_, by default True
    """
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

def AtlasNewROISValues(atlas, netvec, savename):
    """Change the ROI's IDvalue to the network's list values and save a .nii file with the new information.

    Parameters
    ----------
    atlas : str
        the location of a nifti file, with ROI'S enumerated from 1 to n.
    network : list
        len(network)=n
    savename : str
        The name of the new nifti file
    """
    if type(atlas) == str:
        atlas = nib.load(atlas)
    img = atlas.get_fdata()

    base = np.zeros(img.shape)
    for i,x in enumerate(netvec):
        base[img == i+1] = x

    newimg = nib.Nifti1Image(base, atlas.affine, atlas.header)
    makefolder(savename,isfile=True)
    nib.save(newimg, savename)


"""----------------------------------------------------------------------------
Math
----------------------------------------------------------------------------"""

# ---------- \ Transformations / ----------
#asdfasdf

def ztransform(r,inv=False):
    """Fisher's z-Transformation.

    Parameters
    ----------
    r : int, float or numpy array
        _description_
    inv : bool, optional
        Makes the inverse Fisher's z Transformation, by default False

    Returns
    -------
    float or numpy array
        The output type depends on the input.
    """

    if not inv:
        z = np.log((1 + r) / (1 - r)) * (1/2)
    else:
        #z = (np.exp(2 * r) - 1)/(np.exp(2 * r) + 1)
        z = np.tanh(r)
    

    return z

def FisherMean(array, axis=None,ignorenan = False):
    """Calculate the mean value(s), at the fisher's space, of one of it's dimensions

    Parameters
    ----------
    array : array_like
        A (n)-dimension array    
    axis : int, tuple of int, None, optional
        The dimension in which the calculation will be done, by default None
    ignorenan : bool, optional
        ignore NaNs values, by default False

    Returns
    -------
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

    Parameters
    ----------
        vec: array like
            coordenates, len(vec) = 3

        tmtx: array like:
            Transformation matrix, .shape =(4, 4)

        inv: bool, optional
            if True, invertes the transformation matrix.By default False

    Returns
    -------
        list
            new coordenates, len(list) = 3.
    """

    vec = np.array(vec)
    agvec = np.array(np.append(vec, 1))

    if inv:
        tmtx = np.linalg.inv(np.array(tmtx))
        
    newvec = np.dot(tmtx,agvec)

    return newvec[:3]

def SphereMatrix(d):

    mat = np.zeros((d,d,d))
    r = d/2
    for x in range(d):
        for y in range(d):
            for z in range(d):
                
                if np.sqrt((x-r)**2+(y-r)**2+(z-r)**2)<=r:
                    mat[x,y,z]=1
    
    return mat

"""----------------------------------------------------------------------------
Statstics
----------------------------------------------------------------------------"""
def isparametric(vec,p=0.05):
    """
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html

    Parameters
    ----------
    vec : array like
        _description_
    p : float, optional
        _description_, by default 0.05

    Returns
    -------
    bool
        _description_
    """
    return stats.shapiro(vec).pvalue > p

def ttest(a,b,parametric_p=0.5):
    """_summary_

    Parameters
    ----------
        a : array like
            _description_
        b : array like
            _description_
        p : float, optional
            _description_, by default 0.05

    Returns
    -------
        _type_: _description_

    
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
    """
    par = isparametric(a,parametric_p) and isparametric(b,parametric_p)
    if par:
        d = stats.ttest_ind(a,b)
    else:
        d = stats.mannwhitneyu(a,b)

    results = {
        'Parametric': par,
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

def altas_mni(mni,atlas):
   
    if type(atlas) == str:
        atlas = nib.load(atlas)

    atlasimg = atlas.get_fdata()
    
    x,y,z = np.round(NonlinearTransformation3D(np.array(mni),atlas.affine,inv=True))
    x = int(x); y = int(y); z = int(z)
    
    roinumb = int(atlasimg[x,y,z])
    coor = [x,y,z]

    return roinumb, coor

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

"""----------------------------------------------------------------------------
NIRS
----------------------------------------------------------------------------"""

def nirs_shitty_atlas(sd,atlas,channels,savename):
    
    sphere = SphereMatrix(sd)
    
    if type(atlas) == str:
        atlas = nib.load(atlas)
    aimg = atlas.get_fdata()

    nirsimg = np.zeros(aimg.shape)
    
    
    coor = {'x':[],'y':[],'z':[]}
    for i,cn in enumerate(channels):
        _,c = altas_mni(channels[cn],atlas)
        
        for j,u in enumerate(coor):
            aux = int(sd/2)
            a0 = c[j]-aux
            af = c[j]+aux
            
            coor[u] = [a0,af]

        for x1,x2 in enumerate(range(coor['x'][0],coor['x'][1])):
            for y1,y2 in enumerate(range(coor['y'][0],coor['y'][1])):
                for z1,z2 in enumerate(range(coor['z'][0],coor['z'][1])):
                    if sphere[x1,y1,z1] != 0:
                        nirsimg[x2,y2,z2] = cn
    
    aimg[aimg==0] = 1e5
    aimg[aimg!=1e5] = 0
    nirsimg = nirsimg-aimg
    nirsimg[nirsimg<0] = 0

    img = nib.Nifti1Image(nirsimg, atlas.affine, atlas.header)
    nib.save(img, savename)
    






