import glob 
from matplotlib.image import imread 
from matplotlib.pyplot import imshow, plot , figure, show
import pylab 
from scipy.ndimage.interpolation import rotate, shift 

from scipy.optimize import minimize 

import numpy as np 
from numpy import NaN, zeros , arange
import functools 
#from sympy.geometry.util import convex_hull 
from scipy.spatial.qhull import ConvexHull 
from itertools import combinations 
import os
import math
from scipy.optimize.optimize import fmin
from functools import partial
from constants import Dir
import image_utils

piece_files = sorted( glob.glob("build/piece*single_rot.png") ) 



    






 
def normpdf(x, mu, sigma, radius_max):
    
    u = (x-mu)/abs(sigma)
    y = (1/(np.sqrt(2*np.pi)*abs(sigma)))*np.exp(-u*u/2)
    y[ np.abs(x-mu) > radius_max  ] = 0.
    return y

def gaussian2d(sz, mu, sigma, radius_max = 20. ):
    x = normpdf(x=arange(sz[0]), mu=mu[0], sigma = sigma, radius_max=radius_max)
    y = normpdf(x=arange(sz[1]), mu=mu[1], sigma = sigma, radius_max=radius_max)
    t = np.outer(x,y)
    sc = 1/np.max(t)
    #sc=1.
    return t*sc


class GaussianImMode:
    Add = "Add"
    Sub = "Sub"

def build_simple_img(sz, x0,x1,y0,y1, (gX, gY, w0),  mode ):
    im = zeros(sz)
    
    im[x0:x1,y0:y1] = 1.
        
    if mode ==GaussianImMode.Add: 
        im += gaussian2d(sz, mu=(gX,gY), sigma=w0) * 2
    elif mode ==GaussianImMode.Sub: 
        im -= gaussian2d(sz, mu=(gX,gY), sigma=w0) * 2
    else:
        assert False
        
    
    return im


def build_img_from_p(p, sz, x0,x1,y0,y1, direction, mode=GaussianImMode.Add):
    (d0, w0) = p
    X = (x0+x1)/2.
    Y = (y0+y1)/2.
    
    if direction == Dir.Right:
        gX, gY = x1+np.fabs(d0), Y
    elif direction == Dir.Left:
        gX, gY = x0-np.fabs(d0), Y
    elif direction == Dir.Up:
        gX, gY = X, y1+np.fabs(d0)
    elif direction == Dir.Down:
        gX, gY = X, y0 - np.fabs(d0)
    else:  
        assert False
    im = build_simple_img(sz, x0,x1,y0,y1, (gX, gY, w0), mode=mode )
    
    im = np.clip(im,0,1)
    
    return im


def min_func_x(p, x0,x1,y0,y1, im_norm, direction, mode):
    im = build_img_from_p(p=p, sz=im_norm.shape, x0=x0, x1=x1, y0=y0, y1=y1, direction=direction, mode=mode)
    diff = im_norm - im
    res = np.sum(diff**2)
    return res
       
       
       
    
reses = []    




def fit_piece(fname,fname_idx): 
    print fname 
    im_piece = imread(fname) 
    im_piece = (im_piece[:,:,0]*255.).astype('int') 
    
    
    # Derotate the piece:
    im_rot = image_utils.align_to_axes(im_piece)
    

    # Find the corners:
    pts = image_utils.get_corners_of_piece(im_rot)    
    
    
    im_rot =im_rot.astype('float')
    
    
    
    pylab.close('all')
    pts= np.array(pts)
    print pts
    x0,x1 = np.min(pts[:,0]), np.max(pts[:,0])
    y0,y1 = np.min(pts[:,1]), np.max(pts[:,1])
    print x0,x1
    
    w = ((x1-x0)/2 + (y1-y0)/2)/2. * 0.3
    
    
    p0 = (0,w,) 

    im_norm = im_rot - np.min(im_rot)
    im_norm/= np.max(im_norm) 
        
    rect = zeros(im_norm.shape)
    rect[x0:x1,y0:y1] = 1.
    
    ## Adding:
    dirs = [Dir.Right, Dir.Left, Dir.Up, Dir.Down]
    im_opt_dirs_out = []
    mode = GaussianImMode.Add
    for i, dir in enumerate(dirs): 
        p_dir_out = fmin(partial(min_func_x, x0=x0,x1=x1,y0=y0,y1=y1, im_norm=im_norm, direction=dir,mode=mode), p0)
        im_opt_dirs_out.append(build_img_from_p(p_dir_out, sz=im_rot.shape, x0=x0,x1=x1,y0=y0,y1=y1, direction=dir, mode=mode))
        reses.append( p_dir_out )
    
    
    im_opt_out = im_opt_dirs_out[0] + im_opt_dirs_out[1] + im_opt_dirs_out[2] + im_opt_dirs_out[3] - (3*rect)
    im_opt_out = np.clip(im_opt_out,0.,1.)
    
    
    f,axes = pylab.subplots(2)
    pylab.sca(axes[0])
    imshow(im_opt_out.T)
    pylab.sca(axes[1])
    imshow(im_rot.T)
    #pylab.title("p[1]=%f"%p_right_out[1] )
    pylab.savefig('analysis/%03d_bb_fit1_out.png'%fname_idx)





    
    dir = Dir.Right
    mode = GaussianImMode.Sub
    im_opt_dirs_in = []
    for dir in dirs:
        #i=0
        p_dir_in = fmin(partial(min_func_x, x0=x0,x1=x1,y0=y0,y1=y1, im_norm=im_norm, direction=dir,mode=mode), p0)
        im_opt_dirs_in.append( build_img_from_p(p_dir_in, sz=im_rot.shape, x0=x0,x1=x1,y0=y0,y1=y1, direction=dir,mode=mode) )
        reses.append( ( -p_dir_out[0], -p_dir_out[1]) )
    
    im_opt_in = im_opt_dirs_in[0] + im_opt_dirs_in[1] + im_opt_dirs_in[2] + im_opt_dirs_in[3] - (3*rect)
    im_opt_in = np.clip(im_opt_in,0.,1.)
    
    f,axes = pylab.subplots(2)
    pylab.sca(axes[0])
    imshow(im_opt_in.T)
    pylab.sca(axes[1])
    imshow(im_rot.T)
    pylab.savefig('analysis/%03d_bb_fit1_in.png'%fname_idx)


    f,axes = pylab.subplots(2)
    pylab.sca(axes[0])
    imshow(( im_opt_in + im_opt_out - rect).T)
    pylab.sca(axes[1])
    imshow(im_rot.T)
    pylab.savefig('analysis/%03d_bb_fit_piece.png'%fname_idx)

    

    
    pylab.close('all')
    #pylab.show() 
    
    
    
if not os.path.exists("analysis"):
    os.makedirs("analysis")

for i,piece_file in enumerate(piece_files): 
    print piece_file 
    fit_piece(piece_file,i) 
    #if(i>20):
    #    break

reses = np.array(reses)
pylab.figure()
pylab.scatter( reses[:,0],reses[:,1] )
pylab.xlabel("GaussianDist (d0)")
pylab.ylabel("Strength (m)")
pylab.figure()
#pylab.scatter(reses[:,1] )
#pylab.scatter( reses[:,0],reses[:,2] )
pylab.show()
    #break
    #break
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

def test_fits():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.interpolate as si
    
    points = [[0, 0], [0, 2], [2, 3], [4, 0], [6, 3], [8, 2], [8, 0]];
    points = np.array(points)
    x = points[:,0]
    y = points[:,1]
    
    t = range(len(points))
    ipl_t = np.linspace(0.0, len(points) - 1, 100)
    
    x_tup = si.splrep(t, x, k=3)
    y_tup = si.splrep(t, y, k=3)
    
    x_list = list(x_tup)
    xl = x.tolist()
    x_list[1] = xl + [0.0, 0.0, 0.0, 0.0]
    
    y_list = list(y_tup)
    yl = y.tolist()
    y_list[1] = yl + [0.0, 0.0, 0.0, 0.0]
    
    x_i = si.splev(ipl_t, x_list)
    y_i = si.splev(ipl_t, y_list)
    
    #==============================================================================
    # Plot
    #==============================================================================
    
    fig = plt.figure()
    
    ax = fig.add_subplot(231)
    plt.plot(t, x, '-og')
    plt.plot(ipl_t, x_i, 'r')
    plt.xlim([0.0, max(t)])
    plt.title('Splined x(t)')
    
    ax = fig.add_subplot(232)
    plt.plot(t, y, '-og')
    plt.plot(ipl_t, y_i, 'r')
    plt.xlim([0.0, max(t)])
    plt.title('Splined y(t)')
    
    ax = fig.add_subplot(233)
    plt.plot(x, y, '-og')
    plt.plot(x_i, y_i, 'r')
    plt.xlim([min(x) - 0.3, max(x) + 0.3])
    plt.ylim([min(y) - 0.3, max(y) + 0.3])
    plt.title('Splined f(x(t), y(t))')
    
    ax = fig.add_subplot(234)
    for i in range(7):
        vec = np.zeros(11)
        vec[i] = 1.0
        x_list = list(x_tup)
        x_list[1] = vec.tolist()
        x_i = si.splev(ipl_t, x_list)
        plt.plot(ipl_t, x_i)
    plt.xlim([0.0, max(t)])
    plt.title('Basis splines')
    plt.show()
    
    
    
    
    
    
    
    
    
    #im_size=(30.,30.)
    #im_centre = (15.,15.)
    #ma2 = 0.
    #ma1 = 0.0
    #ma0 = im_centre[0]
    #
    #mb2 = 0.
    #mb1 = 0.
    #mb0 = im_centre[1]
    
    #x = np.arange(0, im_size[0])
    #y = np.arange(0, im_size[1])
    
    #x = x - im_centre[0]
    #y = y - im_centre[1]
    
   
    
    #x_valid_mat = np.tile(( y > ((x)**2 * ma2 + x*ma1 + ma0) ), (im_size[1],1 ) )
    #y_valid_mat = np.tile(( x > ((y)**2 * mb2 + y*mb1 + mb0) ), (im_size[0],1 ) ).T
    #print x_valid_mat.shape
    #print x_valid_mat


    #im = x_valid_mat & y_valid_mat
    #imshow(im.T ) 
    #pylab.show()
    
    #, cmap, norm, aspect, interpolation, alpha, vmin, vmax, origin, extent, shape, filternorm, filterrad, imlim, resample, url, hold)
    #
    
    

#test_fits()
#exit(0)    