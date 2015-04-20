import glob 
from matplotlib.image import imread 
from matplotlib.pyplot import imshow, plot , figure, show, subplots, legend,\
    grid
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
from constants import Dir, EdgeType
import image_utils
from caching import joblib_memory
import pickle
from fitbspline import test_fits

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
        print direction
        assert False
    im = build_simple_img(sz, x0,x1,y0,y1, (gX, gY, w0), mode=mode )
    
    im = np.clip(im,0,1)
    
    return im


def min_func_x(p, x0,x1,y0,y1, im_norm, direction, mode):
    im = build_img_from_p(p=p, sz=im_norm.shape, x0=x0, x1=x1, y0=y0, y1=y1, direction=direction, mode=mode)
    diff = im_norm - im
    res = np.sum(diff**2)
    return res
       


def get_edge_types( im_norm,(x0,x1,y0,y1), img_idx ):  
    rect = zeros(im_norm.shape)
    rect[x0:x1,y0:y1] = 1.

    w = ((x1-x0)/2 + (y1-y0)/2)/2. * 0.3
    p0 = (0,w,) 
    

    do_plot=False
    
    im_opt_dirs_out = []
    im_opt_dirs_in = []
    
    edge_types = []
    for direction in Dir.directions: 
        p_dir_out = fit_knobdule(x0, x1, y0, y1, p0, im_norm, direction, mode=GaussianImMode.Add)
        im_out = build_img_from_p(p_dir_out, sz=im_norm.shape, x0=x0,x1=x1,y0=y0,y1=y1, direction=direction, mode=GaussianImMode.Add)
        im_opt_dirs_out.append(im_out)
        reses.append( p_dir_out )

        p_dir_in = fit_knobdule(x0, x1, y0, y1, p0, im_norm, direction, mode=GaussianImMode.Sub) 
        im_in = build_img_from_p(p_dir_in, sz=im_norm.shape, x0=x0,x1=x1,y0=y0,y1=y1, direction=direction,mode=GaussianImMode.Sub)
        im_opt_dirs_in.append(im_in)
        reses.append( ( -p_dir_in[0], -p_dir_in[1]) )

        im_out_diff = im_out - rect
        im_in_diff = im_in - rect
        
        xs,origin,sum_axis ={
            Dir.Right: (range(im_norm.shape[0]), x1, 1),       
            Dir.Left:  (range(im_norm.shape[0]), x0, 1),
            Dir.Up: (range(im_norm.shape[1]), y1, 0),       
            Dir.Down:  (range(im_norm.shape[1]), y0, 0),
        }[direction]
        
            
        out_sum = np.sum( im_out_diff, sum_axis)
        in_sum = np.sum( im_in_diff, sum_axis)
        out_weighted = np.fabs(xs-origin)**2 * out_sum 
        in_weighted = np.fabs(xs-origin)**2 * in_sum
        
        if do_plot:
            f,axes = subplots(2)
            axes[0].imshow( im_out_diff.T )
            axes[1].imshow( im_in_diff.T * -1 )
        
            figure()
            plot(xs,out_sum, label="out" )
            plot(xs,in_sum, label="in" )
            plot(xs,out_weighted, label="out-weighted" )
            plot(xs,in_weighted, label="in-weighted" )
            grid()
            plot([origin], [0], 'x')
            legend()
            pylab.show()
        
        weight = np.sum(out_weighted) + np.sum(in_weighted)
        weights.append(weight)
            
        if(weight < -5000):
            edge_types.append( EdgeType.Innie)
        elif(weight > 5000):
            edge_types.append( EdgeType.Outtie )
        else:
            edge_types.append( EdgeType.Flattie )

    
    im_opt_out = np.clip(sum(im_opt_dirs_out) - (3*rect),0.,1.)     
    im_opt_in = np.clip(sum(im_opt_dirs_in) - (3*rect) ,0.,1.)
    
    f,axes = pylab.subplots(2,2)
    axes[0][0].imshow(im_norm.T)
    axes[0][1].imshow(( im_opt_in + im_opt_out - rect).T)
    axes[1][0].imshow(im_opt_in.T)
    axes[1][1].imshow(im_opt_out.T)
    pylab.savefig('analysis/%03d_bb_fit1.png'%img_idx)

    return tuple(edge_types)








     
       
    
reses = []    
weights = []



@joblib_memory.cache
def fit_knobdule(x0, x1, y0, y1, p0, im_norm, direction, mode):
    p_dir_out = fmin(partial(min_func_x, x0=x0, x1=x1, y0=y0, y1=y1, im_norm=im_norm, direction=direction, mode=mode), p0)
    return p_dir_out


def fit_piece(fname,fname_idx): 
    print fname 
    im_piece = imread(fname) 
    im_piece = (im_piece[:,:,0]*255.).astype('int') 
    
    
    # Derotate the piece:
    im_rot = image_utils.align_to_axes(im_piece)
    

    # Find the corners:
    pts = np.array(image_utils.get_corners_of_piece(im_rot) )    
    x0,x1 = np.min(pts[:,0]), np.max(pts[:,0])
    y0,y1 = np.min(pts[:,1]), np.max(pts[:,1])

    
    
    im_rot =im_rot.astype('float')
    im_norm = im_rot - np.min(im_rot)
    im_norm/= np.max(im_norm) 
      

    # Find if each edge is an 'innie', and 'outtie', or a 'flattie':
    edge_types = get_edge_types( im_norm,(x0,x1,y0,y1), img_idx=fname_idx )
    
    #print type(edge_types)
    
    #if fname_idx != 2:
    #    return
    
    test_fits(im_norm, (x0,x1,y0,y1), edge_types, piece_idx = fname_idx)
    
    #assert(0)
    pylab.show()
    pylab.close('all')
    
    
    
    print edge_types  
    #pylab.show() 
    
    

    
if not os.path.exists("analysis"):
    os.makedirs("analysis")

for i,piece_file in enumerate(piece_files): 
    print piece_file 
    fit_piece(piece_file,i) 
    #if(i>5):
    #    break

reses = np.array(reses)
weights = np.array(weights)


pylab.figure()
pylab.scatter( reses[:,0],reses[:,1] )
pylab.xlabel("GaussianDist (d0)")
pylab.ylabel("Strength (m)")

pylab.figure()
pylab.plot( np.fabs(reses[:,0]) * reses[:,1], 'o' )


pylab.figure()
pylab.plot(range(len(weights)), weights, 'o' )

with open("weights.pckl","w") as f:
    pickle.dump(weights, f)
pylab.show()    
    
    
    
    
    
    
    
    
    
    
    