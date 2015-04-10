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

piece_files = sorted( glob.glob("build/piece*single_rot.png") ) 


class Dir:
    Right = "Right"
    Left = "Left"
    Up = "Up"
    Down = "Down"
    

def eval_pts_as_rect(pt0,pt1,pt2,pt3):
    #print "Eval'ing"
    #print pt0
    #print pt1
    #print pt2
    #print pt3
    
    A = pt1 - pt0
    B = pt2 - pt1
    C = pt3 - pt2
    D = pt0 - pt3
    
    lenA = np.sqrt( np.dot(A,A) )
    lenB = np.sqrt( np.dot(B,B) )
    lenC = np.sqrt( np.dot(C,C) )
    lenD = np.sqrt( np.dot(D,D) )
    
    angle0 = np.dot(D,A)/(lenD*lenA)
    angle1 = np.dot(A,B)/(lenA*lenB)
    angle2 = np.dot(B,C)/(lenB*lenC)
    angle3 = np.dot(C,D)/(lenC*lenD)
    
    d = np.array( [angle0,angle1,angle2,angle3])
    
    res =  np.dot(d,d)
    #print res
    return res
    
def eval_area_pts(pt0,pt1,pt2,pt3):
    A = pt1 - pt0
    B = pt2 - pt1
    C = pt3 - pt2
    D = pt0 - pt3
    
    lenA = np.sqrt( np.dot(A,A) )
    lenB = np.sqrt( np.dot(B,B) )        
    lenC = np.sqrt( np.dot(C,C) )
    lenD = np.sqrt( np.dot(D,D) )
    
    return ((lenA+lenC)/2.) *  ((lenB+lenD)/2.) 



def get_corners_of_piece(piece):

    imshow(piece.T) 
    #pylab.show() 
    
    points = np.transpose( np.nonzero(piece) ) 
    #print points.shape 
    x = ConvexHull(points) 
    #print x 
    #print x.simplices 
    cvx_hull_pts =  points[x.vertices,0:2] # for simplex in x.simplices] 
    plot(cvx_hull_pts[:,0], cvx_hull_pts[:,1],'rx') 
    #print cvx_hull_pts 
    
    
    squares = combinations( range(cvx_hull_pts.shape[0]), 4 ) 
    
    res = []
    for i,square_idxs in enumerate(squares): 
        #print square_idxs 
        sc = eval_pts_as_rect(
                cvx_hull_pts[square_idxs[0]], 
                cvx_hull_pts[square_idxs[1]],
                cvx_hull_pts[square_idxs[2]],
                cvx_hull_pts[square_idxs[3]],
                )
        area = eval_area_pts(
                cvx_hull_pts[square_idxs[0]], 
                cvx_hull_pts[square_idxs[1]],
                cvx_hull_pts[square_idxs[2]],
                cvx_hull_pts[square_idxs[3]],
                             )
        res.append( (sc,area, square_idxs))
    
    res.sort()
    res = res[0:5]
    
    
    for sc,area,inds in res:
        #inds = res[i][1]
        pylab.plot( cvx_hull_pts[inds,0], cvx_hull_pts[inds,1], 'g-' )
        pylab.plot( [ cvx_hull_pts[inds[-1],0],cvx_hull_pts[inds[0],0] ], [ cvx_hull_pts[inds[-1],1],cvx_hull_pts[inds[0],1] ], 'g-' )
        
    # Lets pick the res with the biggest area:
    res.sort(key=lambda e:e[1], reverse=True)
    
    _,_,inds_max = res[0]
    
    
    pt1 = cvx_hull_pts[inds_max[0],:]
    pt2 = cvx_hull_pts[inds_max[1],:]
    pt3 = cvx_hull_pts[inds_max[2],:]
    pt4 = cvx_hull_pts[inds_max[3],:]
    return [pt1,pt2,pt3,pt4]

def get_rotation_angle( pts ):    
    pt1,pt2,pt3,pt4 = pts
    vA = pt2-pt1
    vB = (pt3-pt2)
    vC = (pt4-pt3)* -1
    vD = (pt1-pt4)* -1
    
    angle1 = np.arctan2( vA[1], vA[0]) * (360./(2*np.pi)) + 360
    angle2 = np.arctan2( vB[1], vB[0]) * (360./(2*np.pi)) + 360 - 90.
    angle3 = np.arctan2( vC[1], vC[0]) * (360./(2*np.pi)) + 360
    angle4 = np.arctan2( vD[1], vD[0]) * (360./(2*np.pi)) + 360 - 90
    
    angle1 = math.fmod(angle1, 180.)
    angle2 = math.fmod(angle2, 180.)
    angle3 = math.fmod(angle3, 180.)
    angle4 = math.fmod(angle4, 180.)
    
    print "%0.2f" % angle1, "%0.2f" % angle2, "%0.2f" % angle3, "%0.2f" % angle4
    median_angle = np.median([angle1,angle2,angle3,angle4] )
    median_angle = math.fmod(median_angle, 90.)
    return median_angle


def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
 
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2) 
 
def normpdf(x, mu, sigma):
    radius_max = 20. 
    u = (x-mu)/abs(sigma)
    
    y = (1/(np.sqrt(2*np.pi)*abs(sigma)))*np.exp(-u*u/2)
    
    y[ np.abs(x-mu) > radius_max  ] = 0.
    
    return y

def gaussian2d(sz, mu, sigma):
    x = normpdf(x=arange(sz[0]), mu=mu[0], sigma = sigma)
    y = normpdf(x=arange(sz[1]), mu=mu[1], sigma = sigma)
    t = np.outer(x,y)
    sc = 1/np.max(t)
    #sc=1.
    return t*sc

def build_simple_img2(sz, X0,X1,Y0,Y1, (d0, w0,m0,  d1, w1, m1,  d2, w2, m2,  d3, w3,m3) ):
    #im =zeros(sz)
    #im[X0:X1,Y0:Y1] = 1.0
    X = (X0+X1)/2.
    Y = (Y0+Y1)/2.
    
    #m = 5. * 1e3
    t0 = gaussian2d(sz, mu=(X1+d0, Y), sigma=w0) * m0 
    t1 = gaussian2d(sz, mu=(X0-d1, Y), sigma=w1) * m1 
    t2 = gaussian2d(sz, mu=(X, Y1+d2), sigma=w2) * m2 
    t3 = gaussian2d(sz, mu=(X, Y0-d3), sigma=w3) * m3 
    ts = [t0,t1,t2,t3]
    ts = [np.clip(t,-1,1,) for t in ts]
    
    tot = ts[0]+ts[1]+ts[2]+ts[3]
    print tot.min(), tot.max()
    assert tot.min()>=-1 and tot.max() <=1.
    tot = np.clip(tot,-1,1)
    
    im = zeros(sz)
    
    im[X0:X1,Y0:Y1] = 1.
    
    im += tot
    im = np.clip(im,0,1)
    
    #figure()
    #imshow(im)
    #show()
    
    return im


def build_simple_img(sz, x0,x1,y0,y1, p ):
    (d0, w0) = p
    X = (x0+x1)/2.
    Y = (y0+y1)/2.
    gX, gY = x1+np.fabs(d0), Y
 
    im = zeros(sz)
    im[x0:x1,y0:y1] = 1.
    im += gaussian2d(sz, mu=(gX,gY), sigma=w0) * 2
    im = np.clip(im,0,1)
    
    return im

def min_func_x(p, x0,x1,y0,y1, im_norm, direction):

    
    im = build_simple_img(im_norm.shape, x0,x1,y0,y1, p )
    diff = im_norm - im
    res = np.sum(diff**2)
    return res
       
    
reses = []    

def fit_piece(fname,fname_idx): 
    print fname 
    im_piece = imread(fname) 
    
    im_piece = (im_piece[:,:,0]*255.).astype('int') 
    
    #pt1,pt2,pt3,pt4 = get_corners_of_piece(im_piece)    
    
    # Derotate the piece:
    pts = get_corners_of_piece(im_piece)    
    median_angle = get_rotation_angle(pts)
    im_rot = rotate( im_piece, -median_angle)
    
    #figure()
    #imshow(im_rot,)
    #pylab.savefig('analysis/%03d_bb_rot.png'%fname_idx)
    
    # Find the corners:
    pts = get_corners_of_piece(im_rot)    
    centre = np.mean(np.array(pts), axis=0)
    corner_size = (30,30)
    f,axes = pylab.subplots(2,2)
    axes = [axes[0][0], axes[0][1],axes[1][0],axes[1][1]]
    
    #print x>centre[0], y>centre[1]
    pts = sorted(pts, key=lambda pt: (pt[0]>centre[0],pt[1]>centre[1]))
    for i,(x,y) in enumerate(pts):
        pylab.sca(axes[i])
        print x>centre[0], y>centre[1]
        
        xmin,xmax = x-corner_size[0], x+corner_size[0]
        ymin,ymax = y-corner_size[1], y+corner_size[1]
        
        xmin = np.clip(xmin, 0, im_rot.shape[0])
        ymin = np.clip(ymin, 0, im_rot.shape[1])
        xmax = np.clip(xmax, 0, im_rot.shape[0]-1)
        ymax = np.clip(ymax, 0, im_rot.shape[1]-1)
        
        corner = im_rot[ xmin:xmax, ymin:ymax ]
        
        imshow(corner.T)
        axes[i].plot([corner_size[0]],[corner_size[1]])
    
    pylab.savefig('analysis/%03d_bb_corners.png'%fname_idx)
    #pylab.show()
    
    
    figure()
    imshow(im_rot,)
    pylab.savefig('analysis/%03d_bb_rot.png'%fname_idx)
    
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
        
        
    min_func = partial(min_func_x, x0=x0,x1=x1,y0=y0,y1=y1, im_norm=im_norm, direction=Dir.Down)
    
    p = fmin(min_func, p0)
    im_opt = build_simple_img(im_rot.shape, x0,x1,y0,y1, (p) )
    
    
    f,axes = pylab.subplots(2)
    pylab.sca(axes[0])
    imshow(im_opt)
    pylab.sca(axes[1])
    imshow(im_rot)
    pylab.title("p[1]=%f"%p[1] )
    pylab.savefig('analysis/%03d_bb_fit1.png'%fname_idx)
    #pylab.show()
    reses.append( p )
    
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
pylab.xlabel("Strength (m)")
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