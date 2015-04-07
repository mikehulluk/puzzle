import glob 
from matplotlib.image import imread 
from matplotlib.pyplot import imshow, plot , figure
import pylab 
from scipy.ndimage.interpolation import rotate, shift 

from scipy.optimize import minimize 

import numpy as np 
from numpy import NaN, zeros 
import functools 
#from sympy.geometry.util import convex_hull 
from scipy.spatial.qhull import ConvexHull 
from itertools import combinations 
import os
import math

piece_files = sorted( glob.glob("build/piece*single_rot.png") ) 



    

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

def fit_piece(fname,fname_idx): 
    print fname 
    piece = imread(fname) 
    
    
    piece = (piece[:,:,0]*255.).astype('int') 
    #print piece.shape 
    
    #Normalise: 
    #piece = piece - np.min(piece) 
    #piece = piece / np.max(piece) 
    
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
    pylab.plot( cvx_hull_pts[inds_max,0], cvx_hull_pts[inds_max,1], 'm-', lw=3 )
    pylab.plot( [ cvx_hull_pts[inds_max[-1],0],cvx_hull_pts[inds_max[0],0] ], [ cvx_hull_pts[inds_max[-1],1],cvx_hull_pts[inds_max[0],1] ], 'm-', lw=3 )
        
        
    
    pylab.savefig('analysis/%03d_bb.png'%fname_idx)
    
    pt1 = cvx_hull_pts[inds_max[0],:]
    pt2 = cvx_hull_pts[inds_max[1],:]
    pt3 = cvx_hull_pts[inds_max[2],:]
    pt4 = cvx_hull_pts[inds_max[3],:]
    
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
    #print 'Combinations:', len(list(squares)) 
    
    
    #n_cvx_hull_pts = cvx_hull_pts.shape[0] 
    #angle_matrix = zeros((n_cvx_hull_pts,n_cvx_hull_pts)) 
    
    im_rot = rotate( piece, -median_angle)
    
    figure()
    imshow(im_rot,)
    pylab.savefig('analysis/%03d_bb_rot.png'%fname_idx)
    
    
    pylab.close('all')
    #pylab.show() 
    
    
    
    
if not os.path.exists("analysis"):
    os.makedirs("analysis")

for i,piece_file in enumerate(piece_files): 
    print piece_file 
    fit_piece(piece_file,i) 
    
    #break