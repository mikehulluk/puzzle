

import numpy as np
from scipy.ndimage.interpolation import rotate
from matplotlib.pyplot import imshow
from scipy.spatial.qhull import ConvexHull
from itertools import combinations
import math
import pylab
from caching import joblib_memory


@joblib_memory.cache
def align_to_axes(im_piece):
    pts = get_corners_of_piece(im_piece)
    median_angle = get_rotation_angle(pts)
    im_rot = rotate(im_piece, -median_angle)
    return im_rot


def eval_pts_as_rect(pt0,pt1,pt2,pt3):

    
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


@joblib_memory.cache
def get_corners_of_piece(piece):
    plot = False
    
    
    if plot:
        imshow(piece.T) 
        
    points = np.transpose( np.nonzero(piece) ) 
    x = ConvexHull(points) 
    cvx_hull_pts =  points[x.vertices,0:2]  
    
    if plot:
        pylab.plot(cvx_hull_pts[:,0], cvx_hull_pts[:,1],'rx') 
 

    squares = combinations( range(cvx_hull_pts.shape[0]), 4 ) 
    
    res = []
    for i,square_idxs in enumerate(squares): 
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
    
    if plot:
        for sc,area,inds in res:
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








def extract_corners(im_rot):
    plot = False
    
    pts = get_corners_of_piece(im_rot)
    
    centre = np.mean(np.array(pts), axis=0)
    corner_size = (30,30)
    if plot:
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
        
        if plot:
            imshow(corner.T)
            axes[i].plot([corner_size[0]],[corner_size[1]])
    
    #pylab.savefig('analysis/%03d_bb_corners.png'%fname_idx)