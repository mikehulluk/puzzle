'''
Created on 2 Apr 2015

@author: mjh2
'''
import glob
from matplotlib.image import imread
from matplotlib.pyplot import imshow, hist, figure, plot
import pylab
from scipy.ndimage.interpolation import rotate

from scipy.optimize import minimize

import numpy as np
from numpy import NaN, pi, zeros_like
import functools
from scipy.spatial._plotutils import convex_hull_plot_2d
from scipy.spatial import ConvexHull
from matplotlib.figure import Figure



import sklearn
import sklearn.cluster



import scipy


# def gen_rect_image(shape, (rx0,ry0,rx1,ry1) ):
#     im = np.zeros(shape)
#     im[rx0:rx1, ry0:ry1] = 1.0
#     return im
#     
# 
# def test_img(im, (rot, rx0,ry0,rx1,ry1) ):
#     
#     rx0,ry0,rx1,ry1 = int(rx0),int(ry0),int(rx1),int(ry1)
#     
#     if rx1-rx0 < 10:
#         return NaN
#     if ry1-ry0 < 10:
#         return NaN
#     
#     print (rot, rx0,ry0,rx1,ry1)
#     piece_rot = rotate(im, rot)
#     
#     test_rect = gen_rect_image(piece_rot.shape, (rx0,ry0,rx1,ry1) )
#     
#     #pylab.imshow(test_rect) # , cmap, norm, aspect, interpolation, alpha, vmin, vmax, origin, extent, shape, filternorm, filterrad, imlim, resample, url, hold)
#     
#     #pylab.show()
#     
#     
#     return  - np.sum( np.fabs(piece_rot - test_rect ) )
    
from scipy.stats.distributions import  vonmises
    


def fit_piece(fname):
    
    
    print fname
    
    # Load and normalise:
    piece = (imread(fname) * 255).astype('int')
    piece = piece[:,:,0]
    piece = piece - np.min(piece)
    piece = piece / np.max(piece)
    print piece.shape
    piece_orig = piece[:]
    
    hist( piece.flatten() )
    figure()
    imshow(piece.T)
    
    
    
    # Find the non-zero indices, and build a convex hull:
    points = np.transpose( np.nonzero(piece) )
    centroid = np.mean(points, axis=0)
    
    hull = ConvexHull(points)
    
    
    
    
    import matplotlib.pyplot as plt
    plt.plot(points[:,0], points[:,1], 'o')
    plt.plot([centroid[0]], [centroid[1]] , 'ro')
    
    for simplex in hull.simplices:
        plt.plot(points[simplex,0], points[simplex,1], 'rx-')
    
    
    
    
    
    
    hull_vertices = [ points[vtx] for vtx in hull.vertices ]
    
    #print hull_vertices
    ##Merge anything within 30 pixels:
    #fitter = sklearn.cluster.DBSCAN(eps=30, min_samples=2, metric='euclidean', algorithm='auto', leaf_size=30, p=None, random_state=None)
    #X = np.array( hull_vertices )
    #merged_pts = fitter.fit(X)
    #print merged_pts.__dict__
    #print merged_pts.labels_
    #labels = set( merged_pts.labels_ )
    #pts_per_label = [ ]
    #for pt in merged_pts:
    #    plt.plot( [pt[0]],[pt[1]], 'rx-')
    
    
    
    
        
    angles = []
    horizontal = (1,0)
    for pt  in hull_vertices:
        #pt = points[vertex]
        plt.plot([pt[0]],[pt[1]],'gx')
        a = horizontal
        b = pt - centroid
        angle = np.arctan2(b[1],b[0])
        
        r = 50.
        x = r * np.cos(angle)
        y = r    * np.sin(angle)
        plt.plot([centroid[0], centroid[0]+x], [centroid[1], centroid[1]+y], 'rx-')

        angles.append(angle)
        
    
    print angles
    
    # Now, lets scale each angle, by the size of its neighbours:
    # angles run between -pi and +pi
    new_angles = []
    for i,angle in enumerate(angles):
        sf = 0
        for j,angle_n in enumerate(angles):
            if i==j:
                continue
            angle_dist = (angle - angle_n)
            if angle_dist > np.pi:
                angle_dist -=np.pi
            angle_dist =np.fabs(angle_dist)
            
            sf += 1/((angle_dist)**1.8)
            
            
        new_angles.append( (angle, 1/sf) )
        
    figure()
    hist(angles, bins=200)
    
    figure()
    data = np.array(new_angles)
    print data.shape
    plot(data[:,0], data[:,1], 'x')
    
    #pylab.close('all')
    figure()
    # Angles are hopefully within a factor of  pi/4
    ang_idx = np.linspace(-pi, pi*2, 1000)
    ang_sum = zeros_like(ang_idx)
    angles = new_angles
    new_angles = []
    for angle,val in angles:
        if angle<0:
            angle +=  pi
        for i in range(4):
            if angle > pi/2.:
                angle -= pi/2
        print angle, val
        
        plot([angle],[val], 'x')
        #ang_vals = scipy.stats.norm(ang_idx-angle)# * val
        #pylab.plot(ang_idx,ang_vals)
        new_angles.append((angle,val) )
        plt_pts = scipy.stats.norm.pdf(ang_idx, loc=angle, scale=0.1) * val
        #plt_pts = scipy.stats.vonmises.pdf(ang_idx, loc=angle, scale=0.1 ) * val
         
        plot(ang_idx, plt_pts)
        ang_sum += plt_pts
        
    plot(ang_idx, ang_sum)    
    
    max_angle = np.argmax(ang_sum)
    
    angle = ang_idx[max_angle]
    
    
    im_unrot = rotate(piece,-angle * (360/ (2*pi)) +45) #, axes, reshape, output, order, mode, cval, prefilter)
    
    pylab.close('all')
    
    figure()
    imshow(im_unrot)
    figure()
    imshow(piece_orig)
    
    
    
    
    
    
    pylab.show()
    
    

    
    
    
    
    
    #func = functools.partial( test_img, piece)
    #x0 = (0, 0,0, piece.shape[0], piece.shape[1] )
    #res = minimize(func, x0,) # args, method, jac, hess, hessp, bounds, constraints, tol, callback, options)
    
    #print res.x
    
    #imshow(piece)
    #pylab.show()
    
    
    
    
def main():    
    piece_files = sorted( glob.glob("build/piece*single_rot.png") )[11:]
    for piece_file in piece_files:
        print piece_file
        fit_piece(piece_file)
        
        #break
    
    
main()
pylab.show()