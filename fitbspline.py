from matplotlib.pyplot import imshow, show, plot, title, subplots
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from constants import Dir
from constants import EdgeType
from matplotlib.pylab import figure
import itertools
from numpy import zeros_like, hstack, zeros, CLIP, clip
import pylab
from scipy.ndimage.measurements import label
from scipy.signal.signaltools import convolve2d
#from scipy.optimize.optimize import fmin
from scipy.optimize import minimize
from skimage.morphology._skeletonize import medial_axis
from numpy.linalg.linalg import norm
from itertools import chain
from scipy.misc.pilutil import imresize
import skimage
from logilab.common import optparser
import skimage.draw
import functools
from caching import joblib_memory
from scipy.optimize._basinhopping import basinhopping



def evaluate_bspline(points, npoints):
    points = np.array(points)
    x = points[:,0]
    y = points[:,1]

    t = range(len(points))
    ipl_t = np.linspace(0.0, len(points) - 1, npoints)

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

    return (x_i, y_i)
   


def get_initial_ctrl_points((X0,X1,Y0,Y1), direction, edge_type, w):
    X = (X0+X1)/2.
    Y = (Y0+Y1)/2.
    
    print direction
    startpt, endpt = {
        Dir.Up:    ([X1,Y1],[X0,Y1]),
        Dir.Down:  ([X0,Y0],[X1,Y0]),
        Dir.Right: ([X1,Y0],[X1,Y1]),
        Dir.Left:  ([X0,Y1],[X0,Y0]),
        }[direction]
        
    if edge_type == EdgeType.Flattie:
        startpt = np.array(startpt)
        endpt = np.array(endpt)
        pts =  [startpt, startpt + 0.3 * (endpt-startpt), startpt +  0.7 * (endpt-startpt) , endpt]
        pts = [ pt.tolist() for pt in pts] 
        return pts
    w_out_mul = 1. if edge_type ==EdgeType.Outtie else -1.0     
    
    
    #points = [[0,0],[1,1],[2,2],[3,3]]
    if direction == Dir.Up:
        points = [
                  startpt,
                  [(X+startpt[0])/2., Y1],
                  [X+w/2, Y1],
                  [X+w, Y1+w*1.5*w_out_mul],
                  [X,   Y1+w*2.0*w_out_mul],
                  [X-w, Y1+w*1.5*w_out_mul],
                  [X-w/2, Y1],
                  [(X+endpt[0])/2., Y1],
                  endpt,
                  ]
    if direction == Dir.Down:
        points = [
                  startpt,
                  [(X+startpt[0])/2., Y0],
                  [X-w/2, Y0],
                  [X-w, Y0-w*1.5*w_out_mul],
                  [X,   Y0-w*2.0*w_out_mul],
                  [X+w, Y0-w*1.5*w_out_mul],
                  [X+w/2, Y0],
                  [(X+endpt[0])/2, Y0],
                  endpt,
                  ]
    if direction == Dir.Right:
        points = [
                  startpt,
                  [X1, (Y+startpt[1])/2.],
                  [X1, Y-w/2],
                  [X1+w*1.5*w_out_mul, Y-w],
                  [X1+w*2.0*w_out_mul, Y],
                  [X1+w*1.5*w_out_mul, Y+w],
                  [X1, Y+w/2],
                  [X1, (Y+endpt[1])/2.],
                  endpt,
                  ]
    if direction == Dir.Left:
        points = [
                  startpt,
                  [X0, (Y+startpt[1])/2.],
                  [X0, Y+w/2],
                  [X0-w*1.5*w_out_mul, Y+w],
                  [X0-w*2.0*w_out_mul, Y],
                  [X0-w*1.5*w_out_mul, Y-w],
                  [X0, Y-w/2],
                  [X0, (Y+endpt[1])/2.],
                  endpt,
                  ]
    
    return points




class PieceSplineTemplate(object):
    def __init__(self, (X0,X1,Y0,Y1), edge_types, im_norm):
        self.edge_types = edge_types
        self.im_norm = im_norm
        self.X0 = X0
        self.X1 = X1
        self.Y0 = Y0
        self.Y1 = Y1

        npts = {EdgeType.Flattie:4, EdgeType.Innie:9, EdgeType.Outtie:9}
        self.pts_per_direction = [ npts[edge_type] for edge_type in edge_types ]
        
        
        self.call_cnt = 0
    
    def get_intial_ctrl_pts(self):
        w0 = 10
        return [get_initial_ctrl_points((self.X0,self.X1,self.Y0,self.Y1), direction=direction, edge_type=edge_type, w=w0) for (direction, edge_type) in zip( Dir.directions, self.edge_types) ]
        
    @classmethod
    def control_points_to_curves(self, ctrl_pts, npoints=300):
        return [ evaluate_bspline(points, npoints=npoints) for points in ctrl_pts]


    def plot_ctrl_points(self, ctrl_pts, plot_image=True, ax=None, plot_spline=False):
        curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts)
        
        if ax is None:
            ax = pylab.gca()
        
        if plot_image:
            ax.imshow(self.im_norm.T)
        
        for i in range(4):
            ctrlx, ctrly =  zip(*ctrl_pts[i])
            ax.plot(ctrlx, ctrly, '-mo', markersize=2)
            ax.plot(curves[i][0], curves[i][1], '-g', markersize=2)
            
        if plot_spline:
            curves = self.control_points_to_curves(ctrl_pts=ctrl_pts)
            for c in curves:
                pts = np.array(c).T
                ax.plot( pts[:,0], pts[:,1], 'k-')    
                
            
            


    def ctrl_points_to_optvector(self, ctrl_pts):
        
        dir_right_start = ctrl_pts[0][0]
        dir_right_end =   ctrl_pts[0][-1]
        
        dir_up_start = ctrl_pts[1][0]
        dir_up_end =   ctrl_pts[1][-1]

        dir_left_start = ctrl_pts[2][0]
        dir_left_end = ctrl_pts[2][-1]   


        dir_down_start = ctrl_pts[3][0]
        dir_down_end =   ctrl_pts[3][-1]
        
        
        p0 = dir_right_end
        p1 = dir_up_end
        p2 = dir_left_end
        p3 = dir_down_end 
        
        
        do_check = True
        if do_check:

            
            # Double check values:
            tp0 = np.array(p0) - np.array(dir_up_start)
            tp1 = np.array(p1) - np.array(dir_left_start)
            tp2 = np.array(p2) - np.array(dir_down_start)
            tp3 = np.array(p3) - np.array(dir_right_start)
            
            assert np.linalg.norm(tp0) < 1.
            assert np.linalg.norm(tp1) < 1. 
            assert np.linalg.norm(tp2) < 1. 
            assert np.linalg.norm(tp3) < 1.  
            

        
        p = [p0,p1,p2,p3, ] + \
                ctrl_pts[0][1:-1] + \
                ctrl_pts[1][1:-1] + \
                ctrl_pts[2][1:-1] + \
                ctrl_pts[3][1:-1]
        if do_check:
            print p   
            expected_len = sum(self.pts_per_direction) - 4
            print "Expected len:", expected_len
            assert(len(p) == expected_len)
            
        p = list(itertools.chain(*p))
        return p

    def optvector_to_ctrl_points(self, p):
        pts_x, pts_y = p[::2], p[1::2]
        pts = zip(pts_x,pts_y)
        p0,p1,p2,p3 = pts[:4]
        
        dir_offset = 4
        right_offset = dir_offset  + 0
        up_offset = right_offset + self.pts_per_direction[0] -2
        
        left_offset = up_offset + self.pts_per_direction[1] -2
        
        down_offset = left_offset + self.pts_per_direction[2] -2
        
        vecs = [
                [p3] + pts[right_offset:right_offset+self.pts_per_direction[0] -2] + [p0], #Right
                [p0] + pts[up_offset:up_offset+self.pts_per_direction[1] -2] + [p1], #Up
                [p1] + pts[left_offset:left_offset+self.pts_per_direction[2] -2] + [p2], #Left
                [p2] + pts[down_offset:down_offset+self.pts_per_direction[3] -2] + [p3], #Down
                ]
        return vecs
    
    
    def image_distance(self, p, piece_idx, upsample):
        do_plot=False
        
        X = (self.X1 - self.X0) /2
        Y = (self.Y1 - self.Y0) /2
        ctrl_pts = self.optvector_to_ctrl_points(p)
        
        
        
        sz = self.im_norm.shape
        sz_upsample = (sz[0]*upsample, sz[1]*upsample)
        
        
        pts = []
        curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts, npoints=100 * upsample)
        for i in range(4):
            for x,y in zip( *curves[i] ):
                pts.append( (x,y) )
        
        xvals = np.linspace(0, sz[0], num=sz_upsample[0])
        yvals = np.linspace(0, sz[1], num=sz_upsample[1])
        
        indices = set()
        
        pts_upscale = []
        for x,y in pts:
            xbin = np.digitize([x],xvals)[0]
            ybin = np.digitize([y],yvals)[0]
            
            #im_new[xbin,ybin] = 1.
            indices.add((xbin,ybin))
        
            p = (xbin,ybin)
            if not pts_upscale or pts_upscale[-1] != p: 
                pts_upscale.append(p)
        
        poly_coords = np.array(pts_upscale).astype('float')
        
        [rr,cc] = skimage.draw.polygon( poly_coords[:,0] ,poly_coords[:,1], sz_upsample )
        
        op = zeros(sz_upsample)
        op[rr,cc] = 1.0
        
        #print op.shape
        
        
        # Downsampled image:
        im_new = imresize(op, sz,interp='bilinear').astype('float') / 255.


        im_diff = (im_new - self.im_norm) **2
        

        #print " >> ", im_diff[100,100]
        diff = np.sum(im_diff)
        
        print self.call_cnt, diff
        
        
        
        if self.call_cnt % 50 == 0:
            f,(ax1,ax2) = pylab.subplots(1, 2, squeeze=True)
            f.suptitle('Res: %0.2f'% diff )
            pylab.sca(ax1)
            ax1.set_title("diff"); pylab.imshow(im_diff.T, vmin=0,vmax=1.);pylab.colorbar(ticks=[0.,1.])
             
            self.plot_ctrl_points(ctrl_pts=ctrl_pts, plot_image=True, ax=ax2, plot_spline=True)
            
            
            pylab.savefig("analysis/fitting_%03d_upsample%02d_iteration_%05d.png"%(piece_idx, upsample, self.call_cnt) )
            pylab.close("all")
            #show()
        self.call_cnt += 1
        
        return diff
    
    
    
                    
        
        
    
         
#@joblib_memory.cache
def test_fits(im_norm, (X0,X1,Y0,Y1), edge_types, piece_idx):
    
    
    
    
    tmpl = PieceSplineTemplate((X0,X1,Y0,Y1), edge_types, im_norm)
    
    ctrl_pts_0 = tmpl.get_intial_ctrl_pts()
    
    curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts_0)
    
    figure()
    
    tmpl.plot_ctrl_points(ctrl_pts=ctrl_pts_0) 
    
    figure()
    p0 = tmpl.ctrl_points_to_optvector(ctrl_pts_0)
    p0_dash = tmpl.optvector_to_ctrl_points(p0)
    tmpl.plot_ctrl_points(ctrl_pts=p0_dash)
    
    pylab.close('all')
    
    
    
    # Rough fit:
    method = 'Powell'
    fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 1) 
    p0 = minimize(fit_func, p0, method='Nelder-Mead', tol=0.1).x
    #fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 1) 
    #p0 = minimize(fit_func, p0, method=method, tol=0.1).x    
    #fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 1) 
    #p0 = minimize(fit_func, p0, method='Anneal', tol=0.1).x    
    
    fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 2) 
    minimizer_kwargs = {"method": "Nelder-Mead"}
    p0 = basinhopping(fit_func, p0, minimizer_kwargs=minimizer_kwargs, niter=200).x
    
    #print ret
    #assert(0)
    
    #fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 3) 
    #p0 = minimize(fit_func, p0, method='Nelder-Mead', tol=0.5).x    
    ## Refine the fit:
    #fit_func = functools.partial( tmpl.image_distance, piece_idx=piece_idx, upsample = 5) 
    #p = minimize(fit_func, p0, method='Nelder-Mead', tol=0.1).x

    p = p0

    print p
    
    p_dash = tmpl.optvector_to_ctrl_points(p)
    
    
    # Save the final fit:
    figure()
    f,ax = subplots(1,1,squeeze=True) 
    tmpl.plot_ctrl_points(ctrl_pts=p_dash, ax = ax, plot_spline=True)
    pylab.savefig("analysis/fitting_%02d.png" % piece_idx)
    
    #pylab.show()
    
    pylab.close('all')
    
    #pylab.show()
    
   
    #plt.show()
    
    
    
    
    
    
    
    