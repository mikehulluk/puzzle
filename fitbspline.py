from matplotlib.pyplot import imshow, show, plot, title
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

def evaluate_bspline(points):
    points = np.array(points)
    x = points[:,0]
    y = points[:,1]

    t = range(len(points))
    ipl_t = np.linspace(0.0, len(points) - 1, 300)

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
                  [X+w*2, Y1],
                  [X+w/2, Y1],
                  [X+w, Y1+w*2*w_out_mul],
                  [X-w, Y1+w*2*w_out_mul],
                  [X-w/2, Y1],
                  [X-w*2, Y1],
                  endpt,
                  ]
    if direction == Dir.Down:
        points = [
                  startpt,
                  [X-w*2, Y0],
                  [X-w/2, Y0],
                  [X-w, Y0-w*2*w_out_mul],
                  [X+w, Y0-w*2*w_out_mul],
                  [X+w/2, Y0],
                  [X+w*2, Y0],
                  endpt,
                  ]
    if direction == Dir.Right:
        points = [
                  startpt,
                  [X1, Y-w*2],
                  [X1, Y-w/2],
                  [X1+w*2*w_out_mul, Y-w],
                  [X1+w*2*w_out_mul, Y+w],
                  [X1, Y+w/2],
                  [X1, Y+w*2],
                  endpt,
                  ]
    if direction == Dir.Left:
        points = [
                  startpt,
                  [X0, Y+w*2],
                  [X0, Y+w/2],
                  [X0-w*2*w_out_mul, Y+w],
                  [X0-w*2*w_out_mul, Y-w],
                  [X0, Y-w/2],
                  [X0, Y-w*2],
                  endpt,
                  ]
    
    return points



def build_distsquare_im(sz, pt):
    xs = np.arange(0, sz[0])
    ys = np.arange(0, sz[1])
    
    dxs = xs - pt[0]
    dys = ys - pt[1]
    
    dXsq = np.tile(dxs**2, (sz[1],1)).T
    dYsq = np.tile(dys**2, (sz[0],1))
    
    #print dXsq.shape
    #print dYsq.shape
    
    #figure(); imshow(dXsq.T); pylab.colorbar()
    #figure(); imshow(dYsq.T); pylab.colorbar()
    #pylab.show()
    
    
    dist_sq = dXsq+ dYsq
    
    #figure(); imshow(dist_sq.T); pylab.colorbar();
    #pylab.plot([pt[0]],[pt[1]],'x')
    #pylab.show()

    return dist_sq

class PieceSplineTemplate(object):
    def __init__(self, (X0,X1,Y0,Y1), edge_types, im_norm):
        self.edge_types = edge_types
        self.im_norm = im_norm
        self.X0 = X0
        self.X1 = X1
        self.Y0 = Y0
        self.Y1 = Y1

        npts = {EdgeType.Flattie:4, EdgeType.Innie:8, EdgeType.Outtie:8}
        self.pts_per_direction = [ npts[edge_type] for edge_type in edge_types ]
        
        
        self.call_cnt = 0
    
    def get_intial_ctrl_pts(self):
        w0 = 10
        return [get_initial_ctrl_points((self.X0,self.X1,self.Y0,self.Y1), direction=direction, edge_type=edge_type, w=w0) for (direction, edge_type) in zip( Dir.directions, self.edge_types) ]
        
    @classmethod
    def control_points_to_curves(self, ctrl_pts_s ):
        return [ evaluate_bspline(points) for points in ctrl_pts_s]


    def plot_ctrl_points(self, ctrl_pts, plot_image=True):
        curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts)
        
        if plot_image:
            imshow(self.im_norm.T)
        
        for i in range(4):
            ctrlx, ctrly =  zip(*ctrl_pts[i])
            plt.plot(ctrlx, ctrly, '-mo')
            plt.plot(curves[i][0], curves[i][1], '-g')


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
    
    
    def image_distance(self, p):
        X = (self.X1 - self.X0) /2
        Y = (self.Y1 - self.Y0) /2
        ctrl_pts = self.optvector_to_ctrl_points(p)
        
        
        im = zeros(self.im_norm.shape)
        
        
        curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts)
        assert len(curves) == 4
        pts = []
        for i in range(4):
            for x,y in zip( *curves[i] ):
                pts.append( (x,y) )
        #print pts
        
        
        
        
        pts_internal = []
        for i in range(4):
            xys = list(zip( *curves[i] ))
            n_xys = len(xys)
            for j,(x,y) in enumerate(xys):
                if j != n_xys -1:
                    _pt1 = np.array( xys[j+1] )
                    _pt0 = np.array( xys[j] )
                    dir0 =  _pt1 - _pt0
                    dir0 /= norm(dir0)
                    dir_perp = np.array( [-dir0[1], dir0[0]] )
                    in_pt = _pt0 + dir_perp
                    pts_internal.append( (in_pt[0],in_pt[1]) )
        
        
        pts_internal =[ (int(x), int(y)) for x,y in pts_internal ]
        
        #print pts_internal
        
        #figure();
        #ptsX,ptsY = zip(*pts)
        #ptsINX,ptsINY = zip(*pts_internal)
        #figure();
        #plot(ptsX,ptsY, 'g')
        #plot(ptsINX,ptsINY, 'orange')
        
        
        #show()
        
        
        # Find the distance of each pixel to the puzzle piece:
        im_dist = None
        for pt in pts:
            im = build_distsquare_im(sz=self.im_norm.shape, pt=pt)
            if im_dist is None:
                im_dist = im
            else:
                im_dist = np.minimum(im_dist, im)
        im_dist = np.sqrt(im_dist)
        im_dist = 1./(1+im_dist)
        
        
        # Make a thick outline:
        im_outline = zeros(im_dist.shape)
        im_outline_inside = zeros(im_dist.shape)
        for (x,y) in pts:
            im_outline[ int(x), int(y)] = 1.
        for (x,y) in pts_internal:
            im_outline_inside[ int(x), int(y)] = 1.
        
        im_outline_thick = clip(im_outline + im_outline_inside, 0,1)
        
        #figure(); title("IM_OUTLINE");imshow(im_outline_thick.T); pylab.colorbar();
        
        
        labelled_im, num_comps = label(1-im_outline_thick, ) 
        #figure(); imshow(labelled_im.T); pylab.colorbar();
        
        internal_label = np.where(labelled_im==2, 1,0)
        #figure(); title("internal_label"); imshow(internal_label.T); pylab.colorbar();
        
        int2 = clip(im_outline_thick - im_outline, 0,1)
        int2 = np.where(int2>0.5, 1,0 )
        #figure(); title("int2"); imshow(int2.T); pylab.colorbar();
        
        
        
        im_new = zeros(im_dist.shape)
        im_new[np.where(internal_label)] = 1
        im_new[np.where(int2)] = 2
        
        #figure(); title("IM_internal"); imshow(im_new.T); pylab.colorbar();
        
        
        
        im_final = np.clip( im_new + im_dist, 0., 1.)
        #figure(); title("im_final"); imshow(im_final.T); pylab.colorbar();
        
            
        
        im_diff = (im_final - self.im_norm) **2
        
        #figure(); title("diff"); imshow(diff.T); pylab.colorbar();
        #show()
        #imshow(diff)
        #pylab.show()
        print " >> ", im_diff[100:100]
        diff = np.sum(im_diff)
        
        print self.call_cnt, diff
        
        
        
        if self.call_cnt % 10 == 0:
            figure(); title("diff"); imshow(im_diff.T, vmin=0,vmax=1.); cbar = pylab.colorbar(ticks=[0.,1.]) # ; cbar.vmin=0.0, vmax=1.0
            pylab.savefig("analysis/fitting_%04d.png"%self.call_cnt)
            pylab.close("all")
            #show()
        self.call_cnt += 1
        
        return diff
        
        
        
            
        
        
    
         

def test_fits(im_norm, (X0,X1,Y0,Y1), edge_types):

     
    
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
    
    tmpl.image_distance( p0 )
    
     
    
    #p = minimize(tmpl.image_distance, p0, method='Nelder-Mead') #, maxiter=1000,maxfun=1000)
    
    p =  np.array([ 107.16361691,  115.12643946,   51.98507228,  114.63592808,
         51.77483912,   39.95369715,  108.007557  ,   40.08943   ,
        108.44392364,   54.44782343,   99.78542787,   71.16922936,
        127.09163972,   69.52290702,  122.6773119 ,   88.97755043,
        108.06189382,   80.79077168,  107.13091062,   97.45627773,
         98.49957312,  111.20334618,   80.27707309,  113.27985432,
         86.76490985,  136.13309015,   70.23518408,  138.15139454,
         78.87377298,  114.78391742,   62.17470862,  110.72425835,
         52.32592034,   93.77664198,   52.30943143,   65.18547157,
         53.43415438,   40.71936311,   92.72767322,   40.00852822])
           
    print p
    
    p_dash = tmpl.optvector_to_ctrl_points(p)
    tmpl.plot_ctrl_points(ctrl_pts=p_dash)
    #tmpl.plot_ctrl_points(ctrl_pts=p0_dash, plot_image=False)
    
    print
    for (t,tt) in zip(p,p0 ):
        print t,tt
    #imshow()
    pylab.show()
    
    assert(0)
    #fig = plt.figure()
    #imshow(im_norm.T) 
    #for i in range(4):
    #    ctrlx, ctrly =  zip(*ctrl_pts_0[i])
    #    plt.plot(ctrlx, ctrly, '-mo')
    #    plt.plot(curves[i][0], curves[i][1], '-g')
    
    
    #for i, (direction, edge_type) in enumerate(zip( Dir.directions, edge_types)):
    #      
    #    points  = get_initial_ctrl_points((X0,X1,Y0,Y1), direction=direction, edge_type=edge_type, w=10)
    #    (x_ctrl, y_ctrl) = zip(*points)
    #    (x_i, y_i)  = evaluate_bspline(points)
    #
    #    
    #    #==============================================================================
    #    # Plot
    #    #==============================================================================
#
#        plt.plot(x_ctrl, y_ctrl, '-og')
#        plt.plot(x_i, y_i, 'r')
        #plt.title('Splined f(x(t), y(t))')
    
        #break
   
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
    
    
#if __name__ == '__main__':
#    test_fits()
#exit(0)    