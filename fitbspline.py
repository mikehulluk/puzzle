from matplotlib.pyplot import imshow
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from constants import Dir
from constants import EdgeType
from matplotlib.pylab import figure
import itertools


def evaluate_bspline(points):
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

    return (x_i, y_i)
   


def get_initial_ctrl_points((X0,X1,Y0,Y1), direction, edge_type, w):
    X = (X0+X1)/2.
    Y = (Y0+Y1)/2.
    
    print direction
    startpt, endpt = {
        Dir.Up:    ([X0,Y1],[X1,Y1]),
        Dir.Down:  ([X0,Y0],[X1,Y0]),
        Dir.Right: ([X1,Y0],[X1,Y1]),
        Dir.Left:  ([X0,Y0],[X0,Y1]),
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
                  [X-w*2, Y1],
                  [X-w/2, Y1],
                  [X-w, Y1+w*2*w_out_mul],
                  [X+w, Y1+w*2*w_out_mul],
                  [X+w/2, Y1],
                  [X+w*2, Y1],
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
                  [X0, Y-w*2],
                  [X0, Y-w/2],
                  [X0-w*2*w_out_mul, Y-w],
                  [X0-w*2*w_out_mul, Y+w],
                  [X0, Y+w/2],
                  [X0, Y+w*2],
                  endpt,
                  ]
    
    return points



class PieceSplineTemplate(object):
    def __init__(self, (X0,X1,Y0,Y1), edge_types):
        self.edge_types = edge_types
        self.X0 = X0
        self.X1 = X1
        self.Y0 = Y0
        self.Y1 = Y1

        npts = {EdgeType.Flattie:4, EdgeType.Innie:8, EdgeType.Outtie:8}
        self.pts_per_direction = [ npts[edge_type] for edge_type in edge_types ]
    
    def get_intial_ctrl_pts(self):
        w0 = 10
        return [get_initial_ctrl_points((self.X0,self.X1,self.Y0,self.Y1), direction=direction, edge_type=edge_type, w=w0) for (direction, edge_type) in zip( Dir.directions, self.edge_types) ]
        
    @classmethod
    def control_points_to_curves(self, ctrl_pts_s ):
        return [ evaluate_bspline(points) for points in ctrl_pts_s]


    def plot_ctrl_points(self, ctrl_pts):
        curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts)
        fig = plt.figure()
        
        for i in range(4):
            ctrlx, ctrly =  zip(*ctrl_pts[i])
            plt.plot(ctrlx, ctrly, '-mo')
            plt.plot(curves[i][0], curves[i][1], '-g')


    def ctrl_points_to_optvector(self, ctrl_pts):
        #directions = [Right, Left, Up, Down]
        
        dir_right_start = ctrl_pts[0][0]
        dir_right_end =   ctrl_pts[0][-1]
        
        dir_down_start = ctrl_pts[3][0]
        dir_down_end =   ctrl_pts[3][-1]
        
        dir_left_start = ctrl_pts[1][0]
        dir_left_end =   ctrl_pts[1][-1]
        
        dir_up_start = ctrl_pts[2][0]
        dir_up_end =   ctrl_pts[2][-1]
        
        
        print "dir_right_start", dir_right_start
        print "dir_right_end", dir_right_end  
        
        print "dir_down_start", dir_down_start  
        print "dir_down_end",dir_down_end   
        
        print "dir_left_start",dir_left_start  
        print "dir_left_end",dir_left_end  
        
        print "dir_up_start",dir_up_start  
        print "dir_up_end",dir_up_end
        print 

        p0 = dir_up_end
        p1 = dir_down_end
        p2 = dir_down_start
        p3 = dir_up_start
        
        # Double check values:
        tp0 = np.array(p0) - np.array(dir_right_end)
        tp1 = np.array(p1) - np.array(dir_right_start)
        tp2 = np.array(p2) - np.array(dir_left_start)
        tp3 = np.array(p3) - np.array(dir_left_end)
        
        assert np.linalg.norm(tp0) < 1.
        assert np.linalg.norm(tp1) < 1. 
        assert np.linalg.norm(tp2) < 1. 
        assert np.linalg.norm(tp3) < 1.  
        
        print 'c0', len(ctrl_pts[0]),  ctrl_pts[0]
        print 'c1', len(ctrl_pts[1]),  ctrl_pts[1]
        print 'c2', len(ctrl_pts[2]),  ctrl_pts[2]
        print 'c3', len(ctrl_pts[3]),  ctrl_pts[3]
        
        #print 'c0', len(ctrl_pts[0][1:-1]),  ctrl_pts[0][1:-1]
        #print 'c1', len(ctrl_pts[1][1:-1]),  ctrl_pts[1][1:-1]
        #print 'c2', len(ctrl_pts[2][1:-1]),  ctrl_pts[2][1:-1]
        #print 'c3', len(ctrl_pts[3][1:-1]),  ctrl_pts[3][1:-1]
        
        p = [p0,p1,p2,p3, ] + \
                ctrl_pts[0][1:-1] + \
                ctrl_pts[1][1:-1] + \
                ctrl_pts[2][1:-1] + \
                ctrl_pts[3][1:-1]
        print p   
        expected_len = sum(self.pts_per_direction) - 4
        print "Expected len:", expected_len
        assert(len(p) == expected_len)
        
        p = list(itertools.chain(*p))
        return p

    def optvector_to_ctrl_points(self, p):
        #directions = [Right, Left, Up, Down]
        print "Undoing"
        print p
        pts_x, pts_y = p[::2], p[1::2]
        pts = zip(pts_x,pts_y)
        print  len(pts), pts
        assert len(pts) ==  sum(self.pts_per_direction) - 4
        p0,p1,p2,p3 = pts[:4]
        
        dir_offset = 4
        right_offset = dir_offset  + 0
        left_offset = right_offset + self.pts_per_direction[0] -2
        up_offset = left_offset + self.pts_per_direction[1] -2
        down_offset = up_offset + self.pts_per_direction[2] -2
        
        vecs = [
                [p1] + pts[right_offset:right_offset+self.pts_per_direction[0] -2] + [p0], #Right
                [p2] + pts[left_offset:left_offset+self.pts_per_direction[1] -2] + [p3], #Left
                [p3] + pts[up_offset:up_offset+self.pts_per_direction[2] -2] + [p0], #Up
                [p2] + pts[down_offset:down_offset+self.pts_per_direction[3] -2] + [p1], #Down
                ]
        print 
        
        print vecs[0]
        print vecs[1]
        print vecs[2]
        print vecs[3]
        return vecs
        
         

def test_fits(im_norm, (X0,X1,Y0,Y1), edge_types):

     
    
    tmpl = PieceSplineTemplate((X0,X1,Y0,Y1), edge_types)
    
    ctrl_pts_0 = tmpl.get_intial_ctrl_pts()
    
    curves = PieceSplineTemplate.control_points_to_curves(ctrl_pts_0)
    
    figure()
    imshow(im_norm.T)
    tmpl.plot_ctrl_points(ctrl_pts=ctrl_pts_0) 
    
    figure()
    p0 = tmpl.ctrl_points_to_optvector(ctrl_pts_0)
    p0_dash = tmpl.optvector_to_ctrl_points(p0)
    tmpl.plot_ctrl_points(ctrl_pts=p0_dash)
    
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