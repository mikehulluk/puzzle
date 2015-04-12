from matplotlib.pyplot import imshow
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from constants import Dir
from fit_puzzle2 import EdgeType


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
   


def get_ctrl_points((X0,X1,Y0,Y1), direction, edge_type, w):
    X = (X0+X1)/2
    Y = (Y0+Y1)/2
    points = [
              [X0,Y1],
              [X-w*2, Y1],
              [X-w/2, Y1],
              
              [X-w, Y1+w*2],
              [X+w, Y1+w*2],
              
              [X+w/2, Y1],
              [X+w*2, Y1],
              [X1,Y1],
              ]
    return points


def test_fits(im_norm, (X0,X1,Y0,Y1), edge_types):

    
  
    
    w = 10
    points  = get_ctrl_points((X0,X1,Y0,Y1), direction=Dir.Up, edge_type=EdgeType.Outtie, w=10)
    (x_ctrl, y_ctrl) = zip(*points)
    (x_i, y_i)  = evaluate_bspline(points)

    
    #==============================================================================
    # Plot
    #==============================================================================
    
    fig = plt.figure()
    
    imshow(im_norm.T)
    plt.plot(x_ctrl, y_ctrl, '-og')
    plt.plot(x_i, y_i, 'r')
    plt.title('Splined f(x(t), y(t))')
    
    
    
    #ax = fig.add_subplot(234)
    #for i in range(7):
    #    vec = np.zeros(11)
    #    vec[i] = 1.0
    #    x_list = list(x_tup)
    #    x_list[1] = vec.tolist()
    #    x_i = si.splev(ipl_t, x_list)
    #    plt.plot(ipl_t, x_i)
    #plt.xlim([0.0, max(t)])
    #plt.title('Basis splines')
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