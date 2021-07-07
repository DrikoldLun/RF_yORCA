import numpy as np
from numpy import sin, cos
import obspy
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

def vectorangle(a,b):
    cosangle = np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))
    angle = np.arccos(cosangle)
    return angle

def faultvector(strike,dip,rake):
    strike, dip, rake = map(np.deg2rad,[strike,dip,rake])
    vector = {}
    vector['strike'] = np.array([sin(strike),cos(strike),0])
    vector['dip'] = np.array([cos(dip)*cos(strike),-cos(dip)*sin(strike),-sin(dip)])
    vector['rake'] = rake_vector = np.array([cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike),cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike),sin(rake)*sin(dip)])
    # normal vector should be the cross product of dip vector and strike vector
    vector['normal'] = np.cross(vector['dip'],vector['strike'])
    #vector['normal'] = np.array([sin(dip)*cos(strike),-sin(dip)*sin(strike),cos(dip)])
    vector['normal'] = vector['normal'] / np.linalg.norm(vector['normal'])
    #vector['normal'] = np.array([sin(dip)*cos(strike),-sin(dip)*sin(strike),cos(dip)])
    # p_vector = normal_vector - rake_vector
    vector['p'] = vector['normal'] - vector['rake']
    #vector['p'] = np.array([sin(dip)*cos(strike) - cos(rake)*sin(strike) + sin(rake)*cos(dip)*cos(strike),-sin(dip)*sin(strike) - cos(rake)*cos(strike) - sin(rake)*cos(dip)*sin(strike),cos(dip) - sin(rake)*sin(dip)])
    vector['p'] = vector['p'] / np.linalg.norm(vector['p'])
    #p_vector = np.array([sin(dip)*cos(strike) - cos(rake)*sin(strike) + sin(rake)*cos(dip)*cos(strike),-sin(dip)*sin(strike) - cos(rake)*cos(strike) - sin(rake)*cos(dip)*sin(strike),cos(dip) - sin(rake)*sin(dip)])
    # t_vector = normal_vector + rake_vector
    vector['t'] = vector['normal'] + vector['rake']
    #vector['t'] = np.array([sin(dip)*cos(strike) + cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike),-sin(dip)*sin(strike) + cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike),cos(dip) + sin(rake)*sin(dip)])
    vector['t'] = vector['t'] / np.linalg.norm(vector['t'])
    vector['null'] = np.cross(vector['rake'],vector['normal'])
    vector['null'] = vector['null'] / np.linalg.norm(vector['null'])
    #t_vector = np.array([sin(dip)*cos(strike) + cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike),-sin(dip)*sin(strike) + cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike),cos(dip) + sin(rake)*sin(dip)])
    for index in ['null','t','p']:
        if vector[index][2] > 0:
            vector[index] = -1 * vector[index]
    return vector

def thetaphi(rayvec,vector):
    theta = vectorangle(vector['normal'],rayvec)
    phi2 = vectorangle(vector['rake'],rayvec)
    phi = np.arccos(cos(phi2)/sin(theta))
    return theta, phi

# x1 - E, x2 - N, x3 - Z
def evtvec(azi,srcdepth,distance_in_degree,phase='S'):
    arrivals = model.get_travel_times(source_depth_in_km=srcdepth,
                                      distance_in_degree=distance_in_degree,
                                      phase_list=[phase])
    takeoff_angle = np.deg2rad(arrivals[0].takeoff_angle)
    azi = np.deg2rad(azi)
    rayvec = np.array([sin(azi)*sin(takeoff_angle), cos(azi)*sin(takeoff_angle), -cos(takeoff_angle)]) #? upgoing or downgoing
    rayvec = rayvec / np.linalg.norm(rayvec)
    return rayvec

def radiationfactor(azi,srcdepth,distance_in_degree,strike,dip,rake,phase='P',isnorm=False):
    vector = faultvector(strike,dip,rake)
    rayvec = evtvec(azi,srcdepth,distance_in_degree,phase=phase)
    theta, phi = thetaphi(rayvec,vector)
    #print(np.rad2deg(np.array([theta, phi])))
    if phase == 'P':
        F = sin(2*theta)*cos(phi)
    elif phase == 'S':
        if not isnorm:
            F = cos(2*theta)*cos(phi)
        else:
            F = (cos(2*theta)*cos(phi))**2 + (cos(theta)*sin(phi))**2
    return F

if __name__ == '__main__':
    #strike, dip, rake = 347, 85, -177
    #strike1, dip1, rake1 = 257, 87, -5
    #strike, dip, rake = 100, 15, -50
    #strike1, dip1, rake1 = 239, 79, -99
    strike, dip, rake = 32, 43, 103
    strike1, dip1, rake1 = 194, 48, 78
    #strike, dip, rake = 124, 66, -160
    #strike1, dip1, rake1 = 26, 72, -25
    #strike, dip, rake = 268, 16, 61
    #strike1, dip1, rake1 = 118, 76, 98
    '''
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    vectorp = faultvector(strike, dip, rake)['p']
    vectorn = faultvector(strike, dip, rake)['normal']
    vectort = faultvector(strike, dip, rake)['t']
    #vector1 = faultvector(strike1, dip1, rake1)['p']
    #print(np.dot(vector,vector1))
    ax.plot([0,vectorp[0]],[0,vectorp[1]],[0,vectorp[2]],label='p')
    ax.plot([0,vectorn[0]],[0,vectorn[1]],[0,vectorn[2]],label='normal')
    ax.plot([0,vectort[0]],[0,vectort[1]],[0,vectort[2]],label='t')
    ax.set(xlabel='E',ylabel='N',zlabel='Z')
    ax.legend()
    plt.axis('equal')
    plt.show()
    '''
    #print(radiationfactor(50,30,50,strike,dip,rake,phase='P'))
    #print(radiationfactor(50,30,50,strike1,dip1,rake1,phase='P'))
