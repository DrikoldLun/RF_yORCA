import numpy as np
from numpy import sin, cos, arcsin, arccos, arctan2, deg2rad, rad2deg, pi, mod
#from seispy import distaz
from scipy import interpolate
import math

def sind(deg):
    rad = np.radians(deg)
    return np.sin(rad)

def cosd(deg):
    rad = np.radians(deg)
    return np.cos(rad)

def tand(deg):
    rad = np.radians(deg)
    return np.tan(rad)

def cotd(deg):
    rad = np.radians(deg)
    return np.cos(rad) / np.sin(rad)

def asind(x):
    rad = np.arcsin(x)
    return np.degrees(rad)

def acosd(x):
    rad = np.arccos(x)
    return np.degrees(rad)

def atand(x):
    rad = np.arctan(x)
    return np.degrees(rad)

def km2deg(km):
    radius = 6371
    circum = 2*pi*radius
    conv = circum / 360
    deg = km / conv
    return deg

def deg2km(deg):
    radius = 6371
    circum = 2*pi*radius
    conv = circum / 360
    km = deg * conv
    return km

def rad2deg(rad):
    deg = rad*(360/(2*pi))
    return deg

def skm2sdeg(skm):
    sdeg = skm * deg2km(1)
    return sdeg

def sdeg2skm(sdeg):
    skm = sdeg / deg2km(1)
    return skm

def srad2skm(srad):
    sdeg = srad * ((2*pi)/360)
    return sdeg / deg2km(1)

def skm2srad(skm):
    sdeg = skm * deg2km(1)
    return rad2deg(sdeg)


def rot3D(bazi, inc):
    """
    :param bazi:
    :param inc:
    :return:
    M = [cos(inc)     -sin(inc)*sin(bazi)    -sin(inc)*cos(bazi);
        sin(inc)      cos(inc)*sin(bazi)     cos(inc)*cos(bazi);
        0              -cos(bazi)             sin(bazi)];
    """

    if isinstance(inc, float) or isinstance(inc, int):
        value31 = 0
    elif isinstance(inc, np.ndarray):
        value31 = np.repeat(0, len(inc))
    else:
        raise TypeError('Input args sould be in \'float\', \'int\', or \'numpy.ndarray\'')

    inc = inc / 180 * pi
    bazi = bazi / 180 * pi

    M = np.array([[np.cos(inc), -np.sin(inc)*np.sin(bazi), -np.sin(inc)*np.cos(bazi)],
                  [np.sin(inc), np.cos(inc)*np.sin(bazi), np.cos(inc)*np.cos(bazi)],
                  [value31, -np.cos(bazi), np.sin(bazi)]])
    return M


def rotateSeisZENtoLQT(Z, E, N, bazi, inc):
    M = rot3D(bazi, inc)
    ZEN = np.array([Z, E, N])
    LQT = np.dot(M, ZEN)
    return LQT[0, :], LQT[1, :], LQT[2, :]


def spherical2cartesian(lon, lat, dep):
    cola = 90. - lat
    r = 6371 - dep
    x = r * sind(cola) * cosd(lon)
    y = r * sind(cola) * sind(lon)
    z = r * cosd(cola)
    return x, y, z


def rotateSeisENtoTR(E, N, BAZ):
    angle = mod(BAZ+180, 360)
    R = N*cosd(angle) + E*sind(angle)
    T = E*cosd(angle) - N*sind(angle)
    return T, R

#def rssq(x):
#    return np.sqrt(np.sum(x**2)/len(x))

def rms(data):
    return np.sqrt(np.mean(np.power(data,2)))

def snr(data,fs,t_boundary,time_before,time_after):
    tb = t_boundary-time_before
    te = t_boundary+time_after
    nb = int(tb*fs-1)
    nc = int(t_boundary*fs-1)
    ne = int(te*fs-1)

    n_win = data[nb:nc]
    s_win = data[nc+1:ne+1]

    return max(np.abs(s_win))/rms(n_win)

'''
def snr(x, y):
    spow = rssq(x)**2
    npow = rssq(y)**2
    return 10 * np.log10(spow / npow)
'''

def latlon_from(lat1, lon1, azimuth, gcarc_dist):
    lat2 = asind ((sind (lat1) * cosd (gcarc_dist)) + (cosd (lat1) * sind (gcarc_dist) * cosd (azimuth)))
    if isinstance(gcarc_dist, np.ndarray):
        lon2 = np.zeros_like(lat2)
        for n in range(len(gcarc_dist)):
            if cosd (gcarc_dist[n]) >= (cosd (90 - lat1) * cosd (90 - lat2[n])):
                lon2[n] = lon1 + asind (sind (gcarc_dist[n]) * sind (azimuth) / cosd (lat2[n]))
            else:
                lon2[n] = lon1 + asind (sind (gcarc_dist[n]) * sind (azimuth) / cosd (lat2[n])) + 180
    elif isinstance(azimuth, np.ndarray):
        lon2 = np.zeros_like(lat2)
        for n in range(len(azimuth)):
            if cosd (gcarc_dist) >= (cosd (90 - lat1) * cosd (90 - lat2[n])):
                lon2[n] = lon1 + asind (sind (gcarc_dist) * sind (azimuth[n]) / cosd (lat2[n]))
            else:
                lon2[n] = lon1 + asind (sind (gcarc_dist) * sind (azimuth[n]) / cosd (lat2[n])) + 180
    else:
        if ( cosd (gcarc_dist) >= (cosd (90 - lat1) * cosd (90 - lat2))):
            lon2 = lon1 + asind (sind (gcarc_dist) * sind (azimuth) / cosd (lat2))
        else:
            lon2 = lon1 + asind (sind (gcarc_dist) * sind (azimuth) / cosd (lat2)) + 180
    return lat2, lon2

'''
def geoproject(lat_p, lon_p, lat1, lon1, lat2, lon2):
    azi = distaz(lat1, lon1, lat2, lon2).baz
    dis_center = distaz(lat1, lon1, lat_p, lon_p).delta
    azi_center = distaz(lat1, lon1, lat_p, lon_p).baz
    dis_along = atand(tand(dis_center))*cosd(azi-azi_center)
    (lat, lon) = latlon_from(lat1, lon1, azi, dis_along)
    return lat, lon
'''

def extrema(x, opt='max'):
    if opt == 'max':
        idx = np.intersect1d(np.where(np.diff(x) > 0)[0]+1, np.where(np.diff(x) < 0)[0])
    elif opt == 'min':
        idx = np.intersect1d(np.where(np.diff(x) < 0)[0]+1, np.where(np.diff(x) > 0)[0])
    else:
        raise ImportError('Wrong Options!!!')
    return idx

def slantstack(seis, timeaxis, rayp_range, tau_range, ref_dis, dis):
    EV_num = seis.shape[0]
    tmp = np.zeros([EV_num, tau_range.shape[0]])
    amp = np.zeros([rayp_range.shape[0], tau_range.shape[0]])
    for j in range(rayp_range.shape[0]):
        for i in range(EV_num):
            seis[i, :] = seis[i, :] / np.max(np.abs(seis[i, :]))
            tps = tau_range - rayp_range[j] * (dis[i] - ref_dis)
            tmp[i, :] = interpolate.interp1d(timeaxis, seis[i, :], fill_value='extrapolate')(tps)
        amp[j, :] = np.mean(tmp, axis=0)
    return amp

def distance(lat1, lon1, lat2, lon2):
    radius = 6371
    lat1, lon1, lat2, lon2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat*0.5)**2 + cos(lat1)*cos(lat2)*sin(dlon/2)**2
    dis = rad2deg(2*arcsin(np.sqrt(a)))

    y = sin(lon2-lon1)*cos(lat2)
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1)
    brng = rad2deg(arctan2(y,x))
    azi = (brng+360.0)%360.0

    return dis, azi

def latlon_from(lat1,lon1,azimuth,distance):
    #distance must be an array
    #gcarc_dist = km2deg(distance)
    gcarc_dist = distance
    lat1, lon1, azimuth, gcarc_dist = map(deg2rad, [lat1, lon1, azimuth, gcarc_dist])
    lat2 = arcsin((sin(lat1)*cos(gcarc_dist))+(cos(lat1)*sin(gcarc_dist)*cos(azimuth)))
    lon2 = np.zeros(len(gcarc_dist))
    for n in range(len(gcarc_dist)):
        if cos(gcarc_dist[n]) >= cos(pi/2-lat1)*cos(pi/2-lat2[n]):
            lon2[n] = lon1 + arcsin(sin(gcarc_dist[n])*sin(azimuth)/cos(lat2[n]))
        else:
            lon2[n] = lon1 + arcsin(sin(gcarc_dist[n])*sin(azimuth)/cos(lat2[n]))+pi
    return rad2deg(lat2), rad2deg(lon2)
