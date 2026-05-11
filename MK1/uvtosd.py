from math import *

def uvTosd(u, v, uk, vk):
    # Convert u, v to the oblique aspect
    dv = vk - v

    # Latitude
    sarg = sin(u) * sin(uk) + cos(u) * cos(dv)*cos(uk)
    s = asin(sarg)
    
    #Longitude
    num = sin(dv) * cos(u)
    denom = -sin(u) * cos(uk) + cos(u) * cos(dv) * sin(uk)

    d = atan2(num, denom)

    return s, d