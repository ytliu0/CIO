import math

from mod_functions import mod2pi, mod2pi_omgDf
from Eo_Vondrak_IAU2000A_spline import Eo_Vondrak_IAU2000A_spline
from Eo_Vondrak_longT import Eo_Vondrak_longT

def ERA_from_UT1(jd0_ut1, jd1_ut1):
    """
    Calculate ERA at UT1 jd_ut1 = jd0_ut1 + jd1_ut1 from the definition of UT1.
    Return ERA in radian in the range [-pi, pi).
    """
    D0 = math.floor(jd0_ut1 - 2451545 + jd1_ut1)
    fday = (jd0_ut1 - math.floor(jd0_ut1)) + (jd1_ut1 - math.floor(jd1_ut1))
    fday -= math.floor(fday)
    ERA = mod2pi_omgDf(0.01720217957524373, D0, 0) + fday*6.300387486754831 - 1.38822409435583
    return mod2pi(ERA)

def GAST_from_Eo(jd0_ut1, jd1_ut1, Eo):
    """
    Calculate GAST at UT1 jd_ut1 = jd0_ut1 + jd1_ut1 from ERA and Eo.
    Return ERA in radian in the range [-pi, pi).
    """
    D0 = math.floor(jd0_ut1 - 2451545 + jd1_ut1)
    fday = (jd0_ut1 - math.floor(jd0_ut1)) + (jd1_ut1 - math.floor(jd1_ut1))
    fday -= math.floor(fday)
    ERA = mod2pi_omgDf(0.01720217957524373, D0, 0) + fday*6.300387486754831 - 1.38822409435583
    return mod2pi(ERA - Eo)

def GAST_Vondrak_IAU2000A_spline(jd0_ut1, jd1_ut1, jd0_tt, jd1_tt):
    """
    Calculate GAST at UT1 jd_ut1 = jd0_ut1 + jd1_ut1 from ERA and Eo calculated by the spline formula at TT jd_tt = jd0_tt + jd1_tt.
    Return GAST in radian in the range [-pi, pi).
    """
    Eo = Eo_Vondrak_IAU2000A_spline(jd0_tt, jd1_tt)
    return GAST_from_Eo(jd0_ut1, jd1_ut1, Eo)

def GAST_Vondrak_longT(jd0_ut1, jd1_ut1, jd0_tt, jd1_tt):
    """
    Calculate GAST at UT1 jd_ut1 = jd0_ut1 + jd1_ut1 from ERA and Eo calculated by the long time fitting formula at TT jd_tt = jd0_tt + jd1_tt.
    Return GAST in radian in the range [-pi, pi).
    """
    Eo = Eo_Vondrak_longT(jd0_tt, jd1_tt)
    return GAST_from_Eo(jd0_ut1, jd1_ut1, Eo)
