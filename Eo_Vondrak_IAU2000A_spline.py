import numpy as np
from fundamental_arguments import fundamental_arguments
from Dpsi_cos_epsilonA import Dpsi_cos_epsilonA

def Eo_Vondrak_IAU2000A_spline(jd0, jd1):
    """
    Calculate the equation of origin Eo compatible with Vondrak et al/IAU2000A precession-nutation model at TT Julian date jd = jd0 + jd1 using a spline fitting formula. 

    JD must be in the range so that |T| = |(jd-2451545)/36525| <= 60.

    The interval [-60,60] are divided into sub-intervals [-60, -40], [-40, -20], [-20, -5], [-5, 5], [5, 20], [20, 40] and [40, 60]. Each sub-interval has a fitting formula. The inner boundary points -40, -20, -5, 5, 20 and 40 are the knots in the regression spline. The function Eo(T) and its first derivative Eo'(T) are continuous, but the second and higher derivatives of Eo are discontinuous at the knots.

    Following SOFA, Julian date is specified by two parts jd0 and jd1 in any way users may find convenient. For example, JD(TT)=2450123.7 could be expressed in any of these ways, among others.
                jd0             jd1
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)

    Accuracy of the spline formula:
    time period      estimated max error    rms error
    --------------------------------------------------
    -60 < T < -40         2.97 mas           0.596 mas
    -40 < T < -20         2.32 mas           0.467 mas
    -20 < T < -5          2.04 mas           0.408 mas
    -5 < T < 5            2.07 mas           0.396 mas
    5 < T < 20            2.37 mas           0.408 mas
    20 < T < 40           2.74 mas           0.491 mas
    40 < T < 60           3.63 mas           0.824 mas

    The fitting formula and code were developed by Yuk Tung Liu in June 2025.

    Eo is returned in radians.
    """
    jd_int = np.floor(jd0 + jd1)
    fday = (jd0 - np.floor(jd0)) + (jd1 - np.floor(jd1))
    fday -= np.floor(fday)
    T = ((jd_int - 2451545) + fday)/36525.0

    if abs(T) > 60:
        raise RuntimeError('Requested time is out of range.')
    
    F = fundamental_arguments(jd_int, fday)
    return Eop_Vondrak_IAU2000A_spline(T, F[4]) - Dpsi_cos_epsilonA(T, F)

def Eop_Vondrak_IAU2000A_spline(T, Omg):
    """
    Calculate Eo + Dpsi cos(epsilon_A) using the spline formula
    """
    T0, cpoly, csin, ph = set_Eop_coefficients(T)
    Tp = T - T0
    Eop = cpoly[0] + Tp*(cpoly[1] + Tp*(cpoly[2] + Tp*(cpoly[3] + Tp*cpoly[4])))
    ang = np.array([Omg, 2*Omg]) + ph
    Eop += sum(csin*np.sin(ang))
    return Eop

def set_Eop_coefficients(T):
    """
    Set the coefficients of the spline formula for Eo + Dpsi cos(epsilon_A)
    """
    csin = np.array([1.278687035263072e-08, 2.991955490317251e-10])
    ph = np.array([-3.141431849335106, -3.129942218845127])
    if abs(T) <= 5:
        T0 = 0
        cpoly = np.array([-7.029051838429728e-08, -0.02236036588274203, -6.744772398120004e-06, 3.326168239200108e-11, 1.260687080703534e-10])
    elif T >= -20 and T < -5:
        T0 = -12.5
        cpoly = np.array([0.2784536399120333, -0.02219271446784961, -6.627806792518711e-06, -6.300029177548349e-09, 1.280006947605919e-10])
    elif T >= -40 and T < -20:
        T0 = -30
        cpoly = np.array([0.6648424195980085, -0.02196935373260266, -6.052550663498622e-06, -1.580284734141516e-08, 1.422757542796084e-10])
    elif T < -40:
        T0 = -50
        cpoly = np.array([1.101957675015858, -0.02175077321573791, -4.76323677159108e-06, -2.7019621176594e-08, 1.327902426205912e-10])
    elif T > 5 and T < 20:
        T0 = 12.5
        cpoly = np.array([-0.28055535693713, -0.02252797753551305, -6.623865056562917e-06, 6.487274305837375e-09, 1.322936992702824e-10])
    elif T >= 20 and T < 40:
        T0 = 30
        cpoly = np.array([-0.6767759851418548, -0.02275091013887167, -6.027979162032721e-06, 1.644498248728261e-08, 1.50119330008187e-10])
    else:
        T0 = 50
        cpoly = np.array([-1.134049896561772, -0.02296751467215479, -4.684438743393382e-06, 2.807675328930574e-08, 1.332196864951549e-10])
    return T0, cpoly, csin, ph
