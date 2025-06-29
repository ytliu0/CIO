import numpy as np
import math
from s_Vondrak_longT import calc_sA_Vondrak_fit
from mod_functions import mod2pi
from Eo_Vondrak_IAU2000A_spline import Eo_Vondrak_IAU2000A_spline

def Eo_Vondrak_longT(jd0, jd1):
    """
    Calculates the equation of origin Eo compatible with Vondrak et al precession model at TT Julian date jd = jd0 + jd1 using a fitting formula valid for long period from J2000. Eo is returned in radians.

    JD must be in the range so that |T| = |(jd-2451545)/36525| <= 2000

    Following SOFA, Julian date is specified by two parts jd0 and jd1 in any way
    users may find convenient. For example, JD(TT)=2450123.7 could be expressed
    in any of these ways, amomg others.

                jd0             jd1
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)

    The fitting formulas include IAU 2000A nutation for |T| <= 59.8, but ignore
    nutation for |T| >= 60. For 59.8 < |T| < 60, a weighted average of the two
    sets of fitting formulas is calculated to ensure continuity.

    Accuracy of the fitting formula:
    For |T| <=59.8, the accuracy is the same as Eo_Vondrak_IAU2000A_spline():
    time period      estimated max error    rms error
    --------------------------------------------------
    -59.8 < T < -40       2.97 mas           0.596 mas
    -40 < T < -20         2.32 mas           0.467 mas
    -20 < T < -5          2.04 mas           0.408 mas
    -5 < T < 5            2.07 mas           0.396 mas
    5 < T < 20            2.37 mas           0.408 mas
    20 < T < 40           2.74 mas           0.491 mas
    40 < T < 59.8         3.63 mas           0.824 mas

    For T > 59.8, the estimated maximum error is 51" and the rms error is 18".
    For T < -59.8, the estimated maximum error is 69" and the rms error is 22".

    The fitting formulas and code were developed by Yuk Tung Liu in June 2025.

    Return Eo in radians
    """
    T = ((jd0 - 2451545) + jd1)/36525
    if abs(T) > 2000:
        raise RuntimeError('Request time is out of range')
    
    if abs(T) <= 59.8: return Eo_Vondrak_IAU2000A_spline(jd0, jd1)
    if abs(T) >= 60: return Eo_Vondrak_from_s(T)

    # T is near a boundary. Calculate s by taking a weighted average of two fitting formulas
    r = 0.2; Tb = 59.9;
    x = (T + Tb)/r if T < 0 else (T - Tb)/r
    w = np.sin(0.5*np.pi*(x + 0.5))**2
    if T < 0:
        Eo1 = Eo_Vondrak_IAU2000A_spline(jd0, jd1)
        Eo2 = Eo_Vondrak_from_s(T)
    else:
        Eo1 = Eo_Vondrak_from_s(T)
        Eo2 = Eo_Vondrak_IAU2000A_spline(jd0, jd1)
    return w*Eo1 + (1-w)*Eo2

def Eo_Vondrak_from_s(T):
    """
    Calculate precession contribution to Eo at TT jd = jd0+jd1
    Vondrak et al's precession model is used
    """
    s = calc_sA_Vondrak_fit(T)
    pb = PB_Vondrak(T)
    X = pb[2][0]; Y = pb[2][1]; a = 1.0/(1.0 + pb[2][2])
    RST = [1-a*X*X, -a*X*Y, -X]
    p = pb[0][0]*RST[0] + pb[0][1]*RST[1] + pb[0][2]*RST[2]
    q = pb[1][0]*RST[0] + pb[1][1]*RST[1] + pb[1][2]*RST[2]
    return mod2pi(s - math.atan2(q, p))

def PB_Vondrak(T):
    """
    Calculate the PB matrix at TT Julian century T.
    Precession matrix calculated according to J. Vondrak, N. Capitaine, P. Wallace, A&A 534, A22 (2011), DOI: 10.1051/0004-6361/201117274. P is valid for |T| < 2000. The periodic coefficients are taken from Tables 4 and 6 of the paper. The matrix is computed using Eq. (20) in the paper. 
    """
    omega = np.array([0.01559490024120026, 0.0244719973015758, 0.02151775790129995, 0.01169573974755144, 0.02602271819084525, 0.01674533688817117, 0.0397997422384214, 0.02291460724719032, 0.03095165175950535, 0.01427996660722633, 0.03680403764749055, 0.008807750966790847, 0.02007407446383254, 0.0489420883874403, 0.03110487775831478, 0.01994662002279234, 0.04609144151393476, 0.01282282715750936])
    cPsiA = np.array([-0.1076593062579846, 0.05932495062847037, -0.007703729840835942, 0.01203357586861691, 0.000728786082003343, -6.609012098588148e-05, 0.001888045891520004, 0.009848668946298234, 0.001763501537747769, -0.004347554865592219, -0.004494201976897112, 0.000179723665294558, -0.002897646374457124, 0.0003213481408001133, 0, 0, 0, 0])
    sPsiA = np.array([-0.01572365411244583, -0.01924576393436911, 0.03441793111567203, -0.009229382101760265, 0.0007099369818066644, 0.00630563269451746, 0.008375146833970948, 0.001453733482001713, -0.005900793277074788, -0.002285254065278213, -0.002141335465978059, -0.0004177599299066708, -0.001494779621447613, -0.002049868015261339, 0, 0, 0, 0])
    cOmgA = np.array([0.00614611792998422, 0.008253100851149026, -0.01440165141619654, 0.003363590350788535, -7.138615291626988e-05, -0.002504786979418468, -0.00172978832643207, -0.0006280861013429611, 0.001241749955604002, 0.0009224361511874661, 0.0004610771596491818, -0.001613979006196489, 0.0006367428132294327, 0.0004010956619564596, 0, 0, 0, 0])
    sOmgA = np.array([-0.04155568953790275, 0.0257426196723017, -0.002959273392809311, 0.004475809265755418, 1.822441292043207e-05, -0.0001972760876678778, 0.000389971927172294, 0.003913904086152674, 0.0004058488092230152, -0.001787289168266385, -0.0009302656497305446, -2.067134029104406e-05, -0.0013107116813526, 5.625225752812272e-05, 0, 0, 0, 0])
    cChiA = np.array([-0.06673908312554792, 0.06550733801292973, -0.007055149797375992, 0.005111848628877972, 0, -0.0005444464620177098, 0.0009830562551572195, 0.009386235733694169, 0, -0.003177877146985308, -0.004324046613805478, 0, 0, -0.001615990759958801, 0.001587849478343136, -0.002398762740975183, 0.002838548328494804, 0.0005357813386138708])
    sChiA = np.array([-0.01069967856443793, -0.02029794993715239, 0.03266650186037179, -0.0041544791939612, 0, 0.004640389727239152, 0.008287602553739408, 0.0007486759753624905, 0, -0.00118062300801947, -0.001970956729830991, 0, 0, -0.002165451504436122, -0.005086043543188153, -0.001461733557390353, 0.0002004643484864111, 0.000690981600754813])
    
    psiA = 0.04107992866630529 + T*(0.02444817476355586 + T*(-3.592047589119096e-08 + 1.401111538406559e-12*T))
    omgA = 0.4086163677095374 + T*(-2.150908863572772e-06 + T*(7.078279744199225e-12 + 7.320686584753994e-13*T));
    chiA = -9.530113429264049e-05 + T*(3.830798934518299e-07 + T*(7.13645738593237e-11 - 2.957363454768169e-13*T));

    cosAng = np.cos(omega*T)
    sinAng = np.sin(omega*T)
    psiA += sum(cPsiA*cosAng + sPsiA*sinAng)
    omgA += sum(cOmgA*cosAng + sOmgA*sinAng)
    chiA += sum(cChiA*cosAng + sChiA*sinAng)
    cEps = 0.9174821430652418; sEps = 0.397776969112606;
    sPsi = math.sin(psiA); cPsi = math.cos(psiA);
    sOmg = math.sin(omgA); cOmg = math.cos(omgA);
    sChi = math.sin(chiA); cChi = math.cos(chiA);
    
    p = [[0,0,0],[0,0,0],[0,0,0]]
    p[0][0] = cChi*cPsi + sChi*cOmg*sPsi
    p[0][1] = (-cChi*sPsi + sChi*cOmg*cPsi)*cEps + sChi*sOmg*sEps
    p[0][2] = (-cChi*sPsi + sChi*cOmg*cPsi)*sEps - sChi*sOmg*cEps
    p[1][0] = -sChi*cPsi + cChi*cOmg*sPsi
    p[1][1] = (sChi*sPsi + cChi*cOmg*cPsi)*cEps + cChi*sOmg*sEps
    p[1][2] = (sChi*sPsi + cChi*cOmg*cPsi)*sEps - cChi*sOmg*cEps
    p[2][0] = sOmg*sPsi
    p[2][1] = sOmg*cPsi*cEps - cOmg*sEps
    p[2][2] = sOmg*cPsi*sEps + cOmg*cEps
    # frame bias matrix
    b = [[0.9999999999999942, -7.078279744199226e-8, 8.05614893899716e-8],
         [7.078279477859602e-8, 0.999999999999997, 3.306041454222148e-8],
         [-8.056149173008023e-8, -3.30604088398539e-8, 0.9999999999999962]]
    pb = [[0,0,0],[0,0,0],[0,0,0]]
    pb[0][0] = p[0][0]*b[0][0] + p[0][1]*b[1][0] + p[0][2]*b[2][0]
    pb[0][1] = p[0][0]*b[0][1] + p[0][1]*b[1][1] + p[0][2]*b[2][1]
    pb[0][2] = p[0][0]*b[0][2] + p[0][1]*b[1][2] + p[0][2]*b[2][2]
    pb[1][0] = p[1][0]*b[0][0] + p[1][1]*b[1][0] + p[1][2]*b[2][0]
    pb[1][1] = p[1][0]*b[0][1] + p[1][1]*b[1][1] + p[1][2]*b[2][1]
    pb[1][2] = p[1][0]*b[0][2] + p[1][1]*b[1][2] + p[1][2]*b[2][2]
    pb[2][0] = p[2][0]*b[0][0] + p[2][1]*b[1][0] + p[2][2]*b[2][0]
    pb[2][1] = p[2][0]*b[0][1] + p[2][1]*b[1][1] + p[2][2]*b[2][1]
    pb[2][2] = p[2][0]*b[0][2] + p[2][1]*b[1][2] + p[2][2]*b[2][2]
    return pb
