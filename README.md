# CIO Locator and Equation of Origin

CIO locator *s* and equation of origin *Eo* play an important role in the transformation between Geocentric Celestial Reference System (GCRS) and Terrestrial Intermediate Reference System (TIRS), which are important in astronomical calculations. The formulas currently available are based on the precession and nutation models recommended by IAU in 2000 and 2006. However, the precession models are intended for high-accuracy applications over a limited time span. For calculations involving more than few hundred years from 2000, the precession model developed by [Vondrák, Capitaine and Wallace](https://www.aanda.org/articles/aa/full_html/2011/10/aa17274-11/aa17274-11.html) in 2011 is more appropriate. 

This repository provides python functions to calculate the CIO locator *s* and equation of origin *Eo* compatible with the Vondrák et al/IAU2000A precession-nutation model in -4000-8000 to milliarcsecond (mas) accuracy based on semi-analytic formulas I developed. Although they are not as accurate as the microarcsecond (&mu;as) accuracy of the Lagrange interpolation method described [here](http://ytliu.epizy.com/eclipse/Eo.html#sect_Lagrange), it's still useful, especially for times more than hundreds of years from 2000. 

As explained in [Capitaine &amp; Wallace 2006](https://ui.adsabs.harvard.edu/link_gateway/2006A&A...450..855C/doi:10.1051/0004-6361:20054550), the accuracy in the location of the celestial intermediate pole (CIP) is limited by the uncertainties in the parameters of precession models and nutation amplitudes. The accuracy that can be achieved for predicting the CIP location is of order 2 mas after one century, even though the models provide &mu;as precision for any given values of the precession-nutation parameters. It is therefore useful to consider simple semi-analytic expressions of *s* and *Eo* to mas accuracy. The semi-analytic formulas used here are explained on [this page](http://ytliu.epizy.com/eclipse/Eo.html#sect_Semi_analytic). 

In addition to the formulas of *s* and *Eo* for -4000-8000, I also developed another set of expressions for *s* and *Eo* that can be used over ±200 millennia time span. These formulas ignore nutation in the calculation and are therefore not as accurate. 

The Jupyter notebook `examples.ipynb` shows examples of using the functions in this repository.

The following is a list of main functions provided in this repository. Most functions take the two-part Julian date `jd0` and `jd1` as input arguments. Note that `jd0` and `jd1` must be numbers, not arrays.

- `s_Vondrak_IAU2000A_spline(jd0, jd1)` in `s_Vondrak_IAU2000A_spline.py`: Calculate the CIO locate *s* compatible with the Vondrák et al/IAU2000A precession-nutation model at TT Julian date jd = jd0 + jd1 using a spline fitting formula. This formula covers the time span from -4000 to 8000.

- `Eo_Vondrak_IAU2000A_spline(jd0, jd1)` in `Eo_Vondrak_IAU2000A_spline.py`: Calculate the equation of origin *Eo* compatible with the Vondrák et al/IAU2000A precession-nutation model at TT Julian date jd = jd0 + jd1 using a spline fitting formula. This formula covers the time span from -4000 to 8000.

- `s_Vondrak_longT(jd0, jd1)` in `s_Vondrak_longT.py`: Calculate *s* compatible with the Vondrák et al/IAU2000A model at TT Julian date jd = jd0 + jd1. This function covers ±200 millennia time span. It returns the same values as `s_Vondrak_IAU2000A_spline(jd0, jd1)` in -4000-8000, but ignores nutation outside that time interval.

- `Eo_Vondrak_longT(jd0, jd1)` in `Eo_Vondrak_longT.py`: Calculate *Eo* compatible with the Vondrák et al/IAU2000A model at TT Julian date jd = jd0 + jd1. This function covers ±200 millennia time span. It returns the same values as `Eo_Vondrak_IAU2000A_spline(jd0, jd1)` in -4000-8000, but ignores nutation outside that time interval.

- `GAST_from_Eo(jd0_ut1, jd1_ut1, Eo)` in `ERA_GAST.py`: Calculate the Greenwich apparent sidereal time (GAST) at UT1 Julian date `jd_ut1 = jd0_ut1 + jd1_ut1` from *Eo*. It simply subtracts *Eo* from the Earth rotation angle (ERA) computed using the equation defining UT1.


Following [SOFA](http://www.iausofa.org/), Julian date is specified by two parts jd0 and jd1 in any way users may find convenient. For example, JD = 2450123.7 could be expressed in any of these ways, amomg others.

                jd0             jd1
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method) 

