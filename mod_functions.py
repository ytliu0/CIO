import math

# restrict x to the range [-pi, pi) by subtracting integer multiples of 2 pi.
def mod2pi(x): 
  return x -2*math.pi*math.floor(0.5*x/math.pi + 0.5)

def quotient_remainder(x, k):
  """
  Given x and k, calculate q = floor(x/k + 0.5) and r = x - q*k.
  x and k are assumed to be integers.
  """
  q = math.floor(x/k + 0.5)
  r = x - q*k
  return q,r

def mod2pi_omgDf(omg, D, f):
  """
  Calculate mod(omg*(D + f), 2*pi), where D is an integer and omg < 2*pi.
  This calculation takes into account the integer nature of D to prevent the
  loss of precision caused by truncation error.
  """
  tpi = 2*math.pi
  x = omg*f; ph = omg*D; omg1 = omg; qD=D; rD=0;
  while abs(ph) > tpi:
    p = abs(tpi/omg1) + 0.5
    k = math.floor(p)
    qD, rD = quotient_remainder(qD, k)
    x += omg1*rD
    omg1 *= (k-p+0.5)
    ph = omg1*qD
  return mod2pi(x + ph)