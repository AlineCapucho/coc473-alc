import numpy as np
import matrixLib as mtx
from math import ceil

def __getNlFunction(c, theta1, theta2):
  c2 = c[0]
  c3 = c[1]
  c4 = c[2]
  f1 = 2*(c3**2) + (c2**2) + 6*(c4**2) - 1
  f2 = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - theta1
  f3 = 60*(c3**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2) + 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3)*c4 + 3*c2 - theta2
  result = [-f1,-f2,-f3]
  return result


def __getJacobian(c):
  c2 = c[0]
  c3 = c[1]
  c4 = c[2]

  df1c2 = 2*c2 #derivada parcial da função 1 em relação a c2
  df1c3 = 4*c3 #derivada parcial da função 1 em relação a c3
  df1c4 = 12*c4 #derivada parcial da função 1 em relação a c4

  df2c2 = 12*c3*c2 + 36*c3*c4
  df2c3 = 24*(c3**2) + 6*(c2**2) + 36*c2*c4+ 108*(c4**2)
  df2c4 = 36*c3*c2 + 216*c3*c4

  df3c2 = 120*(c3**2)*c2 + 576*(c3**2)*c4 + 504*(c4**2)*c2 + 1296*(c4**3) + 72*(c2**2)*c4 + 3
  df3c3 = 240*(c3**3) + 120*c3*(c2**2) + 1152*c3*c2*c4 + 4464*c3*(c4**2)
  df3c4 = 576*(c3**2)*c2 + 4464*(c3**2)*c4 + 504*c4*(c2**2)+ 3888*(c4**2)*c2 + 13392*(c4**3) + 24*(c2**3)
  
  result = [[df1c2,df1c3,df1c4],[df2c2,df2c3,df2c4],[df3c2,df3c3,df3c4]]
  
  return result


def __getApprJacobian(jacAppr, y, deltaX):
  deltaX = -deltaX

  term1 = np.matmul(jacAppr, deltaX) #jacAppr*deltaX
  term2 = y+term1 # y - jacAppr*deltaX
  term3 = np.outer(term2, deltaX) # (y - jacAppr*deltaX)*(deltaX)T
  term4 = np.matmul(deltaX, deltaX) # (deltaX)T*deltaX
  term5 = np.divide(term3,term4)

  result = jacAppr - term5

  return result


def __getRootFunction(c,x):
  c1 = float(c[0])
  c2 = float(c[1])
  c3 = float(c[2])
  c4 = float(c[3])
  f = c1**(c2*x) + c3*(x**c4)

  return f


def __getRootDiffFunction(c,x):
  c1 = float(c[0])
  c2 = float(c[1])
  c3 = float(c[2])
  c4 = float(c[3])

  p1 = ((c1*c2)**(c2*x))*np.log(c1)
  p2 = c3*c4*(x**(c4-1))

  df = p1 + p2

  return df


def newton(theta1, theta2, tolm):
  nInter = 1000
  count = 0
  x0 = [1, 0, 0]
  x1 = [0]*3
  s = [0]*3
  tol = [0]*3
  
  while(count < nInter):
    func = np.array(__getNlFunction(x0, theta1, theta2))
    jac = np.array(__getJacobian(x0))
    s = mtx.luDecMod(jac, func, 0)
    for index in range(len(x0)):
      x1[index] = x0[index] + s[index]
      if(x1[index] != 0):
        tol[index] = abs((x1[index]-x0[index]))/x1[index]
      else:
        tol[index] = 0

    aux = list(filter(lambda tolValue: (tolm >= tolValue), tol))
    if len(aux) == len(tol):
      return x1

    x0 = x1.copy()
    
    count += 1
  raise Exception('Convergence unreachable.')


def broyden(theta1, theta2, tolm):
  nInter = 1000
  count = 0
  x0 = np.array([1, 0, 0])
  x1 = np.array([0]*3)
  y = np.array([0]*3)
  deltaX = np.array([0]*3)
  func0 = np.array(__getNlFunction(x0, theta1, theta2))
  jacAppr = np.array(__getJacobian(x0))
  
  while(count < nInter):
    deltaX = np.matmul(-np.linalg.inv(jacAppr),func0)
    x1 = x0 - deltaX

    if(mtx.__eucNorm(x1)!=0):
      tol = mtx.__eucNorm(deltaX)/mtx.__eucNorm(x1)
    else:
      tol = 0
    if(tolm >= tol):
      return x1

    func1 = -np.array(__getNlFunction(x1, theta1, theta2))
    y = func1 - func0
    jacAppr = __getApprJacobian(jacAppr, y, deltaX)

    x0 =  x1.copy()
    func0 = func1.copy()

    count += 1

  raise Exception('Convergence unreachable.')


def rootBisect(c, a, b, tolm):
  x = (a+b)/2
  while(abs(b-a)>tolm and abs(a-b)>tolm):
    f = __getRootFunction(c, x)
    if(f > 0):
      if(abs(a) > abs(b)):
        b = x
      else:
        a = x
    elif(f < 0):
      if(abs(a) > abs(b)):
        a = x
      else:
        b = x
    else:
      return [x]
    
    x = (a+b)/2
  
  return [x]


def rootNewton(c, a, b, tolm):
  x0 = (a+b)/2
  x1 = 0
  nInter = 1000
  count = 0
  while(count < nInter):
    x1 = x0 - (__getRootFunction(c,x0)/__getRootDiffFunction(c,x0))

    tol = abs(x1 - x0)
    if(tol <= tolm):
      return [x1]
    x0 = x1
    count += 1

  raise Exception('Convergence unreachable.')


def __getIntFunction(c, x):
  c1 = float(c[0])
  c2 = float(c[1])
  c3 = float(c[2])
  c4 = float(c[3])

  return (c1**(c2*x))+ c3*(x**c4)


def integralPol(c, a, b, numP):
  if b < a:
    a, b = b, a
  L = b - a

  delta = L/(numP-1)

  wiValues = {
    1: {
      1: L,
    },
    2: {
      1: L/2,
      2: L/2,
    },
    3: {
      1: L/6,
      2: (2*L)/3,
      3: L/6,
    },
    4: {
      1: L/8,
      2: (3*L)/8,
      3: (3*L)/8,
      4: L/8,
    },
    5: {
      1: (7*L)/90,
      2: (16*L)/45,
      3: (2*L)/15,
      4: (16*L)/45,
      5: (7*L)/90,
    },
    6: {
      1: (19*L)/288,
      2: (75*L)/288,
      3: (50*L)/288,
      4: (50*L)/288,
      5: (75*L)/288,
      6: (19*L)/288,
    },
    7: {
      1: (41*L)/(140*6),
      2: (216*L)/(140*6),
      3: (27*L)/(140*6),
      4: (272*L)/(140*6),
      5: (27*L)/(140*6),
      6: (216*L)/(140*6),
      7: (41*L)/(140*6),
    },
    8: {
      1: (751*L)/(17280),
      2: (3577*L)/(17280),
      3: (1323*L)/(17280),
      4: (2989*L)/(17280),
      5: (2989*L)/(17280),
      6: (1323*L)/(17280),
      7: (3577*L)/(17280),
      8: (751*L)/(17280),
    },
    9: {
      1: (989*L)/(28350),
      2: (5888*L)/(28350),
      3: (928*L)/(28350),
      4: (10496*L)/(28350),
      5: (4540*L)/(28350),
      6: (10496*L)/(28350),
      7: (928*L)/(28350),
      8: (5888*L)/(28350),
      9: (989*L)/(28350),
    },
    10: {
      1: (2857*L)/(89600),
      2: (15741*L)/(89600),
      3: (1080*L)/(89600),
      4: (19344*L)/(89600),
      5: (5778*L)/(89600),
      6: (5778*L)/(89600),
      7: (19344*L)/(89600),
      8: (1080*L)/(89600),
      9: (15741*L)/(89600),
      10: (2857*L)/(89600),
    },
  }

  xiValues = {
    1: {
      1: (a+b)/2,
    },
    2: {
      1: a,
      2: b,
    },
    3: {
      1: a,
      2: (a+b)/2,
      3: b,
    },
    4: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: b,
    },
    5: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: b,
    },
    6:{
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: a + (4*delta),
      6: b,
    },
    7: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: a + (4*delta),
      6: a + (5*delta),
      7: b,
    },
    8: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: a + (4*delta),
      6: a + (5*delta),
      7: a + (6*delta),
      8: b,
    },
    9: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: a + (4*delta),
      6: a + (5*delta),
      7: a + (6*delta),
      8: a + (7*delta),
      9: b,
    },
    10: {
      1: a,
      2: a + delta,
      3: a + (2*delta),
      4: a + (3*delta),
      5: a + (4*delta),
      6: a + (5*delta),
      7: a + (6*delta),
      8: a + (7*delta),
      9: a + (8*delta),
      10: b,
    },
  }

  result = 0

  for i in range(1, numP+1):
    x = xiValues[numP][i]
    wi = wiValues[numP][i]
    f = __getIntFunction(c,x)

    result += (wi*f)

  return [result]


def integralGauss(c, a, b, numP):
  if a > b:
    a, b = b, a
  L = b - a
  
  wiValues = {
    1: {
      1: 0,
    },
    2: {
      1: 1,
      2: 1,
    },
    3: {
      1: 0.888888888888889,
      2: 0.555555555555556,
      3: 0.555555555555556,
    },
    4: {
      1: 0.652145154862546,
      2: 0.652145154862546,
      3: 0.347854845137454,
      4: 0.347854845137454,
    },
    5: {
      1: 0.568888888888889,
      2: 0.478628670499367,
      3: 0.478628670499367,
      4: 0.236926885056189,
      5: 0.236926885056189,
    },
    6: {
      1: 0.360761573048139,
      2: 0.360761573048139,
      3: 0.467913934572691, 
      4: 0.467913934572691, 
      5: 0.171324492379170,
      6: 0.171324492379170,
    },
    7: {
      1: 0.417959183673469,
      2: 0.381830050505119,
      3: 0.381830050505119,
      4: 0.279705391489277,
      5: 0.279705391489277,
      6: 0.129484966168870,
      7: 0.129484966168870,
    },
    8: {
      1: 0.3626837833783620,
      2: 0.3626837833783620,
      3: 0.3137066458778873,
      4: 0.3137066458778873,
      5: 0.2223810344533745,
      6: 0.2223810344533745,
      7: 0.1012285362903763,
      8: 0.1012285362903763,
    },
    9: {
      1: 0.3302393550012598,
      2: 0.1806481606948574,
      3: 0.1806481606948574,
      4: 0.0812743883615744,
      5: 0.0812743883615744,
      6: 0.3123470770400029,
      7: 0.3123470770400029,
      8: 0.2606106964029354,
      9: 0.2606106964029354,
    },
    10: {
      1: 0.2955242247147529,
      2: 0.2955242247147529,
      3: 0.2692667193099963,
      4: 0.2692667193099963,
      5: 0.2190863625159820,
      6: 0.2190863625159820,
      7: 0.1494513491505806,
      8: 0.1494513491505806,
      9: 0.0666713443086881,
      10: 0.0666713443086881,
    },
  }

  xiValues = {
    1: {
      1: 2.000000000000000,
    },
    2: {
      1: -0.5773502691896257,
      2: 0.5773502691896257,
    },
    3: {
      1: 0.0000000000000000,
      2: -0.7745966692414834,
      3: 0.7745966692414834,
    },
    4: {
      1: -0.3399810435848563,
      2: 0.3399810435848563,
      3: -0.8611363115940526,
      4: 0.8611363115940526,
    },
    5: {
      1: 0.0000000000000000,
      2: -0.5384693101056831,
      3: 0.5384693101056831,
      4: -0.9061798459386640,
      5: 0.9061798459386640,
    },
    6:{
      1: 0.6612093864662645,
      2: -0.6612093864662645,
      3: -0.2386191860831969,
      4: 0.2386191860831969,
      5: -0.9324695142031521,
      6: 0.9324695142031521,
    },
    7: {
      1: 0.0000000000000000,
      2: 0.4058451513773972,
      3: -0.4058451513773972,
      4: -0.7415311855993945,
      5: 0.7415311855993945,
      6: -0.9491079123427585,
      7: 0.9491079123427585,
    },
    8: {
      1: -0.1834346424956498,
      2: 0.1834346424956498,
      3: -0.5255324099163290,
      4: 0.5255324099163290,
      5: -0.7966664774136267,
      6: 0.7966664774136267,
      7: -0.9602898564975363,
      8: 0.9602898564975363,
    },
    9: {
      1: 0.0000000000000000,
      2: -0.8360311073266358,
      3: 0.8360311073266358,
      4: -0.9681602395076261,
      5: 0.9681602395076261,
      6: -0.3242534234038089,
      7: 0.3242534234038089,
      8: -0.6133714327005904,
      9: 0.6133714327005904,
    },
    10: {
      1: -0.1488743389816312,
      2: 0.1488743389816312,
      3: -0.4333953941292472,
      4: 0.4333953941292472,
      5: -0.6794095682990244,
      6: 0.6794095682990244,
      7: -0.8650633666889845,
      8: 0.8650633666889845,
      9: -0.9739065285171717,
      10: 0.9739065285171717,
    },
  }

  result = 0

  for i in range(1, numP+1):
    x = (a + b + (xiValues[numP][i])*L)/2
    wi = wiValues[numP][i]
    f = __getIntFunction(c,x)

    result += (wi*f)

  result = (result*L)/2

  return [result]


def diffStepFwd(c, x, deltaX):
  if (deltaX == 0):
    return ["Unable to complete. DeltaX must be different than 0."]

  fDelta = __getIntFunction(c, x+deltaX)
  f = __getIntFunction(c, x)

  result = (fDelta - f)/deltaX

  return [result]


def diffStepBack(c, x, deltaX):
  if (deltaX == 0):
    return ["Unable to complete. DeltaX must be different than 0."]

  fDelta = __getIntFunction(c, x-deltaX)
  f = __getIntFunction(c, x)

  result = (f - fDelta)/deltaX

  return [result]


def diffCentral(c, x, deltaX):
  if (deltaX == 0):
    return ["Unable to complete. DeltaX must be different than 0."]

  fDelta1 = __getIntFunction(c, x+deltaX)
  fDelta2 = __getIntFunction(c, x-deltaX)

  result = (fDelta1 - fDelta2)/(2*deltaX)

  return [result]


def diffRe(c, deltaX_1, deltaX_2, tolm):
  pass


def diffEq():
  pass