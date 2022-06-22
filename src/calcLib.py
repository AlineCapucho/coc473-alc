import numpy as np
import matrixLib as mtx
from math import sqrt, sin, cos, atan, pi

def __getNlFunction(c, theta1, theta2):
  c2 = c[0]
  c3 = c[1]
  c4 = c[2]
  f1 = 2*(c3**2) + (c2**2) + 6*(c4**2) - 1
  f2 = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - theta1
  f3 = 60*(c3**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2) + 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3)*c4 + 3*c2 - theta2
  result = [-f1,-f2,-f3]
  return result


def __getNlFunctionEx(c):
  x1 = c[0]
  x2 = c[1]
  f1 = x1 + 2*x2 - 2
  f2 = (x1**2) + 4*(x2**2) - 4
  result = [f1,f2]
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


def __getJacobianEx(c):
  x1 = c[0]
  x2 = c[1]

  f1x1 = 1
  f1x2 = 2
  f2x1 = 2*x1
  f2x2 = 8*x2

  result = [[f1x1,f1x2],[f2x1,f2x2]]

  return result


def __getApprJacobian(jacAppr, y, deltaX):
  deltaX = -deltaX

  term1 = np.matmul(jacAppr, deltaX) #jacAppr*deltaX
  term2 = y+term1 # y - jacAppr*deltaX
  term3 = np.outer(term2, deltaX) # (y - jacAppr*deltaX)*(deltaX)T
  # term4 = np.outer(deltaX, deltaX) # (deltaX)T*deltaX
  term4 = np.matmul(deltaX, deltaX) # (deltaX)T*deltaX
  # term4 = np.multiply(deltaX, deltaX) # (deltaX)T*deltaX
  # term4 = np.dot(deltaX, deltaX) # (deltaX)T*deltaX
  term5 = np.divide(term3,term4)

  result = jacAppr - term5

  return result

def newton(theta1, theta2, tolm):
  nInter = 1000
  count = 0
  x0 = [1, 0, 4]
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
  x0 = np.array([2, 3])
  x1 = np.array([0]*2)
  y = np.array([0]*2)
  deltaX = np.array([0]*2)
  func0 = np.array(__getNlFunctionEx(x0))
  jacAppr = np.array(__getJacobianEx(x0))
  
  while(count < nInter):
    deltaX = np.matmul(-np.linalg.inv(jacAppr),func0)
    x1 = x0 + deltaX

    if(mtx.__eucNorm(x1)!=0):
      tol = mtx.__eucNorm(deltaX)/mtx.__eucNorm(x1)
    else:
      tol = 0
    if(tolm >= tol):
      print(x1)
      return x1

    func1 = np.array(__getNlFunctionEx(x1))
    y = func1 - func0
    jacAppr = __getApprJacobian(jacAppr, y, deltaX)

    x0 = x1
    func0 = func1

    count += 1

  raise Exception('Convergence unreachable.')


def root():
  pass


def integral():
  pass


def diffDf():
  pass


def diffRe():
  pass


def diffEq():
  pass