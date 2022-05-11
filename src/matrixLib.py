import numpy as np
from math import sqrt

def fwdSub(matrix, vector):
  y = np.zeros([vector.size], dtype=float)
  y[0] = vector[0]
  for i in range (1, vector.size, 1):
    y[i] = (vector[i] - sum((y[j]*matrix[i][j]) for j in range(0, i)))
  return y

def retSub(matrix, vector):
  x = np.zeros([vector.size], dtype=float)
  x[vector.size-1] = vector[vector.size-1]/matrix[vector.size-1][vector.size-1]
  for i in range (vector.size - 2, -1, -1):
    x[i] = (vector[i] - sum((x[j]*matrix[i][j]) for j in range(vector.size-1, i, -1)))/matrix[i][i]
  return x

def eucNorm(vector):
  eucNorm = 0
  for i in range(0, vector.size, 1):
      eucNorm += pow(vector[i], 2)
  return(sqrt(eucNorm))

def infNorm(vector):
  infNorm = 0
  for i in range(0, vector.size, 1):
      if abs(vector[i]) > infNorm:
        infNorm = abs(vector[i])
  return(infNorm)

def luDec(matrix, vector):
  for k in range (0, matrix[0].size - 1, 1):
    for i in range (k+1, matrix[0].size, 1):
      matrix[i][k] = (matrix[i][k])/(matrix[k][k])
    for j in range (k+1, matrix[0].size, 1):
      for l in range (k+1, matrix[0].size, 1):
        matrix[l][j] = matrix[l][j]-(matrix[l][k]*matrix[k][j])

  y = fwdSub(matrix, vector)
  x = retSub(matrix, y)
  return x

def choleskyDec(matrix, vector):
  lMatrix = np.zeros([matrix[0].size,matrix[0].size], dtype=float)
  for i in range(0, matrix[0].size, 1):
    lMatrix[i][i] = sqrt(matrix[i][i] - sum((pow(lMatrix[i][k], 2)) for k in range(0, i)))
    for j in range (i+1, matrix[0].size, 1):
      lMatrix[j][i] = (1/(lMatrix[i][i]))*(matrix[i][j] - sum((lMatrix[i][k]*lMatrix[j][k]) for k in range(0, i)))

  y = np.zeros([vector.size], dtype=float)
  y[0] = vector[0]/lMatrix[0][0]
  for i in range (1, vector.size, 1):
    y[i] = (vector[i] - sum((y[j]*lMatrix[i][j]) for j in range(0, i)))/lMatrix[i][i]

  x = np.zeros([y.size], dtype=float)
  x[y.size-1] = y[y.size-1]/lMatrix[y.size-1][y.size-1]
  for i in range (y.size - 2, -1, -1):
    x[i] = (y[i] - sum((x[j]*lMatrix[j][i]) for j in range(y.size-1, i, -1)))/lMatrix[i][i]

  return x

def itrJacobi(matrix, vector, tolm):
  x0 = np.zeros([vector.size], dtype=float)
  x = np.zeros([vector.size], dtype=float)
  cont = 0
  res = 1
  while(res > tolm):
    for i in range(0, vector.size, 1):
      x[i] = (vector[i] - sum((matrix[i][j]*x0[j] if j!=i else 0) for j in range (0, vector.size, 1)))/matrix[i][i]
    res = eucNorm(np.subtract(x, x0))/eucNorm(x)
    x0 = np.copy(x)
    cont += 1
  return x

def itrGaussSeidel(matrix, vector, tolm):
  x0 = np.zeros([vector.size], dtype=float)
  x = np.zeros([vector.size], dtype=float)
  cont = 0
  res = 1
  while(res > tolm):
    for i in range(0, vector.size, 1):
      x[i] = (vector[i] - sum((matrix[i][j]*x[j] if j!=i else 0) for j in range (0, vector.size, 1)))/matrix[i][i]
    res = eucNorm(np.subtract(x, x0))/eucNorm(x)
    x0 = np.copy(x)
    cont += 1
  return x


def powerMet(matrix, vector):
  print('')

def jacobiMet(matrix, vector):
  print('')

def interpolation():
  print('')

def regression():
 print('')


