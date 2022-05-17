import numpy as np
from math import sqrt, sin, cos, atan, pi

np.seterr(all='raise')

def __fwdSub(matrix, vector):
  """"Calculates the foward substitution in LU decomposition"""
  y = np.zeros([vector.size], dtype=float)
  y[0] = vector[0]
  for i in range (1, vector.size, 1):
    y[i] = (vector[i] - sum((y[j]*matrix[i][j]) for j in range(0, i)))
  return y

def __retSub(matrix, vector):
  """"Calculates the retro substitution in LU decomposition"""
  x = np.zeros([vector.size], dtype=float)
  x[vector.size-1] = vector[vector.size-1]/matrix[vector.size-1][vector.size-1]
  for i in range (vector.size - 2, -1, -1):
    x[i] = (vector[i] - sum((x[j]*matrix[i][j]) for j in range(vector.size-1, i, -1)))/matrix[i][i]
  return x

def __detTrg(matrix):
  """Calculates the determinant for an upper or lower triangular matrix"""
  det = 0
  for i in range(0, matrix[0].size, 1):
    det += matrix[i][i]
  return det

def __eigenDet(eigen):
  det = 1
  for i in range(0, eigen.size, 1):
    det *= eigen
  return det

def __diagDom(matrix):
  """"Calculates wether a given matrix is diagonally dominant"""
  sumL = 0
  sumC = 0
  for i in range(0, matrix[0].size, 1):
    for j in range(0, matrix[0].size, 1):
      if i!=j:
        sumL += matrix[i][j]
        sumC += matrix[j][i]
    if sumL > matrix[i][i] or sumC > matrix[i][i]:
      return False
    sumL = 0
    sumC = 0
  return True

def __symmetric(matrix):
  if( matrix.any() != np.transpose(matrix).any):
    return False
  return True

def __eucNorm(vector):
  """Calculates the euclidian norm of a given matrix"""
  eucNorm = 0
  for i in range(0, vector.size, 1):
      eucNorm += pow(vector[i], 2)
  return(sqrt(eucNorm))

def __infNorm(vector):
  """Calculates the infinite norm of a given matrix"""
  infNorm = 0
  for i in range(0, vector.size, 1):
      if abs(vector[i]) > infNorm:
        infNorm = abs(vector[i])
  return(infNorm)

def __diagonal(matrix, tolm):
  for i in range(0, matrix[0].size, 1):
    for j in range(0, matrix[0].size, 1):
      if j!=i and matrix[i][j]>tolm:
        return False
  return True

def luDec(matrix, vector, detCod):
  """Solves a given linear equation system through the LU Decomposition
  Optional: calculates the determinant of the A matrix of the system"""
  try:
    for k in range (0, matrix[0].size - 1, 1):
      for i in range (k+1, matrix[0].size, 1):
        matrix[i][k] = (matrix[i][k])/(matrix[k][k])
      for j in range (k+1, matrix[0].size, 1):
        for l in range (k+1, matrix[0].size, 1):
          matrix[l][j] = matrix[l][j]-(matrix[l][k]*matrix[k][j])
  except:
    return ['Unable to complete LU decomposition. Matrix is singular.']
  else:
    y = __fwdSub(matrix, vector)
    x = __retSub(matrix, y)
    result = [
      "X vector: " + str(x)
    ]
    if detCod>0:
      det = __detTrg(matrix)
      result.append(
        "Determinant: " + str(det)
      )

  return result
  
def choleskyDec(matrix, vector):
  """Solves a given linear equation system through the Cholesky Decomposition"""
  try:
    lMatrix = np.zeros([matrix[0].size,matrix[0].size], dtype=float)
    for i in range(0, matrix[0].size, 1):
      lMatrix[i][i] = sqrt(matrix[i][i] - sum((pow(lMatrix[i][k], 2)) for k in range(0, i)))
      for j in range (i+1, matrix[0].size, 1):
        lMatrix[j][i] = (1/(lMatrix[i][i]))*(matrix[i][j] - sum((lMatrix[i][k]*lMatrix[j][k]) for k in range(0, i)))
  except:
    return ['Unable to complete Cholesky Decomposition. Matrix is not symmetric positive definite.']
  else:
    y = np.zeros([vector.size], dtype=float)
    y[0] = vector[0]/lMatrix[0][0]
    for i in range (1, vector.size, 1):
      y[i] = (vector[i] - sum((y[j]*lMatrix[i][j]) for j in range(0, i)))/lMatrix[i][i]

    x = np.zeros([y.size], dtype=float)
    x[y.size-1] = y[y.size-1]/lMatrix[y.size-1][y.size-1]
    for i in range (y.size - 2, -1, -1):
      x[i] = (y[i] - sum((x[j]*lMatrix[j][i]) for j in range(y.size-1, i, -1)))/lMatrix[i][i]

    result = [
      "X vector: " + str(x)
    ]

    return result

def itrJacobi(matrix, vector, tolm):
  """Solves a given linear equation system through the Jacobi Iteration"""
  try:
    result = []
    if (__diagDom(matrix) == False):
      result.append('Matrix is not diagonal dominant. Convergence not guaranteed.')

    x0 = np.zeros([vector.size], dtype=float)
    x = np.zeros([vector.size], dtype=float)
    cont = 0
    res = 1
    resHist = []

    while(res > tolm):
      for i in range(0, vector.size, 1):
        x[i] = (vector[i] - sum((matrix[i][j]*x0[j] if j!=i else 0) for j in range (0, vector.size, 1)))/matrix[i][i]
      res = __eucNorm(np.subtract(x, x0))/__eucNorm(x)
      resHist.append(res)
      x0 = np.copy(x)
      cont += 1
  
    result.append(
      'X vector: ' + str(x)
    )
    result.append(
      'Number of iterations: ' + str(cont)
    )
    result.append(
      'Error History: '
    )
    result.append(
      resHist
    )

  except:
    result.append('Unable to complete operation. Could not converge.')
  
  return result

def itrGaussSeidel(matrix, vector, tolm):
  """Solves a given linear equation system through the Gauss-Seidel Iteration"""
  try:
    result = []
    if (__diagDom(matrix) == False):
      result.append('Matrix is not diagonal dominant. Convergence not guaranteed.')

    x0 = np.zeros([vector.size], dtype=float)
    x = np.zeros([vector.size], dtype=float)
    cont = 0
    res = 1
    resHist = []
    
    while(res > tolm):
      for i in range(0, vector.size, 1):
        x[i] = (vector[i] - sum((matrix[i][j]*x[j] if j!=i else 0) for j in range (0, vector.size, 1)))/matrix[i][i]
      res = __eucNorm(np.subtract(x, x0))/__eucNorm(x)
      resHist.append(res)
      x0 = np.copy(x)
      cont += 1
    
    result.append(
      'X vector: ' + str(x)
    )
    result.append(
      'Number of iterations: ' + str(cont)
    )
    result.append(
      'Error History: '
    )
    result.append(
      resHist
    )

  except:
    result.append(
      'Unable to complete operation. Could not converge.'
    )
  
  return result

def powerMet(matrix, tolm, detCod):
  x = np.full(matrix[0].size, 1)
  x[0] = 1
  ups0 = __infNorm(x)
  cont = 0
  res = 1
  result = []
  
  while(res > tolm):
    y = np.matmul(matrix, x)
    ups1 = __infNorm(y)
    res = abs(ups1 - ups0)/ups1
    if res>tolm:
      cont += 1
      x = y/ups1
      ups0 = ups1
  
  result = [
    'X vector:' + str(x),
    'Number of iterations: ' + str(cont),
    'Residue: ' + str(res)
  ]

  if detCod>0:
    det = __eigenDet(x)
    result.append(
      'Determinant: ' + str(det)
    )

  return result

def jacobiMet(matrix, tolm, detCod):
  if(__symmetric(matrix) == False):
    return['Could not complete Jacobi Method. Matrix is not symmetric.']
  else:
    pMatrix = np.zeros([matrix[0].size, matrix[0].size], dtype=float)
    aMatrix = np.copy(matrix)
    cont = 0
    max = 0
    line = 0
    column = 0
    o = 0

    while(__diagonal(pMatrix, tolm) == False):
      cont += 1

      for i in range(0, aMatrix[0].size, 1):
        for j in range(0, aMatrix[0].size, 1):
          if j!=i and abs(aMatrix[i][j])>max:
            max = abs(aMatrix[i][j])
            line = i
            column = j

      for i in range(0, aMatrix[0].size, 1):
        for j in range(0, aMatrix[0].size, 1):
          if j!=i and abs(aMatrix[i][j])>max:
            max = abs(aMatrix[i][j])
            line = i
            column = j
        pMatrix[i][i] = 1
      
      if(aMatrix[line][line] != aMatrix[column][column]):
        o = (atan*(aMatrix[line][column]/(aMatrix[line][line] - aMatrix[column][column])))/2
      else:
        o = pi/4

      pMatrix[line][column] = sin(o)
      pMatrix[column][line] = -sin(o)
      pMatrix[line][line] = cos(o)
      pMatrix[column][column] = cos(o)

      apmult = np.matmul()
      if(__diagonal(pMatrix, tolm) == False):
        aMatrix = pMatrix
  
  result = []

  if detCod>0:
    det = __eigenDet(x)
    result.append(
      'Determinant: ' + str(det)
    )

  return result

def interpolation():
  print('')

def regression():
 print('')


