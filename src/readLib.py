import numpy as np
import matrixLib as mtx

def runOperation(readFile, writeFile, task):
  if task == 1:
    results = readTask1(readFile)
  elif task == 2:
    results = readTask2(readFile)
  elif task == 3:
    results = readTask3(readFile)
  write(writeFile, results)

def readTask1(file):
  with open(file, 'r') as f:
    order = int(f.readline())
    icod = int(f.readline())
    idet = float(f.readline())
    aMatrix = np.zeros([order,order], dtype=float)
    bVector = np.zeros([order], dtype=float)
    for i in range(0, order, 1):
      line = (f.readline()).split()
      for j in range(0, order+1, 1):
        if j!=order:
          aMatrix[i][j] = float(line[j])
        else:
          bVector[i]=float(line[j])
    if icod==3 or icod==4:
      tolm = float(f.readline())

    if(icod==1):
      result = mtx.luDec(aMatrix, bVector, idet)
    elif(icod==2):
      result = mtx.choleskyDec(aMatrix, bVector)
    elif(icod==3):
      result = mtx.itrJacobi(aMatrix, bVector, tolm)
    elif(icod==4):
      result = mtx.itrGaussSeidel(aMatrix, bVector, tolm)
    if(idet>0 and icod in [2,3,4]):
      result.append("Warning! Could not calculate determinant.")
  
  return result

def readTask2(file):
  with open(file, 'r') as f:
    order = int(f.readline())
    icod = int(f.readline())
    idet = float(f.readline())
    aMatrix = np.zeros([order,order], dtype=float)
    for i in range(0, order, 1):
      line = (f.readline()).split()
      for j in range(0, order, 1):
        if j!=order:
          aMatrix[i][j] = float(line[j])
    tolm = float(f.readline())

    if(icod==1):
      result = mtx.powerMet(aMatrix, tolm)
    elif(icod==2):
      result = mtx.jacobiMet(aMatrix, tolm)

  return result

def readTask3(file):
  with open(file, 'r') as f:
    icod = f.readline()
    n = f.readline()
    x = f.readline()

    if(icod==1):
      result = mtx.interpolation()
    elif(icod==2):
      result = mtx.regression()

  return result

def write(file, values):
  with open(file, 'w') as f:
    for value in values:
      f.write(str(value)+"\n")
