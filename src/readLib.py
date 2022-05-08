import numpy as np
import matrixLib as mtx

def readFileT1(file):
  with open(file, 'r') as f:
    order = int(f.readline())
    icod = int(f.readline())
    idet = int(f.readline())
    aMatrix = np.zeros([order,order], dtype=float)
    bVector = np.zeros([order], dtype=float)
    for i in range(0, order, 1):
      line = (f.readline()).split()
      for j in range(0, order+1, 1):
        if j!=order:
          aMatrix[i][j] = int(line[j])
        else:
          bVector[i]=int(line[j])
    if icod==3 or icod==4:
      tolm = int(f.readline())

    if(icod==1):
      mtx.luDec(aMatrix, bVector)
    elif(icod==2):
      mtx.choleskyDec(aMatrix, bVector)
    elif(icod==3):
      mtx.itrJacobi(aMatrix, bVector, tolm)
    elif(icod==4):
      mtx.itrGaussSeidel(aMatrix, bVector, tolm)

def readFileT2(file):
  with open(file, 'r') as f:
    order = f.readline()
    icod = f.readline()
    idet = f.readline()
    tolm = f.readline()

    if(icod==1):
      mtx.powerMet()
    elif(icod==2):
      mtx.jacobiMet()


def readFileT3(file):
  with open(file, 'r') as f:
    icod = f.readline()
    n = f.readline()
    x = f.readline()

    if(icod==1):
      mtx.interpolation()
    elif(icod==2):
      mtx.regression()
