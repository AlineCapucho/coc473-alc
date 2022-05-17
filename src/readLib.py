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
    with open(file, "r") as f:
        order = int(f.readline())
        icod = int(f.readline())
        idet = float(f.readline())
        aMatrix = np.zeros([order, order], dtype=float)
        bVector = np.zeros([order], dtype=float)
        for i in range(0, order, 1):
            line = (f.readline()).split()
            for j in range(0, order + 1, 1):
                if j != order:
                    aMatrix[i][j] = float(line[j])
                else:
                    bVector[i] = float(line[j])
        if icod == 3 or icod == 4:
            tolm = float(f.readline())

        if icod == 1:
            result = mtx.luDec(aMatrix, bVector, idet)
        elif icod == 2:
            result = mtx.choleskyDec(aMatrix, bVector)
        elif icod == 3:
            result = mtx.itrJacobi(aMatrix, bVector, tolm)
        elif icod == 4:
            result = mtx.itrGaussSeidel(aMatrix, bVector, tolm)
        if idet > 0 and icod in [2, 3, 4]:
            result.append("Warning! Could not calculate determinant.")

    return result


def readTask2(file):
    with open(file, "r") as f:
        order = int(f.readline())
        icod = int(f.readline())
        idet = float(f.readline())
        aMatrix = np.zeros([order, order], dtype=float)
        for i in range(0, order, 1):
            line = (f.readline()).split()
            for j in range(0, order, 1):
                if j != order:
                    aMatrix[i][j] = float(line[j])
        tolm = float(f.readline())

        if icod == 1:
            result = mtx.powerMet(aMatrix, tolm)
            if idet > 0:
                result.append("Warning! Could not calculate determinant.")
        elif icod == 2:
            result = mtx.jacobiMet(aMatrix, tolm, idet)

    return result


def readTask3(file):
    with open(file, "r") as f:
        print("starting t3 task")
        icod = float(f.readline())
        n = int(f.readline())
        x = float(f.readline())
        matrix = np.zeros([n - 1, n - 1])
        for i in range(0, n - 1, 1):
            line = (f.readline()).split(" ")
            matrix[i][0] = float(line[0])
            matrix[i][1] = float(line[1])

        if icod == 1:
            return mtx.interpolation(x, matrix[0], matrix[1])
        elif icod == 2:
            return mtx.regression(x, matrix[0], matrix[1])

    return result


def write(file, values):
    with open(file, "w") as f:
        for value in values:
            f.write(str(value) + "\n")
