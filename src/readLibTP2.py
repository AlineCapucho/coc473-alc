from xml.etree.ElementTree import C14NWriterTarget
import numpy as np
import calcLib as calc

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
        icod = int(f.readline())

        tethas = (f.readline()).split()
        theta1 = float(tethas[0])
        theta2 = float(tethas[1])

        tolm = float(f.readline())

        if icod == 1:
            result = calc.newton(theta1, theta2, tolm)
        elif icod == 2:
            result = calc.broyden(theta1, theta2, tolm)

    return result

def readTask2(file):
    with open(file, "r") as f:
        icod = int(f.readline())
        if(icod != 4):
            icod2 = int(f.readline())

        c = (f.readline()).split()

        if(icod == 1):
            limits = (f.readline()).split()
            a = float(limits[0])
            b = float(limits[1])
            tolm = float(f.readline())

        if(icod == 2):
            limits = (f.readline()).split()
            a = float(limits[0])
            b = float(limits[1])
            numP = int(f.readline())

        if(icod == 3):
            x = float(f.readline())
            deltaX = float(f.readline())
        
        if(icod == 4):
            x = float(f.readline())
            deltaX1 = float(f.readline())
            deltaX2 = float(f.readline())
        
        if icod == 1:
            if icod2 == 1:
                result = calc.rootBisect(c, a, b, tolm)
            if icod2 == 2:
                result = calc.rootNewton(c, a, b, tolm)
        elif icod == 2:
            if icod2 == 1:
                result = calc.integralPol(c, a, b, numP)
            if icod2 == 2:
                result = calc.integralGauss(c, a, b, numP)
        elif icod == 3:
            if icod2 == 1:
                result = calc.diffStepFwd(c, x, deltaX)
            if icod2 == 2:
                result = calc.diffStepBack(c, x, deltaX)
            if icod2 == 3:
                result = calc.diffCentral(c, x, deltaX)
        elif icod == 4:
            result = calc.diffRe(c, x, deltaX1, deltaX2)
            
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
