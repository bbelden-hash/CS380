import numpy as np
import sympy as sp
import time

def read_file(file):

    try:
        with open(file, "r") as file:
            fileName = file.readlines()
        return fileName
    except FileNotFoundError:
        print("error: file is not valid or non-existent")
        return None

def get_mean(file):

    sumX = 0
    sumY = 0
    sumZ = 0
    num = 0

    header = True

    for line in file:

        if header:
            if line.strip() == "end_header":
                header = False
            continue

        coordinates = line.split()
        if float(coordinates[0]) == 3:
            break
        else:
            x = float(coordinates[0])
            y = float(coordinates[1])
            z = float(coordinates[2])

            sumX += x 
            sumY += y 
            sumZ += z

            num += 1
    
    meanX = sumX / num
    meanY = sumY / num
    meanZ = sumZ / num

    theMeans = np.array([meanX, meanY, meanZ], dtype = np.float64)
    return theMeans

# create a function that takes a point cloud, and subtracts the mean from every point -> translating the point cloud so its center is at the origin
# points' relative position and shape stay the same
def move_point_cloud(file, theMeans, size):

    # a matrix that represents the transform towards the origin 
    transform = np.array([
        [1, 0, 0, -theMeans[0]],
        [0, 1, 0, -theMeans[1]],
        [0, 0, 1, -theMeans[2]],
        [0, 0, 0, 1]
    ])

    points = np.zeros((4, size), dtype = np.float64)
    num = 0

    header = True

    for line in file:

        if header:
            if line.strip() == "end_header":
                header = False
            continue

        coordinates = line.split()
        if float(coordinates[0]) == 3:
            break
        else:
            x = float(coordinates[0])
            y = float(coordinates[1])
            z = float(coordinates[2])

            points[0, num] = x
            points[1, num] = y
            points[2, num] = z
            points[3, num] = 1

            num += 1
                
    centered = np.matmul(transform, points)
    return centered

# how points are spread out in space
# if x tends to increase when y increases, covariance is positive / if one tends to increase while other decreases, covariance is negative
# my covariance function for learning purposes --> pretty slow
def my_covariance(centeredCloud):

    xValues = centeredCloud[0].copy()
    yValues = centeredCloud[1].copy()
    zValues = centeredCloud[2].copy()

    num = 0

    C_xx = 0
    C_xy = 0
    C_xz = 0
    C_yx = 0
    C_yy = 0
    C_yz = 0
    C_zx = 0
    C_zy = 0
    C_zz = 0

    for x in range(len(xValues)):
        C_xx += xValues[x] ** 2
        num += 1

    C_xx = C_xx / num

    for x in range(len(xValues)):
        for y in range(len(yValues)):
            if x == y:
                C_xy += xValues[x] * yValues[y]

    C_xy = C_xy / num

    for x in range(len(xValues)):
        for z in range(len(zValues)):
            if x == z:
                C_xz += xValues[x] * zValues[z]
    
    C_xz = C_xz / num

    for i in range(num):
        C_yx += yValues[i] * xValues[i]
        C_yy += yValues[i] ** 2
        C_yz += yValues[i] * zValues[i]

    C_yx = C_yx / num
    C_yy = C_yy / num
    C_yz = C_yz / num

    for i in range(num):
        C_zx += zValues[i] * xValues[i]
        C_zy += zValues[i] * yValues[i]
        C_zz += zValues[i] ** 2

    C_zx = C_zx / num
    C_zy = C_zy / num
    C_zz = C_zz / num 

    covarianceMatrix = np.array([
        [C_xx, C_xy, C_xz],
        [C_yx, C_yy, C_yz],
        [C_zx, C_zy, C_zz]
    ])

    return covarianceMatrix

# NumPy covariance function for speed
def numpy_covariance(centeredCloud, n):

    coordinates = centeredCloud[:n, :]
    covarianceMatrix = coordinates @ coordinates.T / coordinates.shape[1]

    return covarianceMatrix

# for a 3x3 matrix ...
def my_eigenvalues(coMatrix):

    a = sp.symbol('a')
    
    C = sp.Matrix(coMatrix)
    characteristicMatrix = C - a * sp.eye(3)
    detCo = characteristicMatrix.det()

    eigenvalues = sp.solve(detCo, a)

    return eigenvalues

def numpy_eigen(coMatrix):

    eigenvalues, eigenvectors = np.linalg.eigh(coMatrix)
    return eigenvalues, eigenvectors

def spectralPCA(filename):

    # read the file
    myFile = read_file(filename)

    # get the mean
    mean = get_mean(myFile)

    # move point cloud to origin
    cloudToOrigin = move_point_cloud(myFile, mean, TBD)

    # get the covariance matrix
    covariance = numpy_covariance(cloudToOrigin, 3)

    # get the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = numpy_eigen(covariance)

    # get eigenvector corresponding to smallest eigenvalue
    normal = eigenvectors[:, 0]

    #return mean and normal
    return mean, normal

# PCA w/ spectral analysis
meanSpectral, normalSpectral = spectralPCA(<filename>)







