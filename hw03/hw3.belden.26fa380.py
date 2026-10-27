import numpy as np
import sympy as sp
import time

def dot(V_1, V_2):

    sum = 0
    for i in range(len(V_1)):
        product = V_1[i] * V_2[i]
        sum += product

    return sum

def build_vector_w_distance(P_1, P_2):

    v = np.zeros(3)
    for i in range(len(P_1)):
        coordinate = P_2[i] - P_1[i]
        v[i] = coordinate

    return v

def my_interesting_build_planar_cloud(P: np.ndarray, N: np.ndarray, numOfPoints):
    """ 
        build a random planar point cloud 
        probability of randoming finding {numOfPoints} valid points on plane is highly unlikely with this technique,
        although it is technically mathematically sound
    """

    vectorsOnPlane = 0
    while vectorsOnPlane != 2:
        testP = np.random.uniform(-10, 10, size = 3)

        testV = build_vector_w_distance(P, testP)
        orthoWithDot = dot(N, testV)

        if abs(orthoWithDot) < 1e-6:
            if vectorsOnPlane == 0:
                v_1 = np.array([testV[0], testV[1], testV[2]])
                vectorsOnPlane += 1
            elif vectorsOnPlane == 1:
                v_2 = np.array([testV[0], testV[1], testV[2]])
                vectorsOnPlane += 1

    pointsOnPlane = np.zeros((numOfPoints, 3))
    for i in range(numOfPoints):

        # random coefficients
        a = np.random.uniform(-1, 1)
        b = np.random.uniform(-1, 1)

        # construct a point on the plane
        p = P + a * v_1 + b * v_2

        pointsOnPlane[i, 0] = p[0]
        pointsOnPlane[i, 1] = p[1]
        pointsOnPlane[i, 2] = p[2]

        time.sleep(.1)

    return pointsOnPlane

def create_ply(file, planePoints):

    vertex = planePoints.shape[0]

    with open(file, "w") as ply:
        ply.write("ply\n")
        ply.write("format ascii 1.0\n")
        ply.write(f"element vertex {vertex}\n")
        ply.write("property float x\n")
        ply.write("property float y\n")
        ply.write("property float z\n")
        ply.write("element face 0\n")
        ply.write("end_header\n")

        for i in range(vertex):
            x = planePoints[i, 0]
            y = planePoints[i, 1]
            z = planePoints[i, 2]

            ply.write(f"{x} {y} {z}\n")

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

# PCA w/ spectral analysis w/ random planar cloud
P = np.array([1, 2, 3])
N = np.array([0, 0, 1])
planarPointCloud = my_interesting_build_planar_cloud(P, N, 100)
create_ply("planar_plane.ply", planarPointCloud)
meanSpectral, normalSpectral = spectralPCA("planar_plane.ply")









