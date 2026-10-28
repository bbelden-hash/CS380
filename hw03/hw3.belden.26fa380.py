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

def two_norm(V_1):

    squaredSum = 0
    for i in range(len(V_1)):
       squaredSum += V_1[i] ** 2
    
    norm = squaredSum ** .5
    return norm

# sloppy cross product, bad stuff
def my_cross(V_1, V_2):

    i = (V_1[1] * V_2[2]) - (V_1[2] * V_2[1])
    j = (V_1[0] * V_2[2]) - (V_1[2] * V_2[0])
    k = (V_1[0] * V_2[1]) - (V_1[1] * V_2[0])

    j = -j

    cross = np.array([i, j, k])
    return cross

def my_interesting_build_planar_cloud(P, N, numOfPoints, timeout=10.0):
    """ 
        build a random planar point cloud 
        probability of randoming finding {numOfPoints} valid points on plane is highly unlikely with this technique,
        although it is technically mathematically sound
    """
    start = time.perf_counter()
    vectorsOnPlane = 0

    while vectorsOnPlane != 2:

        if time.perf_counter() - start >= timeout:
            raise TimeoutError(f"could not find two independent vectors within {timeout} secs")

        testP = np.random.uniform(-10, 10, size=3)
        testV = build_vector_w_distance(P, testP)

        if two_norm(testV) < 1e-12:
            continue

        if abs(dot(N, testV)) < 1e-6:

            if vectorsOnPlane == 0:
                v_1 = testV.copy()
                vectorsOnPlane = 1

            else:
                cross = my_cross(v_1, testV)

                if two_norm(cross) > (1e-10 * two_norm(v_1) * two_norm(testV)):
                    v_2 = testV.copy()
                    vectorsOnPlane = 2

    pointsOnPlane = np.zeros((numOfPoints, 3))

    for i in range(numOfPoints):
        a = np.random.uniform(-1, 1)
        b = np.random.uniform(-1, 1)

        pointsOnPlane[i] = P + a * v_1 + b * v_2

    return pointsOnPlane

def build_planar_cloud(P: np.ndarray, N: np.ndarray, numOfPoints: int) -> np.ndarray:

    P = np.asarray(P, dtype=float)
    N = np.asarray(N, dtype=float)

    # validate inputs
    if P.shape != (3,) or N.shape != (3,):
        raise ValueError("P and N must each have shape (3,)")

    if not np.all(np.isfinite(P)) or not np.all(np.isfinite(N)):
        raise ValueError("P and N must contain finite values")

    if not isinstance(numOfPoints, (int, np.integer)) or numOfPoints < 0:
        raise ValueError("numOfPoints must be a nonnegative integer")

    norm_N = np.linalg.norm(N) # length of normal vector

    if norm_N == 0:
        raise ValueError("The normal vector N cannot be zero")

    N = N / norm_N # unit

    # pick an axis that is not parallel to N
    axis = np.zeros(3) # [0.0, 0.0, 0.0]
    axis[np.argmin(np.abs(N))] = 1.0 # smallest absolute component in the normal vector 

    # first unit vector in the plane
    v_1 = np.cross(N, axis)
    v_1 = v_1 / np.linalg.norm(v_1)

    # second unit vector in the plane, perpendicular to v_1
    v_2 = np.cross(N, v_1)

    # generate random coefficients for all points at once
    coefficients = np.random.uniform(-10, 10, size=(numOfPoints, 2))

    pointsOnPlane = (P + coefficients[:, 0, None] * v_1 + coefficients[:, 1, None] * v_2)

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

def read_cloud(file):

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

def spectralPCA(filename, size):

    # read the file
    myFile = read_cloud(filename)

    # get the mean
    mean = get_mean(myFile)

    # move point cloud to origin
    cloudToOrigin = move_point_cloud(myFile, mean, size)

    # get the covariance matrix
    covariance = numpy_covariance(cloudToOrigin, 3)

    # get the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = numpy_eigen(covariance)

    # get eigenvector corresponding to smallest eigenvalue
    normal = eigenvectors[:, 0]

    #return mean and normal
    return mean, normal

# PCA w/ spectral analysis w/ random planar cloud
P = np.array([4.0, 2.0, 0.0])
N = np.array([7.5, 2.1, 9.4])
numOfPoints = 100

start1 = time.perf_counter()
planarPointCloud1 = build_planar_cloud(P, N, numOfPoints)
end1 = time.perf_counter()

create_ply("planar_plane.ply", planarPointCloud1)
meanSpectral, normalSpectral = spectralPCA("planar_plane.ply", 100)
print("mean from spectral w/ coplanar test:", meanSpectral)
print("normal from spectral w/ coplanar test:", normalSpectral, "\n")

try:
    start2 = time.perf_counter()
    planarPointCloud2 = my_interesting_build_planar_cloud(P, N, numOfPoints)
    end2 = time.perf_counter()

    print("generating planar point cloud w/ random until dot = 0 between norm and vector:", end2 - start2, "sec\n")

except TimeoutError as error:
    print("random method failed", error)

print("generating planar point cloud w/ numpy and cross:", end1 - start1, "sec\n")

"""
The reconstructed plane represents the best-fit plane that PCA estimates from your 3D point cloud. 
It tells you the orientation and position of the flat surface that best describes the points.

In your code, you already know the original plane. 
But imagine you didn't know its normal vector and only had the 100 XYZ points. 
PCA lets you estimate the plane from those points alone.

In short: you center the points at the origin to calculate the normal 
then use the saved mean to locate the reconstructed plane in its original position

"""