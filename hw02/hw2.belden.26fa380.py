import numpy as np 
import math
import random
import time

def readFile(file):
    try:
        with open(file, "r") as file:
            fileName = file.readlines()
        return fileName
    except FileNotFoundError:
        print("error: file is not valid or non-existent")
        return None

def getPoints(file, start, offset):
    
    header = True
    vertex = 0
    col = 0

    points = np.zeros((3, offset), dtype = np.float64)

    for lines in file:

        if header:
            if lines.strip() == "end_header":
                header = False
            continue

        if start <= vertex < (start + offset):

            elements = lines.split()

            for n in range(3):
                component = float(elements[n])
                points[n, col] = component

            col += 1

            if col == offset:
                break

        vertex += 1
    return points

def getATOM(file, start, offset, name):

    # for ATOM:
    # col 1-4: "ATOM"
    # col 7-11: Atom serial number
    # col 13-16: Atom name
    # col 17: Alternate location indicator
    # col 18-20: Residue name
    # col 22: Chain identifier
    # col 23-26: Residue sequence number
    # col 27: Code for insertions of residues
    # col 31-38: X orthogonal A coordinate
    # col 39-46: Y orthogonal A coordinate
    # col 47-54: Z orthogonal A coordinate
    # col 55-60: Occupancy
    # col 61-66: Temperature factor
    # col 73-76: Segment identifier
    # col 77-78: Element symbol
    # col 79-80: Charge

    # for HETATM:
    # col 1-6: "HETATM"
    # col 7-80: same as ATOM records

    noATOM = 0
    noHETATM = 0
    currATOM = 0
    col = 0

    atoms = np.zeros((15, offset), dtype = object)

    for line in file:

        if "PROTEIN ATOMS" in line:
            noATOM = int(line.split(":")[1].strip()) # parsing until line has 'PROTEIN ATOMS' in which line is split -> index 0 everything before ':', index 1 everything after ':'
            
            if noATOM == 0:
                print(f"situation: nothing to parse for 'ATOM; number of 'ATOM' = {noATOM} in {name}")
            else:
                print(f"PROTEIN ATOMS in {name}: {noATOM}")
        
        if "HETEROGEN ATOMS" in line:
            noHETATM = int(line.split(":")[1].strip())
            
            if noHETATM == 0 and noATOM == 0:
                print(f"situation: nothing to parse for {name}; ATOM and HETATM = 0")
                return None
            elif noHETATM == 0:
                print(f"situation: nothing to parse for 'HETATM'; number of 'HETATM' = {noHETATM} in {name}")
            else:
                print(f"HETEROGEN ATOMS in {name}: {noHETATM}")

        matter = line

        if matter[0:6].strip() == "ATOM":

            if start <= currATOM < (start + offset):

                serialNumber = int(matter[6:11].strip())
                atomName = matter[12:16].strip()
                altLoc = matter[16:17].strip()
                resName = matter[17:20].strip()
                chainId = matter[21:22].strip()
                resSeq = int(matter[22:26].strip())
                iCode = matter[26:27].strip()
                x = float(matter[30:38].strip())
                y = float(matter[38:46].strip())
                z = float(matter[46:54].strip())
                occupancy = float(matter[54:60].strip())
                tempFactor = float(matter[60:66].strip())
                segmentID = matter[72:76].strip()
                element = matter[76:78].strip()
                charge = matter[78:80].strip()

                for n in range(15):

                    if n == 0:
                        atoms[n, col] = serialNumber
                    elif n == 1:
                        atoms[n, col] = atomName
                    elif n == 2:
                        atoms[n, col] = altLoc
                    elif n == 3:
                        atoms[n, col] = resName
                    elif n == 4:
                        atoms[n, col] = chainId
                    elif n == 5:
                        atoms[n, col] = resSeq
                    elif n == 6:
                        atoms[n, col] = iCode
                    elif n == 7:
                        atoms[n, col] = x 
                    elif n == 8:
                        atoms[n, col] = y 
                    elif n == 9:
                        atoms[n, col] = z 
                    elif n == 10:
                        atoms[n, col] = occupancy
                    elif n == 11:
                        atoms[n, col] = tempFactor
                    elif n == 12:
                        atoms[n, col] = segmentID
                    elif n == 13:
                        atoms[n, col] = element
                    elif n == 14:
                        atoms[n, col] = charge

                col += 1
                if col == offset:
                    break
            
            currATOM += 1

    return atoms

def meanBunny(file):

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

    means = np.array([meanX, meanY, meanZ], dtype = np.float64)
    return means

def stdDev(variance):

    xVar = variance[0]
    yVar = variance[1]
    zVar = variance[2]

    xStd = math.sqrt(xVar)
    yStd = math.sqrt(yVar)
    zStd = math.sqrt(zVar)

    stdDevs = np.array([xStd, yStd, zStd], dtype = np.float64)
    return stdDevs

def variance(file, mean):

    forX = 0
    forY = 0
    forZ = 0
    num = 0

    xMean = mean[0]
    yMean = mean[1]
    zMean = mean[2]

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

            tempX = xMean - x 
            tempY = yMean - y 
            tempZ = zMean - z 

            squaredX = tempX ** 2
            squaredY = tempY ** 2
            squaredZ = tempZ ** 2

            forX += squaredX
            forY += squaredY
            forZ += squaredZ 

            num += 1

    varX = forX / num 
    varY = forY / num 
    varZ = forZ / num 

    variance = np.array([varX, varY, varZ], dtype = np.float64)
    return variance

def matMul(A, B):

    rowsA = 0
    for a in A[:, 0]:
        rowsA += 1

    colsB = 0
    for b in B[0, :]:
        colsB += 1

    colsA = 0
    for a in A[0, :]:
        colsA += 1
    
    rowsB = 0
    for b in B[:, 0]:
        rowsB += 1

    if colsA == rowsB:
        C = np.zeros((rowsA, colsB), dtype = np.float64)
    else:
        print("error: to be able to matrix multiply; columns of A and rows of B must equal, unable to proceed")
        return None

    # n is which row of A we are currently in, for a 3x2 ... (A[n, :])
    #[ 1 2 ] n = 0
    #[ 3 4 ] n = 1
    #[ 7 8 ] n = 3
    # ...
    for n in range(len(A)):
        # m is a row of B working with, B[0, :] =
        # [ 5 6 ] is has len(2), so m = 0, 1
        for m in range(len(B[0, :])):
            sum = 0

            # for 3x2, len(2), ...
            # allows for traveral down row for A and down column for B with corresponding row and column numbers
            for i in range(len(A[n, :])):
                sum += A[n, i] * B[i, m]
           
            C[n, m] = sum
    
    return C 

def npMatmul2(A, B):

    C = np.matmul(A, B)
    return C

def npMatmul3(A, B, C):

    D = np.matmul(A, B)
    E = np.matmul(C, D)
    return E

# create a function that takes a point cloud and calculates a matrix that transforms the point cloud into a normalized form 
# move the point cloud so its center is at the origin, scale it so its size is normalized
def movePointCloud_Normalize(file, theMeans, size):

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
                
    movingPoints = matMul(transform, points)

    with open("mean_normalization.X", "w") as dest:
        for col in (movingPoints.T):
            dest.write(f"{col[0]} {col[1]} {col[2]}\n")

    return transform

def compareAndContrast():

    A = np.random.rand(32, 32)
    B = np.random.rand(64, 64)
    C = np.random.rand(128, 128)
    D = np.random.rand(512, 512)

    start1a = time.time()
    matMul(A, A)
    end1a = time.time()
    print("My 32x32 matmul time:", end1a - start1a, "sec")

    start1b = time.time()
    npMatmul2(A, A)
    end1b = time.time()
    print("NumPy 32x32 matmul time:", end1b - start1b, "sec")

    start2a = time.time()
    matMul(B, B)
    end2a = time.time()
    print("My 64x64 matmul time:", end2a - start2a, "sec")

    start2b = time.time()
    npMatmul2(B, B)
    end2b = time.time()
    print("NumPy 64x64 matmul time:", end2b - start2b, "sec")

    start3a = time.time()
    matMul(C, C)
    end3a = time.time()
    print("My 128x128 matmul time:", end3a - start3a, "sec")

    start3b = time.time()
    npMatmul2(C, C)
    end3b = time.time()
    print("NumPy 128x128 matmul time:", end3b - start3b, "sec")

    start4a = time.time()
    matMul(D, D)
    end4a = time.time()
    print("My 512x512 matmul time:", end4a - start4a, "sec")

    start4b = time.time()
    npMatmul2(D, D)
    end4b = time.time()
    print("NumPy 512x512 matmul time:", end4b - start4b, "sec")

    return None

def logFib(F, fibNo):

    # [1 1][F_n  ]     =  [F_n+1] 
    # [1 0][F_n-1]        [F_n  ]
    # F_n+1 = F_n + F_n-1
    #  ^ moves us forward one fibonacci number

    I = np.array([
        [1, 0],
        [0, 1]
    ])

    F_0 = 0
    F_1 = 1

    currentState = np.array([
        [F_0], # F_n
        [F_1]  # F_n-1
    ])

    # one multiplication of fibTransform and currentState gives [F_2, F_1]^T ... -> we want to move our current state forward one fibonacci position
    # first output needs to be (F_n + F_n-1) = F_n+1, so row one of transform needs to be [1 1]
    # second output needs to be F_n, so row two of transform needs to be [1, 0]
    # exponential multiplcation: 2^16 you could multiple 2 15 times (2*2*2*2*2*2*...) or 2^2 = 4, 4^2 = 16, 16^2 = 256, 256^2 = 2^16 -> sqrt(16) = 4 -> 4 operations, we can do same with the matrices
    
    # one multiplication of currentState and fibtransform moves us forward one step, multiply fibtransform twice?
    # [1 1][1 1] = [2 1] F_3 = 2, F_2 = 1, F_1 = 1 [F_3 F_2]
    # [1 0][1 0]   [1 1]                           [F_2 F_1]

    # fibTransform^n = [F_n+1 F_n]
    #                  [F_n F_n-1]

    # fibTransform^8 = [F_9 F_8]
    #                  [F_8 F_7]

    if fibNo == 0:
        return I 
    
    halfFib = logFib(F, fibNo // 2) # recursive call in which the stack dives deep until the base case is reached and halfFib == I

    if (fibNo % 2) == 0:
        return npMatmul2(halfFib, halfFib)
    elif (fibNo % 2) != 0:
        return npMatmul3(F, halfFib, halfFib)

def savePDB(file):

    with open(file, 'r') as source, open('1GCN.X', 'w') as destination:
        destination.write("Glucagon, 1GCN.pdb; Point Cloud ->\n\n")

        for line in source:

            sourceLine = line.split()
            if sourceLine[0] == "ATOM":

                for word in sourceLine:
                    destination.write(word)
                    destination.write(" ")
            
                destination.write("\n")

def printPDB(subMatrix):
    print("\n")
    print("ATOM sub-matrix:")

    # Calculate column widths using the transpose of the submatrix
    col_widths = [max(len(str(item)) for item in col) + 2 for col in subMatrix.T]

    # Print each row without adding an extra index counter on the left
    for row in subMatrix:
        print("".join(f"{str(val):<{col_widths[i]}}" for i, val in enumerate(row)))
       
def pusherman1(file, subMatrix):
    with open(file, "w") as file:
        
        file.write("Random Points, Stanford Bunny\n\n")
        for row in subMatrix:
            for coordinate in row:
                file.write(str(coordinate) + " ")
            file.write("\n")


# function calls

# Data Matrix
# Stanford bunny, extract x, y, z coordinates from bun_zipper.ply file, point location within a base + offset specified by user
# place data into a (3, n) numpy matrix -> first row-x, second row-y, third row-z / first col-p1, second col-p2, ..., n col-p(n)
jeff_the_bunny = readFile("bun_zipper.ply")
extracting_jeff_poor_jon = getPoints(jeff_the_bunny, 30000, 5)
pusherman1("bunny.X", extracting_jeff_poor_jon)

# Write a function to read the atoms of a PDB format into a column-major data matrix,
# using the ATOM and HETATM lines (lec09)
# Use this function to extract the 5 points of 1GCN.pdb starting at 100th point and write this submatrix to ’1GCN.X' 
billy_the_glucagon = readFile("1GCN.pdb")
savePDB("1GCN.pdb")
grabbing_billys_ATOMS = getATOM(billy_the_glucagon, 99, 5, "1GCN.pdb")
printPDB(grabbing_billys_ATOMS)

# Mean Normalization
# finding the mean of a point cloud by adding up all of the x coordinates (x1 + x2 + x3 + ... + x(n)), y coordinates (y1 + y2 + ... y(n)), and z coordinates (z1 + z2 + ... z(n)) ...
# divide each sum by the total number of x, y, z triples in the point cloud
# calculate the variance from a point cloud by subtracting each x, y, z coordinate by their corresponding mean, squaring this value, adding all corresponding component values up, and dividing by the total number of coordinates
# mean normalize the points by moving the 'center' of the point cloud to the origin ... this is done by transforming the points by a 4x4 matrix consisting of the means of the x, y, and z coordinates
angry_jeff = meanBunny(jeff_the_bunny)
jeffs_variance = variance(jeff_the_bunny, angry_jeff)
jeff_likes_to_deviate = stdDev(jeffs_variance)
jeff_moves_towards_origin = movePointCloud_Normalize(jeff_the_bunny, angry_jeff, 35947)
print("jeff the bunny is a peculiar little fellow, his mean is:", angry_jeff)
print("jeff's variance is:", jeffs_variance)
print("jeff's standard deviation is:", jeff_likes_to_deviate)
print("jeff is on the move towards the origin, how is he going to move?:\n", jeff_moves_towards_origin)
print("\n")

# Matrix Multiplication
# the tortoise and the hare
# compareAndContrast()
# print("\n")

# Logarithmic Fibonacci Solution, this was a process but fun ...
fibTransform = np.array([
        [1, 1],
        [1, 0]
    ])
fibNo = 1000
start = time.time()
fibMatrix = logFib(fibTransform, fibNo)
end = time.time()
fib = fibMatrix[0, 1]
print("fibonacci logarithmic, F_", fibNo, ":", fib)
print("Time taken to compute F_", fibNo, "is", end - start, "sec")
















                
                


            

            


                





