import numpy as np 
import math

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
jon_the_bunny = readFile("bun_zipper.ply")
extracting_jon_poor_jon = getPoints(jon_the_bunny, 30000, 5)
pusherman1("bunny.X", extracting_jon_poor_jon)

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
angry_jon = meanBunny(jon_the_bunny)
jons_variance = variance(jon_the_bunny, angry_jon)
jon_likes_to_deviate = stdDev(jons_variance)
print("jon the bunny is a peculiar little fellow, his mean is:", angry_jon)
print("jon's variance is:", jons_variance)
print("jon's standard deviation is:", jon_likes_to_deviate)














                
                


            

            


                





