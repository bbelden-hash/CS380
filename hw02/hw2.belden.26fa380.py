import numpy as np 

def readFile(file):
    try:
        with open(file, "r") as file:
            bunny = file.readlines()
        return bunny
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

            print("LINE:",vertex)
            print("TEXT:", repr(lines))
            print("ELEMENTS", elements)
            print()

            for n in range(3):
                component = float(elements[n])
                points[n, col] = component

            col += 1

            if col == offset:
                break

        vertex += 1
    return points

def getATOM(file, start, offset):

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

    atoms = np.zeros((80, offset), dtype = object)

    for line in file:

        if "PROTEIN ATOMS" in line:
            int(noATOM) = line.split(":")[1].strip() # parsing until line has 'PROTEIN ATOMS' in which line is split -> index 0 everything before ':', index 1 everything after ':'
            print("number of ATOM: ", noATOM)
        
        if "HETEROGEN ATOMS" in line:
            int(noHETATM) = line.split(":")[1].strip()
            print("number of HETATM: ", noHETATM)

        











def pusherman1(file, subMatrix):
    with open(file, "w") as file:
        
        file.write("Random Points, Stanford Bunny\n\n")
        for row in subMatrix:
            for coordinate in row:
                file.write(str(coordinate) + " ")
            file.write("\n")


# function calls

# Stanford bunny, extract x, y, z coordinates from bun_zipper.ply file, point location within a base + offset specified by user
# place data into a (3, n) numpy matrix -> first row-x, second row-y, third row-z / first col-p1, second col-p2, ..., n col-p(n)
jon_the_bunny = readFile("bun_zipper.ply")
extracting_jon_poor_jon = getPoints(jon_the_bunny, 30000, 5)
pusherman1("bunny.X", extracting_jon_poor_jon)

# Write a function to read the atoms of a PDB format into a column-major data matrix,
# using the ATOM and HETATM lines (lec09)
# Use this function to extract the 5 points of 1GCN.pdb starting at 100th point and write this submatrix to ’1GCN.X'

# REMARK   3   PROTEIN ATOMS            : 246                                     
# REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
# REMARK   3   HETEROGEN ATOMS          : 0                                       
# REMARK   3   SOLVENT ATOMS            : 0          
billy_the_glucagon = readFile("1GCN.pdb")
grabbing_billys_ATOMS = grabATOM(billy_the_glucagon, 100, 5)








                
                


            

            


                





