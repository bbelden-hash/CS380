import numpy as np 

def readFile(file):
    try:
        with open(file, "r") as file:
            bunny = file.readlines()
        return bunny
    except FileNotFoundError:
        print("error: file is not valid or non-existent")
        return None

def getPointsPLY(file, start, offset):
    
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

def pusherman(file, subMatrix):
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
extracting_jon_poor_jon = getPointsPLY(jon_the_bunny, 30000, 5)
pusherman("bunny.X", extracting_jon_poor_jon)



                
                


            

            


                





