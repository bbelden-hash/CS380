import numpy as np

# all examples for 3x3
def main(A: np.ndarray):
   """
    A = LU
    L = ones on the diagonal w/ multiplicites under it
    U = like a row echelon form without 1s on the diagonal

    got to be square, QR for rectangles
    """
    if (A.shape[0] != A.shape[1]):
        print("error: LU decomp - got to be square\n")
        return -1

    n = A.shape[0]

    L = np.zeros((n, n))
    for i in range(len(A)):
        for j in range(len(A)):
            if i == j:
                L[i, j] = 1
    
    U = A.copy()

    """
        pivots - U[0, 0], U[1, 1], ..., U[n, n]

        col 1 - > a[1, 1]
        m_21 = a_21 / a_11, m_31 = a_31 / a_11, m_41 = a_41 / a_11, m_51 = a_51 / a_11, ..., 
        start_n = 2
        m_n1 = a_n1 / a_11, m_(n+1)1 = a_(n+1)1 / a_11, m_(n+2)1 = m_(n+2)1 / a_11, ...
        m_n(col) = a_n(col) / a_n(col)(col)

        col 2 -> a[2, 2]
        m_32 = a_32 / a_22, m_42 = a_52 - a_22, ...
        start_n = 3
        m_n2 = a_n2 / a_22, m_(n+1)2 = a_(n+1)2 / a_22, ...
        m_n(col) = a_n(col) / a_n(col)(col)
    """
    num_rows = len(A)
    num_cols = len(A[0])

    for col in range(num_cols):
        for row in range(num_rows):
            if col == (num_cols - 1):
                return L, U 

            if row == col:
                # pivot = U[row, col]

                i = row
                k = row + 1
                pivot_row = U[i] 
                while k < U.shape[0]:

                    L[k][col] = U[k][col] / U[col][col]
                    multiplier = L[k][col]

                    m_row = pivot_row * multiplier
                    new_row = U[k] - m_row

                    U[k] = new_row

                    k += 1


    return L, U

# idk if this works yet but would be cool -> I think the logic is sound
# see LUDecomp.pdf for visual example w/ iterations

                



    


