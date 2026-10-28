import numpy as np

def main(A: np.ndarray):
    """
    -> matrix a transformation of space

        A = UEV^T
        V = col of V tell you directions before the transformation the stretch or compress happens scaled by according to corresponding column of E
        E = stretching or compressing factor along col direction, singular values along its diagonal
        U = final directions of the transformed principal axes, can rotate or reflect the principal axes into their final positions

        E -> sqrt(eigenvalues) of A^TA along diagonal from largest to smallest, left to right
        V -> each col an eigenvector from each eigenvalue of A^TA from largest to smallest, left to right
        U -> u_i = 1/sqrt(eigenvalues_i)(A(eigenvector_i))

        V_T -> if you want to reach x in standard basis coordinate system, you have to move V_Tx along V coordinate system
    """
    A = np.asarray(A, dtype=float)
    m, n = A.shape
    k = min(m, n)

    E = np.zeros((m, n))
    V = np.zeros((n, n))
    U = np.zeros((m, m))

    A_TA = A.T @ A

    # check ||Ax||^2 = (x^T)(A^T)(A)(x)
    x = np.random.rand(n)
    left = np.linalg.norm(A @ x) ** 2
    right = x.T @ A.T @ A @ x
    if not np.allclose(left, right):
        print("error: w/ tester vector x; ||Ax||^2 != (x^T)(A^T)(A)(x) in SVD\n")
        return -1
    
    eigenvalues, V_eigen = np.linalg.eigh(A_TA)

    # reverse the order of eigenvalues and eigenvectors of A^TA from largest to smallest to smallest to largest
    eigenvalues = eigenvalues[::-1]
    V_eigen = V_eigen[:, ::-1]
    # eigenvalues = np.maximum(eigenvalues, 0)

    # get singular values, sqrt(eigenvalue)
    singular_values = np.zeros(A.shape[1])
    for i in range(len(eigenvalues)):
        singular_values[i] = eigenvalues[i] ** .5

    # build E, main diagonal of a rectangular matrix has min(number of rows, number of columns) positions
    for i in range(min(E.shape)):
        for j in range(min(E.shape)):
            if i == j:
                E[i][j] = singular_values[i]

    # build V
    for j in range(V.shape[1]):
        V[:, j] = V_eigen[:, j]

    epsilon = (np.finfo(A.dtype).eps * max(m, n) * singular_values[0])
    vector_epsilon = 1e-10

    # build U
    for i in range(k):
        if singular_values[i] > epsilon:
            u_i = (1 / singular_values[i]) * A @ V_eigen[:, i]
            U[:, i] = u_i

    # fill any missing columns of U
    for j in range(m):

        # if this column already contains a vector, skip it
        if np.linalg.norm(U[:, j]) > vector_epsilon:
            continue

        # try standard basis vectors until one works
        for basis_index in range(m):
            candidate = np.zeros(m)
            candidate[basis_index] = 1.0 # creates [1, 0, 0], [0, 1, 0], or [0, 0, 1]

            # remove components pointing in directions already occupied by columns of U
            for i in range(m):
                if np.linalg.norm(U[:, i]) > vector_epsilon:
                    candidate -= np.dot(U[:, i], candidate) * U[:, i] # make candidate orthogonal to the vectors already in U

            # check the candidate after all projections
            norm = np.linalg.norm(candidate)

            if norm > vector_epsilon:
                U[:, j] = candidate / norm
                break

    if not np.allclose(U.T @ U, np.eye(m)):
        print("error: sub-matrix U cols are not orthonormal in SVD\n")
        return -1

    if not np.allclose(V.T @ V, np.eye(n)):
        print("error: sub-matrix V cols are not orthonormal in SVD\n")
        return -1

    if not np.allclose(A, U @ E @ V.T):
        print("error: SVD reconstruction failed\n")
        return -1
    
    return U, E, V.T

# assumes all eigenvalues are positive
# orthogonality checks need to pass
# A = UEV^T passes
# robust replacement - np.linalg.svd()








