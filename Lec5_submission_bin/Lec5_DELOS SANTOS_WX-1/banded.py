from numpy import copy

def banded(Aa,va,up,down):

    #copy inputs and determine size of system
    A = copy(Aa)
    v = copy(va)
    N = len(v)

    #gaussian elim.
    for m in range(N):
        
        div = A[up,m] #normalization factor
        v[m] /= div   #update the vector first
        for k in range(1,down+1):
            if m+k<N:
                v[m+k] -= A[up+k,m]*v[m]

        for i in range(up): #normalize the pivot row of A and subtract from lower ones
            j = m + up - i
            if j<N:
                A[i,j] /= div
                for k in range(1,down+1):
                    A[i+k,j] -= A[up+k,m]*A[i,j]

    for m in range(N-2,-1,-1): #back-subst.
        for i in range(up):
            j = m + up - i
            if j<N:
                v[m] -= A[i,j]*v[j]

    return v
        