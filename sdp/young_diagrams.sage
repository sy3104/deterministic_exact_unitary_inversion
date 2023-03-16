# -------------------------------------------------------------------------
#
# File: young_diagrams.sage
#
# Description :
# Codes to produce auxiliary values related to Young diagrams and
# representation theory of the special unitary group SU(d) and symmetric
# group S_n. Calculated values are saved in "young_diagrams.m".
#
# -------------------------------------------------------------------------

import numpy as np

N_max = 7   #Maximal number of boxes 

num_partitions_max = len(Partitions(N_max).list())

f = open("young_diagrams.m", "w")

for n in range(1,N_max+1):
    partitions_list = Partitions(n).list()
    prev_partitions_list = Partitions(n-1).list()
    next_partitions_list = Partitions(n+1).list()
    permutation = [n] + list(range(1,n))
    f.write("num_young_diagrams{"+str(n)+"} = "+str(len(partitions_list))+";\n\n")
    for index_mu in range(len(partitions_list)):
        mu = partitions_list[index_mu]
        
        # Calculate the dimension of irrep of S_n and SU(d) for 'd'
        d_mu = mu.dimension()
        m_mu = ""
        for (i,j) in mu.cells():
            m_mu = m_mu + " * (d+(" + str(j-i) + "))"
        m_mu = m_mu + " / " + str(factorial(n)/d_mu)
        m_mu = m_mu[3:]
        f.write("depth{"+str(n)+"}{"+str(index_mu+1)+"} = "+str(len(mu))+";\n")
        f.write("dim{"+str(n)+"}{"+str(index_mu+1)+"} = "+str(d_mu)+";\n")
        f.write("mult{"+str(n)+"}{"+str(index_mu+1)+"} = "+m_mu+";\n")
        
        # Calculate the representation matrix for the permutation in the irrep mu
        sgr = SymmetricGroupRepresentation(mu, "orthogonal")
        permutation_matrix = numerical_approx(sgr(permutation))
        permutation_matrix = matlab.sage2matlab_matrix_string(permutation_matrix)
        f.write("perm{"+str(n)+"}{"+str(index_mu+1)+"} = "+permutation_matrix+";\n")
        
        # Calculate the matrix X^alpha_mu
        st_list_mu = StandardTableaux(mu).list()
        st_list_mu = list(reversed(st_list_mu))
        
        for index_alpha in range(len(prev_partitions_list)):
            alpha = prev_partitions_list[index_alpha]
            st_list_alpha = StandardTableaux(alpha).list()
            st_list_alpha = list(reversed(st_list_alpha))
            X = np.zeros((len(st_list_alpha), len(st_list_mu)))
            for a in range(len(st_list_alpha)):
                st_a = st_list_alpha[a]
                for st_nu in st_a.up_list():
                    if st_nu.shape() == mu:
                        j = st_list_mu.index(st_nu)
                        X[a,j] = 1
            X = matlab.sage2matlab_matrix_string(matrix(X))
            f.write("X{"+str(n)+"}{"+str(index_alpha+1)+"}{"+str(index_mu+1)+"} = "+X+";\n")
        
        f.write("\n")
            
f.close()