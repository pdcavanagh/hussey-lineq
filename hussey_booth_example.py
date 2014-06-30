#!/Users/pdcavanagh/anaconda/bin/python

#******************************************************************************
# Trial 2
# A = x^-1 * B
# residual calculated: r = (A*E) - B
#
# E: estimate of actual values
# B: APXS bulk composition
# A: oxide wt % compositions of phases
#
# Iterative result using a residual calculation based on Booth
# Objective to minimize grand residual: R^2 = sum from i=1 to q (r(i)^2)
#******************************************************************************
import numpy as np

debug=1
results=1

a = np.matrix([[14,5,5],[10,100,160],[10,100,10]])
b = np.matrix([[11.3,6.8,6.8],
               [49.0,100,118],
               [19.0,55,28]])
c = np.matrix([[0.07, 0.10, 0.20],
               [0.20, 0.50, 0.30],
               [0.20, 0.20, 0.60]])
 
b_mod = np.matrix([[11.86,6.5,7.1],
               [51.4,95.0,112.0],
               [18.1,57.5,26.6]])

x = np.matrix(np.zeros(3)).getT()
e=x

# Basic inverse test case, assuming that the experimental values are accurate
if debug:
    print a.getI()*b
    print np.round(a.getI()*b_mod,2)


thresh=0.001
r = a*e - b_mod
if debug:
    print r

a_prime = a.getT()

# Calculate the residual over all of the phases property measurements
for n,i in enumerate(r.getT()):
    r_sqr=1 # Global residual value
    if results:
        print ('\n***************Phase %s***************\n' % n)
        print i
    i=i.getT()
    count=0 
    while r_sqr>thresh:
        count=count+1
        correct_2 = a.getT()*i
        correct_3 = a*correct_2
        correct_4 = i.getT()*correct_3
        correct_denom = sum(np.square(correct_3))
        correct_5 = correct_4/correct_denom  
        correct_term = np.multiply(correct_5,correct_2)
        e=e-correct_term #update the approximation vector 
        i = a*e - b_mod[:,n]
        
        if debug: 
            print ('correct_2:\n%s' % correct_2)
            print ('correct_3:\n%s' % correct_3)
            print ('correct_4:\n%s' % correct_4)
            print ('correct_denom: %s' % correct_denom)
            print ('correct_4/correct_denom: %s' % correct_5)
            print ('correct_term: %s' % correct_term)

        r_sqr=0 #reset the r_sqr value
        for resid in i:
            if debug:
                print ('resid: %s' % resid) 
                print ('r_sqr: %s' % r_sqr)
            r_sqr = resid*resid + r_sqr
    if results: 
        print ('Global residual value (R^2): {0:5.2f}'.format(float(r_sqr)))
        print ('The number of iterations to achieve <0.01 R^2: %s' % count)
        print ('The estimated phase abundaces are:\n%s' % np.round(e,3))
