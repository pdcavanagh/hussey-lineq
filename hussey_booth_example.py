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

debug=0

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

r_sqr=0 # Global residual value
r = a*e - b[:,0]
if debug:
    print r

a_prime = a.getT()

# Calculate the residual over all of the phases property measurements
for i in r:
    if debug:
        print r_sqr 
    r_sqr = i*i + r_sqr

if debug: print 'Global residual value (R^2): {0:5.2f}'.format(float(r_sqr))

correct_2 = a.getT()*r
if debug: print ('correct_2:\n%s' % correct_2)
correct_3 = a*correct_2
if debug: print 'correct_3:'
if debug: print correct_3
correct_4 = r.getT()*correct_3
if debug: print r
if debug: print 'correct_4:'
if debug: print correct_4

correct_denom = sum(np.square(correct_3))
if debug: print ('correct_denom: %s' % correct_denom)
correct_5 = correct_4/correct_denom  
if debug: print 'correct_4/correct_denom: {}'.format(correct_5)

correct_term = np.multiply(correct_5,correct_2)
if debug: print correct_term
print e-correct_term

