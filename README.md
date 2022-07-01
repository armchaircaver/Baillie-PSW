# Baillie-PSW
Python implementation of the Baillie-PSW primality test

Implementation of the Baillie-PSW primality checking algorithm, from the following papers:

Lucas Pseudoprimes
Robert Baillie and Samuel S. Wagstaff, Jr.
Mathematics of Computation Vol. 35, No. 152 (Oct., 1980), pp. 1391-1417
https://www.jstor.org/stable/2006406

Strengthening the Baillie-PSW primality test,
Robert Baillie, Andrew Fiori and Samuel S. Wagstaff, Jr.
Math. Comp. 90 (2021), 1931-1955
https://arxiv.org/pdf/2006.14425v1.pdf,  https://homes.cerias.purdue.edu/~ssw/bfw.pdf

- baillie_psw.py contains the implemetation of the algorithm

- baillie PSW test suite.py  has a suite of tests of the algorithm, using sieve and segmented sieve to generate test numbers

- Compare pseudoprimes.py is a test suite to compare the performance and correctness of the standard and strengthened Lucas algorithms, 
using known pseudoprimes obtained from Pseudoprime Statistics, Tables, and Data (Fermat, Miller-Rabin, Lucas, Fibonacci, Pell, Frobenius, Baillie-PSW), 
Dana Jacobsen, 31 March 2020, http://ntheory.org/pseudoprimes.html 

My implementation of the strengthened algorithm is slightly slower than the standard algorithm, so the standard algorithm is used in baillie_psw
