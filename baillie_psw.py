"""
Implementation of the Baillie-PSW primality checking algorithm, from the following papers:

Lucas Pseudoprimes
Robert Baillie and Samuel S. Wagstaff, Jr.
Mathematics of Computation Vol. 35, No. 152 (Oct., 1980), pp. 1391-1417
https://www.jstor.org/stable/2006406

Strengthening the Baillie-PSW primality test,
Robert Baillie, Andrew Fiori and Samuel S. Wagstaff, Jr.
Math. Comp. 90 (2021), 1931-1955
https://arxiv.org/pdf/2006.14425v1.pdf,  https://homes.cerias.purdue.edu/~ssw/bfw.pdf

"""


# from https://stackoverflow.com/questions/44531084/python-check-if-a-number-is-a-square
def is_square(n):
    if n < 0:
        return False
    if n == 0:
        return True
    x, y = 1, n
    while x + 1 < y:
        mid = (x+y)//2
        if mid**2 < n:
            x = mid
        else:
            y = mid
    return n == x**2 or n == (x+1)**2
#-------------------------------------------------------------------------------

# adapted from https://rosettacode.org/wiki/Jacobi_symbol#Python
def jacobi(a, n):
  a %= n
  result = 1
  while a != 0:
      while a % 2 == 0:
          a //= 2
          n_mod_8 = n % 8
          if n_mod_8 in (3, 5):
              result = -result
      a, n = n, a
      if a % 4 == 3 and n % 4 == 3:
          result = -result
      a %= n
  if n == 1:
      return result
  else:
      return 0
#-------------------------------------------------------------------------------
def miller_rabin_base_2(n):
  # adapted from https://rosettacode.org/wiki/Miller%E2%80%93Rabin_primality_test#Python
  s = 0
  d = n-1
  while d%2==0:
      d>>=1
      s+=1

  def trial_composite(a):
    if pow(a, d, n) == 1:
      return False
    for i in range(s):
      if pow(a, 2**i * d, n) == n-1:
        return False
    return True  

  return not trial_composite(2)
#-------------------------------------------------------------------------------

def D_chooser(n):
  #Choose a D value suitable for the Baillie-PSW test
  D = 5
  j = jacobi(D, n)

  while j > 0:
    D += 2 if D > 0 else -2
    D *= -1

    if D==-15 :
      # check for a square
      if is_square(n):
        # The value of D isn't 0, but we are just communicating
        # that we have found a square
        return (0,0) 

    j = jacobi(D, n)
  return (D,j)
#-------------------------------------------------------------------------------

def div2mod(x,n):
  # divide by 2 modulo n
  # assumes n is odd
  if x & 1:
    return ((x+n)>>1)%n
  return (x>>1)%n
#-------------------------------------------------------------------------------

def U_V_subscript(k, n, P, D):
  U=1
  V=P
  digits = bin(k)[2:]

  for digit in digits[1:]:
    U, V = (U*V) % n,  div2mod(V*V + D*U*U, n)

    if digit == '1':
      U,V = div2mod(P*U + V, n), div2mod(D*U + P*V, n)
  return U, V
#-------------------------------------------------------------------------------

def lucas_pp(n, D, P, Q):                                                                                                                                                                                                                         
  assert n & 1
  U, V = U_V_subscript(n+1, n, P, D)
  return U==0
#-------------------------------------------------------------------------------

def lucas_spp(n, D, P, Q):
  # Lucas strong probable prime test
  # https://arxiv.org/pdf/2006.14425v1.pdf
  # This is a bit slower than lucas_pp, so is not used

  assert n & 1
  
  d = n+1
  s = 0
  while (d & 1) == 0 :
    s+=1
    d >>= 1

  
  U, V = U_V_subscript(d, n, P, D)
  if U==0:
      return True

  Q = pow(Q,d,n)

  for r in range(s):
    if V==0:
        return True  
    V = ( V*V - 2*Q)%n
    Q=pow(Q,2,n)
  
  return False

#-------------------------------------------------------------------------------
def baillie_psw(n):

  if n <= 1: return False
  
  for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
            53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]:
    if n % p == 0:
      return n==p

  if not miller_rabin_base_2(n):
     return False

  
  D,j = D_chooser(n)
  
  if j==0:
    return False # see https://www.jstor.org/stable/2006406

  # even numbers and squares have been eliminated by this point
  return lucas_pp(n, D, 1, (1-D)//4) # slightly faster than lucas_spp
#-------------------------------------------------------------------------------

