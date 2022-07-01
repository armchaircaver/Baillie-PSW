from baillie_psw import is_square, D_chooser, lucas_pp, lucas_spp
from gmpy2 import is_prime
from time import perf_counter


def readfile(filename):
  f = open(filename, 'r')
  numbers = []
  for line in f.readlines():
      fields = line.split(',')
      numbers.append(int(fields[0]))
  return numbers

# pseudoprimes from http://ntheory.org/pseudoprimes.html
spsps = readfile('spsps.txt')
print( 'spsps.txt', len(spsps) )
      
slpsps_baillie = readfile('slpsps-baillie.txt')
print( 'slpsps-baillie.txt', len(slpsps_baillie) )

lpsps_baillie = readfile('lpsps-baillie.txt')
print( 'lpsps-baillie.txt', len(lpsps_baillie) )

#print( set(slpsps_baillie).intersection(set(lpsps_baillie)) )

set_slpsps_baillie = set(slpsps_baillie)
set_lpsps_baillie = set(lpsps_baillie)

for n in spsps:
  if n in set_slpsps_baillie:
    print(n)

for n in spsps:
  if n in set_lpsps_baillie:
    print(n)


def match_lpsps_baillie(MAX):
  begintime= perf_counter()
  for n in range(5,MAX,2):
    D,j = D_chooser(n)
    if j==-1:
      Lpp = lucas_pp(n, D, 1, (1-D)//4)
      prime = is_prime(n)
      expected_result = (Lpp==prime) ^ (n in set_lpsps_baillie)
      if not expected_result :
        print(n, f"prime={prime} Lpp={Lpp}, in lpsps = {n in set_lpsps_baillie}")
        return
    if n%200000==1:
      print("lpsps tested up to ", n, ", ", perf_counter()-begintime, "sec")
  print("lpsps tested up to ", n, ", ", perf_counter()-begintime, "sec")



def match_slpsps_baillie(MAX):
  begintime= perf_counter()
  for n in range(5,MAX,2):
    D,j = D_chooser(n)
    if j==-1:
      Lpp = lucas_spp(n, D, 1, (1-D)//4)
      prime = is_prime(n)
      expected_result = (Lpp==prime) ^ (n in set_slpsps_baillie)
      if not expected_result :
        print(n, f"prime={prime} Lpp={Lpp}, in slpsps = {n in set_slpsps_baillie}")
        return
    if n%200000==1:
      print("slpsps tested up to ", n, ", ", perf_counter()-begintime, "sec")
  print("slpsps tested up to ", n, ", ", perf_counter()-begintime, "sec")


match_slpsps_baillie(10**6)
match_lpsps_baillie(10**6)
    
