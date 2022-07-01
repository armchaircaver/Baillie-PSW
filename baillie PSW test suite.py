
from array import array
from baillie_psw import baillie_psw

def primegen(MAX):
  plist = array('b', [0,0,1]) + array('b',[1,0])*(MAX//2+1)
  for n in range(3,int(MAX**0.5)+1,  2):
    if plist[n]: plist[n*n::2*n]= array('b',[0])*len(plist[n*n::2*n])
  yield 2
  for n in range(3,MAX,2):
    if plist[n]:
      yield n

# from https://tucs.fi/publications/attachment.php?fname=G46.pdf, p 128
for n in (79397009999, 63278892599, 2013745337604001, 894221105778001,
          582920080863121, 443372888629441, 28295303263921,443372888629441 ):
  if baillie_psw(n):
    print(f"FAIL for large carmichael number {n}" )

print(f"large carmichael number tests completed" )
  
#largest known twin prime - a bit too ambitious
"""
p = 2996863034895*(2**1290000)-1
n = p*(p+2)
if baillie_psw(n):
  print(f"FAIL for largest twin prime number" )
else:
  print(f"pass for largest twin prime number" )
"""
          
A217255 =[5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439, 100127, 113573, 115639, 130139, 155819, 158399, 161027, 162133, 176399, 176471, 189419, 192509, 197801, 224369, 230691, 231703, 243629, 253259, 268349, 288919, 313499, 324899]
for n in A217255 :
  if baillie_psw(n):
    print("FAIL", n, baillie_psw(n) )
print("A217255 tests completed")

A217719 = [989, 3239, 5777, 10877, 27971, 29681, 30739, 31631, 39059, 72389, 73919, 75077, 100127, 113573, 125249, 137549, 137801, 153931, 155819, 161027, 162133, 189419, 218321, 231703, 249331, 370229, 429479, 430127, 459191, 473891, 480689, 600059, 621781, 632249, 635627]    
for n in A217719 :
  if baillie_psw(n):
    print("FAIL", n, baillie_psw(n) )
print("A217719 tests completed")

# Strong pseudoprimes to bases 2 and 3.
A072276 = [1373653, 1530787, 1987021, 2284453, 3116107, 5173601, 6787327, 11541307, 13694761, 15978007, 16070429, 16879501, 25326001, 27509653, 27664033, 28527049, 54029741, 61832377, 66096253, 74927161, 80375707, 101649241]
for n in A072276 :
  if baillie_psw(n):
    print("FAIL", n, baillie_psw(n) )
print("A072276 tests completed")


for m in range(2,10000):
  n = m*m
  if baillie_psw(n):
    print("FAIL", n, baillie_psw(n) )
print("small square tests completed")
    
for m in range(10**60, 10**60+2000):
  n = m*m
  if baillie_psw(n):
    print("FAIL", n, baillie_psw(n) )

print("large square tests completed")

for n in primegen(10**6):
  if not baillie_psw(n):
    print("FAIL for prime", n, baillie_psw(n) )
print("prime tests completed")

MAX=10**6
plist = array('b', [0,0,1]) + array('b',[1,0])*(MAX//2+1)
for n in range(3,int(MAX**0.5)+1,  2):
  if plist[n]: plist[n*n::2*n]= array('b',[0])*len(plist[n*n::2*n])

for n in range(MAX):
    if baillie_psw(n) != plist[n] :
      print("FAIL for number", n, baillie_psw(n) )
    if n%1000000==0 and n > 0:
        print(f"test of all numbers up to {n} completed")
print(f"test of all numbers up to {MAX} completed")

delta = 1000000
for exp in range(7,16):
  MIN = 10**exp
  deltalist = array('b',[1])*delta

  for p in primegen(int((MIN+delta)**0.5) + 2):
    start = (MIN//p)*p
    if start < MIN : start += p
    deltalist[start-MIN::p] = array('b',[0])*len(deltalist[start-MIN::p])

  for n in range(MIN, MIN+delta):
    if baillie_psw(n) != deltalist[n-MIN] :
      print("FAIL for number", n, baillie_psw(n) )
  print(f"segmented sieve test of {delta} numbers from 10^{exp} completed")
    
  # search for twin primes
  twincount=0
  for n in range(MIN, MIN+delta-2):
    if deltalist[n-MIN] and deltalist[n-MIN+2] :
      pq = n*(n+2)
      twincount+=1
      if baillie_psw(pq):
        print("FAIL for twin prime product", pq )
        
  print(f"twin prime test of {twincount} numbers from 10^{exp} completed")
    
# large carmichael number
# https://en.wikipedia.org/wiki/Carmichael_number
p=29674495668685510550154174642905332730771991799853043350995075531276838753171770199594238596428121188033664754218345562493168782883
n= p*(313*(p-1) +1)*(353*(p-1)+1)
if baillie_psw(n):
  print("FAIL for large carmichael number" )
else:
  print("pass for large carmichael number" )
  
