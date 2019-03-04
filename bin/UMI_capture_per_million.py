#!/ur/bin/python2.7

import random,time

random.seed(time.time())

with open('/dev/stdin','r') as f:
    UMI=f.read().split("\n")

random.shuffle(UMI)

d={}

for x in UMI[0:1000000]:
    d[x]=""

print(len(d))
