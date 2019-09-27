#!/usr/bin/python2.7

import socket,os,sys,multiprocessing,time
from sys import argv
import psutil
process = psutil.Process(os.getpid())

def myprint(*args):
	import sys
	value=args[0]
	for x in args[1:]:
		value+=" "+str(x)
	x=sys.stderr.write(value+"\n")
	

with open("Concatamer_summary.tab","w") as f:
	f.write("")

def totalsize(o):
    size = 0
    for x in o:
        size += sys.getsizeof(x)
    size += sys.getsizeof(o)
    return size


def get_contamination(data):
	concat='TAGTCGAACTGAAGGTCTCCAGC'
	rev='AGATCGGAAGAGCGGTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
	cc=[]
	rr=[]
	r=0
	n=0
	both=0
	for x in range(0,len(concat)-7):
		cc.append(concat[x:x+8])

	for x in range(0,len(rev)-8):
		rr.append(rev[x:x+8])

	reverse={}
	concatamer={}
	RandC={}
	for x in data:
		b=0
		for s in cc:
			if s in x:
				n+=1
				b+=1
				if s+" "+x in concatamer:
					concatamer[s+" "+x]+=1
				else: 
					concatamer[s+" "+x]=1
				break
		for s in rr:
			if s in x:
				r+=1
				b+=1
				if s+" "+x in reverse:
					reverse[s+" "+x]+=1
				else: 
					reverse[s+" "+x]=1
				break
		if b == 2: 
			both+=1;
			if s+" "+x in RandC: 
				RandC[s+" "+x]+=1
			else: 
				RandC[s+" "+x]=1
	with open("Contamination_sequences.txt","a") as f:
		for x in reverse:
			f.write("reverse\t"+x+"\t"+str(reverse[x])+"\n")
		for x in concatamer:
			f.write("concatamer\t"+x+"\t"+str(concatamer[x])+"\n")
		f.write("\nBoth:\n")
		for x in RandC:
			f.write("Both\t"+x+"\t"+str(RandC[x])+"\n")
	return [r,n,both]


def getlines4(f):
    for i, line in enumerate(f):
        if i % 4 == 1:
            yield line.rstrip()

def get_con(in_file,running,memory_usage):
	running.value=running.value+1
	if (in_file.split(".")[-1]=="gz"):
		import gzip
		myprint("Unzipping & opening",in_file)
		with gzip.open(in_file,'rb') as f:
			data=list(getlines4(f))
	else:
		myprint("Opening",in_file)
		with open(in_file,'r') as f:
			data=list(getlines4(f))

	footprint=totalsize(data)
	memory_usage.value+=footprint

	myprint("parsing seq data:", in_file)
	l=[]
	import math
	t=time.time()
	if len(data)/10==0:
		myvar=1
	else:
		myvar=len(data)/10
	for b in range(0,len(data),int(myvar)):
		l.append(data[b:b+int(myvar)])
	length=len(data)
	del data
	p=multiprocessing.Pool(10)
	counts=p.map(get_contamination,l)
	r=sum([x[0] for x in counts])
	n=sum([x[1] for x in counts])
	both=sum([x[2] for x in counts])
	#both=len([x for x in counts if x == [1,1]])
	n=float(n)/length*100
	r=float(r)/length*100
	both=float(both)/length*100
	with open('Concatamer_summary.tab','a') as f:
		f.write("%s\t%1.2f\t%2.2f\t%2.2f\n"%(in_file.split(".")[0],n,r,both))
	p.terminate()
	memory_usage.value-=footprint
	running.value=running.value-1

running=multiprocessing.Value('i',0)
memory_usage=multiprocessing.Value('i',0)
processes={}
cores=int(argv[1])
n=0
for item in argv:
	if item==argv[0] or item==argv[1]:continue
	myprint("Active Processes:")
	for x in list(processes.keys()):
		if processes[x].is_alive():
			myprint("\t"+processes[x].name)
		else:
			processes[x].terminate()
			del(processes[x])
	while (memory_usage.value+1000000000 > psutil.virtual_memory().available) or (len(processes) >= cores):
		time.sleep(5)
	p=multiprocessing.Process(target=get_con,args=(item,running,memory_usage),name=item)
	processes[n]=p
	n+=1
	p.start()
	time.sleep(1)
	if (running.value > 0):###This should be if, not while!
		myprint("Available Memory:",psutil.virtual_memory().available/1000000000.0)
		myprint("Running Process:")
		for x in list(processes.keys()):
			if processes[x].is_alive():
				myprint("\t"+processes[x].name)
			else:
				processes[x].terminate()
				del(processes[x])
		time.sleep(10)

myprint("Contamination sequences printed")
exit(0)



