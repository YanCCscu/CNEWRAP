#!/usr/bin/env python3
import sys
import subprocess
import time
from multiprocessing.pool import Pool
def runcmd(command, runshell=True, runstdout = None, runstderr = None):
	try:
		current_time = time.strftime("%Y-%m-%d %H:%M:%S",
		time.localtime(time.time()))
		print(current_time, "\n", command, "\n", sep="")
		subprocess.check_call(command, shell=runshell, stdout = runstdout, stderr = runstderr,executable="/bin/bash")
	except:
		sys.exit("Error occured when running command:\n%s" % command)

def runcmd_pipe(command, runshell=True, runstdout = subprocess.PIPE, runstderr = None):
	try:
		proc=subprocess.Popen(command, shell=runshell, stdout = runstdout, text=True, stderr = runstderr,executable="/bin/bash")
		return(proc.stdout)
	except:
		sys.exit("Error occured when running command:\n%s" % command)


def runpool(fun,args,thread_num):
	p = Pool(thread_num)
	res=p.starmap(fun,args)
	p.close()
	p.join()
	return(res)

def subrunpool(fun,args,thread_num):
	subp = Pool(thread_num)
	res=subp.starmap(fun,args)
	subp.close()
	subp.join()
	return(res)
def get_system_id():
	info = {}
	with open("/etc/os-release", "r") as f:
		for line in f:
			if "=" in line:
				key, value = line.strip().split("=", 1)
				info[key] = value.strip('"')
	return info.get("ID", "").lower()  # 'centos', 'ubuntu'

def print_nested_list(mylist, indent=0):
	space = "    " * indent  
	if isinstance(mylist, list):
		print("#----------level %s----------"%indent)
		for item in mylist:
			print_nested_list(item, indent + 1)
	else:
		print(space + str(mylist))
