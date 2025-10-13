#!/usr/bin/env python3
import os,sys
import time
import tempfile
cmdir=os.path.dirname(os.path.abspath(sys.argv[0]))
if os.path.exists(os.path.join(cmdir,"lib/libdrmaa.so.1.0")):
	os.environ['DRMAA_LIBRARY_PATH']=os.path.join(cmdir,"lib/libdrmaa.so.1.0")
else:
	os.environ['DRMAA_LIBRARY_PATH']=os.path.join(cmdir,"../lib/libdrmaa.so.1.0")
#os.system('export DRMAA_LIBRARY_PATH')
import drmaa

def process_bar(num, total,commandline, jobid):
	rate = float(num)/total
	ratenum = int(100*rate)
	message = f'submit job "{commandline}" in jobid: {jobid}'
	r = '[{}{}]{}%'.format('*'*ratenum,' '*(100-ratenum), ratenum)
	sys.stdout.write(f"\r{message}\033[B")
	sys.stdout.flush()
	
	sys.stdout.write(f"\r{r}\033[F")
	sys.stdout.flush()

def process_bar2(num, total,commandline, jobid):
	rate = float(num)/total
	ratenum = int(100*rate)
	message = f'submit job "{commandline}" in jobid: {jobid}'
	r = '[{}{}]{}%'.format('*'*ratenum,' '*(100-ratenum), ratenum)
	sys.stdout.write(f"\r{message}\n")
	sys.stdout.flush()
	
	sys.stdout.write(f"\r{r}")
	sys.stdout.flush()


def submit_job(s, commandline, jobname, sge_para,tmpjobdir):
	"""submit a single job  and return jobid"""
	if isinstance(commandline,list):
		if len(commandline) < 1:
			return None
		icommand="/bin/bash"
		with tempfile.NamedTemporaryFile(mode='w', suffix='.sh',dir=tmpjobdir, delete=False) as temp_script:
			temp_script.write("\n".join(commandline))
			temp_script_path = temp_script.name 
		iargs=[temp_script_path]
	else:
		command_parts = commandline.strip().split()
		if len(command_parts) == 0:
			return None  
		icommand = command_parts[0]
		if len(command_parts) > 1:
			iargs=command_parts[1:]
	
	jt = s.createJobTemplate()
	jt.nativeSpecification = sge_para
	jt.joinFiles = True
	jt.remoteCommand = icommand
	if iargs:
		jt.args = iargs
	jt.jobName = jobname
	jt.outputPath=f':{tmpjobdir}'
	jobid = s.runJob(jt)
	s.deleteJobTemplate(jt)
	return jobid

def subjobs(scripts, maxjobs=100,check_interval=10,jobname="CNEali", sge_para="-l vf=1G,p=1"):
	"""Submit tasks to the SGE cluster and control the number of concurrent tasks """
	totaljobs=len(scripts)
	submitedjobs=0
	tmpjobdir="tmp_%s"%jobname
	os.makedirs(tmpjobdir,exist_ok=True)
	with drmaa.Session() as s:
		running_jobs = set()

		def check_and_submit_jobs():
			nonlocal running_jobs
			nonlocal submitedjobs
			while scripts and len(running_jobs) < int(maxjobs):
				commandline = scripts.pop(0)
				jobid = submit_job(s, commandline, jobname, sge_para,tmpjobdir)
				submitedjobs+=1
				if jobid:
					running_jobs.add(jobid)
					process_bar(submitedjobs,totaljobs,commandline, jobid)
		# submit for the first time, not more than MAX_RUNNING_JOBS
		check_and_submit_jobs()

		while running_jobs:
			# check status every CHECK_INTERVAL
			time.sleep(int(check_interval))
			completed_jobs = []
			
			# check completed jobs
			for jobid in running_jobs:
				try:
					job_info = s.jobStatus(jobid)
					if job_info in [drmaa.JobState.DONE, drmaa.JobState.FAILED]:
						completed_jobs.append(jobid)
				except drmaa.errors.InvalidJobException:
					completed_jobs.append(jobid)  # include those canceled jobs
			
			
			# remove completed jobs
			for jobid in completed_jobs:
				running_jobs.remove(jobid)
			
			# check and submit new jobs
			check_and_submit_jobs()

		print('\n\nAll jobs have completed.')


if __name__=='__main__':
	scripts=[['/data/nfs2/lrh/inputs/sleeper.py 1 afadsd', '/data/nfs2/lrh/inputs/sleeper.py 2 bfadsd'],
		'/data/nfs2/lrh/inputs/sleeper.py 3 cfadsd',
		'/data/nfs2/lrh/inputs/sleeper.py 2 dfadsd',
		'/data/nfs2/lrh/inputs/sleeper.py 1 efadsd']
	subjobs(scripts,maxjobs=2,check_interval=2,jobname="TestPySub", sge_para="-l vf=1G,p=1")
