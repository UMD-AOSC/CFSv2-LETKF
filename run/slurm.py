################################################################################
## slurm.py
## Module for enabling easy launching and monitoring of SLURM HPC jobs from
## within Python.
##
## Travis Sluka
## University of Maryland, 2015
## tsluka@umd.edu
################################################################################

## Setup logging for this module
import logging
log = logging.getLogger(__name__)


## load other modules
import os
import sys
import subprocess as sp
import time


## make sure this module is not being run from the command line
if __name__ == '__main__':
    log.critical("This file should not be run from the command line")
    sys.exit(1)

    
## module wide configurables
maxJobRetries = 5     # number of times to retry running a SLURM job after failure
maxSlurmRetries = 100 # number of times to retry SLURM commands ('sacct', 'squeue', etc.)
sleepDuration = 5     # seconds to wait after a failure
account = None        # account to submit the SLURM job as   


## global parameters
_bannedNodes = []     # list of nodes determined to be bad,
                      # subsequent job launches will not use these nodes

################################################################################

def getJobs(username=os.getenv('USER')):
    """
    Calls 'squeue' to parse a list of SLURM jobs that are running, queued,
    completed, etc. Tries mutliple times in case there is an error with the 
    'squeue' command.
    
    Arguments:
    username -- the name of the user, or the default obtained from the USER
                environment variable

    Returns: a list of jobs for that user
    """
   
    retries = 0
    cmd = 'squeue -u '+username

    ## wrap in a loop and keep trying if it fails
    while retries <= maxSlurmRetries:
        retries += 1
        try:
            ## run the 'squeue' command in the shell and get its output
            proc = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            out = out.strip().split('\n')

            ## make sure we didn't have an error
            if proc.returncode != 0  or ('squeue: error:' in err):
                raise Exception

            ## parse that output
            jobs = []
            for line in out[1:]:
                job = {}
                words = line.split()
                job['id']        = words[0]
                job['partition'] = words[1]
                job['name']      = words[2]
                job['user']      = words[3]                
                job['status']    = words[4]
                job['time']      = words[5]
                job['nodes']     = words[6]
                job['nodelist']  = words[7]
                jobs.append(job)
            return jobs
            
        except:
            ## an error occured, we'll try again after sleeping
            log.warn("Error with running 'slurm.getJobs', retrying")
            time.sleep(retries * sleepDuration)
            continue

    log.critical("Completely unable to run 'slurm.getJobs', giving up")
    raise Exception



################################################################################

def waitForJobs(jobs):
    '''
    Waits for all the jobs with ids given to finish before exiting this 
    function.

    Arguments:
    jobs -- a list of integer ids
    '''
    ## make sure the list of jobs is list of strings
    jobs = [str(j) for j in jobs]

    ## sit in a loop and wait until all jobs given are done
    wait = True
    while wait:
        jobsRunning = getJobs()
        wait = False
        ## for each of the jobs running for the user, make sure
        ## none of those are in the 'jobs' list given to the function.
        ## if they are, sleep and keep waiting.
        for job in jobsRunning:
            if job['id'] in jobs:
                wait = True
                break

        ## if pending jobs were found, sleep for a bit before retrying
        if wait:
            time.sleep(sleepDuration)



################################################################################

def getJobInfo(jobId):
    ## generate the command to run in the shell, depending on the parameters
    ## we want to get values for
    keys = ['jobid','jobname','partition','account','exitcode','ncpus','nodelist','state']
    cmd = 'sacct -P -X -j ' + str(jobId) + ' -o ' + ','.join(keys)

    retries = 0
    ## wrap in a loop and keep trying if it fails
    while retries <= maxSlurmRetries:
        retries += 1
        try:
            ## run the 'squeue' command in the shell and get its output
            proc = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            out = out.strip().split('\n')

            ## make sure we didn't have an error
            if proc.returncode != 0 or len(out) < 2 or ('squeue: error:' in err):
                raise Exception

            ## parse that output
            if len(out) < 2:
                return None
            vals = out[1].split('|')
            ret = {}
            for i in range(len(keys)):
                ret[keys[i]] = vals[i]
            return ret
                
        except:
            ## an error occured, we'll try again after sleeping
            log.warn("Error with running 'slurm.getJobInfo', retrying")
            time.sleep(retries * sleepDuration)
            continue

    log.critical("Completely unable to run 'slurm.getJobInfo', giving up")
    raise Exception


################################################################################

def monitor(jobs):
    '''
    Launches and then monitors the given SLURM jobs. 
    If, while waiting, we detect that a job did not finish, try to restart it
    '''
    watching = jobs     ## list of jobs to wait for, all of them at first
    nextWatching = []   ## the list of jobs to watch after current loop has finished

    ## submit the jobs to the SLURM system
    for j in jobs:
        j.submit()
        j.retries = 0

    ## wait for all the jobs to finish, if a failure is detected try to rerun the job
    cycle = True
    while cycle:
        cycle = False
        jobIdsRunning = [j['id'] for j in getJobs()]
        nextWatching = []
        for j in watching:
            if j.id in jobIdsRunning or (getJobInfo(j.id)['state'] in (
                    'RUNNING', 'PENDING','CONFIGURING','COMPLETING')):
                ## job is still running
                nextWatching.append(j) 
            else:
                info = getJobInfo(j.id)
                ## job finished, make sure it finished correctly
                fail = False
                
                ### check the SLURM job status
                if info['state'] != 'COMPLETED':
                    fail = True
                    log.error('job {0}, {1} completed with state {2}'.format(
                        j.name,j.id,info['state']))

                ## run the user provided check function
                if not fail and j.fnCheck:
                    if not j.fnCheck(j):
                        fail = True

                ## job did not finish correctly, possibly retry
                if fail:
                    log.error(
                        "job {0}, {1} did not finish correctly, rerunning...".format(
                            j.name,j.id))

                    ## add this node to the  banned node list
                    badNodes = info['nodelist']
                    if badNodes not in _bannedNodes:
                        log.warn("banning node(s): "+str(badNodes))
                        _bannedNodes.append(badNodes)

                    ## retry if we can
                    j.retries += 1
                    if j.retries < maxJobRetries:
                        if j.fnRetry:
                            j.fnRetry(j)
                        j.submit()
                        nextWatching.append(j)
                    else:
                        log.critical('Unable to successfully run job '+j.id)
                        raise Exception

        ## if we are still waiting on jobs to finish, sleep for a while
        ## before continuing
        if (len(nextWatching) > 0):
            watching = nextWatching
            cycle = True
            time.sleep(sleepDuration)
        
    

################################################################################

class Job:
    '''Wrapper for a SLURM job submission.
    This wrapper will automatically 

    cmd -- the command the slurm job should run
    runtime -- the max amount of time to run the job e.g. "3:00" for 3 minutes

    Function Handles:
    fnCheck -- a function to check if the job was successful
    fnRetry -- a function to run if a job failed and should be cleaned
               up before retrying
    '''


    
    def __init__(self, cmd, runtime, fnCheck=None, fnRetry=None, output=None,name=None, log=None):
        self.fnCheck   = fnCheck
        self.fnRetry = fnRetry
        self.cmd     = cmd
        self.runtime = runtime
        self.output  = output
        self.id      = None
        self.retries = 0
        self.name    = name

        
    def wait(self):
        waitForJobs([self.id,])

        
    def submit(self):
        '''
        Submits the job to the SLURM queue
        '''

        if account is None:
            raise Exception("slurm.account must be set first")

        ## create the output directory if it does not already exist
        if self.output:
            dirname = os.path.dirname(self.output)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        
        ## create the shell command to run
        shellCmd = 'sbatch -A {0} -t {1}'.format(
            account, self.runtime)
        ## add on optional arguments to the shell command
        if len(_bannedNodes) > 0:
            shellCmd += ' -x ' + ','.join(_bannedNodes)
        if self.output:
            shellCmd += ' -o ' + self.output
        if self.name:
            shellCmd += ' -J ' + self.name
                

        ## add on the main command we want to run
        shellCmd += ' ' + self.cmd

        ## submit the job
        ## wrap in a loop and keep trying if it fails
        retries = 0
        while retries <= maxSlurmRetries:
            retries += 1
            try:
                proc = sp.Popen(shellCmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
                out, err = proc.communicate()
               
                ## make sure slurm didn't have an error
                if proc.returncode != 0 or ('error:' in err):
                    log.debug(str(proc.returncode))
                    log.debug(err)
                    raise Exception

                ##TODO, make sure the job didn't actually get submitted,
                ## check for a job id somehow??

                ## get the job ID
                self.id = out.strip().split()[-1]
                log.debug("SLURM job "+str(self.name)+" submitted as "+str(self.id))

            except:
                ## an error occured, we'll try again after sleeping
                log.warn("Error with submitting SLURM job, retrying")
                time.sleep(retries * sleepDuration)
                continue
            return

        log.critical("Completely unable to run slurm job, giving up")
        raise Exception('critical error')
        
