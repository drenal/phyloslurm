#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""@package docstring

This is a Phylobayes hybrid^2 starter and governor script 
inspired by Tom Williams' script which similar purpose.

What's:
- Phylobayes: phylogenetics program with GTR and CAT support
- hybrid: it can start Phylobayes with srun (Slurm) OR with mpirun OR single threaded
- hybrid again: if the alignments were already started but stopped for some reason, 
                this script restarts them (filenames have to match to use this feature)
- starter and governor: it queries the results each hour and stops Phylobayes if 
                        the chains converged

License: MIT

Authors:
- Lenard Szantho <lenard [\\] drenal [!] eu>
- Tom Williams <tom.a.williams [\\] bristol.ac [!] uk>

Version: v0.5 (2019-05-14)

Todo:
- input error handling
- docstrings
- slurm detection: backfall to mpirun or pure pb
- set subsampling from command line argument

Changelog:
    v0.1 (2019-03-24)
        - basic functionality for the use with Slurm and pb or pb_mpi
        - handles multiple chains
    v0.2 (2019-03-28)
        - detect number of processes from SLURM environment variable
    v0.3 (2019-04-02)
        - bugfix: subtasks may not be started, wait for them and exit if they are not started
    v0.4 (2019-05-13)
        - convenience fixes: min / max sample size, random sleep range
        - bugfix: path is not kept when calling bpcomp and tracecomp
    v0.5 (2019-05-14)
        - splitting python starter script (allowing stand-alone use) and sbatch starter which starts it
        - exit pyhton if chains are stopped for whatever reason
        - class created instead of using global variables

MIT license declaration:

Copyright 2019 Lenard Szantho, Tom Williams

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import os
import re
import time
import subprocess
import shlex
import sys
import argparse
import random

class PBStarter:

    # Convergence crtierion
    # bpcomp: <0.1 good run, <0.3 acceptable run
    BPCOMP_MAXDIFF=0.3
    # tracevomp: <0.1 and >300 good run, <0.3 and >50 acceptable run
    TRCOMP_EFFSIZE=50.0
    TRCOMP_RELDIFF=0.3

    # Stop pb chain if *.trace reaches this size
    MAX_SAMPLE_SIZE=100000
    # Run bpcomp and tracecomp if *.trace reached this size
    MIN_SAMPLE_SIZE=6667

    def detect_slurm(self):
        return True

    def start_pb_jobs(self, infile, pb_args, chains, processes):
        """Starts or restarts the specified number of Phylobayes chains
        
        Arguments:
        infile (string) : alignment file
        pb_args (string) : arguments with which the phylobayes should be started
        chains (int) : number of chains to start (and to check)
        processes (int) : number of processes phylobayes should be started with
        """
        
        # TODO check for SLURM_* environmental variables and decide:
        # whether SLURM is installed -> stand_alone or mpirun execution

        pb_command = "pb"
        if processes > 1:
            pb_command += "_mpi"
        pb_command += " "
        
        for i in range(1, chains+1):
            
            srun_started = False

            pb_out_filename = self.generate_output_filename(infile, pb_args, processes, i)

            print("Following pb chain is going to be started: " + pb_out_filename)
            sys.stdout.flush()
            
            if os.path.exists(pb_out_filename + ".trace"):
                os.system("srun --cpu_bind=v,threads -c 1 -n " + str(processes) + " -o " + pb_out_filename + ".out " + pb_command + pb_out_filename + " &") 
            else:
                os.system("srun --cpu_bind=v,threads -c 1 -n " + str(processes) + " -o " + pb_out_filename + ".out " + pb_command + pb_args + " -d " + infile + " "  + pb_out_filename + " &")
            
            slurm_job_id = os.environ.get("SLURM_JOB_ID")
            if os.environ.get("SLURM_ARRAY_JOB_ID") is not None:
                slurm_job_id = os.environ.get("SLURM_ARRAY_JOB_ID") + "_" + os.environ.get("SLURM_ARRAY_TASK_ID")

            max_wait = 60
            waited = 0
            if not slurm_job_id is None and slurm_job_id != "_":
                while srun_started == False:
                    time.sleep(1)

                    #print(slurm_job_id + "." + str(i-1))
                    #sys.stdout.flush()


                    subtask_id = slurm_job_id + "." + str(i-1)
                    srun_proc = subprocess.Popen(shlex.split("sacct -n -j " + subtask_id ), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    srun_proc_out, srun_proc_error = srun_proc.communicate()

                    waited += 1
                    
                    first_line = ""
                    if len(srun_proc_out.splitlines()) > 0:
                        first_line = srun_proc_out.splitlines()[0]

                    if len(first_line.rstrip()) != 0:
                        srun_started = True
                        continue
                    else:
                        print("Waiting for srun to start for " + subtask_id)
                        sys.stdout.flush()

                    if waited > max_wait:
                        print("Max waiting time reached, srun did not start for " + subtask_id)
                        print("Wasn't able to start every subtask, exiting...")
                        sys.stdout.flush()
                        sys.exit(3)

        return

    def check_convergence(self, infile, pb_args, chains, processes):
        """Checks using bpcomp and tracecomp whether the two or more chains converged
        
        Arguments:
        infile (string) : alignment file
        pb_args (string) : arguments with which the phylobayes was started
        chains (int) : number of chains to start (and to check)
        processes (int) : number of processes phylobayes was started with
        """
        
        # first check whether the phylobayes chains are running (if using Slurm)
        slurm_job_id = os.environ.get("SLURM_JOB_ID")
        if os.environ.get("SLURM_ARRAY_JOB_ID") is not None:
            slurm_job_id = os.environ.get("SLURM_ARRAY_JOB_ID") + "_" + os.environ.get("SLURM_ARRAY_TASK_ID")


        max_wait = 60
        waited = 0
        if not slurm_job_id is None and slurm_job_id != "_":
            every_chain_is_running = True

            for i in range(1, chains+1):
                
                subtask_id = slurm_job_id + "." + str(i-1)
                srun_proc = subprocess.Popen(shlex.split("sacct -n -j " + subtask_id ), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                srun_proc_out, srun_proc_error = srun_proc.communicate()

                first_line = ""
                if len(srun_proc_out.splitlines()):
                    first_line = srun_proc_out.splitlines()[0]

                if len(first_line.rstrip()) != 0 and first_line.find("RUNNING") != -1:
                    #print("Phylobayes job of " + subtask_id + " is running.")
                    continue
                else:
                    print("Phylobayes job of " + subtask_id + " is NOT running.")
                    every_chain_is_running = False
                    sys.stdout.flush()

            if not every_chain_is_running:
                print("Some of the chains looks like not running. Check output logs for more information.")
                print("Exiting...")
                sys.stdout.flush()
                sys.exit(4)

        
        # second check whether enough points have been sampled
        num_samples = []
        for i in range(1, chains+1):
            num_samples.append(sum(1 for line in open(self.generate_output_filename(infile, pb_args, processes, i) + ".trace")) )
            
            # if 100 000 sample size reached, stop Phylobayes
            if num_samples[i-1] > self.MAX_SAMPLE_SIZE:
                stop = os.system("stoppb " + self.generate_output_filename(infile, pb_args, processes, i)  )
        
        # get the smallest smaple size
        min_samples = min(num_samples)
        
        if min_samples >= self.MIN_SAMPLE_SIZE:
            burnin = int(float(min_samples) / 4.0)

            #enough samples. Check bpcomp and tracecomp diagnostics
            
            # use shell=False, otherwise no output is written
            # use universal_newlines to get string output, otherwise it will be binary

            #print("Starting bpcomp, " + self.generate_output_filename(infile, pb_args, processes, 0) + ".bpcomp")

            bp_proc = subprocess.Popen(shlex.split("bpcomp -o " + self.generate_output_filename(infile, pb_args, processes, 0) + ".bpcomp -x " + str(burnin) + " " + " ".join(self.generate_output_filename(infile, pb_args, processes, i) for i in range(1,chains +1)) ) , shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            
            bp_proc_out, bp_proc_error = bp_proc.communicate()
            
            bp_res = 0
            tr_res = 0
            tr_line_number = 0
            
            lines = bp_proc_out.splitlines()
            
            print("Bpcomp:")
            for line in lines:
                print(line.rstrip())
                if line.startswith("maxdiff"):
                    fields = re.split("\s+", line.rstrip())
                    if float(fields[-1]) <= self.BPCOMP_MAXDIFF:
                        print("Converged because maxdiff is " + str(fields[-1]))
                        bp_res = 1
                    else:
                        print("Not converged yet: maxdiff is " + str(fields[-1]))
                        
            #now check tracecomp

            #print("Starting tracecomp, " + self.generate_output_filename(infile,pb_args,processes, 0) + ".tracecomp")
            
            trc_proc = subprocess.Popen(shlex.split("tracecomp -o " + self.generate_output_filename(infile, pb_args, processes, 0) + ".tracecomp -x " + str(burnin) + " " + " ".join(self.generate_output_filename(infile, pb_args, processes, i) for i in range(1,chains+1))), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            
            trc_proc_out, trc_proc_error = trc_proc.communicate()

            t_lines = trc_proc_out.splitlines()
            
            print("Tracecomp:")
            for line in t_lines:
                print(line.rstrip())
                fields = re.split("\s+", line.rstrip())
                if len(fields) > 1:
                    if fields[0] == "name":
                        continue
                    else:
                        tr_line_number += 1
                        if float(fields[1]) >= self.TRCOMP_EFFSIZE:
                            print("Converged : " + fields[0] + " " + fields[1])
                            tr_res += 1
                        else:
                            print("Not converged yet : " + fields[1])
                        if float(fields[2]) <= self.TRCOMP_RELDIFF:
                            tr_res += 1
                            print("Converged : " + fields[0] + " " + fields[2])
                        else:
                            print("Not converged yet : " + fields[2])

            # write out buffered stdout messages 
            sys.stdout.flush()
                            
            if bp_res == 1 and tr_res == tr_line_number:
                for i in range(1, chains+1):
                    stop = os.system("stoppb " + self.generate_output_filename(infile, pb_args, processes, i)  )
                    
                sys.exit(0)
                
            return

        if min_samples >= self.MAX_SAMPLE_SIZE:
            sys.exit(0)

        all_running = False
        for i in range(1, chains+1):
            pb_chain_run = open(self.generate_output_filename(infile, pb_args, processes, i) + ".run")
            if pb_chain_run.read() == "1":
                all_running = True

        if not all_running:
            sys.exit(0)

        # return if we reached this point
        return
        
    def create_tokens_from_filepath(self, filepath):
        """Split the path into parts needed by other functions
        
        Arguments:
        filepath (string): relative or absolute path to file
        
        Returns:
        (string, string, string) : path to the directory, filename with extension, filename without extension    
        """
        dirbase = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        basename = os.path.splitext(filename)[0]
        
        return dirbase, filename, basename 

    def generate_output_filename(self, filepath, pb_args, processes, chain_nr):
        """Generate the output filename passed to phylobayes based on the input filename
        
        Basic return syntax: 
        - alignment file without extension
        - phylobayes parameters escaped 
        - number of processes 
        - "chain" + number of chain
        
        If chain_nr is 0, the _chainX postfix is not included in return
        """
        dirbase, filename, basename = self.create_tokens_from_filepath(filepath)
        
        pb_postfix = re.sub(r"-", "_", re.sub(r"\s+", "", pb_args))
        
        # the first -xxx became _xxx in pb_postfix so no need for "_" there
        base = basename + pb_postfix + "_" + str(processes)

        return_path = ""
        if dirbase:
            return_path = dirbase + "/" + base
        else:
            return_path = base
        
        if chain_nr == 0:
            return return_path

        return_path = return_path + "_chain" + str(chain_nr)
        
        return return_path

def main():
    """Entry method of start_phylobayes.sbatch.py
    
    Parsing input and calling the starter and governor functions.
    """
    parser = argparse.ArgumentParser(prog="start_phylobayes.sbatch.py", 
            description="This program starts and governs two or more Phylobayes chain jobs\n"
                        "either stand-alone or with mpirun or with Slurm's srun\n\n"
                        "Developed by Lenard Szantho <lenard [\\] drenal [!] eu>.", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
            help="path to the alignment. This will be the base of the output files' name.")

    # Optional settings
    parser.add_argument("--pb-args", required=False, default="-cat -gtr", help="Phylobayes arguments to use, always put it between quotation marks! (default: \"-cat -gtr\") The value of this argument will be included in the output files' name. WARNING: Because of the shell, you have to specify this argument the following way: --pb-args=\"-cat -gtr\"")
    parser.add_argument("-c", "--chains", required=False, default=2, type=int, help="number of chains to start (default: 2)")
    parser.add_argument("-p", "--processes", required=False, default=1, type=int, help="number of processes to start (default: 1). \nIf you specify a value greater than 1, pb_mpi will be started instead of pb. If you're using SLURM you can omit this option, because the number of processes can be calculated from $SLURM_NTASKS environment variable. The process number will be included in the output files' name.")

    # Bpcomp and tracecomp settings
    parser.add_argument("-m", "--max-diff-bpcomp", required=False, default=0.3, type=float, help="Required max difference reported by bpcomp to end the Phylobayes chains")
    parser.add_argument("-e", "--effective-size-tracecomp", required=False, default=50.0, type=float, help="Required effective size reported by tracecomp to end Phylobayes chains")
    parser.add_argument("-r", "--relative-diff-tracecomp", required=False, default=0.3, type=float, help="Required relative difference reported by tracecomp to end Phylobayes chains")

    # Scriptspecific behavioural settings
    parser.add_argument("-x", "--max-sample-size", required=False, default=100000, type=int, help="Samaple size above which the Phylobayes chain should be stopped")
    parser.add_argument("-n", "--min-sample-size", required=False, default=6667, type=int, help="Sample size above which bpcomp and tracecomp should be executed periodically")

    parser.add_argument("-y", "--sleep-range-min", required=False, default=3600, type=int, help="Minimum time between the periodical checks")
    parser.add_argument("-z", "--sleep-range-max", required=False, default=3600, type=int, help="Maximum time between the periodical checks.")

    args = parser.parse_args()

    pb_starter = PBStarter()
    
    #
    # error handling TODO
    #
    
    input_filename_path = args.input
    pb_args = args.pb_args
    chains = args.chains
    processes = args.processes
    
    pb_starter.BPCOMP_MAXDIFF = args.max_diff_bpcomp
    pb_starter.TRCOMP_EFFSIZE = args.effective_size_tracecomp
    pb_starter.TRCOMP_RELDIFF = args.relative_diff_tracecomp

    pb_starter.MIN_SAMPLE_SIZE = args.min_sample_size
    pb_starter.MAX_SAMPLE_SIZE = args.max_sample_size
    SLEEP_MIN = args.sleep_range_min
    SLEEP_MAX = args.sleep_range_max

    # if Slurm is used, make sure that the number of processes * chains 
    # doesn't exceed number of tasks allocated for this sbatch job
    slurm_ntasks = os.environ.get("SLURM_NTASKS")
    if not slurm_ntasks is None:
        if processes != chains * int(slurm_ntasks):
            processes = int(int(slurm_ntasks)/chains)

    #print("Parameters: input=" + input_filename_path + " pb_args=" + str(pb_args) + " chains=" + str(chains) + " processes=" + str(processes) + " bpcomp=" + str(pb_starter.BPCOMP_MAXDIFF) + " trcomp_eff=" + str(pb_starter.TRCOMP_EFFSIZE) + " trcomp_reldiff=" + str(pb_starter.TRCOMP_RELDIFF) + " min_sample=" + str(pb_starter.MIN_SAMPLE_SIZE) + " max_sample=" + str(pb_starter.MAX_SAMPLE_SIZE) + " sleep_min=" + str(SLEEP_MIN) + " sleep_max=" + str(SLEEP_MAX) )

    sys.stdout.flush()

    # start or restart jobs
    pb_starter.start_pb_jobs(input_filename_path, pb_args, chains, processes)
    
    # check for convergence once an hour
    while True:
        time.sleep( random.randint(SLEEP_MIN, SLEEP_MAX) )
        pb_starter.check_convergence(input_filename_path, pb_args, chains, processes)


if __name__ == "__main__":
    main()
