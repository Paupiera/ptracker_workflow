Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 2
Job stats:
job           count
----------  -------
test_conda        1
total             1

Select jobs to execute...

[Fri May 31 12:25:41 2024]
rule test_conda:
    output: test.delme
    log: /maps/projects/rasmussen/scratch/ptracker/ptracker/snakemake_output/imout, /maps/projects/rasmussen/scratch/ptracker/ptracker/snakemake_output/imout
    jobid: 0
    reason: Missing output files: test.delme
    resources: mem_mb=1000, mem_mib=954, disk_mb=1000, disk_mib=954, tmpdir=<TBD>, walltime=0, mem_gb=1


            echo "HELLO"
            sleep 10
            snakemake &> test.delme
            
Submitted job 0 with external jobid 'Submitted batch job 302822'.
[Fri May 31 12:26:21 2024]
Error in rule test_conda:
    jobid: 0
    output: test.delme
    log: /maps/projects/rasmussen/scratch/ptracker/ptracker/snakemake_output/imout, /maps/projects/rasmussen/scratch/ptracker/ptracker/snakemake_output/imout (check log file(s) for error details)
    shell:
        
            echo "HELLO"
            sleep 10
            snakemake &> test.delme
            
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    cluster_jobid: Submitted batch job 302822

Error executing rule test_conda on cluster (jobid: 0, external: Submitted batch job 302822, jobscript: /maps/projects/rasmussen/scratch/ptracker/ptracker/.snakemake/tmp.u5fcn5xg/snakejob.test_conda.0.sh). For error details see the cluster log and the log files of the involved rule(s).
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2024-05-31T122539.374864.snakemake.log
