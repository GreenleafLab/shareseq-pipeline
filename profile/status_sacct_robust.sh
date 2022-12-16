#!/usr/bin/env bash
# Adapted from: https://github.com/jdblischak/smk-simple-slurm/
# Modifications:
# - bumped from 5 to 20 attempts, and added 1 second sleep between attempts

# Check status of Slurm job. More robust because runs `sacct` multiple times if
# needed

jobid="$1"

if [[ "$jobid" == Submitted ]]
then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

function get_status(){
  sacct -j "$1" --format State --noheader | head -n 1 | awk '{print $1}'
}

for i in {1..20}
do
  output=`get_status "$jobid"`
  if [[ ! -z $output ]]
  then
    break
  fi
  sleep 1
done

if [[ -z $output ]]
then
  echo sacct failed to return the status for jobid "$jobid" >&2
  echo Maybe you need to use scontrol instead? >&2
  exit 1
fi

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]
then
  echo running
else
  echo failed
fi