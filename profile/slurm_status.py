#!/usr/bin/env python3
# Modified from: https://github.com/Snakemake-Profiles/slurm Copyright (c) 2017 Snakemake-Profiles, MIT license
import json
import os
import re
import requests
import subprocess as sp
import shlex
import sys
import time
import logging

logger = logging.getLogger(__name__)

STATUS_ATTEMPTS = 20
SIDECAR_VARS = os.environ.get("SNAKEMAKE_CLUSTER_SIDECAR_VARS", None)
DEBUG = bool(int(os.environ.get("SNAKEMAKE_SLURM_DEBUG", "0")))

if DEBUG:
    logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)

def get_status_direct(jobid):
    """Get status directly from sacct/scontrol"""
    for i in range(STATUS_ATTEMPTS):
        try:
            sacct_res = sp.check_output(shlex.split(f"sacct -P -b -j {jobid} -n"))
            res = {x.split("|")[0]: x.split("|")[1] for x in sacct_res.decode().strip().split("\n")}
            break
        except sp.CalledProcessError as e:
            logger.error("sacct process error")
            logger.error(e)
        except IndexError as e:
            logger.error(e)
            pass
        # Try getting job with scontrol instead in case sacct is misconfigured
        try:
            sctrl_res = sp.check_output(shlex.split(f"scontrol -o show job {jobid}"))
            m = re.search(r"JobState=(\w+)", sctrl_res.decode())
            res = {jobid: m.group(1)}
            break
        except sp.CalledProcessError as e:
            logger.error("scontrol process error")
            logger.error(e)
            if i >= STATUS_ATTEMPTS - 1:
                print("failed")
                exit(0)
            else:
                time.sleep(1)

    return res[jobid] or ""


jobid = sys.argv[1]

logger.debug("slurm-status.py: direct query")
status = get_status_direct(jobid)

logger.debug("job status: %s", repr(status))

if status == "BOOT_FAIL":
    print("failed")
elif status == "OUT_OF_MEMORY":
    print("failed")
elif status.startswith("CANCELLED"):
    print("failed")
elif status == "COMPLETED":
    print("success")
elif status == "DEADLINE":
    print("failed")
elif status == "FAILED":
    print("failed")
elif status == "NODE_FAIL":
    print("failed")
elif status == "PREEMPTED":
    print("failed")
elif status == "TIMEOUT":
    print("failed")
elif status == "SUSPENDED":
    print("running")
else:
    print("running")
