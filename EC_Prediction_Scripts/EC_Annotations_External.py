#!/usr/bin/env python

import os
import os.path
import sys
import time
import subprocess
import pexpect.pxssh

Input_Dir = "/scratch/j/jparkins/mobolaji/NOD_Mouse/Single_End/Mapped_Reads_ECs_Sam/EC_Mapping/Input"
Output_Dir = "/scratch/j/jparkins/mobolaji/NOD_Mouse/Single_End/Mapped_Reads_ECs_Sam/EC_Mapping/Output"

Split_Dir = os.path.join(Input_Dir, "Split")
Python = "/home/j/jparkins/mobolaji/python"
Script_Dir = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts"

Done = False

print "Commands:"
print " ".join([Python, os.path.join(Script_Dir, "0_Preprocess_Input.py"), Input_Dir, Split_Dir])
print " ".join([Python, os.path.join(Script_Dir, "1-1_Detect_Submission.py"), Split_Dir, Output_Dir])
print " ".join([Python, os.path.join(Script_Dir, "2-1_PRIAM_Submission.py"), Split_Dir, Output_Dir])
print " ".join([Python, os.path.join(Script_Dir, "3-1_BLAST-Diamond_Submission.py"), Split_Dir, Output_Dir]) + "\n"


try:
    ssh = pexpect.pxssh.pxssh()
    hostname = "login.scinet.utoronto.ca"
    username = sys.argv[1]
    password = sys.argv[2]
    ssh.login(hostname, username, password)
    ssh.sendline('gpc')   # run a command
    ssh.prompt()             # match the prompt
    ssh.sendline(" ".join([Python, os.path.join(Script_Dir, "0_Preprocess_Input.py"), Input_Dir, Split_Dir]))
    ssh.prompt()
    print  ssh.before
    print "0_Preprocess_Input complete"
    ssh.sendline(" ".join([Python, os.path.join(Script_Dir, "1-1_Detect_Submission.py"), Split_Dir, Output_Dir]))
    ssh.prompt()
    print  ssh.before
    print "1-1_Detect_Submission complete"
    ssh.sendline(" ".join([Python, os.path.join(Script_Dir, "2-1_PRIAM_Submission.py"), Split_Dir, Output_Dir]))
    ssh.prompt()
    print  ssh.before
    print "2-1_PRIAM_Submission complete"
    ssh.sendline(" ".join([Python, os.path.join(Script_Dir, "3-1_BLAST-Diamond_Submission.py"), Split_Dir, Output_Dir]))
    ssh.prompt()
    print  ssh.before
    print "3-1_BLAST-Diamond_Submission complete"
except pexpect.pxssh.ExceptionPxssh as e:
    print "pxssh failed on login. " + e

while not Done:
    try:
        ssh = pexpect.pxssh.pxssh()
        hostname = "login.scinet.utoronto.ca"
        username = sys.argv[1]
        password = sys.argv[2]
        ssh.login(hostname, username, password)
        ssh.sendline('gpc')   # run a command
        ssh.prompt()             # match the prompt
        ssh.sendline("ls")
        ssh.prompt()
        data = ssh.before
        if "Detect_Watcher.o" in data and "PRIAM_Watcher.o" in data and "BLAST_Watcher.o" in data:
            ssh.sendline(" ".join([Python, os.path.join(Script_Dir, "4_EC_Consolidation.py"), Input_Dir, Split_Dir, Output_Dir]))
            ssh.prompt()
            Done = True
        else:
            ssh.sendline("qstat")
            ssh.prompt()
            print ssh.before
    except pexpect.pxssh.ExceptionPxssh as e:
        print "pxssh failed on login. " + e
    time.sleep(600)