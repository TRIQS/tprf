import sys
import time
import subprocess

import numpy as np

filename = sys.argv[-1]
it = float(sys.argv[1])

bash_command = "free -m | grep Mem >> %s"%filename

subprocess.call(['bash', '-c', "free - m | grep total > %s"%filename])

while True:
    time.sleep(it)
    #print("%.2f"%t)
    subprocess.call(['bash', '-c', bash_command])
