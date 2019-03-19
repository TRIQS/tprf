import sys
import time
import subprocess

import numpy as np

filename = sys.argv.pop()

duration, it = [float(ele) for ele in sys.argv[1:]]

bash_command = "free -m | grep Mem >> %s"%filename

subprocess.call(['bash', '-c', "free - m | grep total > %s"%filename])

for t in np.arange(0, duration, it):
    time.sleep(it)
    print("%.2f"%t)
    subprocess.call(['bash', '-c', bash_command])
