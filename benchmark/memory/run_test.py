import subprocess
import itertools

for norb, nk, nw in itertools.product([1, 2, 4, 6, 10], [128, 256], [5000, 10000]):

    print(norb, nk, nw)

    filename = './data/norb_%s_nk_%s_nw_%s.txt'%(norb, nk, nw)

    observe_memory = subprocess.Popen(['python', './observe_memory.py', '0.01', filename])
    run_gf = subprocess.Popen(['python', './calculate_gf.py', str(norb), str(nk), str(nw)])

    run_gf.wait()
    observe_memory.kill()
