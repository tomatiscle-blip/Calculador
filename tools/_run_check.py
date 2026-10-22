import subprocess, sys
p = subprocess.Popen([sys.executable, '-m', 'app'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
try:
    out, err = p.communicate(timeout=12)
    print('exit', p.returncode)
    print('STDOUT:', out[:1000])
    print('STDERR:', err[:2000])
except subprocess.TimeoutExpired:
    p.kill()
    out, err = p.communicate()
    print('timeout killed; exit', p.returncode)
    print('STDERR:', err[:2000])
