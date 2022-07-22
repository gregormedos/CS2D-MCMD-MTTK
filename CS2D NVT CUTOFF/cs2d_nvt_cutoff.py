from random import randint
import os
import shutil
import subprocess


path = os.getcwd()
print(path)

print("###################################")
print("#    MC AND MD SIMULATION         #")
print("#    NVT ENSEMBLE                 #")
print("#    2D CORE-SOFTENED DISCS       #")
print("###################################")

Ran = int()
Temp = float(input('temperature: '))
GostPoints = int(input('number of density points: '))
GostStart = float(input('starting density: '))
GostStep = float(input('size of density step: '))

print("-----------------------------------")
print("--- PARAMETERS INITIALIZED --------")
print("-----------------------------------")

os.chdir(path + '/' + '{:.3f}'.format(Temp))
for i in range(GostPoints) :
    Ran = randint(1, 32000)
    Gost = GostStart + i * GostStep
    os.mkdir(path + '/' + '{:.3f}'.format(Temp) + '/' + '{:.3f}'.format(Gost))
    os.chdir(path + '/' + '{:.3f}'.format(Temp) + '/' + '{:.3f}'.format(Gost))
    with open(file='md_input', mode='w') as f :
        f.write("'ensemble parameters'\n")
        f.write("'number of particles (N)'                   100\n")
        f.write("'density'                                   " + '{:.3f}'.format(Gost) + 'd0' + '\n')
        f.write("'temperature'                               " + '{:.3f}'.format(Temp) + 'd0' + '\n')
        f.write("'simulation parameters'\n")
        f.write("'LJ potential parameters'                   1.0d0           1.0d0\n")
        f.write("'CS potenital parameters'                   5.0d0           1.0d0           0.7d0\n")
        f.write("'random seed'                               " + str(Ran) + '\n')
        f.write("'number of series of sampling (MAX 20)'     20\n")
        f.write("'equilibration length'                      100000\n")
        f.write("'MC parameters'\n")
        f.write("'number of cycles (N trial moves)'          100000\n")
        f.write("'trial moves per sample (same as N)'        100\n")
        f.write("'MAX random displacement by sigma'          0.1d0\n")
        f.write("'MD parameters'\n")
        f.write("'number of time steps'                      100000\n")
        f.write("'time step'                                 0.001d0\n")
        f.write("'temperature time constant'                 0.1d0\n")
        f.write("'time steps per frame'                      100\n")
        f.write("'correlation function parameters'\n")
        f.write("'g(r) radial interval'                      0.02d0\n")
        f.write("'cycles/steps per g(r) sampling'            10\n")
        f.write("'cycles/steps per msd(n)/msd(t) sampling'   10\n")
    shutil.copy(path + '/' + 'cs2d_nvt_cutoff.exe', path + '/' + '{:.3f}'.format(Temp) + '/' + '{:.3f}'.format(Gost))
    print("-----------------------------------")
    print("--- CS2D_NVT_CUTOFF INITIALIZED ---")
    print("-----------------------------------")
    Process = subprocess.Popen("./cs2d_nvt_cutoff.exe", shell=False)
    Process.wait()
    print("JOB FINISHED")

print("###################################")
print("#    ENDPROGRAM                   #")
print("###################################")
