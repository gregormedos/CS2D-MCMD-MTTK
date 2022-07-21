from random import randint
import os
import shutil
import subprocess


path = os.getcwd()
print(path)

print("###################################")
print("#    MC AND MD SIMULATION         #")
print("#    NPT ENSEMBLE                 #")
print("#    2D CORE-SOFTENED DISCS       #")
print("###################################")

Ran = int()
Pres = float(input('pressure: '))
TempPoints = int(input('number of temperature points: '))
TempStart = float(input('starting temperature: '))
TempStep = float(input('size of temperature step: '))

print("-----------------------------------")
print("--- PARAMETERS INITIALIZED --------")
print("-----------------------------------")

os.chdir(path + '\\' + '{:.3f}'.format(Pres))
for i in range(TempPoints) :
    Ran = randint(1, 32000)
    Temp = TempStart + i * TempStep
    os.mkdir(path + '\\' + '{:.3f}'.format(Pres) + '\\' + '{:.3f}'.format(Temp))
    os.chdir(path + '\\' + '{:.3f}'.format(Pres) + '\\' + '{:.3f}'.format(Temp))
    with open(file='md_input', mode='w') as f :
        f.write("'ensemble parameters'\n")
        f.write("'number of particles (N)'                    100\n")
        f.write("'pressure'                                   " + '{:.3f}'.format(Pres) + 'd0' + '\n')
        f.write("'temperature'                                " + '{:.3f}'.format(Temp) + 'd0' + '\n')
        f.write("'simulation parameters'\n")
        f.write("'starting density'                           0.5d0\n")
        f.write("'LJ potential parameters'                    1.0d0           1.0d0\n")
        f.write("'CS potenital parameters'                    5.0d0           1.0d0           0.7d0\n")
        f.write("'random seed'                                " + str(Ran) + '\n')
        f.write("'number of series of sampling (MAX 20)'      20\n")
        f.write("'equilibration length'                       100000\n")
        f.write("'MC parameters'\n")
        f.write("'number of cycles (N trial moves)'           100000\n")
        f.write("'trial moves per sample (same as N)'         100\n")
        f.write("'trial moves per volume change'              20\n")
        f.write("'MAX random displacement by sigma'           0.1d0\n")
        f.write("'MAX random scaling of volume'               0.2d0\n")
        f.write("'MD parameters'\n")
        f.write("'number of time steps'                       100000\n")
        f.write("'time step'                                  0.001d0\n")
        f.write("'temperature time constant'                  0.1d0\n")
        f.write("'pressure time constant'                     0.5d0\n")
        f.write("'time steps per frame'                       100\n")
        f.write("'correlation function parameters'\n")
        f.write("'g(r) radial interval'                       0.02d0\n")
        f.write("'cycles/steps per g(r) sampling'             10\n")
    shutil.copy(path + '\\' + 'cs2d_npt_cutoff.exe', path + '\\' + '{:.3f}'.format(Pres) + '\\' + '{:.3f}'.format(Temp))
    print("-----------------------------------")
    print("--- CS2D_NPT_CUTOFF INITIALIZED ---")
    print("-----------------------------------")
    Process = subprocess.Popen("cs2d_npt_cutoff.exe", shell=False)
    Process.wait()
    print("JOB FINISHED")

print("###################################")
print("#    ENDPROGRAM                   #")
print("###################################")
