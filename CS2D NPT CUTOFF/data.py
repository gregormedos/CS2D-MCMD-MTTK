import os


path = os.getcwd()
print(path)
temperature = [1.000]
print(temperature)
pritiski = list()
print(pritiski)
izoterme = temperature
izobare = list()

metode = ['mc', 'md']


def data_izoterme(metoda) :
    with open(path + '\\!' + metoda + '_data_izoterme', 'w') as fw :
        for temperatura in izoterme :
            for pritisk in pritiski :
                with open(path + '\\' + '{:.3f}'.format(temperatura) + '\\' + '{:.3f}'.format(pritisk) + '\\!' + metoda + '_data', 'r') as fr :
                    line = fr.readline()
                    fw.write(line)
    with open(path + '\\' + metoda + '_data_izoterme', 'w') as fw :
        for temperatura in izoterme :
            for pritisk in pritiski :
                Temp = '{:16.7E}'.format(temperatura)
                Tlak = '{:16.7E}'.format(pritisk)
                fw.write(Temp)
                fw.write(Tlak)
                with open(path + '\\' + '{:.3f}'.format(temperatura) + '\\' + '{:.3f}'.format(pritisk) + '\\' + metoda + '_data', 'r') as fr :
                    for line in fr :
                        line = line.strip('\n')
                        fw.write(line)
                fw.write('\n')


for metoda in metode :
    data_izoterme(metoda)
    #data_izobare(metoda)


def grafi_izoterme(metoda) :
    with open(path + '\\!' + metoda + '_data_izoterme', 'r') as fr :
        for j in range(len(izoterme)) :
            with open(path + '\\!' + metoda + '_data_izoterma' + str(j+1), 'w') as fw :
                for _ in range(len(pritiski)) :
                    line = fr.readline()
                    fw.write(line)
    with open(path + '\\' + metoda + '_data_izoterme', 'r') as fr :
        for j in range(len(izoterme)) :
            with open(path + '\\' + metoda + '_data_izoterma' + str(j+1), 'w') as fw :
                for _ in range(len(pritiski)) :
                    line = fr.readline()
                    fw.write(line)


for metoda in metode :
    grafi_izoterme(metoda)
    #grafi_izobare(metoda)