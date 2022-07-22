import os


path = os.getcwd()
print(path)
temperature = [0.120]
print(temperature)
gostote = [round(0.100*(i+1), 3) for i in range(8)]
print(gostote)
izoterme = temperature
izohore = list()

metode = ['mc', 'md']


def data_izoterme(metoda) :
    with open(path + '/!' + metoda + '_data_izoterme', 'w') as fw :
        for temperatura in izoterme :
            for gostota in gostote :
                with open(path + '/' + '{:.3f}'.format(temperatura) + '/' + '{:.3f}'.format(gostota) + '/!' + metoda + '_data', 'r') as fr :
                    line = fr.readline()
                    fw.write(line)
    with open(path + '/' + metoda + '_data_izoterme', 'w') as fw :
        for temperatura in izoterme :
            for gostota in gostote :
                Temp = '{:16.7E}'.format(temperatura)
                Gost = '{:16.7E}'.format(gostota)
                fw.write(Temp)
                fw.write(Gost)
                with open(path + '/' + '{:.3f}'.format(temperatura) + '/' + '{:.3f}'.format(gostota) + '/' + metoda + '_data', 'r') as fr :
                    for line in fr :
                        line = line.strip('\n')
                        fw.write(line)
                fw.write('\n')


for metoda in metode :
    data_izoterme(metoda)
    #data_izohore(metoda)


def grafi_izoterme(metoda) :
    with open(path + '/!' + metoda + '_data_izoterme', 'r') as fr :
        for j in range(len(izoterme)) :
            with open(path + '/!' + metoda + '_data_izoterma' + str(j+1), 'w') as fw :
                for _ in range(len(gostote)) :
                    line = fr.readline()
                    fw.write(line)
    with open(path + '/' + metoda + '_data_izoterme', 'r') as fr :
        for j in range(len(izoterme)) :
            with open(path + '/' + metoda + '_data_izoterma' + str(j+1), 'w') as fw :
                for _ in range(len(gostote)) :
                    line = fr.readline()
                    fw.write(line)


for metoda in metode :
    grafi_izoterme(metoda)
    #grafi_izohore(metoda)