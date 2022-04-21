import os
from helperClasses import Read
from scipy.spatial import distance

class Statistics:

    def __init__(self):
        self.readTable = []
        self.entropyTable = []
        self.clusteringTable = []
        self.correctTable = []
        self.filteredPosition = []

    def createReadTable(self):
        file = open('Documents\\5V_allPairs.txt', 'r')
        lines = file.readlines()
        file.close()

        rawTable = []
        count = 0
        print("Reading all reads from csv file")
        for line in lines:
            parts = line.split(' ')
            read = Read(count, parts[7], int(parts[3]))
            rawTable.append(read)
            count = count + 1
        rawTable.sort(key=lambda x: x.startPos, reverse=False)

        # Algorithm to merge duplicate reads
        print("Merging all duplicates")
        currentPos = -1
        sameStartPos = 0
        readTablePos = 0
        self.readTable = []
        for read in rawTable:
            if read.startPos == currentPos:
                duplicate = False
                for i in range(readTablePos - sameStartPos, readTablePos):
                    if self.readTable[i].read == read.read:
                        self.readTable[i].count += 1
                        duplicate = True
                        break
                if not duplicate:
                    read.id = readTablePos
                    self.readTable.append(read)
                    sameStartPos = sameStartPos + 1
                    readTablePos = readTablePos + 1
            else:
                currentPos = read.startPos
                sameStartPos = 1
                read.id = readTablePos
                self.readTable.append(read)
                readTablePos = readTablePos + 1


    def createEntropyTable(self):
        print("Creating Entropy Table")
        file = open('Documents/EntropyFile.txt', 'r')
        lines = file.readlines()
        file.close()

        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 2:
                parts[1] = parts[1].replace(" ", "")
                parts[1] = parts[1].replace("\n", "")
                #read = Read(count, parts[2], int(parts[0]), int(parts[1]), parts[3])
                read = Read(count, parts[2], int(parts[0]), int(parts[1]))
                self.entropyTable.append(read)
                count = count + 1

        for read in self.entropyTable:
            self.correctTable.append([read, [], ""])

    def createFilteredPositionList(self):
        file = open('Documents\\filteredPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.filteredPosition.append(int(numbers[i]))

    def shortestLengtOfFilteredRead(self):
        length = 100
        for i in range(0, len(self.entropyTable)):
            if len(self.entropyTable[i].filteredRead) < length:
                length = len(self.entropyTable[i].filteredRead)
        print(length)

    def createClusteringTable(self, creationMode):

        if creationMode == 0:
            print("Creating Clustering Table")
            for file in os.listdir("Data"):
                if file.startswith("ClusteringFile"):
                    f = open("Data\\" + file, 'r')
                    lines = f.readlines()
                    f.close()
                    count = 0
                    for line in lines:
                        parts = line.split(' ')
                        if len(parts[2]) > 1:
                            self.correctTable[count][1].append(parts[2])
                        count += 1

            print("Collected all files")
            pos = 0
            for elem in self.correctTable:
                for i in range(0, len(elem[0].read)):
                    if elem[0].read[i] != ".":
                        nucleotideCount = [0, 0, 0, 0, 0]
                        nucleotides = ["A", "C", "T", "G", "-"]
                        for h in range(0, len(elem[1])):
                            if len(elem[1][h]) - 1 > i:
                                if elem[1][h][i] == "A":
                                    nucleotideCount[0] += 1
                                elif elem[1][h][i] == "C":
                                    nucleotideCount[1] += 1
                                elif elem[1][h][i] == "T":
                                    nucleotideCount[2] += 1
                                elif elem[1][h][i] == "G":
                                    nucleotideCount[3] += 1
                                elif elem[1][h][i] == "-":
                                    nucleotideCount[4] += 1
                                else:
                                    continue
                        maxValue = -1
                        index = -1
                        for t in range(0, len(nucleotideCount)):
                            if nucleotideCount[t] >= maxValue:
                                maxValue = nucleotideCount[t]
                                index = t
                        if index != -1 and maxValue > 0:
                            elem[2] += nucleotides[index]
                        else:
                            elem[2] += elem[0].read[i]
                    else:
                        elem[2] += "."
                pos += 1
            print("Corrected reads")


            file = open('Documents\\ClusteringFile.txt', 'w')

            for i in range(0, len(self.correctTable)):
                file.write(
                    str(self.correctTable[i][0].startPos) + " " + str(self.correctTable[i][0].count) + " " + self.correctTable[i][
                        2] + "\n")
            file.close()

        file = open('Documents\\ClusteringFile.txt', 'r')
        lines = file.readlines()
        file.close()

        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) != 3:
                continue
            parts[2] = parts[2].replace("\n", "")
            parts[2] = parts[2].replace(" ", "")
            read = Read(count, parts[2], int(parts[0]), int(parts[1]))
            self.clusteringTable.append(read)
            count = count + 1


    def testEnd(self, haplotypes):

        print("Applying similarity test")
        count = 0
        entropyDistribution = [0, 0, 0, 0, 0]
        clusteringDistribution = [0, 0, 0, 0, 0]
        list = []
        arrow = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
        for i in range(0, 10000):
            list.append([i, 0])

        for i in range(0, len(self.clusteringTable)):
            val1 = [0, 0, 0, 0, 0]
            val1[0] = self.similarity(self.entropyTable[i], haplotypes[0], self.entropyTable[i].read)
            val1[1] = self.similarity(self.entropyTable[i], haplotypes[1], self.entropyTable[i].read)
            val1[2] = self.similarity(self.entropyTable[i], haplotypes[2], self.entropyTable[i].read)
            val1[3] = self.similarity(self.entropyTable[i], haplotypes[3], self.entropyTable[i].read)
            val1[4] = self.similarity(self.entropyTable[i], haplotypes[4], self.entropyTable[i].read)
            min1 = 1000
            index1 = -1
            for k in range(0, len(val1)):
                if val1[k] < min1:
                    min1 = val1[k]
                    index1 = k
            entropyDistribution[index1] += self.entropyTable[i].count
            val2 = [0, 0, 0, 0, 0]
            val2[0] = self.similarity(self.entropyTable[i], haplotypes[0], self.clusteringTable[i].read)
            val2[1] = self.similarity(self.entropyTable[i], haplotypes[1], self.clusteringTable[i].read)
            val2[2] = self.similarity(self.entropyTable[i], haplotypes[2], self.clusteringTable[i].read)
            val2[3] = self.similarity(self.entropyTable[i], haplotypes[3], self.clusteringTable[i].read)
            val2[4] = self.similarity(self.entropyTable[i], haplotypes[4], self.clusteringTable[i].read)
            min2 = 1000
            index2 = -1
            if min1 < min2:
                list[self.clusteringTable[i].startPos][1] += 1
            for k in range(0, len(val2)):
                if val2[k] < min2:
                    min2 = val2[k]
                    index2 = k
            if index1 != index2:
                count += 1
                #print(index1, " *", index2, " #", min1, " _", min2, " !", len(self.correctTable[i][1]), " ", self.clusteringTable[i].startPos, " ", i, " ")
            arrow[index1][index2] += 1
            clusteringDistribution[index2] += self.entropyTable[i].count

        #print(list)
        print(count)
        #print(entropyDistribution)
        #print(clusteringDistribution)
        print(arrow)

    def getStatistics(self, haplotypes, mode):
        array = []
        for i in range(0, 100):
            array.append([0, 0, 0, 0, 0])
        if mode == 1:
            for i in range(0, len(self.readTable)):
                if 790 <= self.readTable[i].startPos and self.readTable[i].startPos + len(
                        self.readTable[i].read) <= 5096:
                    distance = [0, 0, 0, 0, 0]
                    for t in range(0, len(haplotypes)):
                        distance[t] = self.similarity(self.readTable[i], haplotypes[t], self.readTable[i].read)
                    minValue, index = self.leastValue(distance)
                    if minValue < 100:
                        array[minValue][index] += 1
        elif mode == 2:
            for i in range(0, len(self.entropyTable)):
                if 790 <= self.entropyTable[i].startPos and self.entropyTable[i].startPos + len(
                        self.entropyTable[i].read) <= 5096:
                    distance = [0, 0, 0, 0, 0]
                    for t in range(0, len(haplotypes)):

                        distance[t] = self.similarity(self.entropyTable[i], haplotypes[t], self.entropyTable[i].read)
                    minValue, index = self.leastValue(distance)
                    if minValue < 100:
                        array[minValue][index] += 1
        else:
            count = 0
            for i in range(0, len(self.clusteringTable)):
                if 790 <= self.clusteringTable[i].startPos and self.clusteringTable[i].startPos + len(
                        self.clusteringTable[i].read) <= 5096:
                    count += self.clusteringTable[i].count
                    distance = [0, 0, 0, 0, 0]
                    for t in range(0, len(haplotypes)):
                        distance[t] = self.similarity(self.clusteringTable[i], haplotypes[t], self.clusteringTable[i].read)
                    minValue, index = self.leastValue(distance)
                    if minValue < 100:
                        array[minValue][index] += self.clusteringTable[i].count
            print(count)
        print(array)

    def getMismatchDistributionOfClusteringTable(self, haplotypes, matchMode, printMode):
        misMatches = []
        matches = []
        starting = []
        ending = []
        for i in range(0, 10000):
            misMatches.append([0, 0, 0, 0, 0])
            matches.append([0, 0, 0, 0, 0])
            starting.append([0, 0, 0, 0, 0])
            ending.append([0, 0, 0, 0, 0])

        for i in range(0, len(self.clusteringTable)):
            distances = []

            distances.append(self.similarity(self.clusteringTable[i], haplotypes[0], self.clusteringTable[i].read))
            distances.append(self.similarity(self.clusteringTable[i], haplotypes[1], self.clusteringTable[i].read))
            distances.append(self.similarity(self.clusteringTable[i], haplotypes[2], self.clusteringTable[i].read))
            distances.append(self.similarity(self.clusteringTable[i], haplotypes[3], self.clusteringTable[i].read))
            distances.append(self.similarity(self.clusteringTable[i], haplotypes[4], self.clusteringTable[i].read))
            minValue = 10000
            index = -1
            for k in range(0, 5):
                if distances[k] < minValue:
                    minValue = distances[k]
                    index = k
            if self.clusteringTable[i].startPos > len(haplotypes[index]):
                return None
            check = False
            for k in range(0, len(self.clusteringTable[i].read)):
                if self.clusteringTable[i].startPos + k < len(haplotypes[index]):
                    if self.clusteringTable[i].read[k] == haplotypes[index][self.clusteringTable[i].startPos + k - 1]:
                        if matchMode == 0:
                            matches[self.clusteringTable[i].startPos + k - 1][index] += 1
                        elif matchMode == 1:
                            if minValue == 0:
                                matches[self.clusteringTable[i].startPos + k - 1][index] += 1
                    if self.clusteringTable[i].read[k] != haplotypes[index][self.clusteringTable[i].startPos + k - 1]:
                        if self.clusteringTable[i].read[k] != ".":
                            misMatches[self.clusteringTable[i].startPos + k - 1][index] += 1
                else:
                    break


        for i in range (790, 5090):
            if misMatches[i] != [0, 0, 0, 0, 0]:
                if printMode == 0:
                    print(i + 1, " ", misMatches[i], " ", matches[i])
                elif printMode == 1:
                    if misMatches[i][0] > matches[i][0] or misMatches[i][1] > matches[i][1] \
                        or misMatches[i][2] > matches[i][2] or misMatches[i][3] > matches[i][3] \
                        or misMatches[i][4] > matches[i][4] or misMatches[i][0] > 200 or \
                            misMatches[i][1] > 200 or misMatches[i][2] > 200 or misMatches[i][3] > 200 or \
                            misMatches[i][4] > 200:
                        print(i + 1, " ", misMatches[i], " ", matches[i])


    def getMismatchDistributionOfReadTable(self, haplotypes):
        misMatches = []
        matches = []

        for i in range(0, 10000):
            misMatches.append([0, 0, 0, 0, 0])
            matches.append([0, 0, 0, 0, 0])

        for i in range(0, len(self.readTable)):
            distances = []

            distances.append(self.similarity(self.readTable[i], haplotypes[0], self.readTable[i].read))
            distances.append(self.similarity(self.readTable[i], haplotypes[1], self.readTable[i].read))
            distances.append(self.similarity(self.readTable[i], haplotypes[2], self.readTable[i].read))
            distances.append(self.similarity(self.readTable[i], haplotypes[3], self.readTable[i].read))
            distances.append(self.similarity(self.readTable[i], haplotypes[4], self.readTable[i].read))
            minValue = 10000
            index = -1
            for k in range(0, 5):
                if distances[k] < minValue:
                    minValue = distances[k]
                    index = k
            if self.readTable[i].startPos > len(haplotypes[index]):
                return None
            check = False
            for k in range(0, len(self.readTable[i].read)):
                if self.readTable[i].startPos + k < len(haplotypes[index]):
                    if self.readTable[i].read[k] == haplotypes[index][self.readTable[i].startPos + k - 1]:
                        matches[self.readTable[i].startPos + k - 1][index] += 1
                    if self.readTable[i].read[k] != haplotypes[index][self.readTable[i].startPos + k - 1]:
                        if self.readTable[i].read[k] != ".":
                            misMatches[self.readTable[i].startPos + k - 1][index] += 1
                else:
                    break
        for i in range(0, len(misMatches)):
            if misMatches[i] != [0, 0, 0, 0, 0]:
                print(i + 1, " ", misMatches[i], " ", matches[i])

    def similarity(self, read, haplotype, corRead):
        if read.startPos > len(haplotype):
            return None
        distance = 0
        for i in range(0, len(read.read)):
            if read.startPos + i < len(haplotype):
                if corRead[i] != haplotype[read.startPos + i - 1]:
                    if corRead[i] != ".":
                        distance += 1
            else:
                break
        return distance

    def leastValue(self, list):
        minValue = list[0]
        index = 0
        for i in range(1, len(list)):
            if list[i] is None:
                continue
            if list[i] < minValue:
                minValue = list[i]
                index = i
        return minValue, index

    def diversityLevel(self, list):
        diversity = 0
        count = 0
        for i in range(0, len(list)):
            for j in range(0, len(list)):
                if i != j:
                    for k in range(0, len(list[i])):
                        if list[i][k] != list[j][k]:
                            diversity += 1
                    count += 1
        print(count)
        print(len(list[0]))
        print("Average Diversity: ", (diversity / count))
        diversity = diversity / count
        print("Diversity level: ", (diversity / len(list[0])))

    def printHaplotypes(self, haplotypes, mode):
        hap1F = ""
        hap2F = ""
        hap3F = ""
        hap4F = ""
        hap5F = ""

        if mode == 0:
            file = open('Documents\\filteredPositions.txt', 'r')
            line = file.readlines()
            file.close()
            file.close()
            numbers = line[0].split(" ")
            filteredPosition = []
            for i in range(0, len(numbers)):
                filteredPosition.append(int(numbers[i]))

            for i in range(3000, 3400):
                if i in filteredPosition:
                    hap1F += haplotypes[0][i]
                    hap2F += haplotypes[1][i]
                    hap3F += haplotypes[2][i]
                    hap4F += haplotypes[3][i]
                    hap5F += haplotypes[4][i]

        else:
            file = open('Documents\\variantPositions.txt', 'r')
            line = file.readlines()
            file.close()
            file.close()
            numbers = line[0].split(" ")
            variantPosition = []
            for i in range(0, len(numbers)):
                variantPosition.append(int(numbers[i]))

            for i in range(790, 5096):
                if i in variantPosition:
                    hap1F += haplotypes[0][i]
                    hap2F += haplotypes[1][i]
                    hap3F += haplotypes[2][i]
                    hap4F += haplotypes[3][i]
                    hap5F += haplotypes[4][i]

        print(hap1F)
        print(hap2F)
        print(hap3F)
        print(hap4F)
        print(hap5F)

def getHaplotypes():
    file = open('Documents\\5VirusMixReference.fasta', 'r')
    lines = file.readlines()
    hap1 = ""
    for i in range(1, 163):
        hap1 += lines[i]
    hap1 = hap1.replace("\n", "")

    hap2 = ""
    for i in range(164, 326):
        hap2 += lines[i]
    hap2 = hap2.replace("\n", "")

    hap3 = ""
    for i in range(327, 489):
        hap3 += lines[i]
    hap3 = hap3.replace("\n", "")

    hap4 = ""
    for i in range(490, 652):
        hap4 += lines[i]
    hap4 = hap4.replace("\n", "")

    hap5 = ""
    for i in range(653, 815):
        hap5 += lines[i]
    hap5 = hap5.replace("\n", "")

    haplotypes = []
    haplotypes.append(hap1)
    haplotypes.append(hap2)
    haplotypes.append(hap3)
    haplotypes.append(hap4)
    haplotypes.append(hap5)

    return haplotypes

if __name__ == '__main__':
    s = Statistics()
    list = getHaplotypes()
    #s.diversityLevel(list)
    s.printHaplotypes(list, 1)
    #s.createEntropyTable()
    #s.createReadTable()
    #s.createClusteringTable(0)
    #s.createFilteredPositionList()
    #s.shortestLengtOfFilteredRead()
    #s.testEnd(list)
    #s.getMismatchDistributionOfReadTable(list)
    #s.getMismatchDistributionOfClusteringTable(list, 0, 1)
    #s.getMismatchDistribution2(list)
    #s.getStatistics(list, 3)
    #s.getMismatchDistribution(list)
