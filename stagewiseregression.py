import numpy
import numpy as np
import os


class StagewiseRegression:

    def __init__(self, dbg, nlist, roh):
        self.p = np.zeros(0)
        self.h = np.zeros(0)
        self.y = np.zeros(0)
        self.dbg = dbg
        self.nlist = nlist
        self.roh = roh

    def setVariablesToInitialConfiguration(self):
        self.p = np.zeros(0)
        self.y = np.zeros(0)

    def calculateP(self, border1, border2):
        print(len(self.dbg.haplotypes), len(self.nlist))
        self.p = np.zeros((len(self.nlist), len(self.dbg.haplotypes)))
        for i in range(0, len(self.nlist)):
            for j in range(0, len(self.dbg.haplotypes)):
                if self.similarity(self.nlist[i], j, border1, border2):
                    self.p[i, j] = 1
        np.save('Haplotypes\\p.npy', self.p)
        return

    def similarity(self, readId, haploId, border1, border2):
        start = 0
        if self.dbg.pairedEndReadList[readId].startPos < border1:
            start = self.dbg.pairedEndReadList[readId].startPos
        else:
            start = border1
        for i in range(start, len(self.dbg.variantPositions)):
            if self.dbg.pairedEndReadList[readId].endPos < i or border2 - 1 < i:
                break
            if self.dbg.pairedEndReadList[readId].startPos > i or border1 > i:
                continue
            posOfRead = i - self.dbg.pairedEndReadList[readId].startPos
            posOfHaplo = i - border1
            if self.dbg.pairedEndReadList[readId].read[posOfRead] != self.dbg.haplotypes[haploId][posOfHaplo]:
                if self.dbg.pairedEndReadList[readId].read[posOfRead] != ".":
                    return False
        return True

    def calculateY(self):
        self.y = np.zeros(len(self.nlist))

        # Create the partition lists
        partitionSets = []
        a = self.dbg.pairedEndReadList[self.nlist[0]].startPos
        b = self.dbg.pairedEndReadList[self.nlist[0]].endPos
        set = []
        for i in range(0, len(self.nlist)):
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos == a and \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos == b:
                set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos != a or \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos != b:
                partitionSets.append(set)
                set = []
                a = self.dbg.pairedEndReadList[self.nlist[i]].startPos
                b = self.dbg.pairedEndReadList[self.nlist[i]].endPos
                set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)

        pos = 0
        for i in range(0, len(partitionSets)):
            for j in range(0, len(partitionSets[i])):
                self.y[pos] = partitionSets[i][j] / sum(partitionSets[i])
                pos += 1
        np.save('Haplotypes\\y.npy', self.y)

    def removeRedundancyAfterP(self):
        i = 0
        while i < len(self.y):
            if np.count_nonzero(self.p[i]) == 0:
                self.p = np.delete(self.p, i, 0)
                self.y = np.delete(self.y, i, 0)
            else:
                i = i + 1

        print("Sum2 ", sum(self.y))
        print("ATT ", len(self.y))
        self.y = self.y / sum(self.y)

        np.save('Haplotypes\\p.npy', self.p)
        np.save('Haplotypes\\y.npy', self.y)

    def removeRedundancyAfterH(self):
        for i in range(0, len(self.h)):
            if self.h[i] == 0:
                self.p = np.delete(self.p, i, 1)
                self.h = np.delete(self.h, i, 1)

    def calculateH(self):
        #self.removeRedundancyAfterP()
        eps = 0.00001
        self.h = np.zeros(len(self.dbg.haplotypes))
        constraint = 0
        print(self.y)
        print(self.p)
        if np.all((self.p == 0)) and np.all((self.y == 0)):
            self.p = np.load('Haplotypes\\p.npy')
            self.y = np.load('Haplotypes\\y.npy')

        m = np.ones((len(self.dbg.haplotypes), len(self.dbg.haplotypes)))
        np.fill_diagonal(m, 0)
        count = 0
        while constraint < 1:
            firstDerivative = 2 * (np.transpose(self.p)).dot(self.p.dot(self.h) - self.y)
            secondDerivative = 2 * self.roh * (m.dot(self.h))
            maxValue = -1000000
            index = -1
            for i in range(0, len(self.h)):
                value = -(firstDerivative[i] + secondDerivative[i])
                if value > maxValue:
                    maxValue = value
                    index = i
            self.h[index] += eps
            constraint = sum(self.h)
            count += 1
        #for i in range(0, len(self.h)):
        #    print(i, round(self.h[i] / sum(self.h), 3))
        #print(sum(abs(self.y - self.p.dot(self.h))))
        return 0

    def pipe(self, border1, border2):
        self.calculateY()
        print("Calculated y")
        self.calculateP(border1, border2)
        print("Calculated P")
        self.calculateH()
        if border1 < 100:
            string1 = "0" + str(border1)
        else:
            string1 = str(border1)
        if border2 < 100:
            string2 = "0" + str(border2)
        else:
            string2 = str(border2)
        if os.path.exists("Haplotypes\\LSR_"+string1+"_"+string2+".txt"):
            os.remove("Haplotypes\\LSR_'+border1+'_'+border2+'.txt")
        file = open('Haplotypes\\LSR_' + string1 + '_' + string2 + '.txt', 'w')
        count = 0
        for i in range(0, len(self.dbg.haplotypes)):
            if self.h[i] >= 0.001:
                file.write(self.dbg.haplotypes[i] + "\n")
                count += 1
        file.close()
        return count
