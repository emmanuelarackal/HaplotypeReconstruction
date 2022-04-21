import math
import random
import numpy
import numpy as np
import os


class LagrangianRegression:

    def __init__(self, dbg, nlist, roh, lmbda, k, theta, delta, accuracy):
        self.p = np.zeros(0)
        self.h = np.zeros(0)
        self.y = np.zeros(0)
        self.dbg = dbg
        self.nlist = nlist
        self.roh = roh
        self.lmbda = lmbda
        self.k = k
        self.theta = theta
        self.delta = delta
        self.accuracy = accuracy

    def calculateP(self, border1, border2):
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
        for i in range(0, len(self.y)):
            if np.count_nonzero(self.p[i]) == 0:
                numpy.delete(self.p, i, 0)
                numpy.delete(self.y, i, 0)

    def removeRedundancyAfterH(self):
        for i in range(0, len(self.h)):
            if self.h[i] == 0:
                self.p = np.delete(self.p, i, 1)
                self.h = np.delete(self.h, i, 1)

    def calculateGradient(self):
        m = numpy.ones((len(self.h), len(self.h)))
        np.fill_diagonal(m, 0)
        ones = np.ones(len(self.h))
        firstDerivative = 2 * (np.transpose(self.p)).dot(self.p.dot(self.h) - self.y)
        secondDerivative = 2 * self.roh * (m.dot(self.h))
        thirdDerivative = - self.lmbda * ones
        forthDerivative = self.k * (sum(self.h) - 1) * ones
        derivative = firstDerivative + secondDerivative + thirdDerivative + forthDerivative
        return derivative

    def getMaxMyu(self, gradient, h):
        maxValue = -10000
        for i in range(0, len(h)):
            if h[i] == 0:
                localMax = max(0, -gradient[i])
            elif h[i] == 1:
                localMax = max(0, gradient[i])
            else:
                localMax = abs(gradient[i])
            if localMax > maxValue:
                maxValue = localMax
        return maxValue

    def getMaxEigenValue(self):
        # Power Method from https://www.sciencedirect.com/topics/mathematics/power-method
        m = numpy.ones((len(self.h), len(self.h)))
        np.fill_diagonal(m, 0)
        ones = np.ones(len(self.h))
        hessian = 2 * (np.transpose(self.p)).dot(self.p) + 2 * self.roh * m + self.k * ones.dot(ones)
        x = np.ones(len(self.h))
        x = x / numpy.sqrt(numpy.sum(x ** 2))
        prevE = 0
        eigenValue = 10
        while abs(eigenValue - prevE) / eigenValue > 0.001:
            x = hessian.dot(x)
            x = x / numpy.sqrt(numpy.sum(x**2))
            prevE = eigenValue
            eigenValue = (hessian.dot(x)).dot(x) / (x.dot(x))
        return eigenValue


    def calculateH(self):
        self.h = np.zeros(len(self.dbg.haplotypes))
        for i in range(0, len(self.h)):
            self.h[i] = random.uniform(0, 1)
        self.h = self.h/sum(self.h)

        if np.all((self.p == 0)) and np.all((self.y == 0)):
            self.p = np.load('Haplotypes\\p.npy')
            self.y = np.load('Haplotypes\\y.npy')

        gradient = self.calculateGradient()
        maxMyu = self.getMaxMyu(gradient, self.h)
        # calculate accur
        rec = max(maxMyu, abs(sum(self.h) - 1))
        outer = 0
        # Required accuracy
        while rec > self.accuracy and outer < 1000:
            # suggested approximation of L
            l = (self.k + 1) * self.getMaxEigenValue()
            inner = 0
            t = 1
            hLine = np.copy(self.h)
            # Required accuracy
            while maxMyu > self.theta * rec and inner < 1000000:
                hHat = self.h - gradient/l
                # PBox
                for i in range(0, len(hHat)):
                    if hHat[i] < 0:
                        hHat[i] = 0
                    elif hHat[i] > 1:
                        hHat[i] = 1

                tLine = 0.5 * (1 + math.sqrt(1 + 4 * math.pow(t, 2)))
                self.h = hHat + (t-1) / tLine * (hHat - hLine)
                hLine = np.copy(hHat)
                t = tLine
                gradient = self.calculateGradient()
                maxMyu = self.getMaxMyu(gradient, hHat)
                inner += 1
            self.lmbda = self.lmbda - self.k*(sum(self.h)-1)
            rec = min(rec, max(maxMyu, abs(sum(self.h) - 1)))
            self.k = self.k * self.delta
            outer += 1
        self.h = self.h/sum(self.h)
        for i in range(0, len(self.h)):
            print(i, round(self.h[i] / sum(self.h), 3))
        print(sum(abs(self.y - self.p.dot(self.h))))
        print("END")

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
        if os.path.exists("Haplotypes\\LSR_" + string1 + "_" + string2 + ".txt"):
            os.remove("Haplotypes\\LSR_'+border1+'_'+border2+'.txt")
        file = open('Haplotypes\\LSR_' + string1 + '_' + string2 + '.txt', 'w')
        count = 0
        for i in range(0, len(self.dbg.haplotypes)):
            if self.h[i] > 0.001:
                file.write(self.dbg.haplotypes[i] + "\n")
                count += 1
        file.close()
        return count
