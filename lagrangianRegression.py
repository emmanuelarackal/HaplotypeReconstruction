import math
import random
import numpy
import numpy as np
import os


class LagrangianRegression:
    """
    This class is used to solve the regression by using the fast projected
    gradient method for support vector machines of Bloom et al.
    """

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

    def setVariablesToInitialConfiguration(self):
        self.p = np.zeros(0)
        self.y = np.zeros(0)
        self.h = np.zeros(0)

    """
    This method computed the matrix P
    @border1 start position of haplotype/region
    @border2 end position of haplotype/region
    @haploFolder folder where the files about haplotypes needs to be stored resp are stored
    """
    def calculateP(self, border1, border2, haploFolder):
        self.p = np.zeros((len(self.nlist), len(self.dbg.haplotypes)))
        for i in range(0, len(self.nlist)):
            for j in range(0, len(self.dbg.haplotypes)):
                if self.similarity(self.nlist[i], j, border1, border2):
                    self.p[i, j] = 1
        np.save(haploFolder+'\\p.npy', self.p)
        return

    """
    This method computes whether a given read matches a given haplotype in the
    overlapping region
    @readId position of read in pairedEndReadList
    @haploId position of haplotypen in haplotypes list
    @border1 start position of haplotype/region
    @border2 end position of haplotype/region
    """
    def similarity(self, readId, haploId, border1, border2):
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

    """
    This method computes the y vector
    @haploFolder folder where the files about haplotypes needs to be stored resp are stored
    """
    def calculateY(self, haploFolder):
        self.y = np.zeros(len(self.nlist))

        # Create the partition lists
        partitionSets = []
        a = self.dbg.pairedEndReadList[self.nlist[0]].startPos
        b = self.dbg.pairedEndReadList[self.nlist[0]].endPos
        set = []
        count = 0
        # computes the individual sets and stores them in partitionSets
        for i in range(0, len(self.nlist)):
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos == a and \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos == b:
                if np.count_nonzero(self.p[i, :]) != 0:
                    count += 1
                    set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos != a or \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos != b:
                if np.count_nonzero(self.p[i, :]) != 0:
                    count += 1
                    partitionSets.append(set)
                    set = []
                    a = self.dbg.pairedEndReadList[self.nlist[i]].startPos
                    b = self.dbg.pairedEndReadList[self.nlist[i]].endPos
                    set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)
        if set:
            partitionSets.append(set)

        pos = 0
        self.y = np.zeros(count)
        # normalizes a partition set and stores the values in vector y
        for i in range(0, len(partitionSets)):
            for j in range(0, len(partitionSets[i])):
                self.y[pos] = partitionSets[i][j] / sum(partitionSets[i])
                pos += 1

        i = 0
        # if there are situations in which a read in y doesn't match to any
        # haplotype remove them. The corresponding algorithm for y is already
        # included in the calculation of y.
        while i < len(self.p):
            if np.count_nonzero(self.p[i]) == 0:
                self.p = np.delete(self.p, i, 0)
            else:
                i = i + 1
        np.save(haploFolder+'\\y.npy', self.y)
        np.save(haploFolder+'\\p.npy', self.p)

        # stores the positions of the reads in pairedEndReadList
        # which were actually used in optimization in nList.txt.
        # This is relevant for the extension.
        file = open(haploFolder+'\\nList.txt', 'w')
        nText = ""
        for i in range(0, len(self.nlist)):
            if i == len(self.nlist) - 1:
                nText += str(self.nlist[i])
            else:
                nText += str(self.nlist[i]) + " "
        file.write(nText)
        file.close()

    """
    After calculating the first fit h, this method removes
    all entries in h having the value 0. Accordingly, all
    columns in P corresponding to the are also removed.
    """
    def removeRedundancyAfterH(self):
        i = 0
        while i < len(self.h):
            if self.h[i] == 0:
                self.p = np.delete(self.p, i, 1)
                self.h = np.delete(self.h, i, 0)
            else:
                i = i + 1

    """
    Computes the gradient of the lagrangia
    """
    def calculateGradient(self):
        m = numpy.ones((len(self.h), len(self.h)))
        np.fill_diagonal(m, 0)
        ones = np.ones(len(self.h))
        gradient = -2 * (np.transpose(self.p)).dot(self.y - self.p.dot(self.h)) + 2 * self.roh * (m.dot(self.h)) - self.lmbda * ones + self.k * (sum(self.h) - 1) * ones
        return gradient

    """
    returns max myu according to the rules mentioned in paper
    """
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

    """
    apporximates largest eigenvalue of Hessian using power method
    """
    def getMaxEigenValue(self):
        m = numpy.ones((len(self.h), len(self.h)))
        np.fill_diagonal(m, 0)
        onesMatrix = np.ones((len(self.h), len(self.h)))
        hessian = 2 * (np.transpose(self.p)).dot(self.p) + 2 * self.roh * m + self.k * onesMatrix
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

    """
    computes fit vector h
    @haploFolder folder where the files about haplotypes needs to be stored resp are stored
    """
    def calculateH(self, haploFolder):
        if np.count_nonzero(self.h) == 0:
            self.h = np.zeros(len(self.dbg.haplotypes))
            # random initialization of h
            for i in range(0, len(self.h)):
                self.h[i] = random.uniform(0, 1)
            self.h = self.h / sum(self.h)

        if np.all((self.p == 0)) and np.all((self.y == 0)):
            self.p = np.load(haploFolder+'\\p.npy')
            self.y = np.load(haploFolder+'\\y.npy')

        gradient = self.calculateGradient()
        maxMyu = self.getMaxMyu(gradient, self.h)
        # calculate accur
        rec = max(maxMyu, abs(sum(self.h) - 1))

        outer = 0
        hLine = np.copy(self.h)
        # Required accuracy
        while rec > 0.0001 and outer < 500:
            # suggested approximation of L
            l = (self.k + 1) * self.getMaxEigenValue()
            inner = 0
            t = 1
            # Required accuracy
            while maxMyu > self.theta * rec and inner < 1000:
                hHat = self.h - gradient/l
                # PBox
                for i in range(0, len(hHat)):
                    if hHat[i] < 0:
                        hHat[i] = 0
                    elif hHat[i] > 1:
                        hHat[i] = 1

                tLine = 0.5 * (1 + math.sqrt(1 + 4 * t * t))
                self.h = hHat + (t-1) / tLine * (hHat - hLine)
                hLine = np.copy(hHat)
                t = tLine
                gradient = self.calculateGradient()
                #maxMyu = self.getMaxMyu(gradient, self.h)
                maxMyu = self.getMaxMyu(gradient, hHat)
                inner += 1

            # Update
            self.lmbda = self.lmbda - self.k*(sum(self.h)-1)
            gradient = self.calculateGradient()
            maxMyu = self.getMaxMyu(gradient, self.h)
            rec = min(rec, max(maxMyu, abs(sum(self.h) - 1)))
            self.k = self.k * self.delta
            outer += 1

        # normalize h if necessary
        self.h = self.h/sum(self.h)

    """
    This method computes the error/cost of the regression
    @h fit vector h
    """
    def calculateFit(self, h):
        res = self.y - self.p.dot(h)
        return res.dot(res)

    """
    The pipe is used during read graph generation to reduce local
    haplotypes
    @haploFolder folder where the files about haplotypes needs to be stored resp are stored
    """
    def pipe(self, border1, border2, haploFolder):
        self.calculateP(border1, border2, haploFolder)
        print("Calculated P")
        self.calculateY(haploFolder)
        print("Calculated y")
        self.calculateH(haploFolder)
        if border1 < 100:
            string1 = "0" + str(border1)
        else:
            string1 = str(border1)
        if border2 < 100:
            string2 = "0" + str(border2)
        else:
            string2 = str(border2)
        if os.path.exists(haploFolder+"\\LSR_" + string1 + "_" + string2 + ".txt"):
            os.remove(haploFolder+"\\LSR_"+border1+'_'+border2+".txt")
        file = open(haploFolder+'\\LSR_' + string1 + '_' + string2 + '.txt', 'w')
        count = 0
        for i in range(0, len(self.dbg.haplotypes)):
            if self.h[i] > 0.001:
                file.write(self.dbg.haplotypes[i] + "\n")
                count += 1
        file.close()
        return count
