class Read:

    def __init__(self, id, read, startPos, count=1, filteredRead=None):
        self.id = id
        self.read = read
        self.startPos = startPos
        self.count = count
        self.filteredRead = filteredRead


class FilteredRead:

    def __init__(self, id, read, startPos, endPos, count=1):
        self.id = id
        self.read = read
        self.startPos = startPos
        self.endPos = endPos
        self.count = count


class Window:

    def __init__(self, id, startPos, endPos, startFilteredList, endFilteredList, list=None, classList=None,
                 haploCount=0):
        if list is None:
            list = []
        if classList is None:
            classList = []
        self.id = id
        self.startPos = startPos
        self.endPos = endPos
        self.startFilteredList = startFilteredList
        self.endFilteredList = endFilteredList
        # structure of list = [read inside window, current assigned class, history of assigned classes with protocol]
        self.list = list
        self.classList = classList
        self.haploCount = haploCount


class HaplotypeClass:

    # structure of haplotypeClass = [id of haplotype, nucleotide sequence, list containing the reads assigned to class]
    def __init__(self, id, haplotype=None, list=None, count=0):
        if list is None:
            list = []
        if haplotype is None:
            haplotype = ""
        self.id = id
        self.haplotype = haplotype
        self.list = list
        self.count = count


class Node_OLC:

    def __init__(self, id, readPointer, part, filteredStartPos, filteredEndPos, predecessor=None, successor=None):
        if predecessor is None:
            predecessor = []
        if successor is None:
            successor = []
        self.id = id
        self.readPointer = readPointer
        # if splitted read, 0 means first, 1 means second part
        self.part = part
        self.filteredStartPos = filteredStartPos
        self.filteredEndPos = filteredEndPos
        self.predecessor = predecessor
        self.successor = successor
