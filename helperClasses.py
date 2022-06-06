class Read:
    """
        Represents a read
        @id id of the read
        @read sequence of read
        @startPos start position of read
        @count number of appearance of read in dataset
        @endPos end position of read
        @filteredRead compact sequence of read
        @filteredStart start position of read in compact form
        @dotStart position where the gap starts, if no gap, then -1
        @dotEnd position where the gap stops, if no gap, then -1
    """
    def __init__(self, id, read, startPos, count=1, endPos=None, filteredRead=None, filteredStart=None, dotStart=None, dotEnd=None):
        self.id = id
        self.read = read
        self.startPos = startPos
        self.endPos = endPos
        self.count = count
        self.filteredRead = filteredRead
        self.filteredStart = filteredStart
        self.dotStart = dotStart
        self.dotEnd = dotEnd


class FilteredRead:

    """
    Represents a compact form read
    @id id of the compact read
    @read compact sequence of read
    @startPos start position of read in compact form
    @endPos end position of read in compact form
    @count number of appearance of read in dataset
    @dotStart position where the gap starts, if no gap, then -1
    @dotEnd position where the gap stops, if no gap, then -1
    """
    def __init__(self, id, read, startPos, endPos, count=1, dotStart=None, dotEnd=None):
        self.id = id
        self.read = read
        self.startPos = startPos
        self.endPos = endPos
        self.count = count
        self.dotStart = dotStart
        self.dotEnd = dotEnd


class Window:
    """
        Represents a window in the error correction, clustering approach
        @id id of the window
        @startPos start position of window on genome
        @endPos start position of window on genome
        @startFilteredList start position in compact form of window on genome
        @endFilteredList end position in compact form of window on genome
        @list list contains for all included reads a tuple with the structure: <read object, current assigned class, history>
        where history corresponds to a list of ids of classes to which the read was assigned to in the past
        @classList list containing the current, existing classes
    """
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
        self.list = list
        self.classList = classList
        self.haploCount = haploCount


class HaplotypeClass:

    """
    Represents a class in the error correction, clustering approach
    @id id of the class
    @haplotype sequence of the class
    @list a list containing the reads assigned to the class
    @count sum of the frequencies of the reads assigned to the class
    """
    def __init__(self, id, haplotype=None, list=None, count=0):
        if list is None:
            list = []
        if haplotype is None:
            haplotype = ""
        self.id = id
        self.haplotype = haplotype
        self.list = list
        self.count = count
