import math
import sys
import time
import Pdata
pdata = Pdata.processingData()
from memory_profiler import memory_usage
class ntnode():
    def __int__(self, position=0):
        self.position = position
        self.snode = []
        self.pnode = []
a = ntnode()
a.__int__(0)
def NetGap(pattern):
    count = 0
    for i in range(SeqNum):
        Nettree = [[] for k in range(len(pattern))]
        CreatNettree(Nettree, pattern, sdb[i])
        UpdateNettree(Nettree)
        while Nettree[0] != []:
            count += 1
            cuttree(Nettree[0][0],Nettree,0)
            UpdateNettree(Nettree)
    return count
def cuttree(node,nettree,levl):
    for i in node.snode:
        i.pnode.remove(node)
    for i in node.pnode:
        i.snode.remove(node)
    if node.snode != []:
        cuttree(node.snode[0],nettree,levl+1)
    nettree[levl].remove(node)
def IfExist(list1, list2):
    for i in list1:
        if i not in list2:
            return 0
    return 1
def CreatNettree(Nettree, pattern, seq):
    for j in range(len(pattern)):
        for i in range(len(seq)):
            if IfExist(pattern[j], seq[i]):
                node = ntnode()
                node.position = i
                node.snode = []
                node.pnode = []
                Nettree[j].append(node)
    for i in range(len(Nettree)-1):
        for j in Nettree[i]:
            k = 0
            while k < len(Nettree[i+1]):
                if j.position < Nettree[i + 1][k].position:
                    j.snode.append(Nettree[i + 1][k])
                    Nettree[i + 1][k].pnode.append(j)
                else:
                    if Nettree[i + 1][k].pnode == []:
                        Nettree[i + 1].remove(Nettree[i + 1][k])
                        k = k - 1
                k = k + 1
def UpdateNettree(Nettree):
    i = len(Nettree) - 2
    while i >= 0:
        j = 0
        while j < len(Nettree[i]):
            j += 1
            if Nettree[i][j-1].snode == []:
                if Nettree[i][j-1].pnode != []:
                    for k in Nettree[i][j-1].pnode:
                        k.snode.remove(Nettree[i][j-1])
                Nettree[i].pop(j-1)
                j -= 1
        i -= 1
def ShowNettree(Nettree):
    print('%%%%%%%')
    for i in Nettree:
        for j in i:
            print(j.position, end=' ')
        print('\n')
    print('%%%%%%%%%%%%%%%%%')
CanNum = 0
def GetUnit(itemsets):
    count = 0
    unit = GetIndex(itemsets[0])[:]
    for i in range(1, len(itemsets)):
        index_i = GetIndex(itemsets[i])
        for j in range(SeqNum):
            unit[j] = sorted(list(set(unit[j]) & set(index_i[j])))
    for i in range(SeqNum):
        count += len(unit[i])
    return count, unit
def ProMatching(list1, itemsets):
    a, list2 = GetUnit(itemsets)
    list3 = [[] for k in range (SeqNum)]
    flag = 0
    count = 0
    for i in range(SeqNum):
        for j in range(len(list1[i])):
            if flag == len(list2[i]):
                break
            for k in range(flag, len(list2[i])):
                if list2[i][k] > list1[i][j]:
                    list3[i].append(list2[i][k])
                    count += 1
                    flag = k + 1
                    break
                if k == len(list2[i]) - 1:
                    flag = len(list2[i])
        flag = 0
    return count, list3
def Join_I(FP, ExpSet, ItemS):
    global CanNum
    while ExpSet != []:
        Temp = []
        for p in ExpSet:
            sufP = p[1:]
            for q in ExpSet:
                preQ = q[:len(q)-1]
                if sufP == preQ:
                    pattern = p[:]
                    pattern.append(q[-1])
                    CanNum += 1
                    qq = []
                    qq.append(q[-1])
                    count, ItemS[str(pattern)] = GetUnit(pattern)
                    if count >= int(minsup):
                        FP.append(pattern)
                        compute_utility([pattern], count)
                        Temp.append(pattern)
                    else:
                        del ItemS[str(pattern)]
        ExpSet = Temp[:]
def Gen_I(FP, ItemS):
    global CanNum
    Item = FP[:]
    ExpSet = []
    for i in range (len(Item)):
        for j in range (i+1, len(Item)):
            pre = Item[i]
            suf = Item[j]
            p = []
            p.append(pre[0])
            p.append(suf[0])
            CanNum += 1
            count, ItemS[str(p)] = GetUnit(p)
            if count >= int(minsup):
                compute_utility([p], count)
                FP.append(p)
                ExpSet.append(p)
            else:
                del ItemS[str(p)]
    Join_I(FP, ExpSet, ItemS)
def Mine_ItemS(FP, ItemS):
    global CanNum
    for i in sort_item:
        CanNum += 1
        count = 0
        index_i = GetIndex(i)
        for j in range(SeqNum):
            count += len(index_i[j])
        if count >= int(minsup):
            p = []
            p.append(i)
            compute_utility([p], count)
            FP.append(p)
            ItemS[str(p)] = [[] for k in range(SeqNum)]
            ItemS[str(p)] = index_i
    Gen_I(FP, ItemS)
def Join_S(FP, ExpSet, ItemS):
    global CanNum
    while ExpSet != []:
        Temp = []
        for p in ExpSet:
            for q in ExpSet:
                sufP = p[1:]
                preQ = q[:len(q)-1]
                if sufP == preQ:
                    pattern = p[:]
                    pattern.append(q[-1])
                    CanNum += 1
                    count = NetGap(pattern)
                    if count >= int(minsup):
                        compute_utility(pattern, count)
                        FP.append(pattern)
                        Temp.append(pattern)
        ExpSet = Temp[:]
def Mine_Pattern(FP, ItemS):
    global CanNum
    ExpSet = []
    temp = FP[:]
    for pre in temp:
        for suf in temp:
            pattern = [pre, suf]
            CanNum += 1
            count = NetGap(pattern)
            if count >= int(minsup):
                compute_utility(pattern, count)
                FP.append(pattern)
                ExpSet.append(pattern)
    Join_S(FP, ExpSet, ItemS)
def Miner():
    FP = []
    ItemS = {}
    Mine_ItemS(FP, ItemS)
    Mine_Pattern(FP, ItemS)
    print("Number of AUNP:" + str(len(AUNP)))
    print("Number of frequent patterns:" + str(len(FP)))
    print("Number of candidate patterns:" + str(CanNum))
def GetIndex(item):
    index2 = []
    for seq in sdb:
        index = []
        for i in range(len(seq)):
            if item in seq[i]:
                index.append(i)
        index2.append(index)
    return index2
def ReadFile(readFileName):
    SeqNum = 0
    sdb = []
    sort_item = []
    itemset = []
    seq = []
    with open(readFileName, 'r') as f:
        lines = f.readlines()
        for line in lines:
            SeqNum += 1
            items_array = line.strip().split(' ')
            for item in items_array:
                if item != '-1':
                    itemset.append(item)
                    if item not in sort_item:
                        sort_item.append(item)
                else:
                    seq.append(itemset)
                    itemset = []
            seq.append(itemset)
            itemset = []
            sdb.append(seq)
            seq = []
        sort_item.sort()
    return SeqNum, sdb, sort_item
Utility = {}
def getUtility(utilityFileName):
    global Utility
    global minsup
    maxau = 0
    pdata = Pdata.processingData()
    lines = pdata.read_file(utilityFileName)
    for line in lines:
        str = line.split(' ')
        Utility[str[0]] = int(str[1])
        if int(str[1]) > maxau:
            maxau = int(str[1])
    minsup = math.ceil(int(minau) / maxau)
def compute_utility(pattern, count):
    global AUNP
    au = 0
    len = 0
    for i in pattern:
        for j in i:
            len += 1
            au += Utility[str(j)]
    au = au * count / len
    if au >= int(minau):
        AUNP.append(pattern)
if __name__ == '__main__':
    try:
        readFileName = sys.argv[1]
    except Exception as e:
        print(e)
        exit(0)
    minsup = 0
    AUNP = []
    S = {}
    last_index = readFileName.rindex(".")
    dataFileName = readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"
    SeqNum, sdb, sort_item = ReadFile(readFileName)
    del pdata
    for minau in sys.argv[2:]:
        getUtility(utilityFileName)
        print('NPSEP-AUNP:', readFileName, 'minau=', minau, ':')
        starttime = time.time()
        Miner()
        endtime = time.time()
        print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")