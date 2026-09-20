import math
import sys
import time
import json
from collections import defaultdict, OrderedDict

#from memory_profiler import memory_usage

import Pdata
pdata = Pdata.processingData()
from memory_profiler import memory_usage

CanNum = 0  #候选模式数量

# ============================================================================
# [位置字典改造说明]
# 本文件的支持度计算方法（单项集用交集 Matching_I、多项集用 DFOM 的 Nettree
# 非重叠匹配算法）与原始 DFOM-HUP 完全一致，未做任何改动；唯一的改动是仿照
# HAUP 的“完整路径可追溯位置索引”设计，把位置字典 ItemS 从“只保存标量位置”
# 改为“保存完整路径元组”，即：
#     ItemS['[\'a\']']       -> 每条序列下形如 (p,) 的位置元组列表（单项）
#     ItemS["['a', 'b']"]   -> 每条序列下形如 (p,) 的位置元组列表（单itemset
#                              内多个item通过Matching_I交集得到，仍是长度1，
#                              因为I-连接不产生新的路径节点，语义上和一个
#                              itemset对应"同一个位置"是一致的）
#     ItemS["[['a'], ['a','b']]"]
#                            -> 每条序列下形如 (p0, p1, ...) 的完整路径元组
#                              列表（多个itemset依次出现，路径长度=itemset数，
#                              由 DFOM 在原有Nettree非重叠匹配过程中顺带记录）
#
# 原版DFOM/Matching_I只返回support计数，从不把size>=2的多itemset模式的出现
# 位置写回ItemS（注释里明确写了"无需存储位置"，Join_S/Mine_Pattern里对应的
# ItemS写入/删除语句也都是被注释掉的），因为DFOM每次都是直接基于size=1的
# 单itemset位置从头计算，不依赖之前保存的size>=2模式的位置。本次修改不改变
# 这一计算方式本身，只是让DFOM在原有的Nettree扫描过程中，"顺手"把每一次
# 成功匹配到的完整链路记录下来一并返回，因此：
#   - Matching_I 的交集匹配条件、DFOM 的 Nettree 构建/回溯条件（含
#     "unit[j][i][k] > Nettree[j-1][i][-1] and unit[j][i][k] > aaa" 这套
#     非重叠占用判据）逐字保留，不做任何调整；
#   - 唯一变化是数组里存的元素从标量int换成了元组，比较/占用判据里凡是需要
#     标量的地方，统一取元组的最后一位（_last()），数值上与原来的标量完全
#     等价；
#   - DFOM 新增返回值 list3（完整路径），Mine_Pattern / Join_S 里不再丢弃
#     多itemset模式的位置，而是和单itemset模式一样写入 ItemS，便于后续和
#     HAUP 做位置可追溯性上的公平对比。
# ============================================================================

def _last(x):
    """从一个位置记录里取出用于数值比较/拼接的标量值。
    位置记录现在统一是元组（如 (5,) 或 (5, 9)），取其最后一位；
    兼容 Nettree 内部用作"占用"哨兵值的裸整数 -1（未曾被赋值时的初始状态）。"""
    return x[-1] if isinstance(x, tuple) else x

#计算单项集模式的支持度（I-连接：同一itemset内多个item要求出现在完全相同的位置，
#因此用位置取值的交集即可；交集结果长度恒为1，只是统一包装成 (pos,) 元组，
#与HAUP的Matching_I"同一itemset扩展不产生新路径节点"的语义保持一致）
def Matching_I(list1, list2):
    list3 = [[] for k in range (SeqNum)]
    count = 0
    for i in range (SeqNum):
        set1 = set(_last(t) for t in list1[i])
        set2 = set(_last(t) for t in list2[i])
        common = sorted(set1 & set2)
        list3[i] = [(pos,) for pos in common]
        count += len(list3[i])
    # print(list3)
    # print(count)
    return count, list3

#计算多项集模式的支持度（S-连接：多个itemset依次非重叠出现，用Nettree贪心
#匹配。匹配/占用判据与原版逐字一致，只是位置记录统一为完整路径元组，
#因此可以在计数的同时，把每一次成功匹配到的完整路径顺带记录并返回）
def DFOM(pattern, ItemS):
    # print(str(pattern))
    count = 0
    Nettree = [[[] for i in range(SeqNum)] for k in range(len(pattern))]
    unit = [[] for k in range(len(pattern))]
    for i in range(len(pattern)):
        unit[i] = ItemS[str(pattern[i])]
    # print(unit)
    # 完整路径结果：list3[i] 保存该序列下每一次成功匹配得到的完整位置路径，
    # 每个路径是长度=len(pattern)的元组，第t个位置对应pattern[t]这个itemset
    # 被匹配到的位置
    list3 = [[] for i in range(SeqNum)]
    for i in range(SeqNum):
        for m in range(len(unit[0][i])):
            bbb = 0
            Nettree[0][i].append(unit[0][i][m])
            for j in range(1, len(unit)):
                n = 1
                if unit[j][i] == []:
                    bbb = 1
                    break
                else:
                    for k in range(len(unit[j][i])):
                        if Nettree[j][i] == []:
                            aaa = -1
                        else:
                            aaa = Nettree[j][i][-1]
                        if _last(unit[j][i][k]) > _last(Nettree[j - 1][i][-1]) and _last(unit[j][i][k]) > _last(aaa):
                            Nettree[j][i].append(unit[j][i][k])
                            if j == len(unit) - 1:
                                count += 1
                                # 本次m从unit[0][i][m]出发的这条链路已经完整匹配到
                                # 最后一个itemset，此刻 Nettree[t][i][-1] (t=0..j)
                                # 依次就是这条链路在各个itemset上选中的位置，
                                # 拼接成完整路径记录下来
                                full_path = tuple(_last(Nettree[t][i][-1]) for t in range(len(unit)))
                                list3[i].append(full_path)
                            n = 0
                            break
                    if n:
                        bbb = 1
                        break
            if bbb:
                break
    # print(str(count))
    return count, list3


#挖掘频繁单项集模式（m长度 m>1，1大小）例如：[abc]
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
                    count, ItemS[str(pattern)] = Matching_I(ItemS[str(p)], ItemS[str(qq)])
                    if count >= int(minsup):
                        compute_utility([pattern], count)
                        FP.append(pattern)
                        Temp.append(pattern)
                        # print(str(pattern) + ":" + str(count))
                        # print(ItemS[str(pattern)])
                    else:
                        del ItemS[str(pattern)]
        ExpSet = Temp[:]

#挖掘频繁单项集模式（2长度。1大小）例如：[ab]
def Gen_I(FP, ItemS):
    global CanNum
    Item = FP[:]
    ExpSet = []   #可扩展集
    for i in range (len(Item)):
        for j in range (i+1, len(Item)):
            pre = Item[i]
            suf = Item[j]
            p = []
            p.append(pre[0])
            p.append(suf[0])
            # print(p)
            # ItemS[str(p)] = [[] for k in range(SeqNum)]
            CanNum += 1
            count, ItemS[str(p)] = Matching_I(ItemS[str(pre)], ItemS[str(suf)])
            if count >= int(minsup):
                compute_utility([p], count)
                FP.append(p)
                ExpSet.append(p)
                # print(p)
                # print(ItemS[str(p)])
                # print(str(p) + ":" + str(ItemS[str(p)]) + ': ' + str(count))
                # print(str(p) + ":" + str(count))
            else:
                del ItemS[str(p)]
    # print(FP)
    # print(ExpSet)
    # print(ItemS)
    Join_I(FP, ExpSet, ItemS)

# 挖掘频繁单项（1长度，1大小）例如：[a]
def Mine_ItemS(FP, ItemS): #FP:频繁模式集(size=1)  ItemS:频繁单项集的出现位置字典
    global CanNum
    for i in sort_item:  # a b c d e f
        CanNum += 1
        count = 0
        for j in range(SeqNum):
            count += len(S[i][j])
        # print(i + ':' + str(count))
        if count >= int(minsup):
            p = []
            p.append(i)
            compute_utility([p], count)
            FP.append(p)
            # 初始化为完整路径：每个原始位置包装成长度为1的元组 (pos,)，
            # 与HAUP的Mine_ItemS单项初始化保持一致
            ItemS[str(p)] = [[] for k in range(SeqNum)]
            for j in range(SeqNum):
                ItemS[str(p)][j] = [(pos,) for pos in S[i][j]]
            # print(str(p) + ":" + str(ItemS[str(p)]) + ': ' + str(count))
            # print(str(p) + ":" + str(count))
    # print(FP)
    # del S
    Gen_I(FP, ItemS)

## 挖掘size>2的多项集模式（通过模式连接得到）
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
                    # count, ItemS[str(pattern)] = Matching_S(ItemS[str(p)], ItemS[str(q[-1])])
                    count, itemS_val = DFOM(pattern, ItemS)
                    # print(str(pattern) + ': ' + str(ItemS[str(pattern)]) + ': ' + str(count))
                    if count >= int(minsup):
                        # 与其余分支保持一致：写入完整路径，供下游按需查询/对比
                        ItemS[str(pattern)] = itemS_val
                        compute_utility(pattern, count)
                        FP.append(pattern)
                        Temp.append(pattern)
                        # print(str(pattern) + ': ' + str(count))
                        # print(str(pattern) + ': ' + str(ItemS[str(pattern)]) + ': ' + str(count))
                    # else:
                    #     del ItemS[str(pattern)]
        ExpSet = Temp[:]

# 挖掘size=2的多项集模式（通过两两拼接枚举得到）
def Mine_Pattern(FP, ItemS):
    global CanNum
    ExpSet = []
    temp = FP[:]
    for pre in temp:
        for suf in temp:
            # pattern = []
            # pattern.append(pre)
            # pattern.append(suf)
            pattern = [pre, suf]
            CanNum += 1
            # count, ItemS[str(pattern)] = Matching_S(ItemS[str(pre)], ItemS[str(suf)])
            count, itemS_val = DFOM(pattern, ItemS)
            if count >= int(minsup):
                # 与其余分支保持一致：写入完整路径，供下游按需查询/对比
                ItemS[str(pattern)] = itemS_val
                FP.append(pattern)
                compute_utility(pattern, count)
                ExpSet.append(pattern)
                # print(str(pattern) + ': ' + str(ItemS[str(pattern)]) + ': '+ str(count))
                # print(str(pattern) + ': ' + str(count))
            # else:
            #     del ItemS[str(pattern)]
    # print(FP)
    # print(ExpSet)
    Join_S(FP, ExpSet, ItemS)


def Miner():
    FP = []
    ItemS = {}
    Mine_ItemS(FP, ItemS) #挖掘频繁单项集模式，例如：[a], [ab], [abc]... 出现位置存储在ItemS[str(itemset)]
    # print("Frequent patterns with size=1:" + str(FP))
    # print("Number of frequent patterns with size=1:" + str(len(FP)))
    # print("Number of candidate patterns with size=1:" + str(CanNum))
    # # print(ItemS)
    Mine_Pattern(FP, ItemS) #挖掘频繁多项集模式, (size>=2). 例如：[a][ab], [a][a][ab]...
    # print("Frequent patterns:" + str(FP))
    print("Number of frequent patterns:" + str(len(FP)))
    print("Number of candidate patterns:" + str(CanNum))

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
    # 向上取整
    minsup = math.ceil(int(minau) / maxau)

    # print(Utility)


# @ti.kernel

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

        # print(pattern)
        # print(au)


if __name__ == '__main__':
    try:
        readFileName = sys.argv[1]
        # minsup = int(sys.argv[2])
    except Exception as e:
        print(e)
        exit(0)
        # readFileName = "../data/SDB2.txt"
    minsup = 0
    # minau = 50000
    AUNP = []
    S = {}

    last_index = readFileName.rindex(".")

    dataFileName = readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"
    # print(utilityFileName)
    # getUtility(utilityFileName)
    SeqNum, S, sort_item = pdata.datap(readFileName, S)
    del pdata
    for minau in sys.argv[2:]:
        getUtility(utilityFileName)
        print('DFOM-AUNP:', readFileName, 'minau=', minau, ':')

        starttime = time.time()
        Miner()
        #a = memory_usage((Miner))
        endtime = time.time()
        #print('Memory usage: ', max(a) - min(a))
        print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
    # print('AUNP-Miner:', readFileName, 'minau=', minau, ':')
    # starttime = time.time()
    # # Miner()
    # a = memory_usage((Miner))
    # endtime = time.time()
    # print('Memory usage: ', max(a) - min(a))
    # print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
