import math
import sys
import time
import json
from collections import defaultdict, OrderedDict
import Pdata

from memory_profiler import memory_usage

# ============================================================================
# [位置字典改造说明]
# 本文件的挖掘策略（S-连接 / I-连接 的生成方式、广度优先的候选扩展顺序）与原始
# HUP-B 完全一致，未做任何改动；唯一的改动是仿照 HAUP 的“完整路径可追溯位置索引”
# 设计，把 P（位置字典）从“只保存每个模式在每条序列中最后一次匹配到的标量位置”
# 改为“保存每个模式从长度1开始累积的完整位置路径（元组）”，即
#     P['a']            -> 每条序列下形如 (p,) 的位置元组列表
#     P['a -1 b']        -> 每条序列下形如 (p_a, p_b) 的位置元组列表
#     P['a b']           -> I-连接不新增路径长度，仍为 (p,) （与 HAUP 的
#                            Matching_I 语义一致：同一itemset内扩展不产生新的
#                            路径节点）
# 这样修改后，HUP-B 与 HAUP 在“位置如何被记录/追溯”这一点上完全对齐，两者的
# 差异就只剩下候选模式的连接与扩展策略本身，可以用于公平地对比“模式连接策略”
# 的有效性。
#
# 匹配判据本身（S-连接的 “>” 判断、I-连接借助 LP 的 “-1 偏移技巧” 判断）
# 全部原样保留，没有做任何“修正”，以确保支持度计数与原版 HUP-B 完全一致，
# 只是把“记录下来的位置”从标量升级为完整路径。
# LP 仍然保持纯标量数组（只用于内部比较，不对外输出），其取值统一改为
# 从新的完整路径元组中取最后一个元素（P[...][t][-1]），语义上等价于原来的
# P[...][t]，因此比较逻辑本身不需要任何调整。
# ============================================================================

# 挖掘（长度为1的频繁-》S-连接-》I-连接）
def Min_FP(SeqNum, S, sort_item, minsup):
    # print("%%",sort_item)
    P = {}  #存储模式的完整出现路径：P[pattern][j] = [(p1,), (p1,p2), ...]（每个位置为一个元组）
    # '[ac][a]': [[(1,) (2,) (4,)], [], [],...]
    LP = {} #存储模式的项集子模式的比较用标量（内部使用，取自完整路径的最后一个位置）
    # [ac][a]
    # LP['[ac][ac]'] = P['[ac]']（取最后一位）
    FP = [] #存储频繁模式
    FP1=[]  #存储长度为1的频繁模式 用于模式增长 a c
    # minsup = 2
    CanNum = 0
    # SeqNum = len(lines_s)
    count = 0
    '''
    S =
    {
    'a': [[0, 1, 3], [1], [0, 2, 4], [0, 1], [1, 3]],
    'b': [[], [0], [2], [0], [4]],
    'c': [[0, 1, 2, 3, 4], [1, 2], [0, 2, 3, 4, 5], [0, 1, 2, 3], [0, 1, 2, 3, 5]],
    'd': [[], [3], [0], [], [0, 3]],
    'e': [[], [1, 3], [4, 5], [3, 4], [5]],
    'f': [[], [0], [1, 4], [4], []]
    }
    '''
    for i in sort_item:  # a b c d e f
        CanNum += 1
        for j in range(SeqNum):
            count += len(S[i][j])
        if count >= minsup:
            FP1.append(i)
            FP.append(i)
            compute_utility(i, count)
            # 初始化为完整路径：每个原始位置包装成长度为1的元组 (pos,)
            P[i] = [[] for k in range(SeqNum)]
            for k in range(SeqNum):
                for pos in S[i][k]:
                    P[i][k].append((pos,))
            LP[i] = [[] for k in range(SeqNum)]
            for k in range(SeqNum):
                for t in range(len(P[i][k])):
                    LP[i][k].append(P[i][k][t][-1] - 1)
        count = 0
    flag = 0
    for fp in FP: # fp: str
        for item in FP1: # item: str
            pattern = fp + " -1 " + item # S连接**************************************
            CanNum += 1
            LP[pattern] = [[] for k in range(SeqNum)]
            # LP[pattern] 与原版一致，取 fp 的（未偏移）标量位置，
            # 只是现在要从完整路径元组里取最后一位
            LP[pattern] = [[pos[-1] for pos in P[fp][k]] for k in range(SeqNum)]
            P[pattern] = [[] for i in range(SeqNum)]
            for i in range(SeqNum):
                for j in range(len(P[fp][i])):
                    last_pos = P[fp][i][j][-1]
                    if flag == len(S[item][i]):
                        # flag = 0
                        break
                    for k in range(flag, len(S[item][i])):
                        if S[item][i][k] > last_pos:
                            # 保存完整路径：在fp的完整路径末尾追加新匹配到的位置
                            P[pattern][i].append(P[fp][i][j] + (S[item][i][k],))
                            count += 1
                            flag = k + 1
                            break
                        if k == len(S[item][i]) - 1:
                            flag = len(S[item][i])
                flag = 0
            if count >= minsup:
                FP.append(pattern)
                compute_utility(pattern, count)
                # 输出模式，出现位置，和支持度
                # print (pattern + ': ', end = "")
                # print (P[pattern])
                # print ('count=' + str(count))
                # print ("**********")
            else:
                del P[pattern]
                del LP[pattern]
            count = 0

            # [a b -1 a b c]   a b c d
            if fp.split()[len(fp.split())-1] < item: # I连接************************************
                CanNum += 1
                pattern = fp + ' ' + item
                LP[pattern] = [[] for i in range(SeqNum)]
                LP[pattern] = LP[fp]
                P[pattern] = [[] for i in range(SeqNum)]
                unit = S[item][:]    # [:]:复制字典中的内容
                pattern_array = fp.strip().split(' ')
                f = -1
                for i in range(len(pattern_array)):
                    if pattern_array[i] == '-1':
                        f = i
                for i in range(f+1, len(pattern_array)):
                    for j in range(len(unit)):
                        unit[j] = sorted(list(set(unit[j]) & set(S[pattern_array[i]][j])))
                for i in range(SeqNum):
                    for j in range(len(LP[fp][i])):
                        if flag == len(unit[i]):
                            break
                        for k in range(flag, len(unit[i])):
                            if unit[i][k] > LP[fp][i][j]:
                                # I-连接不产生新的路径节点（与HAUP的Matching_I语义一致：
                                # 同一itemset内的扩展仍是同一个"位置"），但沿用原版语义，
                                # 用新匹配到的unit值替换路径末位，其余历史路径原样保留
                                P[pattern][i].append(P[fp][i][j][:-1] + (unit[i][k],))
                                count += 1
                                flag = k + 1
                                break
                            if k == len(unit[i]) - 1:
                                flag = len(unit[i])
                    flag = 0
                if count >= minsup:
                    FP.append(pattern)
                    compute_utility(pattern, count)
                    # print(pattern + ': ', end="")
                    # print(P[pattern])
                    # print('count=' + str(count))
                    # print("**********")
                else:
                    del P[pattern]
                    del LP[pattern]
            count = 0
        del P[fp]
        del LP[fp]
    # print ("Frequent patterns:", end = " ")
    # print (FP)
    print ("Number of frequent patterns: " + str(len(FP)))
    print ("Number of candidate patterns: " + str(CanNum))
    print("Number of AUNP:" + str(len(AUNP)))

Utility = {}
def getUtility(utilityFileName):
    global Utility
    global minsup
    pdata1 = Pdata.processingData()
    maxau = 0
    lines = pdata1.read_file(utilityFileName)
    for line in lines:
        str = line.split(' ')
        Utility[str[0]] = int(str[1])
        if int(str[1]) > maxau:
            maxau = int(str[1])
    #向上取整
    minsup = math.ceil(int(minau)/maxau)

    #print(Utility)
# @ti.kernel

def compute_utility(pattern,count):
    global AUNP
    au = 0
    len = 0
    p = pattern.split(" -1 ")
    for i in p:
        j = i.split(" ")
        #print(j)
        for K in j:
            len += 1
            au += Utility[str(K)]
    au = au*count/len
    if au>=int(minau):
        AUNP.append(pattern)

if __name__ == '__main__':
    try:
        readFileName = sys.argv[1]
        # minsup = int(sys.argv[2])
    except Exception as e:
        print(e)
        exit(0)
    # readFileName = "../data/SDB2.txt"
    pdata = Pdata.processingData()
    minsup = 0
    # minau = 50000
    AUNP = []
    S = {}
    last_index = readFileName.rindex(".")

    dataFileName =readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:]+"_utility"+".txt"
    # print(utilityFileName)
    # getUtility(utilityFileName)
    SeqNum, S, sort_item = pdata.datap(readFileName, S)
    del pdata
    # print('NFP-B:', readFileName, 'minsup=', minsup, ':')
    # starttime = time.time()
    # a = memory_usage((Min_FP, (int(SeqNum), S, sort_item, int(minsup))))
    # # Min_FP(SeqNum, S, sort_item, int(minsup)) # lines_s:替换后的序列 S:位置字典 sort_item:字符集[a, b, c, d, e, f]
    # endtime = time.time()
    # print('Memory usage: ', max(a) - min(a))
    # print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
    for minau in sys.argv[2:]:
        getUtility(utilityFileName)
        print('AUNP-B:', readFileName, 'minau=', minau, ':')
        starttime = time.time()
        #a = memory_usage((Min_FP, (int(SeqNum), S, sort_item, int(minsup))))
        Min_FP(SeqNum, S, sort_item, int(minsup)) # lines_s:替换后的序列 S:位置字典 sort_item:字符集[a, b, c, d, e, f]
        endtime = time.time()
        #print('Memory usage: ', max(a) - min(a))
        print ("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
        # print ("Running time: " + str(int(round(endtime)) - int(round(starttime))) + "s")
        # with open(filename, 'w') as f:
        #     for i in range(len(lines_s)):
        #         f.writelines(str(i) + "\t" + lines_s[i])
        #         f.write("\n")
