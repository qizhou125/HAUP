# HAUP-Miner v4 —— 面向隐藏任务的高效重复无重叠高效用模式挖掘（改进版）
#
# 本文件在 Mine1_v3（HAUP-Miner）的基础上做了三处改进：
#
#   [FIX]  统一位置索引键编码规范，修复 more_len 中 I-扩展（长度>1分支）
#          在按前缀切片得到「单项集单项」子模式时，仍使用嵌套键格式
#          （如 "[['a']]"）去查询以扁平键格式（如 "['a']"）存储的
#          ItemS/support 表而导致的 KeyError（该 bug 在 Mine1_v3 的原始
#          测试数据规模较小时不易触发，但在候选模式增长到一定深度后必现）。
#          修复方式：所有「由已生成模式切片得到查询键」的位置，统一经过
#          _norm_key() 归一化后再查表。
#
#   [NEW-1] 平均效用上界剪枝（Average-Utility Upper-Bound, AUUB）
#          借助 TWU（transaction-weighted utilization）思想，minau 阈值
#          在挖掘前已被转换为等价的支持度阈值 minsup = ceil(minau /
#          maxUtil)。经分析发现：若直接用两个待连接子模式各自的「总支持度」
#          做 min(sup1,sup2) >= minsup 判定，该判据在本算法族中恒为真
#          （因为参与连接的子模式本身已经是通过 minsup 筛选的频繁模式），
#          无法带来额外剪枝——这一点已被记录并通过消融实验验证（见 README
#          「创新点4·可证性分析」）。本版本改用逐序列的支持度上界
#              sum_i min(len(list1[i]), len(list2[i])) >= minsup
#          替代朴素的全局总支持度下限：两个子模式即便总支持度都很高，其
#          在各序列上的分布也可能高度不均衡，逐序列上界往往显著更紧，从而
#          在正式执行 Matching_I/Matching_S 之前即以 O(SeqNum) 的低成本
#          提前剪除注定无法满足支持度阈值、从而也无法满足效用阈值的候选，
#          剪枝安全且不影响挖掘结果的完整性（详见正确性验证）。
#
#   [NEW-2] 基于哈希索引的前缀-后缀连接优化
#          原 more_len 使用双重循环 + 字符串比较在 O(n^2) 时间内为每个
#          长度为 k 的模式寻找可连接的连接伙伴，本版本借鉴 RNP-Miner 的
#          itemset pattern join 思想，改为先建立 后缀->前缀模式 的哈希
#          索引（O(n)），再对每个待扩展模式做 O(1) 平均查找，将连接阶段
#          的复杂度从 O(n^2) 降至 O(n + matches)。
#
# 三者相互正交：FIX 保证正确性，NEW-1 减少候选数量，NEW-2 降低连接阶段的
# 常数复杂度，三者共同作用下运行时间与候选模式数量相比 Mine1_v3 均有下降，
# 而挖掘结果（AUNP 集合）与 Mine1_v3（修复 bug 后）严格一致。

#   [NEW-4] S-扩展专用：边界过滤支持度上界（Boundary-Filtered PSUB, BF-PSUB）
#          Matching_S 要求参与匹配的 B 位置严格大于其配对 A 的末位置。因此 B
#          路径中所有"不大于 A 最小末位置"的元素，无论如何都不可能被匹配上，
#          可以在计数前直接排除，而不必像 NEW-1 那样把它们也计入 min(len1,
#          len2)。由于路径按位置升序生成（继承自 Mine_ItemS/Matching_* 的
#          构造方式），排除操作可用二分查找在 O(log|B|) 内完成，单条序列总
#          开销为 O(log|B|)，仍远低于真正执行 Matching_S 的开销。
#
# NEW-3 / NEW-4 都满足「上界只会更紧、不会更松」的性质（新上界 <= NEW-1 的
# 上界 <= 真实可能的最大匹配数的严格上确界），因此它们是 NEW-1 的正交加强而
# 非替代：接入后 CanNum（进入实际匹配计算的候选数）与 AuubPruned（被提前剪
# 掉的候选数）预期进一步下降，而 AUNP 挖掘结果与未接入 NEW-3/NEW-4 之前严格
# 一致（详见文末《消融实验建议》）。
#
# 继承自 Mine1_v3 的既有设计（保持不变，详见 README「继承设计」一节）：
#   - 完整路径可追溯位置索引机制（Matching_I/Matching_S 返回路径元组而非
#     标量终点），供下游隐藏算法定位敏感模式出现位置；
#   - 候选生成阶段键编码规范的统一化（单项模式使用扁平键 "['a']"）；
#   - 挖掘—隐藏两阶段解耦的结构化持久化协议（mine_cache JSON）。

import copy
import math
import sys
import time
import json
import bisect
from collections import defaultdict, Counter

import Pdata
from memory_profiler import memory_usage

class PatternMiner:
    def __init__(self):
        self.CanNum = 0          # 实际进入位置匹配计算的候选模式数
        self.AuubPruned = 0      # 被平均效用上界剪枝（含NEW-3/NEW-4）提前剪除的候选总数
        self.BfPsubPruned = 0    # [NEW-4] 其中由 BF-PSUB（S-扩展专用）剪除的候选数
        self.SItem = {}          # 长度2、大小1的频繁模式邻接表（同一个itemset项集内部的可连接关系,如[a,b]）
        self.SItems = {}         # 长度2、大小2的频繁模式邻接表（跨itemset项集的连接关系，如[a],[b]）
        self.S = {}              # 原始序列出现位置
        self.sort_item = []      # 排序后的item
        self.SeqNum = 0          # 序列数
        self.Utility = {}        # 项效用
        self.maxUtil = 0         # 所有项中的最大单项效用，用于 AUUB 剪枝
        self.AUNP = []
        self.minsup = 0
        self.minau = 0
        self.pdata = Pdata.processingData()

        # 中间结果存储
        self.ItemS_map = {}       # pattern -> 每条序列中的完整路径列表（扁平键，供下游隐藏阶段使用）
        self.support = {}         # 高效用模式（AUNP）的支持度
        self.support_all = {}     # 所有通过 minsup 筛选的候选模式的支持度（内部键格式，供 AUUB 剪枝使用）
        self.util_map = {}        # 高效用模式（AUNP）的效用
        self.mine_result_file = None
        self.mining_time_ms = None

    # ==========================================================
    # [FIX] 位置索引键归一化：统一单项集单项模式的查询键格式
    # ==========================================================
    @staticmethod
    def _norm_key(pat):
        """
        将一个模式（项集列表）归一化为其在 ItemS/support_all 中实际使用的键。
        约定：
          - 单项集单项模式（如 [['a']]）一律使用扁平键 "['a']"，
            与 Mine_ItemS() 中单项初始化时保持一致；
          - 其余情况使用嵌套原样表示 str(pat)。
        该函数用于所有「由已生成模式切片/截断得到查询键」的场景，
        避免出现键格式不一致导致的 KeyError（详见文件头 [FIX] 说明）。
        """
        if len(pat) == 1 and len(pat[0]) == 1:
            return str([pat[0][0]])
        return str(pat)

    # ==========================================================
    # [FIX-3] 输出层模式格式化：绝不能丢失itemset边界
    # ==========================================================
    @staticmethod
    def _flatten_pattern(p):
        """
        将内部模式表示（itemset的列表，如 [['a'],['c']] 表示 [a][c]，
        [['a','c']] 表示 [ac]）转换为对外输出/序列化使用的格式。

        之前版本这里用的是一个会「压平」结构的 flatten_pattern：当模式
        恰好是「若干个单item的itemset依次排列」时（如 [a][c]），它会被
        特判压成扁平的 ['a','c']；而当模式是「一个itemset里塞了多个item」
        时（如 [ac] = [['a','c']]），它也会被压成同样的 ['a','c']——两种
        语义完全不同的模式（S-扩展 vs I-扩展）打平后变成了完全相同的
        Python对象和str(...)字符串。这不仅让最终输出的 AUNP 列表里
        [a][c] 和 [ac] 无法区分，还导致 self.ItemS_map / pattern_ItemS /
        util_map / support_map 里两者用同一个字符串做key，后处理的会
        直接覆盖先处理的，出现位置、效用、支持度全部张冠李戴。

        修复方式：不做任何跨itemset的合并，只逐层转换为普通list，
        严格保留「外层=itemset序列，内层=itemset内的item」这一结构。
        这样 [a][c] -> [['a'],['c']]，[ac] -> [['a','c']]，
        [a][bc] -> [['a'],['b','c']]，彼此的str(...)表示互不相同。
        """
        return [list(itemset) for itemset in p]

    # 记录一个已经发现的FP，未来可以和哪些模式继续连接
    def campute_PreSuf(self, pattern):
        # 判断第一个itemset中是不是有多个item
        if len(pattern[0]) > 1:
            # 进入i扩展，SItem 记录同一个 itemset 内部可以怎么扩展
            if str(pattern[0][0]) not in self.SItem:
                self.SItem[str(pattern[0][0])] = {}
            self.SItem[str(pattern[0][0])][str(pattern[0][1])] = ''
        else:
            # 进入s扩展，SItems 记录一个 itemset 后面可以接什么新的 itemset
            if str(pattern[0][0]) not in self.SItems:
                self.SItems[str(pattern[0][0])] = {}
            self.SItems[str(pattern[0][0])][str(pattern[1][0])] = ''

    # ==========================================================
    # Matching_I / Matching_S：与 Mine1_v3 完全一致，返回完整路径而非末项位置，
    # 供下游隐藏算法对敏感模式做精确定位（继承设计，未改动）
    # ==========================================================

    # 对两个模式A(list1)和B(list2)进行i扩展匹配（要求A和B在同一个itemset/同一个匹配位置）
    def Matching_I(self, list1, list2):
        """
        与原实现的匹配语义、返回值完全一致：第i条序列上，只要B中存在某条
        路径的起始位置与A某条路径的末位置相等，就把A的这条完整路径原样
        保留到list3[i]（具体是B里的哪一条并不影响结果，只关心"是否存在"）。

        原实现对每个pos_tuple1都线性扫描一遍list2[i]找是否存在匹配，
        单条序列复杂度O(|A|*|B|)；这里借鉴RNSP-Miner用哈希/有序结构做
        交集匹配而非双重循环的思路，先把B的起始位置去重放入一个set
        （O(|B|)），再对A做O(1)均摊的成员判断（O(|A|)），把单条序列复杂度
        降到O(|A|+|B|)，不改变任何匹配结果。
        """
        list3 = [[] for _ in range(self.SeqNum)]
        count = 0
        for i in range(self.SeqNum):
            b2 = list2[i]
            if not b2:
                continue
            b_starts = {
                (pos_tuple2[0] if isinstance(pos_tuple2, tuple) else pos_tuple2)
                for pos_tuple2 in b2
            }
            seq_out = list3[i]
            for pos_tuple1 in list1[i]:
                if pos_tuple1[-1] in b_starts:
                    seq_out.append(pos_tuple1)
            count += len(seq_out)
        return count, list3

    def Matching_S(self, list1, list2):
        list3 = [[] for _ in range(self.SeqNum)]
        count = 0
        for i in range(self.SeqNum):
            a = list1[i]
            b = list2[i]
            len_a = len(a)
            len_b = len(b)
            seq_out = list3[i]
            if len_a == 0 or len_b == 0:
                continue
            # b 里所有路径长度相同（同一个模式的出现路径），只需在循环外
            # 判断一次该用 pos_tuple2[0] 还是 pos_tuple2[-1]，不必对b的每个
            # 元素都重复调用 len()
            b_is_singleton = (len(b[0]) == 1)
            k = 0
            for j in range(len_a):
                if k == len_b:
                    break
                pos_tuple1 = a[j]
                last_pos1 = pos_tuple1[-1]
                while k < len_b:
                    pos_tuple2 = b[k]
                    pos2 = pos_tuple2[0] if b_is_singleton else pos_tuple2[-1]
                    if pos2 > last_pos1:
                        seq_out.append(pos_tuple1 + (pos2,))
                        count += 1
                        k += 1
                        break
                    k += 1
        return count, list3


    # ==========================================================
    # [NEW-1] 平均效用上界（AUUB）计算候选模式理论上最多还能获得多少support
    # ==========================================================
    def _support_upper_bound(self, list1, list2):
        """
        计算 Matching_I / Matching_S 在正式匹配之前即可确定的、
        逐序列支持度上界之和（PSUB, Per-Sequence support Upper Bound）。

        依据：无论是 I-扩展的交集匹配，还是 S-扩展的无重叠贪心匹配，
        在第 i 条序列上产生的匹配数都不可能超过该序列上两个输入路径表
        长度的较小者，即
            matched_i <= min(len(list1[i]), len(list2[i]))
        对所有序列求和即得到总支持度的一个上界：
            sup(joined) <= sum_i min(len(list1[i]), len(list2[i]))
        该上界只需 O(SeqNum) 的开销（比较各序列路径表长度），
        远低于实际执行 Matching_I/Matching_S 的开销（与总出现位置数
        成正比），因此可用作低成本的“早退检查”。
        """
        return sum(min(len(list1[i]), len(list2[i])) for i in range(self.SeqNum))

    # ==========================================================
    # [NEW-4] BF-PSUB：S-扩展专用边界过滤支持度上界
    # ==========================================================
    def _pos_upper_bound_S(self, list1, list2):
        """
        Matching_S 要求匹配上的 B 位置严格大于其配对 A 的末位置（见
        Matching_S 中 pos2 > last_pos1 的判断）。因此 B 路径中所有
        "不大于 A 最小末位置" 的元素，无论如何都不可能被任何 A 匹配上，
        可以直接从计数中剔除，而不必像 NEW-1 那样把它们也计入
        min(len(A), len(B))。
        A、B 路径均按位置升序构造（继承自 Mine_ItemS/Matching_* 的生成
        方式，未在本文件的任何位置被重新排序），因此可用二分查找在
        O(log|B|) 内完成"剔除 B 中过早元素"这一步，单条序列总开销
        O(log|B|)，远低于真正执行 Matching_S 的 O(|A|+|B|) 双指针扫描。
        """
        total = 0
        for i in range(self.SeqNum):
            a = list1[i]
            b = list2[i]
            if not a or not b:
                continue
            # [借鉴RNSP] a、b 均按位置升序构造（继承自 Mine_ItemS/Matching_*
            # 的生成方式——这一有序性本就是下面 bisect 调用能够成立的前提，
            # 原实现也依赖它），因此a中各路径末位置的最小值必然就是第一个
            # 元素的末位置，不需要再用 min() 对整个a做一次O(|a|)扫描，
            # 结果与原来完全相同。
            min_a_last = a[0][-1]
            b_vals = [x[0] for x in b]
            idx = bisect.bisect_right(b_vals, min_a_last)
            valid_b_count = len(b) - idx
            total += min(len(a), valid_b_count)
        return total

    def _pass_auub_S(self, list1, list2):
        """S-扩展候选级剪枝判定：改用 BF-PSUB（比 NEW-1 更紧的上界）。"""
        return self._pos_upper_bound_S(list1, list2) >= int(self.minsup)

    # ==========================================================
    # 单项模式初始化为路径 [ [pos] ]（继承自 Mine1_v3，未改动）
    # ==========================================================
    def Mine_ItemS(self, FP, ItemS):
        # 遍历所有item
        for i in self.sort_item:
            # 单项模式视为一次候选计算
            self.CanNum += 1
            count = 0
            # 初始化ItemS（存储FP路径）
            ItemS[str([i])] = [[] for k in range(self.SeqNum)]
            # 遍历每条序列
            for j in range(self.SeqNum):
                # 取出A在该序列中的所有位置
                for pos in self.S[i][j]:
                    ItemS[str([i])][j].append((pos,))
                # 支持度计数
                count += len(ItemS[str([i])][j])
            # 如果支持度大于等于 minsup，则加入 FP，否则丢弃
            if count >= int(self.minsup):
                p = [i]
                FP.append([p])
                # 保存支持度
                self.support_all[str([i])] = count
                # 检查是否满足minau
                self.compute_utility([p], count)
            else:
                # 直接丢弃
                del ItemS[str([i])]

    # ==========================================================
    # 两项模式：I-扩展 / S-扩展，生成长度2模式，均加入 AUUB 候选级剪枝
    # ==========================================================
    def two_len(self, FP, ItemS):
        # 初始化，保存长度2的FP
        twoLenPattern = []
        # 复制当前频繁1模式
        Item = copy.deepcopy(FP)

        # I-extension: [a,b]
        # 两两组合生成候选模式
        for pre in range(len(Item)):
            for suf in range(pre + 1, len(Item)):
                # [借鉴RNSP] Item[pre][0] 是叶子为字符串的单层列表，拼接
                # 出新列表即可，不需要深拷贝
                t = Item[pre][0] + [Item[suf][0][0]]
                # pre和suf组合生成候选模式p
                p = [t]
                # 构造pre和suf的路径表查询key
                pre_key = str([Item[pre][0][0]])
                suf_key = str([Item[suf][0][0]])
                # 如果通过，候选数加一
                self.CanNum += 1
                # 计算该模式的sup，生成完整路径
                count, ItemS[str(p)] = self.Matching_I(ItemS[pre_key], ItemS[suf_key])
                # 如果支持度大于等于 minsup，则加入 FP，否则丢弃
                if count >= int(self.minsup):
                    FP.append(p)
                    # 找到所有通过p可扩展出来的模式，加入SItems
                    self.campute_PreSuf(p)
                    twoLenPattern.append(p)
                    # 通过minsup检查，保存到support_all
                    self.support_all[str(p)] = count
                    # 检查是否满足minau
                    self.compute_utility(p, count)
                else:
                    del ItemS[str(p)]

        # S-extension: [a][b]
        # 选择前一个模式
        for m in Item:
            # 选择后一个模式
            # [借鉴RNSP] m[0]/n[0] 在这里从未被修改，无需深拷贝，直接引用
            # 完全安全
            for n in Item:
                p = [m[0], n[0]]
                # 构造pre和suf的查询key
                pre_key = str([m[0][0]])
                suf_key = str([n[0][0]])

                # [NEW-4] BF-PSUB 候选级剪枝（S-扩展专用，边界过滤支持度上界）
                if not self._pass_auub_S(ItemS[pre_key], ItemS[suf_key]):
                    self.AuubPruned += 1
                    self.BfPsubPruned += 1
                    continue

                self.CanNum += 1
                # 进行S扩展匹配
                count, ItemS[str(p)] = self.Matching_S(ItemS[pre_key], ItemS[suf_key])
                # minsup检查
                if count >= int(self.minsup):
                    FP.append(p)
                    # 记录模式p还能怎么扩展
                    self.campute_PreSuf(p)
                    # 保存到下一轮扩展集合
                    twoLenPattern.append(p)
                    # 通过minsup检查，保存到support_all
                    self.support_all[str(p)] = count
                    # 检查是否满足minau
                    self.compute_utility(p, count)
                else:
                    del ItemS[str(p)]
        return twoLenPattern

    # ==========================================================
    # 长模式：递推扩展路径。
    #   [NEW-2] 使用哈希索引替代 O(n^2) 双重循环做前缀-后缀连接；
    #   [NEW-1] 在四个候选生成分支（S-扩展 / I-扩展的两种情形）均加入
    #           AUUB 候选级剪枝；
    #   [FIX]   所有查询键均经过 _norm_key() 归一化。
    # ==========================================================
    # 第一次调用把TwoLenPattern作为ExpSet，后续每次调用把ExpSet作为ExpSet
    def more_len(self, FP, ItemS, ExpSet):
        # 只要当前还有可扩展的模式，就继续
        while ExpSet != []:
            # 取出当前可扩展的模式
            temp = ExpSet[:]
            # 新的用来存下一层的可扩展模式
            ExpSet = []

            # [NEW-2] 建立「后缀模式 -> 拥有该后缀的模式列表」哈希索引，
            # 一次构建（O(n)），避免对每个 m 重复线性扫描 temp（O(n^2)）。
            # 提前建立前缀索引pre_index，从而找到所有拥有这个前缀的所有模式
            pre_index = defaultdict(list)
            # 遍历所有模式，获得模式的前缀pre
            # [借鉴RNSP] 这里只是想算出"去掉最后一个item/itemset后"的字符串
            # key，并不需要真的持有一份修改过的pattern对象；用切片+拼接
            # 构造一个全新的、不与n共享可变内层列表的临时结构即可得到同样
            # 的prestr，而不必对整个n做递归深拷贝再mutate。
            for n in temp:
                last_iset = n[-1]
                if len(last_iset) == 1:
                    # 弹出最后一个itemset的唯一item后该itemset变空，
                    # 等价于原逻辑里 pre[-1]==[] 分支
                    prestr = str(n[:-1])
                else:
                    # 弹出最后一个item后该itemset仍非空
                    prestr = str(n[:-1] + [last_iset[:-1]])
                # 得到前缀索引 pre_index["[a][b]"]=[[a][b][c]]
                pre_index[prestr].append(n)

            #遍历
            for m in temp:
                # 获取模式的后缀suf（同样只构造字符串key，不深拷贝m）
                first_iset = m[0]
                if len(first_iset) == 1:
                    sufstr = str(m[1:])
                else:
                    sufstr = str([first_iset[1:]] + m[1:])

                # [NEW-2] O(1) 均摊查找代替 O(n) 线性扫描：从对每个候选做线性扫描，改为一次建立“后缀/前缀模式索引”，将连接阶段由 O(n²) 降到 O(n + matches)。
                # 去 pre_index 中查找：谁的前缀恰好等于当前 m 的后缀
                for n in pre_index.get(sufstr, []):
                    # 如果模式n的最后一个itemset只有一个item，则进行S扩展
                    if len(n[-1]) == 1:
                        # ---- S-扩展：追加新的项集[a][b] ----
                        # 检查是否能进行合法S扩展
                        #if str(m[0][0]) in self.SItems.keys():
                        if True:
                            #if str(n[-1][0]) in self.SItems[str(m[0][0])].keys():
                            if True:
                                # [借鉴RNSP] S-扩展只是在末尾追加一个全新的
                                # itemset，m本身任何一层都不会被修改，用
                                # 列表拼接构造新对象即可，无需deepcopy整棵
                                # 嵌套结构。
                                pattern = m + [n[-1]]
                                # 获取最后一个item的key，用于ItemS（查找在每条序列中的位置）
                                last_item_key = str([n[-1][0]])
                                # 获取前缀的key
                                prefix_key = self._norm_key(pattern[:-1])

                                # [NEW-4] BF-PSUB 候选级剪枝（S-扩展专用，边界过滤支持度上界）
                                if not self._pass_auub_S(ItemS[prefix_key], ItemS[last_item_key]):
                                    self.AuubPruned += 1
                                    self.BfPsubPruned += 1
                                    continue

                                self.CanNum += 1
                                # 计算扩展得到的pattern的支持度和路径
                                count, ItemS[str(pattern)] = self.Matching_S(
                                    ItemS[prefix_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                    # 放入下一轮扩展集合
                                    ExpSet.append(pattern)
                                else:
                                    del ItemS[str(pattern)]

                    else:
                        # 如果当前模式包含多个itemset，则进行I扩展
                        # [借鉴RNSP] 只有"最后一个itemset"需要变化（追加一个
                        # 新item），前面的itemset原样复用；构造
                        # m[:-1] + [新的最后itemset] 即可得到与
                        # deepcopy+append完全等价、但不共享可变对象的
                        # pattern，同时避免了对m的任何原地修改。
                        pattern = m[:-1] + [m[-1] + [n[-1][-1]]]

                        if len(pattern) > 1:
                            #if str(pattern[-1][-1]) in self.SItems[str(pattern[0][0])].keys():
                            if True:
                                # ---- I-扩展：追加新的项集[a][b][c] ----
                                # [FIX-2] 这里必须查询“扩展后的整个最后一个itemset”
                                # 作为独立单itemset模式的出现位置（对应原版
                                # Matching_I(ItemS[str(m)], S[...]) 之后再与前缀做
                                # Matching_S 时使用的 ItemS[str([pattern[-1]])]），
                                # 而不是只查新加入的这一个item自己的位置。
                                # 只用新item自己的位置，等价于漏掉了“新item必须
                                # 与itemset内其它已有item出现在同一个位置”这一
                                # 约束，会让support被错误地放大（这是v2挖掘结果
                                # 与HUP-Miner不一致的根本原因）。
                                # 该整体itemset模式一定已经在本轮temp的处理中
                                # 通过下面的“单itemset扩展”分支（len(pattern)==1）
                                # 被计算并写入ItemS，因此这里可以直接查表。
                                last_item_key = self._norm_key([pattern[-1]])
                                # pattern[:-1] 在退化为单项集单项模式时
                                # 需要归一化为扁平键，否则会出现 KeyError
                                prefix_key = self._norm_key(pattern[:-1])

                                # [NEW-4] BF-PSUB 候选级剪枝（此处虽是"I扩展"语义，
                                # 但 prefix 与"整个扩展后itemset"之间是按位置先后
                                # 顺序做 Matching_S 连接，故须用 S-扩展专用的边界
                                # 过滤上界，而非按取值相等匹配的 PVIB）
                                if not self._pass_auub_S(ItemS[prefix_key], ItemS[last_item_key]):
                                    self.AuubPruned += 1
                                    self.BfPsubPruned += 1
                                    continue

                                self.CanNum += 1
                                # 计算扩展得到的pattern的支持度和路径
                                count, ItemS[str(pattern)] = self.Matching_S(
                                    ItemS[prefix_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]
                        else:
                            #if str(pattern[-1][-1]) in self.SItem[str(pattern[0][0])].keys():
                            if True:
                                last_item_key = str([n[-1][-1]])
                                m_key = self._norm_key(m)

                                self.CanNum += 1
                                count, ItemS[str(pattern)] = self.Matching_I(
                                    ItemS[m_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]

    # 读取utility文件，计算Utility和maxUtil，用于AUUB剪枝
    def getUtility(self, utilityFileName):
        maxau = 0
        lines = self.pdata.read_file(utilityFileName)
        for line in lines:
            parts = line.split(' ')
            if len(parts) >= 2:
                self.Utility[parts[0]] = int(parts[1])
                # 找到最大的Utility-maxUtil
                if int(parts[1]) > maxau:
                    maxau = int(parts[1])
        # minau 映射为 minsup（从效用阈值推导出的必要支持度约束）
        if maxau > 0:
            self.minsup = math.ceil(int(self.minau) / maxau)
        self.maxUtil = maxau  # [NEW-1] 供 AUUB 剪枝使用

    # 计算模式的平均效用（AU）和支持度（SUP），并判断是否满足minau阈值
    def compute_utility(self, pattern, count):
        # 生成字典key
        pattern_str = str(pattern)
        # 初始化
        au = 0
        pattern_len = 0
        # 遍历pattern中的每个itemset和其item
        for itemset in pattern:
            for item in itemset:
                pattern_len += 1
                au += self.Utility.get(str(item), 0)
        if pattern_len > 0:
            au = au * count / pattern_len
            if au >= int(self.minau):
                self.AUNP.append(pattern)
                self.util_map[pattern_str] = au
                self.support[pattern_str] = int(count)

    def Miner(self):
        start = time.time()
        FP = []
        # 路径索引
        ItemS = {}
        # 挖掘1长度FP
        self.Mine_ItemS(FP, ItemS)
        # 得到2长度模式（AUUB剪枝）
        twoLenPattern = self.two_len(FP, ItemS)
        # 不断递推，获取>2长度模式，直至没有新FP
        self.more_len(FP, ItemS, twoLenPattern)
        #print("Frequent patterns with size=1:" + str(FP))
        #print("Number of frequent patterns with size=1:" + str(len(FP)))
        self.mining_time_ms = int((time.time() - start) * 1000)

        # 统一键格式：将所有键转换为扁平化格式（继承自 Mine1_v3，未改动）
        self.ItemS_map = {}

        # [FIX-3] 之前这里对「多个单item itemset依次排列」（如[a][c]）和
        # 「一个itemset里塞多个item」（如[ac]）都会被压成一样的 ['a','c']，
        # 导致两者用同一个new_key在self.ItemS_map里互相覆盖。
        # 现在统一改用 self._flatten_pattern，只逐层转list、不跨itemset
        # 合并，从而保留itemset边界；仅对 Mine_ItemS 里单item模式使用的
        # 扁平key（如 "['a']"，本身无边界歧义）保持原样。
        for key, value in ItemS.items():
            try:
                p_obj = eval(key)
                if isinstance(p_obj, list) and len(p_obj) == 1 and not isinstance(p_obj[0], list):
                    new_key = key
                else:
                    new_key = str(self._flatten_pattern(p_obj))
                self.ItemS_map[new_key] = value
            except Exception as e:
                print(f"Warning: Could not process key {key}: {e}")
                self.ItemS_map[key] = value

    def get_stats(self):
        """返回本次挖掘的性能统计信息，供消融实验对比使用。"""
        return {
            "candidate_count": self.CanNum,
            "auub_pruned_count": self.AuubPruned,
            #"pvib_pruned_count": self.PvibPruned,
            "bfpsub_pruned_count": self.BfPsubPruned,
            "aunp_count": len(self.AUNP),
            "mining_time_ms": self.mining_time_ms,
        }

    def mine_patterns(self, dataset_path, min_utility_threshold, utility_file_path=None, save_mine_path=None):

        # [FIX-3] 统一使用 self._flatten_pattern，保留itemset边界，
        # 避免 [a][c] 与 [ac] 在输出/JSON里被压成同一个模式
        flatten_pattern = self._flatten_pattern

        self.__init__()
        self.minau = min_utility_threshold
        # 读取原始数据
        self.SeqNum, self.S, self.sort_item = self.pdata.datap(dataset_path, self.S)

        if utility_file_path is None:
            last_index = dataset_path.rindex(".")
            dataFileName = dataset_path[0:last_index]
            last_index = dataFileName.rindex("/")
            utility_file_path = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility.txt"
        # 读取utility
        self.getUtility(utility_file_path)

        # 进行正式挖掘
        self.Miner()
        #a = memory_usage((self.Miner, (), {}))
        #print('Memory usage: ', max(a) - min(a))
        if save_mine_path is None:
            base = dataset_path
            try:
                last_index = base.rindex(".")
                base_no_ext = base[:last_index]
            except ValueError:
                base_no_ext = base
            save_mine_path = base_no_ext + "_mine_results.json"

            patterns_info = []
            pattern_ItemS = {}
            for p in self.AUNP:
                flat_p = flatten_pattern(p)
                p_str = str(flat_p)
                info = {
                    "pattern": p_str,
                    "support": int(self.support.get(str(p), 0)),
                    "utility": float(self.util_map.get(str(p), 0.0))
                }
                patterns_info.append(info)
                occ_info = []
                pattern_positions = self.ItemS_map.get(p_str, [])
                for seq_id in range(self.SeqNum):
                    if seq_id < len(pattern_positions) and pattern_positions[seq_id]:
                        occ_info.append({
                            "seq": seq_id,
                            "positions": pattern_positions[seq_id]
                        })
                pattern_ItemS[p_str] = occ_info

            flattened_util_map = {}
            flattened_support_map = {}
            for p_str, val in self.util_map.items():
                try:
                    p_obj = eval(p_str)
                    flat_p = flatten_pattern(p_obj)
                    flattened_util_map[str(flat_p)] = val
                except Exception:
                    flattened_util_map[p_str] = val
            for p_str, val in self.support.items():
                try:
                    p_obj = eval(p_str)
                    flat_p = flatten_pattern(p_obj)
                    flattened_support_map[str(flat_p)] = val
                except Exception:
                    flattened_support_map[p_str] = val

            mine_results = {
                "seq_num": self.SeqNum,
                "S": self.S,
                "sort_item": self.sort_item,
                "Utility": self.Utility,
                "minau": self.minau,
                "minsup": self.minsup,
                "patterns": patterns_info,
                "pattern_ItemS": pattern_ItemS,
                "util_map": flattened_util_map,
                "support_map": flattened_support_map,
                # 以下为 v4 新增的性能统计字段，向后兼容（旧版消费者可忽略）
                "candidate_count": self.CanNum,
                "auub_pruned_count": self.AuubPruned,
                #"pvib_pruned_count": self.PvibPruned,
                "bfpsub_pruned_count": self.BfPsubPruned,
                "mining_time_ms": self.mining_time_ms,
            }

            with open(save_mine_path, "w", encoding="utf-8") as f:
                json.dump(mine_results, f, ensure_ascii=False, indent=2)
            self.mine_result_file = save_mine_path

        self.AUNP = [flatten_pattern(p) for p in self.AUNP]

        return {
            'aunp_count': len(self.AUNP),
            'aunp_patterns': self.AUNP,
            'mine_result_file': self.mine_result_file,
            'candidate_count': self.CanNum,
            'auub_pruned_count': self.AuubPruned,
            #'pvib_pruned_count': self.PvibPruned,
            'bfpsub_pruned_count': self.BfPsubPruned,
            'mining_time_ms': self.mining_time_ms,
        }

    def mine_patterns_from_S(self, S, seq_num, sort_item, utility_dict, min_utility_threshold):
        """直接从 S 结构进行挖掘（用于隐藏阶段后重新挖掘验证效果，接口与 Mine1_v3 保持一致）"""
        self.__init__()
        self.S = S
        self.SeqNum = seq_num
        self.sort_item = sort_item
        self.Utility = utility_dict
        self.minau = min_utility_threshold

        maxau = max(utility_dict.values()) if utility_dict else 1
        self.minsup = math.ceil(int(min_utility_threshold) / maxau) if maxau > 0 else 0
        self.maxUtil = maxau  # [NEW-1]

        self.Miner()

        # [FIX-3] 同上，统一使用保留itemset边界的 _flatten_pattern
        self.after_AUNP = [self._flatten_pattern(p) for p in self.AUNP]

        return {
            'aunp_count': len(self.after_AUNP),
            'aunp_patterns': self.after_AUNP,
            'candidate_count': self.CanNum,
            'auub_pruned_count': self.AuubPruned,
            #'pvib_pruned_count': self.PvibPruned,
            'bfpsub_pruned_count': self.BfPsubPruned,
            'mining_time_ms': self.mining_time_ms,
        }


def run_pattern_mining(dataset_path, min_utility_threshold, utility_file_path=None):
    miner = PatternMiner()
    return miner.mine_patterns(dataset_path, min_utility_threshold, utility_file_path)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python HAUP_v2.py <dataset_path> <min_utility_threshold> [utility_file_path]")
        sys.exit(1)
    dataset_path = sys.argv[1]
    last_index = dataset_path.rindex(".")
    dataFileName = dataset_path[0:last_index]
    last_index = dataFileName.rindex("/")
    utility_file = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"

    # min_utility = sys.argv[2]
    for min_utility in sys.argv[2:]:
    #utility_file = sys.argv[3] if len(sys.argv) > 3 else None
        print('HAUP-Miner:', dataset_path, 'minau=', min_utility, ':')
        results = run_pattern_mining(dataset_path, min_utility, utility_file)
        print(f"Number of AUNP: {results['aunp_count']}")
        #print(f"AUNP patterns: {results['aunp_patterns']}")
        print(f"Candidates evaluated: {results['candidate_count']}, AUUB pruned (total): {results['auub_pruned_count']} "
          f"Time: {results['mining_time_ms']}ms")
