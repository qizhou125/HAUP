import re
class processingData(object):
    def read_file(self, readfilename):
        lines_s = []
        with open(readfilename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                lines_s.append(line.strip())
        return lines_s
    def write_file(self, lines_s, filename):
        with open(filename, 'w') as f:
            for i in range(len(lines_s)):
                f.writelines(str(i) + "\t" + lines_s[i])
                f.write("\n")
    def item_sotrd(self, items):
        items_no_repeat = []
        for item in items:
            if item not in items_no_repeat:
                items_no_repeat.append(item)
        sort_items = list(sorted(items_no_repeat))
        return sort_items, items_no_repeat
    def item_to_dict(self, sort_items, len_lines, S):
        for i in sort_items:
            S[i] = [[] for i in range(len_lines)]
        return S
    def replace_seq(self, lines_s, sort_item):
        flag = 1
        for i in range(len(sort_item)):
            for j in range(len(lines_s)):
                if flag == 1:
                    lines_s[j] = lines_s[j] + " "
                p1 = re.compile(" " + sort_item[i] + " ")
                lines_s[j], number = re.subn(p1, " " + str(i) + "* ", lines_s[j])
            flag = 0
        for i in range(len(lines_s)):
            lines_s[i] = lines_s[i].replace('*', '')
        return lines_s
    def split_array(self, lines):
        lines = lines.replace('  ', ' ')
        items_array = lines.strip().split(' ')
        s_array = [[]]
        i = 0
        for item in items_array:
            if item != '-1':
                s_array[i].append(item)
            else:
                i = i + 1
                s_array.append([])
        return items_array, s_array
    def Statistics_items(self, lines, items):
        lines = lines.replace('  ', ' ')
        items_array = lines.strip().split(' ')
        for item in items_array:
            if item != '-1' and item not in items:
                items.append(item)
        return items
    def General_Sn(self, lines_s, S):
        for i in range(len(lines_s)):
            items_array, s_array = self.split_array(lines_s[i])
            sort_item, items_no_repeat = self.item_sotrd(items_array)
            for item in sort_item:
                if item != '-1':
                    for k in range(len(s_array)):
                        if item in s_array[k]:
                            S[item][i].append(k)
        return S
    def datap(self, readFileName, S):
        lines_s = self.read_file(readFileName)
        items = []
        for lines in lines_s:
            items = self.Statistics_items(lines, items)
        sort_item, items_no_repeat = self.item_sotrd(items)
        self.write_file(sort_item, 'sort_item.txt')
        self.item_to_dict(sort_item, len(lines_s), S)
        self.General_Sn(lines_s, S)
        return len(lines_s), S, sort_item
    def utilityp(self, readFileName, U):
        lines_s = self.read_file(readFileName)
        for lines in lines_s:
            key_value = lines.strip().split(' ')
            U[key_value[0]] = key_value[1]
    def __del__(self):
        pass
