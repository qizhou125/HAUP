class processingData:
    """最小桩实现：仅保证 HAUP_v2.py / HAUP_bitdict.py 能够 `import Pdata` 成功。
    单元测试中我们通过 mine_patterns_from_S() 直接注入 self.S，不会真正调用
    这里的方法；如果调用到，说明测试用法有误，直接抛错更安全。
    """

    def datap(self, dataset_path, S):
        raise NotImplementedError("stub Pdata.datap() should not be called in unit tests")

    def read_file(self, filename):
        raise NotImplementedError("stub Pdata.read_file() should not be called in unit tests")
