"""
设计感知的统计（结果包的后处理，``omicsweave run`` 结束时与 ``omicsweave stats`` 运行）：

- :mod:`.models`：``design.models`` 中的公式模型（合并主效应、交互项，G3）；
- :mod:`.power`：每个（组学, 组织, 比较）的可检出效应（G3）；
- :mod:`.robustness`：逐只去掉动物的稳健性与 Type S / Type M 估计（E5）；
- :mod:`.derived`：声明式派生特征（G8）；
- :mod:`.post`：按顺序运行以上步骤，并给差异结果加 ``estimand_type``、``n_total`` 与可选的 ``lfsr``。
"""
