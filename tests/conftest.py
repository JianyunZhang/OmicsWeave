"""
pytest 全局配置。

``--project-dir <目录>`` 指定实验项目目录（例如 ``local/projects/<课题名>``）。
需要真实实验数据的测试只在该目录中找到对应数据时运行，否则自动跳过。
"""
import matplotlib

from omicsweave.core.project import set_project_dir

# 测试在无界面环境中运行：统一用不开窗口的 Agg 后端。否则后端取决于先导入了哪个绘图模块
# （MFA、整合的绘图模块在导入时设为 Agg，代谢组的没有），旧流程的 pyplot 绘图会改用 Tk，
# 在 Windows 上偶尔因读不到 Tcl 文件而失败。
matplotlib.use("Agg")


def pytest_addoption(parser):
    parser.addoption(
        "--project-dir",
        default=None,
        help="实验项目目录；需要真实实验数据的测试从这里读取数据",
    )


def pytest_configure(config):
    # 必须在收集测试模块之前设置：各测试的 skipUnless 条件在导入时求值
    project_dir = config.getoption("--project-dir")
    if project_dir:
        set_project_dir(project_dir)
