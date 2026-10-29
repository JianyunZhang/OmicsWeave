"""
pytest 全局配置。

``--project-dir <目录>`` 指定实验项目目录（例如 ``local/projects/<课题名>``）。
需要真实实验数据的测试只在该目录中找到对应数据时运行，否则自动跳过。
"""
from omicsweave.core.project import set_project_dir


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
