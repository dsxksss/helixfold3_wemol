import logging
import os
import sys

from py_script_template.cli import get_cli_argument
from .utils import set_logging_default_config, set_progress_value


def main() -> int:
    try:
        set_logging_default_config()
        logging.info("My Python Script Template!")
        # 获取当前脚本的文件路径
        script_path = __file__

        # 获取当前脚本所在的目录
        script_dir = os.path.dirname(script_path)

        set_progress_value(5)
        set_logging_default_config()

        # 获取cli参数配置文件
        config_file = os.path.join(script_dir, "..", "..", "cli_config.toml")
        arguments = get_cli_argument(config_file)
        logging.debug(f"Input Arguments: {arguments}")
        return 0
    except Exception as e:
        sys.stderr.write(f"{e}")
        return 100
