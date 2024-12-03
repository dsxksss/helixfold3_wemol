import argparse
import sys
import toml
from typing import Dict


def get_cli_argument(config_path: str) -> Dict:
    """从TOML配置文件读取并解析命令行参数

    Args:
        config_path (str): TOML配置文件路径

    Returns:
        Dict: 解析后的参数字典

    配置文件格式示例:
    ```toml
    description = "工具描述"
    [arguments]
    参数名称 = {
        required = true,      # 是否必需
        help = "参数说明",    # 帮助信息
        choices = [],        # 可选值列表
        type = "float",     # 参数类型
        default = 0.01      # 默认值
    }
    ```
    """

    # 读取TOML配置文件
    config = toml.load(config_path)

    # 创建参数解析器，并使用TOML配置来设置解析器的描述
    parser = argparse.ArgumentParser(description=config.get("description", "CLI Tool"))

    # 从TOML配置中添加参数
    for arg, details in config.get("arguments", {}).items():
        parser.add_argument(
            f"--{arg}",
            required=details.get("required", False),
            choices=details.get("choices"),
            default=details.get("default"),
            help=details.get("help"),
            type=eval(details.get("type", "str")),  # 默认类型为str
        )

    # 解析CLI参数
    parsed_args = parser.parse_args(sys.argv[1:])

    # 将解析后的参数转换为字典
    args_dict = vars(parsed_args)

    return args_dict


# 调用示例
if __name__ == "__main__":
    config_file = "cli_config.toml"  # 替换为实际的配置文件路径
    arguments = get_cli_argument(config_file)
    print(arguments)
