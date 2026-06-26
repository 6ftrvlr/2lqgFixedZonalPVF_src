#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量生成数值模拟数据保存路径。

功能：
1. 支持从 txt 文件中读取 case 名称；
2. 支持从 csv 文件的 script_name 字段中读取 case 名称；
3. 为每个 case 创建同名文件夹；
4. 在每个 case 文件夹下创建：
   ./data_Saving
   ./data_Saving/fig
   ./data_Saving/raw
"""

import argparse
import csv
from pathlib import Path


def read_case_names_from_txt(case_file: Path):
    """
    从普通文本文件中读取 case 名称。

    支持格式：
        case_001
        case_002
        Re100_Ma02

    空行会被忽略。
    以 # 开头的行会被视为注释。
    """
    case_names = []

    with case_file.open("r", encoding="utf-8") as f:
        for line in f:
            name = line.strip()

            if not name:
                continue

            if name.startswith("#"):
                continue

            case_names.append(name)

    return case_names


def read_case_names_from_csv(csv_file: Path, field_name: str = "script_name"):
    """
    从 csv 文件中读取 case 名称。

    默认读取字段：
        script_name

    示例 csv：

        script_name,Re,Ma
        case_001,100,0.2
        case_002,200,0.5
    """
    case_names = []

    with csv_file.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise RuntimeError(f"CSV file is empty or invalid: {csv_file}")

        if field_name not in reader.fieldnames:
            raise RuntimeError(
                f"Field '{field_name}' not found in CSV file: {csv_file}\n"
                f"Available fields: {reader.fieldnames}"
            )

        for row in reader:
            name = row.get(field_name, "").strip()

            if not name:
                continue

            case_names.append(name)

    return case_names


def remove_duplicates_keep_order(items):
    """
    去除重复 case name，同时保持原始顺序。
    """
    seen = set()
    result = []

    for item in items:
        if item in seen:
            continue

        seen.add(item)
        result.append(item)

    return result


def create_case_dirs(base_dir: Path, case_name: str):
    """
    为单个 case 创建目录结构。
    """

    case_dir = base_dir / case_name

    dirs_to_create = [
        case_dir / "data_Saving",
        case_dir / "data_Saving" / "fig",
        case_dir / "data_Saving" / "raw",
    ]

    for d in dirs_to_create:
        d.mkdir(parents=True, exist_ok=True)

    return case_dir


def main():
    parser = argparse.ArgumentParser(
        description="Batch create directory structures for numerical simulation cases."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "-f",
        "--case-file",
        help="存放 case 名称的普通文本文件，例如 cases.txt",
    )

    input_group.add_argument(
        "-c",
        "--csv-file",
        help="存放 case 信息的 csv 文件，将读取 script_name 字段作为 case name",
    )

    parser.add_argument(
        "--csv-field",
        default="script_name",
        help="csv 中用于读取 case name 的字段名，默认 script_name",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="case 文件夹的输出根目录，默认当前目录",
    )

    parser.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="允许重复 case name。默认会自动去重。",
    )

    args = parser.parse_args()

    base_dir = Path(args.output_dir).resolve()
    base_dir.mkdir(parents=True, exist_ok=True)

    if args.case_file is not None:
        case_file = Path(args.case_file).resolve()

        if not case_file.exists():
            raise FileNotFoundError(f"Case file not found: {case_file}")

        case_names = read_case_names_from_txt(case_file)

    elif args.csv_file is not None:
        csv_file = Path(args.csv_file).resolve()

        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        case_names = read_case_names_from_csv(
            csv_file=csv_file,
            field_name=args.csv_field,
        )

    else:
        raise RuntimeError("Either --case-file or --csv-file must be provided.")

    if not args.allow_duplicates:
        case_names = remove_duplicates_keep_order(case_names)

    if not case_names:
        raise RuntimeError("No valid case names found.")

    print(f"Output directory: {base_dir}")
    print(f"Number of cases : {len(case_names)}")
    print("-" * 50)

    for case_name in case_names:
        case_dir = create_case_dirs(base_dir, case_name)
        print(f"Created: {case_dir}")

    print("-" * 50)
    print("All case directories have been created.")


if __name__ == "__main__":
    main()
