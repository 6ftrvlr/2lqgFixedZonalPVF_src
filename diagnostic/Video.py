#!/usr/bin/env python3
# make_video.py

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def find_images(input_dir: Path):
    """
    查找 t_[nums].png 格式的图片，并按 nums 数值排序。
    例如：
        t_1.png
        t_2.png
        t_10.png
    """
    pattern = re.compile(r"^t_(\d+)\.png$", re.IGNORECASE)

    items = []
    for p in input_dir.iterdir():
        if not p.is_file():
            continue

        m = pattern.match(p.name)
        if m:
            num = int(m.group(1))
            items.append((num, p))

    items.sort(key=lambda x: x[0])
    return [p for _, p in items]


def run_ffmpeg(
    images,
    output: Path,
    fps: int,
    codec: str,
    crf: int,
    preset: str,
):
    """
    为了兼容 t_1.png、t_2.png、t_10.png 这种不一定连续的文件名，
    这里创建临时目录，把图片链接为 frame_000001.png 这种连续序列。
    """

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)

        for i, img in enumerate(images, start=1):
            target = tmp_dir / f"frame_{i:06d}.png"

            try:
                os.link(img, target)  # 优先硬链接，速度快，不占额外空间
            except OSError:
                shutil.copy2(img, target)  # 跨磁盘或系统不支持硬链接时退化为复制

        input_pattern = str(tmp_dir / "frame_%06d.png")

        if codec == "h264":
            video_codec = "libx264"
            extra_codec_args = [
                "-pix_fmt", "yuv420p",
            ]

        elif codec == "h265":
            video_codec = "libx265"
            extra_codec_args = [
                "-pix_fmt", "yuv420p",
                "-tag:v", "hvc1",  # 提高 macOS / QuickTime 兼容性
            ]

        else:
            raise ValueError(f"Unsupported codec: {codec}")

        cmd = [
            "ffmpeg",
            "-y",
            "-framerate", str(fps),
            "-i", input_pattern,

            # 防止某些图片宽高为奇数时编码失败
            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",

            "-c:v", video_codec,
            "-preset", preset,
            "-crf", str(crf),

            *extra_codec_args,

            # 让 MP4 更适合在线播放/预览
            "-movflags", "+faststart",

            str(output),
        ]

        print("Running command:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Create compressed video from t_[nums].png images."
    )

    parser.add_argument(
        "-i", "--input-dir",
        default=".",
        help="输入图片目录，默认当前目录",
    )

    parser.add_argument(
        "-o", "--output",
        default="output.mp4",
        help="输出视频文件名，默认 output.mp4",
    )

    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="视频帧率，默认 30",
    )

    parser.add_argument(
        "--codec",
        choices=["h264", "h265"],
        default="h264",
        help="编码格式：h264 兼容性最好；h265 通常体积更小",
    )

    parser.add_argument(
        "--crf",
        type=int,
        default=None,
        help=(
            "压缩质量参数。数值越大，体积越小，画质越低。"
            "h264 推荐 23-30；h265 推荐 26-34。"
        ),
    )

    parser.add_argument(
        "--preset",
        default="slow",
        choices=[
            "ultrafast", "superfast", "veryfast", "faster",
            "fast", "medium", "slow", "slower", "veryslow"
        ],
        help="编码速度/压缩率权衡。越慢通常体积越小，默认 slow",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output = Path(args.output).resolve()

    images = find_images(input_dir)

    if not images:
        raise RuntimeError(f"No images like t_[nums].png found in: {input_dir}")

    if args.crf is None:
        crf = 26 if args.codec == "h264" else 30
    else:
        crf = args.crf

    print(f"Found {len(images)} images.")
    print(f"First image: {images[0].name}")
    print(f"Last image : {images[-1].name}")
    print(f"Output     : {output}")
    print(f"FPS        : {args.fps}")
    print(f"Codec      : {args.codec}")
    print(f"CRF        : {crf}")
    print(f"Preset     : {args.preset}")

    run_ffmpeg(
        images=images,
        output=output,
        fps=args.fps,
        codec=args.codec,
        crf=crf,
        preset=args.preset,
    )

    print("Done.")


if __name__ == "__main__":
    main()
