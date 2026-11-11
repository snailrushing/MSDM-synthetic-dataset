
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = BASE
OUT_DIR = os.path.join(BASE, "data_curves")
FONT_PATH = os.path.join(os.path.dirname(BASE), "times_simsun.ttf")


FILES = [
    "tanker_seastate_6_angle_0_86.csv",
    "tanker_seastate_6_angle_45_86.csv",
    "tanker_seastate_6_angle_90_86.csv",
    "tanker_seastate_6_angle_135_86.csv",
    "tanker_seastate_6_angle_180_86.csv",
    "tanker_seastate_5_angle_0_86.csv",
    "tanker_seastate_5_angle_45_86.csv",
    "tanker_seastate_5_angle_90_86.csv",
    "tanker_seastate_5_angle_135_86.csv",
    "tanker_seastate_5_angle_180_86.csv",
    "tanker_seastate_4_angle_0_86.csv",
    "tanker_seastate_4_angle_45_86.csv",
    "tanker_seastate_4_angle_90_86.csv",
    "tanker_seastate_4_angle_135_86.csv",
    "tanker_seastate_4_angle_180_86.csv",
]


font_manager.fontManager.addfont(FONT_PATH)
FONT_FAMILY = font_manager.FontProperties(fname=FONT_PATH).get_name()
FONT_SIZE = 16  

plt.rcParams.update(
    {
        "font.family": FONT_FAMILY,
        "font.size": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "axes.linewidth": 1.0,
    }
)


PLOTS = [
    ("heave_m", "Heave (m)", "#41B5C2"),
    ("heave_vel_mps", "Heave velocity (m·s⁻¹)", "#722169"),
    ("roll_deg", "Roll (°)", "#EB6D98"),
    ("pitch_deg", "Pitch (°)", "#5E8034"),
]


def plot_one(csv_path, png_path):
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    idx = np.arange(len(data))
    cm = 1 / 2.54

    fig_w = 15.0  
    sub_h = 2.0  
    n_sub = len(PLOTS)
    gap = 0.4 
    left_m = 1.5  
    right_m = 0.3  
    bottom_m = 1.6  
    top_m = 1.2  

    fig_h = bottom_m + n_sub * sub_h + (n_sub - 1) * gap + top_m
    fig, axes = plt.subplots(n_sub, 1, sharex=True, figsize=(fig_w * cm, fig_h * cm))

    left = left_m / fig_w
    width = (fig_w - left_m - right_m) / fig_w
    ax_h = sub_h / fig_h
    ax_gap = gap / fig_h
    top = 1.0 - top_m / fig_h

    for ax, (col, label, color) in zip(axes, PLOTS):
        ax_bottom = top - ax_h
        ax.set_position([left, ax_bottom, width, ax_h])
        top = ax_bottom - ax_gap
        ax.plot(idx, data[col], color=color, linewidth=1.0, label=label)

        ax.tick_params(axis="x", bottom=False, labelbottom=False)

    axes[-1].tick_params(axis="x", bottom=True, labelbottom=True)
    axes[-1].set_xticks(np.arange(0, len(data) + 1, 1500))  
    axes[-1].set_xlabel("Index")


    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    axes[0].legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=n_sub,
        frameon=False,
        columnspacing=1.0,
        handlelength=0.7,
        handletextpad=0.5,
        fontsize=14,  
    )

    fig.savefig(png_path, dpi=600)
    plt.close(fig)


def main():
    src_dir = os.path.join(DATA_DIR, "test")
    dst_dir = os.path.join(OUT_DIR, "test_english")
    os.makedirs(dst_dir, exist_ok=True)
    for name in FILES:
        plot_one(os.path.join(src_dir, name), os.path.join(dst_dir, name[:-4] + ".png"))
        print(f"{name[:-4]}.png", flush=True)


if __name__ == "__main__":
    main()
