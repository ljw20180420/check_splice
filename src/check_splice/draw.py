import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from numpy.typing import ArrayLike


def _white2red() -> LinearSegmentedColormap:
    return LinearSegmentedColormap(
        name="white2red",
        segmentdata={
            "red": [(0.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
            "green": [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)],
            "blue": [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)],
        },
    )


def _heatmap(mat: ArrayLike, extent: list[int]) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.matshow(mat, vmin=0, cmap=_white2red(), extent=extent)
    fig.colorbar(im, fraction=0.046, pad=0.04)
    ax.xaxis.set_label_position("top")
    fig.tight_layout()

    return fig, ax
