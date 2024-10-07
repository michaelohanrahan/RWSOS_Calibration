"""Plotting utility."""
import matplotlib.pyplot as plt


def set_ax_props(
    ax: plt.axes,
    ax_props: dict,
    title: str | None =  None,
    title_props: dict = {},
):
    """_summary_."""
    ax.tick_params(axis="both", which="both", **ax_props)
    
    if title is not None:
        ax.set_title(title, **title_props)
    return ax