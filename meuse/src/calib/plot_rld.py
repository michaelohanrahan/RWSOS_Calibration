import matplotlib.pyplot as plt


def main(
    data,
    maxima,
    minima,
):
    """_summary_"""
    fig = plt.figure()
    fig.set_size_inches(12,5)
    ax = fig.add_subplot(111)
    ax.set_position([0.05, 0.07, 0.91, 0.85])
    
    ax.plot(
        range(0,len(data)),
        data,
        **{"linewidth":1.0, "color": "#000000"}
    )

    ax.plot(
        maxima[0],
        maxima[1],
        **{"marker":"+", "markeredgecolor":"#FF0000", "linestyle":'None'}
    )

    ax.plot(
        minima[0],
        minima[1],
        **{"marker":"o", "markerfacecolor":"None", "markeredgecolor":"#0000FF", "linestyle":'None'}
    )

    ax.set_xlim([0, len(data)])
    ax.set_ylim([0,max(data)*1.15])
    ax.tick_params(axis="both", which="both", **{"labelsize": 7})

    plt.show()

    return fig


if __name__ == "__main__":
    pass