import matplotlib.pyplot as plt
import xarray as xr


def main():
    fig = plt.figure()
    fig.set_size_inches(8,5)
    ax = fig.add_subplot(111)
    ax.set_position([0.05, 0.07, 0.92, 0.85])

    ds = xr.open_dataarray(r"c:\temp\test.nc")
    ds = ds.load()
    ds.isel(lon=1,lat=1).plot(ax=ax)
    pass


if __name__ == "__main__":
    main()