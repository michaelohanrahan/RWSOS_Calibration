"""Small time module of puget."""
import datetime


def shift_time(
    time: str,
    delta: int | float,
    format: str = "%Y-%m-%dT%H:%M:%S",
    units: str = "hours", 
):
    """_summary_."""
    time_obj = datetime.datetime.strptime(time, format)
    time_obj = time_obj - datetime.timedelta(
        **{units: delta}
    )
    # return back to a string
    return time_obj.strftime(format=format)


def check_360_day(
    time: str,
    calendar: str,
    format: str = "%Y-%m-%dT%H:%M:%S",   
):
    """_summary_."""
    if calendar != "360_day":
        return time
    
    # continue if 360 calendar 
    time_obj = datetime.datetime.strptime(time, format)
    if time_obj.day != 31:
        return time
    
    out_time = shift_time(
        time,
        delta=1,
        format=format,
        units="days",
    )

    return out_time

