import platform

def syscheck():
    if platform.system() == 'Windows':
        return 'p:'
    elif platform.system() == 'Linux':
        return '/p/'