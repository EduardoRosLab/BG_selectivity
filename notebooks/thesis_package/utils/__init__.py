import numpy as np
from copy import deepcopy
from time import time as system_time
from math import floor
from .custom_exceptions import ParamError


def split_in_groups(an_iterable, proportions):
    proportions = proportions[:]
    groups = [[]]

    pos = 0.0
    step = sum(proportions) / len(an_iterable)
    for e in an_iterable:
        if proportions and pos > proportions[0]:
            proportions.pop(0)
            groups.append([])
            pos = 0
        pos += step
        groups[-1].append(e)

    for group in groups:
        if len(group) == 0: raise ParamError("Unknown variable", "A channel has zero neurons.")

    return groups


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def copy(self):
        new_dotdict = dotdict()
        for key, value in self.items():
            if type(value) is dotdict:
                new_value = value.copy()
            else:
                new_value = deepcopy(value)
            new_dotdict[key] = new_value
        return new_dotdict


from time import time
def verbose(verbose=True):

    def decorator(f):
        if not verbose: return f

        def wrapped_func(*args, **kwargs):
            print("Entering {}...".format(f.__name__))
            running_time = time()
            f(*args, **kwargs)
            running_time = time() - running_time
            print("Exited {} after {} seconds.\n".format(f.__name__, running_time))

        return wrapped_func

    return decorator


def max_dif(signal1, signal2, t=(0, 150)):
    #return np.max(np.abs(signal1[t[0]:t[1]] - signal2[t[0]:t[1]]))
    return np.percentile(np.abs(signal1[t[0]:t[1]] - signal2[t[0]:t[1]]), 95.0) #Más estable con percentil 95 que con el máximo.

def transient_selectivity(signal1, signal2, t_stable=(2000, 2500), t_transit=(1500, 1700)):
    _S1 = np.mean(signal1[t_stable[0]:t_stable[1]])
    _S2 = np.mean(signal2[t_stable[0]:t_stable[1]])
    dif = max_dif(signal1, signal2, t_transit)
    return 1.0 - ((_S1 - _S2) / dif)

def steadystate_selectivity(signal, t_pre=(0, 50), t_stable=(300, 400)):
    SP_X = np.mean(signal[t_pre[0]:t_pre[1]])
    S_X = np.mean(signal[t_stable[0]:t_stable[1]])
    return 100*(1.0 - (S_X / SP_X))


def log_progress(sequence, every=None, size=None, name='Items'):
    # https://github.com/kuk/log-progress

    from ipywidgets import IntProgress, HTML, VBox
    from IPython.display import display

    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = int(size / 200)     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    full_init_time = system_time()
    elapsed_time = system_time()
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{name}: {index} / ?'.format(
                        name=name,
                        index=index
                    )
                else:
                    time_per_case = (elapsed_time-full_init_time) / index

                    rem_time = time_per_case * (size-index+1)
                    hours = floor(rem_time / 3600.0)
                    minutes = floor((rem_time - hours*3600.0) / 60.0)
                    seconds = floor(rem_time - minutes*60.0 - hours*3600)
                    rem_time = "{hours}:{minutes:02}:{seconds:02}".format(hours=hours, minutes=minutes, seconds=seconds)

                    ela_time = elapsed_time - full_init_time
                    hours = floor(ela_time / 3600.0)
                    minutes = floor((ela_time - hours*3600.0) / 60.0)
                    seconds = floor(ela_time - minutes*60.0 - hours*3600)
                    ela_time = "{hours}:{minutes:02}:{seconds:02}".format(hours=hours, minutes=minutes, seconds=seconds)

                    progress.value = index
                    label.value = u'{name}: {index} / {size} - Ela: {ela_time}\nRem: {rem_time}'.format(
                        name=name,
                        index=index,
                        size=size,
                        ela_time=ela_time,
                        rem_time=rem_time
                    )
            yield record
            elapsed_time = system_time()
    except:
        progress.bar_style = 'danger'
        raise
    else:
        ela_time = system_time() - full_init_time
        hours = floor(ela_time / 3600.0)
        minutes = floor((ela_time - hours * 3600.0) / 60.0)
        seconds = floor(ela_time - minutes * 60.0 - hours * 3600)
        ela_time = "{hours}:{minutes:02}:{seconds:02}".format(hours=hours, minutes=minutes, seconds=seconds)

        progress.bar_style = 'success'
        progress.value = index
        label.value = "{name}: {index}. Total time: {ela_time}".format(
            name=name,
            index=str(index or '?'),
            ela_time=ela_time
        )


import numpy

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise(ValueError, "smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
                raise(ValueError, "Input vector needs to be bigger than window size.")
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=numpy.ones(window_len,'d')
        else:
                w=eval('numpy.'+window+'(window_len)')
        y=numpy.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]
