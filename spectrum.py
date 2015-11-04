#!/usr/bin/env python
"""
This plot displays the audio waveform, spectrum, and spectrogram from the 
microphone.

Based on updating_plot.py
"""
import sys

# Major library imports
try:
    import pyaudio
except ImportError:
    sys.exit('You need pyaudio installed to run this demo.')

from numpy import linspace, short, fromstring, hstack, transpose, reshape, zeros
import numpy as np
from scipy import fft
import base64
# Enthought library imports
from chaco.default_colormaps import jet
from enable.api import Component, ComponentEditor
from traits.api import HasTraits, Instance, Button, Int, Float
from traitsui.api import Item, Group, View, Handler, HGroup
from pyface.timer.api import Timer

# Chaco imports
from chaco.api import Plot, ArrayPlotData, HPlotContainer

"""
NUM_SAMPLES = 1024
SAMPLING_RATE = 11025
SPECTROGRAM_LENGTH = 100
"""
C = 299792458
B = 330000000

FREQUENCY_DIVIDER = 12
CUT = 100
DELAY_TIME = 100
NUM_SAMPLES = 8200
SAMPLING_RATE = 40000
SPECTROGRAM_LENGTH = 50

# ============================================================================
# Create the Chaco plot.
# ============================================================================


def _create_plot_component(obj):
    # Setup the spectrum plot
    # frequencies = linspace(0.0, float(SAMPLING_RATE)/2, num=NUM_SAMPLES/2)
    num_samples = obj.num_samples - CUT
    frequencies = linspace(0.0, float(SAMPLING_RATE)/(2*FREQUENCY_DIVIDER), num=num_samples/(2*FREQUENCY_DIVIDER))
    obj.spectrum_data = ArrayPlotData(frequency=frequencies)
    # empty_amplitude = zeros(NUM_SAMPLES/2)
    empty_amplitude = zeros(num_samples/(2*FREQUENCY_DIVIDER))
    obj.spectrum_data.set_data('amplitude', empty_amplitude)

    obj.spectrum_plot = Plot(obj.spectrum_data)
    obj.spectrum_plot.plot(("frequency", "amplitude"), name="Spectrum",
                           color="red")
    obj.spectrum_plot.padding = 50
    obj.spectrum_plot.title = "Spectrum"
    spec_range = obj.spectrum_plot.plots.values()[0][0].value_mapper.range
    spec_range.low = 0.0
    spec_range.high = 10.0
    obj.spectrum_plot.index_axis.title = 'Frequency (Hz)'
    obj.spectrum_plot.value_axis.title = 'Amplitude'

    # Time Series plot
    # times = linspace(0.0, float(NUM_SAMPLES)/SAMPLING_RATE, num=NUM_SAMPLES)
    times = linspace(0.0, float(num_samples)/SAMPLING_RATE, num=num_samples)
    obj.time_data = ArrayPlotData(time=times)
    # empty_amplitude = zeros(NUM_SAMPLES)
    empty_amplitude = zeros(num_samples)
    obj.time_data.set_data('amplitude', empty_amplitude)

    obj.time_plot = Plot(obj.time_data)
    obj.time_plot.plot(("time", "amplitude"), name="Time", color="blue")
    obj.time_plot.padding = 50
    obj.time_plot.title = "Time - White jack"
    obj.time_plot.index_axis.title = 'Time (seconds)'
    obj.time_plot.value_axis.title = 'Amplitude'
    time_range = obj.time_plot.plots.values()[0][0].value_mapper.range
    time_range.low = -0.2
    time_range.high = 0.2

    # Time Series plot1
    # times = linspace(0.0, float(NUM_SAMPLES)/SAMPLING_RATE, num=NUM_SAMPLES)
    obj.time_data1 = ArrayPlotData(time=times)
    # empty_amplitude = zeros(NUM_SAMPLES)
    obj.time_data1.set_data('amplitude', empty_amplitude)

    obj.time_plot1 = Plot(obj.time_data1)
    obj.time_plot1.plot(("time", "amplitude"), name="Time", color="blue")
    obj.time_plot1.padding = 50
    obj.time_plot1.title = "Time - Red jack"
    obj.time_plot1.index_axis.title = 'Time (seconds)'
    obj.time_plot1.value_axis.title = 'Amplitude'
    time_range = obj.time_plot1.plots.values()[0][0].value_mapper.range
    time_range.low = -0.2
    time_range.high = 0.2

    # Spectrogram plot
    # spectrogram_data = zeros((NUM_SAMPLES/2, SPECTROGRAM_LENGTH))
    spectrogram_data = zeros((num_samples/(2*FREQUENCY_DIVIDER), SPECTROGRAM_LENGTH))
    obj.spectrogram_plotdata = ArrayPlotData()
    obj.spectrogram_plotdata.set_data('imagedata', spectrogram_data)
    spectrogram_plot = Plot(obj.spectrogram_plotdata)

    max_time = float(SPECTROGRAM_LENGTH * num_samples) / SAMPLING_RATE
    # max_freq = float(SAMPLING_RATE / 2)
    # boundary = C * obj.num_samples / (4*B * FREQUENCY_DIVIDER)
    boundary = C * num_samples / (4*B * FREQUENCY_DIVIDER)
    spectrogram_plot.img_plot('imagedata',
                              name='Spectrogram',
                              xbounds=(0, max_time),
                              # ybounds=(0, max_freq),
                              # ybounds=(0, max_freq * C * num_samples / (2*B*SAMPLING_RATE)),
                              ybounds=(0, boundary),
                              colormap=jet,
                              )
    range_obj = spectrogram_plot.plots['Spectrogram'][0].value_mapper.range
    range_obj.high = 5
    range_obj.low = 0.0
    spectrogram_plot.title = 'Spectrogram'
    obj.spectrogram_plot = spectrogram_plot

    container = HPlotContainer()
    container.add(obj.spectrum_plot)
    container.add(obj.time_plot)
    container.add(obj.time_plot1)
    container.add(obj.spectrogram_plot)

    return container

_stream = None


def get_normalized_audio():
    global _stream
    if _stream is None:
        pa = pyaudio.PyAudio()
        _stream = pa.open(format=pyaudio.paInt16, channels=2, rate=SAMPLING_RATE,
                          input=True, frames_per_buffer=NUM_SAMPLES)
    audio_data = fromstring(_stream.read(NUM_SAMPLES), dtype=short)
    audio_data = reshape(audio_data, (NUM_SAMPLES, 2))

    return audio_data / 32768.0


def get_audio_data(num_samples):
    audio_data = get_normalized_audio()

    flanks = get_stream_flanks(audio_data[:, 1])

    # I delete the first and the last period
    rising_flanks = map(lambda x: x + CUT/2, flanks[2:-2:2])
    """
    flanks = map(lambda x: x + CUT/2, flanks[2:-2:2])
    if flanks:
        normalized_data = np.mean([audio_data[i:i+num_samples] for i in flanks], axis=0)
    else:
        normalized_data = audio_data[0:num_samples]

    return abs(fft(normalized_data[:, 0]/2))[:num_samples/(2*FREQUENCY_DIVIDER)], normalized_data
    """

    f = lambda x, y: np.pad(x, (0, y), mode='constant')

    if flanks:
        num_samples2 = int(round(np.mean(map(lambda x, y: x-y, flanks[1::2], flanks[0::2])))) - CUT

        length = num_samples if num_samples2 > num_samples else num_samples2
        dif = num_samples - length
        normalized_data = f(np.mean([audio_data[i:i+length] for i in rising_flanks], axis=0), dif)
    else:
        normalized_data = audio_data[0:num_samples]
        length = num_samples
        dif = 0

    return f(abs(fft(normalized_data[:length, 0]/2))[:length/(2*FREQUENCY_DIVIDER)],
             dif/(2*FREQUENCY_DIVIDER)), normalized_data, length/(2*FREQUENCY_DIVIDER)


# HasTraits class that supplies the callable for the timer event.
class TimerController(HasTraits):
    value = Float(0)

    view = View(Group(Item('value', label='Value', style='readonly')))

    def __init__(self):
        super(TimerController, self).__init__()

        self.__measure_clutter = True
        self.__clutter = None

    def onTimer(self, *args):
        spectrum, time, length = get_audio_data(self.num_samples - CUT)

        if self.__measure_clutter:
            self.__measure_clutter = False
            self.__clutter = spectrum

        # last_spectrum = spectrum - self.__clutter

        final_spectrum = spectrum - np.pad(self.__clutter[:length], (0, len(spectrum) - length), mode='constant')

        self.value = final_spectrum.argmax() * float(SAMPLING_RATE)/(2*FREQUENCY_DIVIDER)/((self.num_samples - CUT)/(2*FREQUENCY_DIVIDER)-1)
        # print(len(final_spectrum), final_spectrum.argmax(), float(SAMPLING_RATE)/(2*FREQUENCY_DIVIDER) / ((self.num_samples - CUT)/(2*FREQUENCY_DIVIDER)-1), self.value)

        self.spectrum_data.set_data('amplitude', final_spectrum)
        self.time_data.set_data('amplitude', time[:, 0])
        self.time_data1.set_data('amplitude', time[:, 1])
        spectrogram_data = self.spectrogram_plotdata.get_data('imagedata')
        spectrogram_data = hstack((spectrogram_data[:, 1:],
                                   transpose([final_spectrum])))

        self.spectrogram_plotdata.set_data('imagedata', spectrogram_data)
        self.spectrum_plot.request_redraw()
        return

    def measure_clutter(self):
        self.__measure_clutter = True

    def reset_clutter(self):
        self.__clutter = np.zeros(len(self.__clutter))


# ============================================================================
# Attributes to use for the plot view.
size = (900, 500)
title = "Audio Spectrum"

# ============================================================================
# Demo class that is used by the demo.py application.
# ============================================================================


class DemoHandler(Handler):

    def closed(self, info, is_ok):
        """ Handles a dialog-based user interface being closed by the user.
        Overridden here to stop the timer once the window is destroyed.
        """
        info.object.timer.Stop()
        return


class Demo(HasTraits):

    plot = Instance(Component)

    controller = Instance(TimerController, ())

    timer = Instance(Timer)

    remove_clutter = Button("Remove Clutter")
    reset_clutter = Button("Reset Clutter")

    # value = Int(4)
    value = Int(0)

    traits_view = View(
        Group(
            HGroup(
                Item('remove_clutter', width=.10, height=32, show_label=False),
                Item('reset_clutter', width=.10, height=32, show_label=False)),
            Item("controller", style='custom', show_label=False),
            Item('plot', editor=ComponentEditor(size=size),
                 show_label=False),
            orientation="vertical"),
        resizable=True, title=title,
        width=size[0], height=size[1],
        handler=DemoHandler
    )

    def __init__(self, num_samples, **traits):
        super(Demo, self).__init__(**traits)
        self.controller.num_samples = num_samples
        self.plot = _create_plot_component(self.controller)

    def _remove_clutter_fired(self):
        self.controller.measure_clutter()

    def _reset_clutter_fired(self):
        self.controller.reset_clutter()

    def edit_traits(self, *args, **kws):
        # Start up the timer! We should do this only when the demo actually
        # starts and not when the demo object is created.
        self.timer = Timer(20, self.controller.onTimer)
        return super(Demo, self).edit_traits(*args, **kws)

    def configure_traits(self, *args, **kws):
        # Start up the timer! We should do this only when the demo actually
        # starts and not when the demo object is created.
        self.timer = Timer(20, self.controller.onTimer)
        return super(Demo, self).configure_traits(*args, **kws)


def get_stream_flanks(stream):
    window = (max(stream) - min(stream)) / 4
    window = window if window > 0.5 else 0.5
    flanks = [i for i, value in enumerate(stream) if abs(stream[i-1] - value) > window]
    return sum([[flanks[i-1], val] for i, val in enumerate(flanks) if val - flanks[i - 1] > DELAY_TIME], [])


def get_stream_num_samples_per_period(stream):
    num_samples = 0
    flanks = get_stream_flanks(stream)

    if len(flanks) > 5:
        print(flanks[3] - flanks[2])
        # num_samples = int(round(np.mean(map(lambda x, y: x-y, flanks[1::2], flanks[0::2]))))
        # I remove the first and the last block.
        num_samples = int(round(np.mean(map(lambda x, y: x-y, flanks[3:-2:2], flanks[2:-2:2]))))
    return num_samples


if __name__ == "__main__":
    quantity_samples = 0
    while quantity_samples == 0:
        quantity_samples = get_stream_num_samples_per_period(get_normalized_audio()[:, 1])

    print("working", quantity_samples)

    popup = Demo(quantity_samples)
    try:
        popup.configure_traits()
    finally:
        if _stream is not None:
            _stream.close()
        print("closing")
