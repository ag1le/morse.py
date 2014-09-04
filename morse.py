#!/usr/bin/env python
# morse.py--  morse decoder 
#
# Copyright (C) 2014   Mauri Niininen, AG1LE
#
#
# bmorse.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# morse.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bmorse.py.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import time
import string
import numpy as np
from numpy.lib import stride_tricks
import pyaudio
import math
import cmath
from scipy.io import wavfile
from scipy.signal import butter, filtfilt, periodogram
import matplotlib.pyplot as plt
from optparse import OptionParser
from array import *
from collections import deque



# Global command line options 
verbosity = None
plotter = None
agc = None
fft_scan = None

MORSE_FREQUENCY = 600.0
DIT_MAGIC = 1200  	# Dit length is 1200/WPM msec 
DEFAULT_WPM = 35   	#  WPM = 1.2*samplerate / (twodits/2)
TRACKING_FILTER_SIZE = 10
UPPER_WPM  = 60		# maximum speed
LOWER_WPM  = 5 		# minimum speed 
UPPER_THRESHOLD = 0.5
LOWER_THRESHOLD = 0.5
	
Codebook = {
  '.-'	:'A', '-...':'B', '-.-.':'C', '-..'	:'D', '.'	:'E',
  '..-.':'F', '--.'	:'G', '....':'H', '..'	:'I', '.---':'J',
  '-.-':'K', '.-..' : 'L', '--' :'M', '-.' :'N', '---':'O',
  '.--.' : 'P', '--.-' : 'Q', '.-.':'R', '...':'S', '-'  :'T',
  '..-':'U', '...-' : 'V', '.--':'W', '-..-' : 'X', '-.--' : 'Y',
  '--..' : 'Z', '.----' : '1', '..---' : '2', '...--' : '3', 
  '....-' : '4', '.....' : '5', '-....' : '6', '--...' : '7', 
  '---..' : '8','----.' : '9','-----' : '0',
  '-...-' : '=', '.-.-':'~', '.-...' :'<AS>', '.-.-.' : '<AR>', '...-.-' : '<SK>',
  '-.--.' : '<KN>', '..-.-' : '<INT>', '....--' : '<HM>', '...-.' : '<VE>',
  '.-..-.' : '\\', '.----.' : '\'', '...-..-' : '$', '-.--.' : '(', '-.--.-' : ')', 
  '--..--' : ',', '-....-' : '-', '.-.-.-' : '.', '-..-.' : '/', '---...' : ':', 
  '-.-.-.' : ';', '..--..' : '?', '..--.-' : '_', '.--.-.' : '@', '-.-.--' : '!'
}

""" short time fourier transform of audio signal """
def stft(sig, frameSize, overlapFac=0.5, window=np.hanning):
    win = window(frameSize)
    hopSize = int(frameSize - np.floor(overlapFac * frameSize))
    
    # zeros at beginning (thus center of 1st window should be for sample nr. 0)
    samples = np.append(np.zeros(np.floor(frameSize/2.0)), sig)    
    # cols for windowing
    cols = np.ceil( (len(samples) - frameSize) / float(hopSize)) + 1
    # zeros at end (thus samples can be fully covered by frames)
    samples = np.append(samples, np.zeros(frameSize))
    
    frames = stride_tricks.as_strided(samples, shape=(cols, frameSize), strides=(samples.strides[0]*hopSize, samples.strides[0])).copy()
    frames *= win
    
    return np.fft.rfft(frames)    


# returns a simple rolling average of n most recent values
# Adapted from: http://www.raspberrypi.org/forums/viewtopic.php?f=32&t=69797
class rolling_avg :
   
    def __init__(self, n=10,debug=False):
        "determine lengh of roll at instantiation"
        self.n = n
        self.xqueue = deque('')
        self.debug = debug
        
    def rolling_avg(self,x):
        # if the queue is empty then fill it with values of x
        if(self.xqueue == deque([])):
            for i in range(self.n):
                self.xqueue.append(x)
        self.xqueue.append(x)
        self.xqueue.popleft()
        avg = 0
        for i in self.xqueue:
            avg += i
        avg = avg/float(self.n)
        if self.debug:
            print("Rolling Avg:")
            for i in self.xqueue:
                print(i)
            print("avg: %f" % avg)
        return avg
 
# returns AGC decay values 
def decayavg(average,input, weight):

	if (weight <= 1.0): 
		return input
	else:
		return input * (1.0 / weight) + average * (1.0 - (1.0 / weight))

# clamps output between mn and mx values 
def clamp(x, mn, mx): 
	if x > mx:
		return mx
	if x < mn:
		return mn
	else:
		return x
	

  
class Morse:
	
	# initialize Morse object
	def __init__(self,sig,samplerate):
		self.last = 0
		self.lastmark = 0
		self.mark = 0
		self.space = 0
		self.twodits = 2*DIT_MAGIC*(samplerate/1000)/DEFAULT_WPM # assume 8 KHz sample rate
		self.ticks = 0
		self.sigma = 0.35  	# used in PNN: 0.3 .. 0.4  produces good results
		self.cws = ""		# cw string to collect . and - based on symbols received
		self.ra = rolling_avg(TRACKING_FILTER_SIZE,False)
		self.ra.rolling_avg(self.twodits )
		self.dit_low_limit = 2 * DIT_MAGIC / UPPER_WPM   #  40 msec in # of samples
		self.dit_high_limit = 2 * DIT_MAGIC / LOWER_WPM   # 240 msec in # of samples
		
	def addchar(self,ch):
		self.cws += ch
		#sys.stdout.write(ch)
		#sys.stdout.flush()
	
	def printchar(self,ch):		# character space detected
		self.addchar(ch)
		try:					# try to find sequence from Codebook
			val = Codebook[self.cws]
		except:
			val = '*'			# output '*' when cannot find sequence from Codebook
		sys.stdout.write(val)
		sys.stdout.flush()
		self.cws = ''

	def printword(self,ch):		# word space detected
		self.printchar(ch)		# print last character in word
		sys.stdout.write(' ')	# print word space 
		sys.stdout.flush()
	
	# Probabilistic Neural Network - find best matching symbol from mark,space duration pair
	def pnn(self,m,s):
		
		# Symbols are defined by [mark, space] duration examples  
		# Classes are normalized: dit = 0.1  dah = 0.3 char space =0.3 wordspace = 0.7
		# Adding more timing examples may help in accuracy
		# Class S0        S1         S2        S3        S4       S5          S6 noise    S7 noise S8 noise 
		w = [[0.1,0.1],[0.1, 0.3],[0.1,0.7],[0.3,0.1],[0.3,0.3],[0.3,0.7],[0.00,0.05],[0.000,0.5],[0.0,0.8]]
		
		resval = np.linspace(0,1,num=9)
		for i in range(0,9): # go through examples for each class in w[]
			v = 0.0; 
			# PATTERN layer - calculates PDF function for each class 
			v = v + pow(m-w[i][0],2) + pow(s-w[i][1],2)
			v = math.exp(-v/(2 * pow(self.sigma,2)))
			resval.flat[i] = v
			if verbosity: 
				print "pnn: m%f s%f pnn[%d] %f" % (m,s,i,v)
		# OUTPUT layer - select best match  
		val = np.nanargmax(resval)
		if verbosity: 
			print "pnn: argmax %d" % val
		return val

	# decode symbols S0...S5 into characters 
	def decode(self,m, s):  
		self.ticks += m + s
		ten_dits = 5.0*self.twodits # normalize  dit = 0.1 dash = 0.3
		sym = self.pnn(m/ten_dits,s/ten_dits)    
	
		if verbosity: 
			print "\nticks:%f m:%f \ts:%f \t 2dit:%d \t " % (self.ticks, m, s, self.twodits)
			print "\nSymbol S%d " % sym

		if sym ==0:
			self.addchar(".")
		elif sym == 1:
			self.printchar(".")
		elif sym == 2: 
			self.printword(".")
		elif sym == 3: 
			self.addchar("-")
		elif sym == 4: 
			self.printchar("-")
		elif sym == 5: 
			self.printword("-")
		else:
			sys.stdout.write('')  # not known symbol - noise?

	# update speed tracking from (dit,dash) pair over rolling average
	def update_tracking(self, dit, dash):
		if (dit > self.dit_low_limit and dit < self.dit_high_limit): 
			#print "\ndit:%f dash:%f" %(dit,dash)
			self.twodits = self.ra.rolling_avg((dash + dit) / 2.)	
		if (dash > 3*self.dit_low_limit and dash < 3*self.dit_high_limit):
			#print "\ndit:%f dash:%f" %(dit,dash)
			self.twodits = self.ra.rolling_avg((dash + dit) / 2.)	
	
	# detect KEYDOWN/KEYUP edges, measure timing and decode symbols 		  	
	def edge_recorder(self, v, upper, lower):
		KEYUP = 1
		KEYDOWN = 2
		if (v > upper):
			if (self.last == KEYUP):
				# calculate speed when received dit-dah  or dah-dit sequence 
				if (self.lastmark > 2*self.mark): 
					if verbosity: 
						print "update1: %f %f" % (self.mark, self.lastmark)
					self.update_tracking(self.mark, self.lastmark)
				if (self.mark > 2*self.lastmark): 
					if verbosity:
						print "update2: %f %f" % (self.lastmark, self.mark)
					self.update_tracking(self.lastmark, self.mark)

				# decode received "mark-space" symbol 
				self.decode(self.mark, self.space)
				self.lastmark = self.mark
				self.mark = 0
				self.space = 0
				self.last = KEYDOWN
				return self.twodits
			self.mark +=1 
		elif (v < lower):  
			self.last = KEYUP
			self.space +=1
		return self.twodits
# end Morse class

# decode signal envelope into Morse symbols and then characters
def decode_stream(signal,samplerate):
	# create morse object
	m = Morse(signal,samplerate)

	# assume 10ms signal rise time 
	bfv = (samplerate * .010)   
	# moving average filter to smooth signal envelope - reduce noise spikes
	env = np.resize(np.convolve(signal, np.ones(bfv)/bfv),len(signal))
	mx = np.nanmax(env)
	mn = np.nanmin(env)
	mean = np.mean(env)
		
	# prepare arrays to collect plotting data	
	agcv = np.arange(len(env))
	twodits = np.arange(len(env))
	t = np.linspace(0,1,len(env))

	agcpeak = mx
	i = 0
	while i < len(env): 

		# AGC is useful if signal has rapid amplitude variations due to fading, QSB etc.
		# In computer generated audio the amplitude is not varied 
		# Parameters control attack/decay time: fast attack (5) - slow decay (700) 
		if agc:
			z = env[i] - mean	
			if (z > agcpeak):
				agcpeak = decayavg(agcpeak,z,5)
			else:
				agcpeak = decayavg(agcpeak,z,700)
			agcv[i] = agcpeak
			if agcpeak > 0:
				z /= agcpeak
				z = clamp(z,0.,1.)
			up  = UPPER_THRESHOLD
			down = LOWER_THRESHOLD
		else:
			# calculate signal threshold if no AGC is used
			z = env[i] 
			up   = UPPER_THRESHOLD * (mx - mn)
			down = LOWER_THRESHOLD * (mx - mn)
		
		# capture estimated speed over time for plotting 
		twodits[i] = m.edge_recorder(z,up,down)
		
		i += 1
	
	# plot key variables 
	if plotter:
		ax1=plt.subplot(3,1,1)
		plt.plot(signal,'g-') #,t,up*signal,'r--')
		ax1.set_title("Signal")

		ax2=plt.subplot(3,1,2)	
		plt.plot(1.2*samplerate*2/twodits)  # plot speed estimate in WPM
		ax2.set_title("WPM")	
		
		if agc:
			ax3=plt.subplot(3,1,3)
			plt.plot( agcv,'g-')
			ax3.set_title("AGC")
		plt.show()

def demodulate(x,Fs,freq):
	# demodulate audio signal with known CW frequency 
	t = np.arange(len(x))/ float(Fs)
	y =  x*((1 + np.sin(2*np.pi*freq*t))/2 )	
	
	#calculate envelope and low pass filter this demodulated signal
	#filter bandwidth impacts decoding accuracy significantly 
	#for high SNR signals 50 Hz is better, for low SNR 20Hz is better
	# 25Hz is a compromise - could this be made an adaptive value? 
	Wn = 40./ (Fs/2.)  	# 25 Hz cut-off for lowpass  
	b, a = butter(2, Wn)  	# 2nd order butter filter
	z = filtfilt(b, a, abs(y))
	
	#pass envelope magnitude to decoder 
	decode_stream(z,Fs)


# process audio file by demodulator and envelope detector 
def process(fname):
	Fs, x = wavfile.read(fname)
	a = string.split(fname,".wav")
	b = string.split(a[0],"cw")
	sys.stdout.write(b[1])
	sys.stdout.write(",")
	# find frequency peaks of high volume CW signals 
	if fft_scan:
		f,s = periodogram(x,Fs,'blackman',4096,'linear',False,scaling='spectrum')
		# download peakdetect from # https://gist.github.com/endolith/250860
		from peakdetect import peakdet
		threshold = max(s)*0.4  # only 0.4 ... 1.0 of max value freq peaks included
		maxtab, mintab = peakdet(abs(s[0:len(s)/2-1]), threshold,f[0:len(f)/2-1] )

	if plotter:
		plt.plot(f[0:len(f)/2-1],abs(s[0:len(s)/2-1]),'g-')
		print maxtab
		from matplotlib.pyplot import plot, scatter, show
		scatter(maxtab[:,0], maxtab[:,1], color='blue')
		plt.show()
	
	# process all CW stations with higher than threshold volume
	if fft_scan:
		for freq in maxtab[:,0]:
			print "\nfreq:%5.2f" % freq
			demodulate(x,Fs,freq)
	else:
		demodulate(x,Fs,MORSE_FREQUENCY)	
	
def main(*args, **kwargs):
  
	global verbosity
	global plotter
	global agc
	global fft_scan
	
	parser = OptionParser(usage="%prog [OPTIONS] <audio files>\nDecodes morse code from .WAV audio files")

	parser.add_option("-v", "--verbose",
	  action="store_true",
	  dest="verbose",
	  default=False,
	  help="Prints details about errors and calls.")
	parser.add_option("-p", "--plot",
	  action="store_true",
	  dest="plotter",
	  default=False,
	  help="Plot signal, speed estimate and AGC if used")
	parser.add_option("-a", "--agc",
	  action="store_true",
	  dest="agc",
	  default=False,
	  help="Use automatic gain control")
	parser.add_option("-f", "--fft",
	  action="store_true",
	  dest="fft",
	  default=False,
	  help="Use automatic FFT frequency scan")

	(options, args) = parser.parse_args()
	if options.verbose:
		verbosity = True
	if options.plotter:
		plotter = True
	if options.agc:
		agc = True
	if options.fft:
		fft_scan = True
	if len(args) < 1:
		print 'usage: [OPTIONS] <audio files>' 
		exit(1)

	#process all audio files given as arguments
	print "ID,Prediction"
	for i in range(0,len(args)):
		process(args[i])
		print ""
		
	 
    
if __name__ == "__main__":
	main()
