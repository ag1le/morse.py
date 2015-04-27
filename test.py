#!/usr/bin/env python
import sys
from math import exp, sqrt
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy.signal import butter, filtfilt, periodogram

class BreakIt(Exception): pass
class BreakIt2(Exception): pass
VARFLAG = 1
SYMFLAG = 1
SYMPLOT = 1
Debug = 0 
ELM_STATES = 6
RATE_STATES = 5 
PATHS = 20	# min paths is 7 ..10 -> spdhat=20  11..15 -> spdhat=40 (16..20 ->spdhat 50) 21..25->spdhat=60  26..30 -> 70   31.. -> 80
NDELAY = 200

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

class BayesMorse:
	
	# ltrstate to element state mapping 
	# 	 1  2  3  4  5  6  7  8   9 10 11 12 13 14 15 16
	# 	.^ .~ .w .p -^ -~ -w -p  ^. ^- ~. ~- w. w- p. p-
	# 	K=0 DIT, K=1 DAH, K=2 E-SPC, K=3 CHR-SPC, K=4 WRD-SPC, K=5 PAUSE 
	ltr_to_elm_state = [2, 3, 4, 5,  2, 3, 4, 5,  0, 1, 0, 1,  0, 1, 0, 1]
	dit_dah_states = [ 1, 1, 0, 0, 0, 0 ]


	memdel = [\
		[0, 0, 2, 2, 5, 10],\
		[0, 0, 2, 2, 5, 10],\
		[2, 2, 0, 0, 0, 0],\
		[2, 2, 0, 0, 0, 0],\
		[2, 2, 0, 0, 0, 0],\
		[2, 2, 0, 0, 0, 0]]

	def init(self,sample_dur):
		self.init = 0
		self.sample_duration = sample_dur
		self.char = ''
		self.PATHS = PATHS
		self.ltrstate = [5]*PATHS # initialize to 5  
		self.ltrlast = 0
		self.dur = [1000.]*PATHS 
		self.wpm = [(i/5+2)*10 for i in range(PATHS)] #[40 for i in range(PATHS)] #
		self.pathsv = [5]*PATHS
		self.ykkip  = [.5]*PATHS
		self.pkkip  = [.1]*PATHS
		self.sort   = [0]*PATHS
		self.ltrsav = [5]*ELM_STATES*RATE_STATES*PATHS    # was 5
		self.dursav = [0.]*ELM_STATES*RATE_STATES*PATHS
		self.wpmsav = [20]*ELM_STATES*RATE_STATES*PATHS   # was 20 
		self.ykksv  = [0.]*ELM_STATES*RATE_STATES*PATHS
		self.pkksv  = [0.]*ELM_STATES*RATE_STATES*PATHS
		self.Pold	= [1.]*ELM_STATES*RATE_STATES*PATHS 
		self.Pnew   = [[0. for col in range(PATHS)] for row in range(ELM_STATES*RATE_STATES)]
		self.lkhd   = [[0. for col in range(PATHS)] for row in range(ELM_STATES*RATE_STATES)]
		self.lastrs = -1
		# the following are used in trelis()
		self.ltrsv = [0]*NDELAY
		self.ipnod = [1]*PATHS
		self.lmdsav = [[0. for col in range(NDELAY)] for row in range(PATHS)] 
		self.pthtrl = [[0. for col in range(NDELAY)] for row in range(PATHS)] 
		self.nbuf = -1
		self.ndelst = 0
	
	def normalize(self,lst):
		s = sum(lst)
		return map(lambda x: float(x)/s, lst)
    
	def view_P(self):
		#Pin = [[self.ltrsav[k+j*ELM_STATES*RATE_STATES] for j in range(PATHS)] for k in range(ELM_STATES*RATE_STATES)]
		fig,ax1 = plt.subplots(nrows=1) #, figsize=(6,10))
		ax1.imshow(self.lkhd, extent=[0,PATHS,0,ELM_STATES*RATE_STATES])
		ax1.set_title('lkhd')
		plt.show()

		
	#=====================================
	def xtrans(self, elemtype, dur, wpm): 
		""" Calculates keystate transition probability conditioned on 
		    elementype, current duration and data rate in WPM
		    
		elemtype:  	0=dit, 1=dah, 2=ele-space, 3=chr-space, 4=wrd-space, 5=pause
		dur: 		current elemtype duration in milliseconds
		wpm:		current data rate in Word per minute (WPM)
		
		Returns keystate transition probability
		"""
		# TABLES  CONTAIN DENSITY PARMS FOR EACH ELEMTYPE AND DATA RATE.
		# K=0 DIT, K=1 DAH, K=2 E-SPC, K=3 CHR-SPC, K=4 WRD-SPC, K=5 PAUSE
		elem_length = [1, 3, 1, 3, 7, 14]
		# was  3. 1.5  1.
		aparm = [3.,3.,3.,3., 1.5, .1]
		
		mscale = elem_length[elemtype]
		rscale = 1200. / wpm
		alpha = mscale * aparm[elemtype] 
		b0 = dur / (mscale * rscale)
		b1 = (dur + self.sample_duration) / (mscale * rscale) 
		
		if (b1 <= 1.):
			p1 = 1.- .5 * exp(alpha * (b1 - 1.))
			p0 = 1.- .5 * exp(alpha * (b0 - 1.))
			return p1 / p0

		if ((b0 < 1.) and (b1 > 1.)):
			p1 = .5 * exp(-alpha * (b1 - 1.))      #check if -.5  or + .5 
			p0 = 1. - .5 * exp(alpha * (b0 - 1.)) 
			return p1 / p0
		return exp(-alpha * (b1 - b0))

	#==================================================
	def spdtr(self, rate_state, wpm, nxt_elm, cur_elm): 
		"""
		# 	THIS FUNCTION RETURNS THE DATA RATE (SPEED) TRANSITION 
		# 	PROBABILITY BASED ON THE CURRENT ELEM TYPE. THE ALLOWABLE 
		# 	TRANSITION PROBS ARE STORED IN THE TABLE RTRANS. 

		# 	VARIABLES: 
		# 	rate_state	- DATA RATE STATE TO WHICH PATH IS BEING EXTENDED 
		# 	wpm		- DATA RATE ON CURRENT PATH 
		# 	nxt_elm 	- ELEM TYPE FOR NEXT STATE 
		# 	cur_elm	- ELEM TYPE ON CURRENT PATH 



		#PAGES 103-104 IN THESIS - SYMBOL CONDITIONAL TRANSITION PROBABILITIES 
		#IF SAVED ELEMENT AND NEW ELEMENT ARE THE SAME THEN THERE CAN BE NO SPEED CHANGE: 
		"""		
		mempr  = [ \
			[0, 0, 1, 2, 1, 2],\
			[0, 0, 1, 2, 1, 2],\
			[1, 1, 0, 0, 0, 0],\
			[1, 1, 0, 0, 0, 0],\
			[1, 1, 0, 0, 0, 0],\
			[1, 1, 0, 0, 0, 0]]

		#rtrans[2][5] - symbol conditional speed transition probabilities - Page 104 - Table X 
		#used in spdtr() 	
		# 1st row:	dot, dash, e-sp, w-s  by rate_state
		# 2nd row:	c-sp, pause   by rate_state 
		rtrans =[[.1,  .2, .4, .2, .1],\
				[ .15, .2, .3, .2, .15]]
		#SAVED ELEMENT AND NEW ELEMENT ARE THE SAME 
		if (cur_elm == nxt_elm):
			# DON'T MAKE SPEED CHANGES DURING ELEMENT DURATION
			if (rate_state != 2):
				return 0.	

		#OTHERWISE, OBTAIN SPEED TRANSITION PROB 
		wpm_delta = self.memdel[nxt_elm][cur_elm]
		index     = mempr[nxt_elm][cur_elm]
		
		if (index == 0): 
			return 0.
			
		wpm_change = (rate_state - 2) * wpm_delta
		new_wpm = wpm + wpm_change
		ret_val = rtrans[index-1][rate_state]
		if (Debug):
			print "spd_ret:%f new_wpm:%d wpm_change:%d index %d rate_state %d nxt_elm:%d cur_elm:%d" %(ret_val,new_wpm,wpm_change,index,rate_state,nxt_elm,cur_elm)
			
		if (new_wpm > 80): #if speed rate is > 60 WPM TRANSITION PROBABILITY = 0 
			return 0.	
		if (new_wpm < 5):
			return 0.	  #if speed rate is < 10 WPM TRANSITION PROBABILITY = 0 
		return ret_val


	#===========================================================================		
	def ptrans(self, elem_state, rate_state, ip, ptrx, psum, pint, n): 

		"""
		# 	PTRANS() RETURNS THE PATH CONDITIONAL TRANSITION 
		# 	PROBABILITIES TO EACH ALLOWABLE STATE N. 
		# 	VARIABLES: 
		# 	elem_state-      INPUT CURRENT ELEMENT STATE 
		# 	rate_state-      INPUT CURRENT DATA RATE STATE 
		# 	ltrstate-     INPUT IDENTITY OF CURRENT LTR STATE 
		# 	PTRX-       INPUT KEYSTATE TRANSITION PROBABILITY 

		# 	FUNCTION FUNCTION USED: 
		# 	SPDTR-      RETURNS DATA RATE TRANSITION PROBS,CONDITIONED ON CURRENT SPACE TYPE. 
		"""
		# 	ELEMTR-     ELEMENT TRANSITION PROBABILITY MATRIX 
		#TABLE XII Second Order Markov Symbol Transition Matrix - Page 105 Table XII
		#elemtr[6][16] 		
		### .^  .~  .w  .p   -^  -~  -w  -p    ^.    ^-    ~.    ~-    w.    w-   p.   p-
		# . 
		# -
		# ^
		# ~ 
		# w
		# p
		elemtr = [\
			[.55, .5, .5, .5, .55, .5, .5, .5,   0.,   0.,   0.,   0.,   0.,   0.,  0.,  0.],\
			[.45, .5, .5, .5, .45, .5, .5, .5,   0.,   0.,   0.,   0.,   0.,   0.,  0.,  0.],\
			[ 0., 0., 0., 0., 0.,  0., 0., 0., .581,  .54, .923, .923, .923, .923, .95, .95],\
			[ 0., 0., 0., 0., 0.,  0., 0., 0., .335, .376, .062, .062, .062, .062, .04, .04],\
			[ 0., 0., 0., 0., 0.,  0., 0., 0., .069, .069, .012, .012, .012, .012, .009,.009],\
			[ 0., 0., 0., 0., 0.,  0., 0., 0., .015, .015, .003, .003, .003, .003, .001,.001]]


		wpm = self.wpm[ip]
		ltrstate = self.ltrstate[ip]-1

		# 	IF THE SAVED ELEMENT AND THE ELEMENT OF THE STATE 
		# 	N TO WHICH THE PATH IS BEING EXTENDED ARE THE 
		# 	SAME, THEN THE STATE TRANS PROB IS SIMPLY KEYSTATE TRANS PROB: 
		if (elem_state == self.ltr_to_elm_state[ltrstate]):
			pint[n] = ptrx
			
			# 	IF CURRENT DATA RATE STATE  != 2, THEN RETURN KEYSTATE TRANS PROB
			# See page 104 in thesis 
			if (rate_state != 2):
				pint[n] = ptrx
				psum += pint[n]	 # AG1LE added - debug K=0 state 
				#print "ptrans: rate_state=%d pint[%d]=%f psum=%f"%(rate_state,n,pint[n],psum)
				return psum,pint
		else:
			# 	OTHERWISE: 
			# 	OBTAIN ELEM TRANS PROBS TABLE: 
			pelem = elemtr[elem_state][ltrstate]

			# 	COMPUTE ELEM-CONDITIONAL SPEED TRANSITION PROB: 
			prate = self.spdtr(rate_state, wpm, elem_state, self.ltr_to_elm_state[ltrstate])

			# 	TRANSITION PROBABILITY IS THE PRODUCT: 
			pint[n] = (1. - ptrx) * pelem * prate
			psum += pint[n]
		return psum,pint


	#=====================================
	def trprob(self,ip):
		psum = 0.
		pint = [0. for col in range(ELM_STATES*RATE_STATES)]
		dur = self.dur[ip]
		wpm = self.wpm[ip]

		# LETTER STATE IS ZERO, INITIALIZE TRANSITION PROBABILITIES TO ZERO FOR PATH IP
		if (self.ltrstate[ip] == 0):
			print "trprob: ltrstate[%d]" %(ip)
			for n in range(ELM_STATES*RATE_STATES):
				self.Pnew[n][ip] = 0.
			return 0.
		elemtype = self.ltr_to_elm_state[self.ltrstate[ip]-1]
		
		# COMPUTE KEYSTATE TRANSITION PROBABILITY: 
		ptrx = self.xtrans(elemtype,dur,wpm)

		# 	FOR EACH STATE, COMPUTE STATE TRANSITION PROBABILITY: 
		
		for es in range(ELM_STATES): # 6 element states 0=dit,1=dah, 2=e-spc, 3=chr-s, 4=wrd-s, 5=pause
			for ss in range(RATE_STATES): # 5 speed (rate) states -2 -1 0 1 2 
				n = ss * ELM_STATES + es
				psum,pint = self.ptrans(es, ss, ip, ptrx, psum, pint, n)

		if (psum == 0.0):
			print "\ntrprob: psum = 0"
			return 0.

		# Normalize with psum 
		for n in range(ELM_STATES*RATE_STATES): 
			self.Pnew[n][ip] = pint[n] / psum
		
		return 0.

	#=====================================
	def path(self, ip):
		"""
		# PATH COMPUTES THE LTR STATE, DURATION, AND DATA RATE OF 
		# EACH NEW PATH EXTENDED TO STATE N 

		# VARIABLES: 
		# IP-		SAVED PATH IDENTITY 
		# ltrstate-	LTR STATE OF SAVED PATH 
		# DUR-		DURATION OF ELEMENT ON SAVED PATH 
		# ILRATE-	DATA RATE OF ELEMENT ON SAVED PATH 
		# ltrsav-	NEW LTR STATES FOR EACH PATH EXTENSION 
		# DURSAV-	NEW ELEM DURATIONS FOR EACH PATH EXTENSION 
		# wpmsav-	NEW DATA RATES FOR EACH PATH EXTENSION 
		# J-		NEW PATH IDENTITY 
		"""
		# THE LETTER TRANSITION TABLE, MEMFCN USED TO LABEL NEW PATH EXTENDED TO STATE N 
		# ltrstate to element state mapping 
		# COLS:	 1  2  3  4  5  6  7   8  9 10 11 12 13 14 15 16
		# 	    .^ .~ .w .p -^ -~ -w -p  ^. ^- ~. ~- w. w- p. p-
		# ROWS: K=0 DIT, K=1 DAH, K=2 E-SPC, K=3 CHR-SPC, K=4 WRD-SPC, K=5 PAUSE 
		#
		#         .^ .~  .w  .p   -^ -~  -w  -p  ^.  ^- ~.   ~- w.   w- p.   p-
		memfcn=[[ 9, 11, 13, 15,  9, 11, 13, 15, 9,  0, 11,  0, 13,  0, 15,  0],\
			    [10, 12, 14, 16, 10, 12, 14, 16, 0, 10,  0, 12,  0, 14,  0, 16],\
			    [ 1,  0,  0,  0,  5,  0,  0,  0, 1,  5,  1,  5,  1,  5,  1,  5],\
			    [ 0,  2,  0,  0,  0,  6,  0,  0, 2,  6,  2,  6,  2,  6,  2,  6],\
			    [ 0,  0,  3,  0,  0,  0,  7,  0, 3,  7,  3,  7,  3,  7,  3,  7],\
			    [ 0,  0,  0,  4,  0,  0,  0,  8, 4,  8,  4,  8,  4,  8,  4,  8]] 

	
	#FOR EACH ELEM STATE K, AND EACH SPEED I, COMPUTE: 
		for k in range(ELM_STATES): # 6 element states 0=dit,1=dah, 2=e-spc, 3=chr-s, 4=wrd-s, 5=pause
			for i in range(RATE_STATES): # 5 speed (rate) states -2 -1 0 1 2 
				
	#NEW PATH IDENTITY: 
				j = ip * ELM_STATES*RATE_STATES + i * ELM_STATES + k
				
	#NEW LTR STATE: 
				if (self.ltrstate[ip] == 0):
					self.ltrsav[j] = 0
					#print "ltrstate[%d]=%d"%(ip,self.ltrstate[ip])
					continue
				
				self.ltrsav[j] = (memfcn[k][self.ltrstate[ip]-1])
				#print "path: ltrsav[%3d]=%3d = memfcn[%d][ltrstate[%d]=%d]]=%d-1" %(j,self.ltrsav[j],k,ip,self.ltrstate[ip],memfcn[k][self.ltrstate[ip]])
				if (self.ltrsav[j] == 0):  # check if we got 0
					#print "path327: ltrsav[%d]=%d k=%d ltrstate[ip=%d]=%d"%(j,self.ltrsav[j],k,ip,self.ltrstate[ip])
					continue
	#  NEW DURATION: OBTAIN KEYSTATE OF SAVED PATH AND NEW STATE: 
				ks_s = self.ltr_to_elm_state[self.ltrstate[ip]-1]
				ixl = self.dit_dah_states[ks_s]  
				ixs = self.dit_dah_states[k]	  

	# CALCULATE DURATION - ADD SAMPLE DURATION 5 ms FOR EACH VALID PATH 
				self.dursav[j] = self.dur[ip] * (1 - ixs - ixl + (ixs << 1) * ixl) + self.sample_duration
				#print "dursav[%d]=%f dur[%d]=%f * sum=%d k_new=%d k_saved=%d + sample_dur=%f" %(j,self.dursav[j],ip,self.dur[ip],(1 - ixs - ixl + (ixs << 1) * ixl),k,ks_s,sample_duration)
	# 	NEW DATA RATE - 
				self.wpmsav[j] = self.wpm[ip] + (i - 2) * self.memdel[k][ks_s]
				#print "wpmsav[%d]=%d = wpm[%d]=%d +(i-2)=%d * memdel[%d][%d]=%d" %(j,self.wpmsav[j],ip,self.wpm[ip],(i-2),k,ks_s,self.memdel[k][ks_s])
		return 0
	
	#=====================================
	def model(self,ip, ielm, ixs):
		# 	THIS FUNCTION COMPUTES THE PARAMETERS OF THE 
		# 	OBSERVATION STATE TRANSITION MATRIX PHI AND THE 
		# 	MEASUREMENT MATRIX

		# 	VARIABLES: 
		# 		DUR-	INPUT ELEMENT DURATION 
		# 		WPM-	INPUT SAVED RATE 
		# 		IELM-	INPUT ELEMENT TYPE 
		# 		ISR-	INPUT RATE OF NEW STATE 
		# 		IXS-	INPUT KEYSTATE OF NEW STATE 
		# 		PHI-	OUTPUT STATE TRANSITION MATRIX ENTRY FOR SIGNAL AMPLITUDE STATE 
		# 		QA-	OUTPUT COVARIANCE FOR AMPLITUDE STATE 

		dur = self.dur[ip]
		wpm = self.wpm[ip]

		if (Debug):
			print "model: wpm=%d" %wpm	
		# 	COMPUTE PHI AND AMPLITUDE STATE VARIANCE (Q): 
		r1 = 1200. / wpm
		bauds = dur / r1
		if (bauds >= 14.):
			bauds = 14.
		# element type 'dit' or 'dah'
		if (ielm < 2):  # was < 2
			qa = 1.e-4
			phi = 1.
			return qa,phi
		# else element type el-spc, chr-spc, wrd=spc  or pause
		# next state is also 'dit' or 'dah'
		if (ixs != 0):
			phi = 1.
			qa = exp((bauds - 14.) * .6) * .15
			qa += bauds * .01 * exp((1. - bauds) * .2)
			return qa,phi
		#next state is el-spc, chr-spc, wrd=spc  or pause
		xsamp = r1 * 22.4
		d1 = (-2. / xsamp)
		phi = pow(10.0, d1)

		if (bauds >= 14.):
			phi = 1.

		qa = 0.
		return qa,phi



	#====================================================
	def kalfil(self, z, ip, rn, ixs, kelem, jnode, pinr):

		#   THIS FUNCTION COMPUTES THE ARRAY OF KALMAN FILTER 
		#   RECURSIONS USED TO DETERMINE THE LIKELIHOODS. 

		#   VARIABLES: 
		#       Z -	INPUT MEASUREMENT 
		#       IP -	INPUT PATH IDENTITY 
		#       RN -	INPUT NOISE POWER ESTIMATE 
		#       IXS -	INPUT KEYSTAT OF NEW NODE 
		#       KELEM -	INPUT ELEM STATE OF NEW NODE 
		#       ISRATE 	INPUT SPEED STATE OF NEW NODE 
		#       DUR - 	INPUT CURRENT DURATION OF ELEMENT ON IP 
		#       WPM 	INPUT SPEED STATE ON PATH IP 
		#       LKHDJ - OUTPUT CALCULATED LIKELIHOOD VALUE 

		#   FUNCTIONS USED 
		#       MODEL - OBTAINS THE SIGNAL-STATE-DEPENDENT LINEAR 
		#       	     MODEL FOR THE KALMAN FILTER RECURSIONS 

		#   IF TRANSITION PROBABILITY IS VERY SMALL, DON'T 
		#   BOTHER WITH LIKELIHOOD CALCULATION: 

		if (pinr <= 0.0001):
			return 0.

	#   OBTAIN STATE-DEPENDENT MODEL PARAMETERS: 
		qa,phi = self.model(ip, kelem, ixs)

	# 	COMPUTE MEASUREMENT COEFFICIENT: 
		hz = float(ixs)
		
	# 	GET PREVIOUS ESTIMATES FOR PATH IP 

		ykk = self.ykkip[ip]
		pkk = self.pkkip[ip]

	#  IMPLEMENT KALMAN FILTER FOR THIS TRANSITION 

		ypred = phi * ykk
		ppred = phi * pkk * phi + qa
		pz = hz * ppred + rn
		pzinv = 1. / pz
		g = ppred * hz * pzinv
		pest = (1. - g * hz) * ppred
		zr = z - hz * ypred

		self.ykksv[jnode] = ypred + g * zr
		self.pkksv[jnode] = pest
		if (self.ykksv[jnode] <= .01): #was 0.01
			self.ykksv[jnode] = .01    #was 0.01

	# Computing 2nd power 
		a = .5*pzinv*(zr * zr)
		if (a > 1000.):
			return 0.
		return (1. / sqrt(pz)) * exp(-a)



	#==========================
	def likhd(self, z, rn, ip):
		""" 
		# 	THIS FUNCTION CALCULATES,FOR EACH PATH 
		# 	EXTENSION TO STATE N, THE LIKELIHOOD OF THAT 
		# 	TRANSITION GIVEN THE MEASUREMENTZ. IT USES 
		# 	AN ARRAY OF LINEAR (KALMAN) FILTERS TO DO SO. 

		# 	VARIABLES: 
		# 	Z- 	INPUT MEASUREMENT 
		# 	RN-	INPUT NOISE POWER ESTIMATE 
		# 	IP-	INPUT SAVED PATH IDENTITY 
		# 	LAMBDA-	INPUT SAVED LTR STATE IDENTITY 
		# 	DUR-	INPUT SAVED DURATION OF ELEMENT ON PATH IP 
		# 	ILRATE-	INPUT SAVED DATA RATE (SPEED) 
		# 	P-	INPUT TRANSITION PROBABILITIES 
		# 	LKHD-	OUTPUT COMPUTED LIKELIHOODS FOR EACH TRANS 

		#  FUNCTIONS USED: 
		# 	KALFIL-KALMAN FILTER FOR EACH NEW PATH 
		"""
		if (self.ltrstate[ip] == 0):
			return 0
		
		#OBTAIN SAVED KEYSTATE: 		
		kelem = self.ltr_to_elm_state[self.ltrstate[ip]-1]

		#FOR EACH ELEMENT STATE: 
		for k in range(ELM_STATES):
			for i in range(RATE_STATES):
				#OBTAIN KEYSTATE, RATE STATE, STATE N, NEW NODE: 
				ixs = self.dit_dah_states[k]
				n = i  * ELM_STATES + k
				j = ip * ELM_STATES*RATE_STATES + n
				#COMPUTE AND STORE LIKELIHOOD: 
				self.lkhd[n][ip] = self.kalfil(z, ip, rn, ixs, kelem, j, self.Pnew[n][ip])  
		return 0

	#=====================================
	def probp(self, isave):
		"""
			PROBP COMPUTES THE POSTERIOR PROBABILITY OF EACH NEW PATH 
			VARIABLES: 
			POLD-		INPUT: SAVED PROBS OF PRIOR PATHS 
			PNEW-		INPUT TRANSISTION PROBABILITIES 
			LKHD-		INPUT LIKELIHOODS OF EACH TRANSTION 
			PSUM-		NORMALIZING CONSTANT (SUM OF P(J)) 
			OUTPUT:		COMPUTED POSTERIOR PROBS OF NEW PATHS 
		"""
		psav = [0.]*ELM_STATES*RATE_STATES*PATHS
		psum = 0.

		#FOR EACH SAVED PATH, EACH TRANSITION
		for i in range(isave):
			for n in range(ELM_STATES*RATE_STATES):
				#COMPUTE IDENTITY OF NEW PATH: 
				j = i * ELM_STATES*RATE_STATES + n
				#PRODUCT OF PROBS, ADD TO PSUM 
				#NOTE: Pold[i] index is isave - using savep() stored values from previous sample
				psav[j] = self.Pold[i] * self.Pnew[n][i] * self.lkhd[n][i]
				psum += psav[j]
		#NORMALIZE TO GET PROBABILITIES SAVE: 
		if (psum == 0.0):
			#print "\nprobp: psum = 0"
			return 

		for j in range(isave * ELM_STATES*RATE_STATES):
			self.Pold[j] = psav[j] / psum
		return 

	#=====================================
	def sprob(self,isave):
		'''
			SPROB COMPUTES THE POSTERIOR PROBS OF THE ELEMENT 
			STATES, DATA RATE STATES, AND KEYSTATES BY SUMMING 
			OVER THE APPROPRIATE PATHS. 

			VARIABLE: 
			POLD	INPUT PATH PROBABILITIES 
			ISAVE- 	NUMBER OF PATHS SAVED 
			PSELEM-	OUTPUT ELEMENT PROB 
			KHAT-	OUTPUT ESTIMATED ELEMENT STATE
			SPDHAT-	OUTPUT SPEED ESTIMATE (DATA RATE WPM) 
			PX- 	OUTPUT KEYSTATE PROBABILITY 
		'''
		#INITIALIZE:
		spdhat = 0.
		px = 0.
		pselem = [0.]*ELM_STATES

		#FOR EACH STATE EXTENSION OF PATH M: 
		#OBTAIN ELEMENT STATE PROBS,KEYSTATE PROBS,SPEED EST: 
		for k in range(ELM_STATES):
			pselem[k] = 0.
			for i in range(RATE_STATES):
				n = i * ELM_STATES + k
				for m in range(isave):
					j = m * ELM_STATES*RATE_STATES + n
					pselem[k] += self.Pold[j]
					spdhat += self.wpmsav[j] * self.Pold[j]
					# capture 'dit' and 'dah' probs in px 
					if (k < 2):
						px += self.Pold[j]
					#print "pselem[%d]=%f spdhat=%f wpmsav[%d]=%d Pold=%f"%(k,pselem[k],spdhat,j,self.wpmsav[j],self.Pold[j])
		pelm = 0.
		for k in range(ELM_STATES):
			# IF WANT TO PRINT ELEMENT PROBABILITIES BY SAMPLE ENABLE VARFLAG
			if (VARFLAG):
				sys.stdout.write("\t%4.2f" % (pselem[k]))
			if (pselem[k] >= pelm):
				pelm = pselem[k]
				khat = k
		if (VARFLAG):
			sys.stdout.write( "\t%4.2f\t%2d %4.3f " % (spdhat,khat,pelm))
		return pselem,pelm, spdhat, khat, px


	#=====================================
	def savep(self,isave):
		# 	THIS FUNCTION PERFORMS THE ALGORITM TO SAVE 
		# 	THE PATHS WITH HIGHEST POSTERIOR PROBABILITY. 
		# 	IT WILL SAVE A MINIMUM OF 7 PATHS (ONE FOR EACH *
		# 	ELEMENT STATE AND ONE ADDITIONAL NODE), AND 
		# 	A MAXIMUM OF "PATHS" PATHS.  WITHIN THESE LIMITS, IT 
		# 	SAVED ONLY ENOUGH TO MAKE THE TOTAL SAVED PROBABILITY 
		# 	EQUAL TO POPT. 
		# 	ADDITIONALLY, IT RE-SORTS THE LTRSTATE,DUR,AND WPM 
		# 	ARRAYS TO CORRESPOND TO THE SAVED NODES. 


		# 	VARIABLES: 
		# 		POLD	INPUT PROBABILITY ARRAY OF NEW NODES 
		# 		PATHSV-	OUTPUT ARRAY OF THE PREVIOUS NODES TO 
		# 				WHICH THE SAVED NODES ARE CONNECTED. 
		# 		ISAVE-	INPUT: NO. OF PREVIOUS NODES SAVED 
		#	OUPUT: NO. OF NODES SAVED AT CURRENT STAGE 
		# 		IMAX-	INDEX OF HIGHEST PROBABILITY NODE 
		# 		LTRSAV-	INPUT ARRAY OF LTR STATES AT EACH NEW NODE 
		# 		DURSAV-	INPUT ARRAY OF SAVED DURATIONS 
		# 		WPMSAV-	INPUT ARRAY OF SAVED RATES 
		# 		LTRSTATE-OUTPUT ARRAY OF SAVED LTR STATES, SORTED 
		# 				 ACCORDING TO PROBABILITY 
		# 		DUR-	OUTPUT ARRAY OF SORTED DURATIONS 
		# 		WPM-	OUTPUT ARRAY OF SORTED RATES 
		imax = -1
		popt = .90 # was 0.9f 
		psav = [0]*PATHS  # save max probabilities found from Pold[] into psav[]
		iconv = [0]*PATHS
		ipsav = 0
		jsav = 0
		isavm1 = 0
		nplus1 = 0
		
		# 	SELECT SIX HIGHEST PROB ELEMENT STATE NODES: 
		nsav = 0  #was 0 
		psum = 0.
		for k in range(ELM_STATES):		
			pmax = 0.
			for ip in range(isave):
				for i in range(RATE_STATES):
					j = ip * ELM_STATES*RATE_STATES + i * ELM_STATES + k
					if (self.Pold[j] >= pmax): #was >=
						pmax = self.Pold[j]
						jsav = j
						ipsav = ip
			if (pmax > 0.00000):
				psum += pmax
				psav[nsav] = pmax
				self.pathsv[nsav] = ipsav
				self.sort[nsav] = jsav
				nsav +=1


		# 	SELECT ENOUGH ADDITIONAL NODES TO MAKE TOTAL 
		# 	PROBABILITY SAVED EQUAL TO POPT, OR A MAX OF 'PATHS': 
		while True:
			pmax = 0.
			for ip in range(isave):  
				try:
					for states in range(ELM_STATES*RATE_STATES): 
						j = ip * ELM_STATES*RATE_STATES + states
						for i in range(nsav): 
							if (j == self.sort[i]):
								raise BreakIt
							if (self.Pold[j] > pmax):
								pmax = self.Pold[j]
								jsav = j
								ipsav = ip
				except BreakIt:
					pass
					#print "savep657: psum=%f pmax=%d nsav=%d isave=%d" %(psum,pmax,nsav,isave)
			psum += pmax
			psav[nsav] = pmax
			self.pathsv[nsav] = ipsav
			self.sort[nsav] = jsav
			nsav += 1
			if (psum >= popt): 
				break 
			if (nsav > PATHS-1):
				break
		
	# 	NEW ISAVE EQUALS NO. OF NODES SAVED: 
		isave = nsav

	# 	SORT THE SAVED ARRAYS TO OBTAIN THE ARRAYS 
	# 	TO BE USED FOR THE NEXT ITERATION: 
		if (psum ==0.0):
			print "error: savep line 670: psum = 0"
			return isave,imax
			
		for i in range(isave): 
			self.Pold[i] = psav[i] / psum
			self.ltrstate[i] = self.ltrsav[self.sort[i]]
			self.dur[i] = self.dursav[self.sort[i]]
			self.wpm[i] = self.wpmsav[self.sort[i]]
			self.ykkip[i] = self.ykksv[self.sort[i]]
			self.pkkip[i] = self.pkksv[self.sort[i]]
			

		for i in range(isave):
			iconv[i] = 1

		# FIND IF CHANGE IN SPEED, DURATION OR LETTERSTATE - MARK WITH ZERO IF NO CHANGE
		for np in range(isave-1):
			if (iconv[np] != 0):
				nplus1 = np +1
				for k in range(nplus1,isave):
					try:
						if (iconv[k] == 0):
							raise BreakIt2
						if (self.wpm[k] != self.wpm[np]):
							#sys.stdout.write("w")
							raise BreakIt2
						if (self.dur[k] != self.dur[np]):
							#sys.stdout.write("d")
							raise BreakIt2
						if (self.ltrstate[k] != self.ltrstate[np]):
							#sys.stdout.write("l")
							raise BreakIt2
						iconv[k] = 0
					except BreakIt2:
						pass
						#print "savep724: ltrstate[%d]=%d ltrstate[%d]=%d iconv=%d" %(k,self.ltrstate[k],np,self.ltrstate[np],iconv[k])
		
		# SAVE SORTED VALUES FOR NEXT SAMPLE 
		psum = 0.
		np = 0
		for i in range(1,isave):
			if (iconv[np] != 0):  # was iconv[i] 
				self.Pold[np] = self.Pold[i]
				psum += self.Pold[np]
				#print "Pold[%d]=%f Pold[%d]=%f Psum=%f" %(np,self.Pold[np],i,self.Pold[i],psum)
				self.sort[np] = self.sort[i]
				self.ltrstate[np] = self.ltrstate[i]
				self.dur[np] = self.dur[i]
				self.wpm[np] = self.wpm[i]
				self.ykkip[np] = self.ykkip[i]
				self.pkkip[np] = self.pkkip[i]
				self.pathsv[np] = self.pathsv[i]			
				np +=1
		
	# 	ALSO OBTAIN HIGHEST PROBABILITY NODE: 
		if (psum ==0.0):
			print "error: savep line 724: psum = 0"
			return isave,imax

		isave = np		

		pmax = 0.
		self.Pold = self.normalize(self.Pold)
		imax = self.Pold.index(max(self.Pold))
		
		#psum = sum(self.Pold)
		#for i in range(isave): 
		#	self.Pold[i] /= psum
		#	#if (VARFLAG):
		#	#	sys.stdout.write(" p[%d]%7.5f" %(i,self.Pold[i]))
		#	if (self.Pold[i] > pmax):
		#		pmax = self.Pold[i]
		#		imax = i
		return isave, imax





	#=====================================
	def translate_ltr(self, ltr):
		def traverse_dit():
			self.char += '.'
			return
		def traverse_dah():
			self.char += '-'
			return
		def tree_leaf():
			s = self.char #Codebook[self.char]
			sys.stdout.write(s)
			self.char = ''
			return
		def tree_leaf_s():
			s = self.char #Codebook[self.char]
			sys.stdout.write(s)
			self.char = ''
			sys.stdout.write(" ")
			return
		def do_nothing():
			return
		

		options = {
		15:	traverse_dit, # p.
		13:	traverse_dit, # w.
		11:	traverse_dit, # ~.
		9:  traverse_dit, # ^.
		16:	traverse_dah, # p-
		14:	traverse_dah, # w-	
		12:	traverse_dah, # ~- 
		10:	traverse_dah, # ^-
		6:	tree_leaf,    # -~
		2:	tree_leaf,    # .~
		3:	tree_leaf_s,    # .w
		7:	tree_leaf_s,    # -w
		4:	tree_leaf_s,    # .p
		8:	tree_leaf_s,    # -p
		1:	do_nothing,   # .^
		5:	do_nothing    # -^
		}
	
		# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
		#.^ .~ .w .p -^ -~ -w -p ^. ^- ~. ~- w. w- p. p-
		arr = ['','.^', '.~', '.w','.p', '-^', '-~', '-w', '-p', '^.', '^-', '~.', '~-', 'w.', 'w-', 'p.', 'p-']
		
		if (ltr != self.ltrlast):
			
			if (SYMFLAG):
				sys.stdout.write("\t%2d %s " % (ltr,arr[int(ltr)]))
				sys.stdout.flush()
			try:
				options[ltr]()
			except:
				sys.stdout.write('*')
		
		if (0):
			sys.stdout.write("\n%d %s" % (ltr,arr[int(ltr)]))
			sys.stdout.flush()
			#sys.stdout.write("%2d %s " % (ltr,arr[int(ltr)]))
			
		sys.stdout.flush()
		self.ltrlast = ltr
		return ltr
		
	#=====================================
	def trelis(self, isave, imax):
		# Local variables 
		#int i, k, ip, ieq, ltr, ndel, retstat
		#static int isavg, init=0
		#static float xsavg, xmmax, xnmax
		#int ndlavg
		#static float xdlavg

		#    THIS FUNCTION STORES THE SAVED NODES AT EACH 
		#    STAGE AND FORMS THE TREE OF SAVED PATHS LINKING 
		#    THE NODES. DECODING IS ACCOMPLISHED BY FINDING 
		#    THE CONVERGENT PATH IF IT OCCURS WITHIN A MAXIMUM 
		#    DELAY SET BY THE PARAMETER NDELAY. IF CONVERGENCE 
		#    TO A SINGLE PATH DOES NOT OCCURS, THEN DECODING IS 
		#    DONE BY READING THE LETTER ON THE PATH WITH HIGHEST 
		#    PROBABILITY 

		
		retstat = 0
	# 	STORE PATHSV AND CORRESPONDING LTRSTATE IN THE 
	# 	TRELLIS USING A CIRCULAR BUFFER OF LENGTH NDELAY : 
		self.nbuf += 1
		if (self.nbuf == NDELAY):
			self.nbuf = 0

		for i in range(isave): 
			self.pthtrl[i][self.nbuf] = self.pathsv[i]		
			self.lmdsav[i][self.nbuf] = self.ltrstate[i]
			#print "833: ltrstate[%d]=%d"%(i,self.ltrstate[i])

	# 	PERFORM DYNAMIC PROGRAM ROUTINE TO FIND CONVERGENT PATH: 
		k = 0
		for i in range(isave):
			self.ipnod[i] = i


	#L190:
		condition = True
		while condition: 
			k += 1
			if (k != NDELAY):

			# 	IF IP EQUALS INDEX OF HIGHEST PROBABILITY NODE, STORE NODE TO IMAX 
				for ip in range(isave): #(ip = 1 ip <= *isave ++ip) {
					i = self.nbuf - k + 1
					if (i < 0):
						i = NDELAY + i
					#print "trelis: init ipnod[%d]%d [%d]=%d"%(ip,self.ipnod[ip],i,self.pthtrl[self.ipnod[ip]][i])
					self.ipnod[ip] = self.pthtrl[self.ipnod[ip]][i]
					if (ip == imax):
						imax = self.ipnod[ip]

			# 	IF ALL NODES ARE EQUAL,THEN PATHS CONVERGE: 

				for ieq in range(2,isave): #(ieq = 2 ieq <= *isave ++ieq) {
					if (self.ipnod[0] != self.ipnod[ieq - 1]):	
						continue
					else: 
						condition = False

		# 	PATHS CONVERGE SET NDEL: 
			ndel = k + 1

		# 	IF POINT OF CONVERGENCE IS SAME AS IT WAS ON 
		# 	LAST CALL, THEN NO NEED TO RE-DECODE SAME NODE: 
			if (ndel == self.ndelst + 1):
				self.ndelst = ndel
				return retstat

		# 	IF POINT OF CONVERGENCE OCCURS AT SAME DELAY AS LAST CALL, THEN TRANSLATE: 
			if (ndel == self.ndelst):
				i = self.nbuf - ndel + 1
				if (i < 0):
					i = NDELAY + i
				ltr = self.lmdsav[self.ipnod[0]][i]
				#print "trelis: ipnod[0]:%d i:%d ltr=%d" %(self.ipnod[0],i,ltr)
				self.ndelst = ndel
				retstat = self.translate_ltr(ltr)
				return retstat,imax


		# 	OTHERWISE,POINT OF CONVERGENCE HAS OCCURED 
		# 	EARLIER ON THIS CALL, SO NEED TO TRANSLATE 
		# 	EVERYTHING ON THE CONVERGENT PATH FROM 
		# 	PREVIOUS POINT OF CONVERGENCE TO THIS POINT: 
		#L350:
			kd = 0
			ip = self.ipnod[0]
			for k in range(ndel,self.ndelst): #(k = ndel k <= ndelst ++k) {
				kd +=1
				i = self.nbuf - k + 1
				if (i < 0):
					i = NDELAY + i
				self.ltrsv[kd - 1] = self.lmdsav[ip][i]
				ip = self.pthtrl[ip][i]


		# 	REVERSE ORDER OF DECODED LETTERS, SINCE THEY 
		# 	WERE OBTAINED FROM THE TRELLIS IN REVERSE 
		# 	TRANSLATE EACH: 

			for i in range(kd): #(i = 1 i <= kd ++i) {
				ltr = self.ltrsv[kd-i]
				retstat = self.translate_ltr(ltr)
				print "reverse order"
			self.ndelst = ndel
			return retstat,imax

	#L700:
	# 	PATHS HAVE NOT CONVERGED AT MAXIMUM ALLOWABLE 
	# 	DELAY, SO TRANSLATE WHAT IS ON HIGHEST 
	# 	PROBABILITY PATH: 
		print "L700"
		ndel = NDELAY
		i = self.nbuf - NDELAY + 1
		if (i <= 0):
			i = NDELAY + i
		ltr = self.lmdsav[imax][i]
		retstat = self.translate_ltr(ltr)

	# 	PRUNE AWAY NODES WHICH ARE NOT ON THIS PATH: 
		for k in range(isave): #(k = 1 k <= *isave ++k) {
			if (self.ipnod[k] != imax):
				self.ltrstate[k] = 0 #was 0 
				print "prune away nodes"
	#L800:
		self.ndelst = ndel
		return retstat,imax


	# returns AGC decay values 
	def decayavg(self,average,input, weight):

		if (weight <= 1.0): 
			return input
		else:
			return input * (1.0 / weight) + average * (1.0 - (1.0 / weight))
		
	def noise(self,z,rn,mn):
		a = rn if (z > 2*mn ) else z
		rn = self.decayavg(rn,a,60)
		return rn


def decode_stream(signal,samplerate):
	x = []
	y = []
	y1 = []
	y2 = []
	y3 = []
	y4 = []
	y5 = []
	y6 = []
	rn = 0.05
	_isave = PATHS
	imax = 0
	pselem =[0.]*6
	pelm = 0.
	spdhat = 0.
	px = 0.

	decimate = samplerate/200
	sample_dur = 1000.*decimate/float(samplerate)

	m = BayesMorse()    
	m.init(sample_dur)  

	z = signal[1::decimate]/max(signal)
	mn = abs(min(z))
	print "length=%d min=%f samplerate=%d decimate=%d sample=%f"%(len(z), mn, samplerate, decimate, sample_dur)
	
	for j in range(len(z)):
		rn = m.noise(z[j],rn,mn)
		if (VARFLAG):
			sys.stdout.write("\n%5d %4.2f\t%4.2f"%(j,z[j],rn))
	#process 
		for ip in range(_isave):
			m.trprob(ip)
			m.path(ip)
			m.likhd(z[j], rn, ip)
		m.probp(_isave)
		pselem,pelm, spdhat, khat, px = m.sprob(_isave)
		_isave,imax = m.savep(_isave)
		rs,imax = m.trelis(_isave,imax)
	#end process
		if (SYMPLOT and rs != m.lastrs):
			plt.text(j,1.5,rs)
		m.lastrs = rs
		y.append(z[j])			# green 			- signal
		y1.append(pselem[0])	# red				- dit 
		y2.append(pselem[1])	# blue 				- dah 
		y3.append(pselem[2])	# magenta(violet) 	- el-space
		y4.append(pselem[3])	# cyan(light blue)	- chr-space
		y5.append(rs/10.)	    # yellow			- wrd-space
		y6.append(khat)			# blue * 			- element state 
		x.append(j)
		
		#print percentage done 
		#sys.stdout.write("\r%4.2f %%" % (100.*float(j)/len(z)))
		#sys.stdout.flush()
		
	plt.plot(x,y,'g',x,y1,'r',x,y2,'b',x,y3,'m',x,y4,'c',x,y5,'y',x,y6,'b*')
	plt.show()	



def demodulate(x,Fs,freq):
    # demodulate audio signal with known CW frequency 
    t = np.arange(len(x))/ float(Fs)
    y =  x*((1 + np.sin(2*np.pi*freq*t))/2 )	

    #calculate envelope and low pass filter this demodulated signal
    #filter bandwidth impacts decoding accuracy significantly 
    #for high SNR signals 40 Hz is better, for low SNR 20Hz is better
    # 25Hz is a compromise - could this be made an adaptive value? 
    wn = 30./ (Fs/2.)  	# 30 Hz cut-off for lowpass  
    b, a = butter(3, wn)  	# 3rd order butterworth filter
    z = filtfilt(b, a, abs(y))
    #pass envelope magnitude to decoder 
    decode_stream(abs(z),Fs)

### MAIN PROGRAM ####

Fs, x = wavfile.read("test20db.wav") #test20db.wav cw002.wav
demodulate(x,Fs,600)
