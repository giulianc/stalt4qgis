#! /usr/bin/env python
# -*- coding: utf-8 -*-

# general libraries
import sys
import math as math

# QT modules
from PyQt4 import QtCore, QtGui

# MatPlotLib libraries
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
import matplotlib.pyplot as plt

prefix_icons = "/home/giuliano/0_Configurazione/qGis/plugins-dev/SphericalTrigonometry/"

"""
/***************************************************************************
spherTrig.py

A visual tool for Spherical Trigonometry Education
"""
vers = 'E-0.00'
build_date = '2014-06-06'
author = 'giuliano curti'
mail = 'giulianc51 at gmail dot com'
copyright = '2013-2018 giuliano curti'
license = 'GPL v2 (http://www.gnu.org/licenses/gpl-2.0.html)'
credits = ''
"""
 ***************************************************************************/

purtroppo non riesco ad usare in modo interattivo MatPlotLib
pertanto la procedura è strutturata in modulo, ognuno dei quali
illustra un particolare algoritmo di trigonometria sferica

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation version 2                                *
 *                                                                         *
 ***************************************************************************/

------------- to do --------------
dbAddCircle()
	usare ma fare attenzione che alcuni cerchi hanno 1 parametro ed altri due
dbGetPointRectCds())
	cambiare codice 'int' in 'int2c'
vettori
	attenzione: cambiata la chiave, ora è AB al posto di A,B
	aggiungere i tipi
		nor(mali)
			da punto a vettore
			a due vettori
cancellazione punto
	eliminare gli antipodali
valutare salvataggio/lettura di file
editing
	vertici
	vettori
	archi
	great circle
inquiry
	archi
	angoli
	aree?
controllare
	archi fra polari (risolto?)
	triangoli polari
normale da punto a vettore
	controllo che ci siano almeno 2(+1) punti
tangente da punto a vettore
	controllo che ci siano almeno 2(+1) punti
gestione cerchi
	arc e circle si possono fondere
	memorizzare la normale forse semplifica
	controllo che ci siano almeno 2(+1) punti
gestione archi
	controllo che ci siano almeno 2(+1) punti
todhunter
	Todhunter_22 convenzione lunghezza archi
	Todhunter 28 principle of duality (p.14)
	Todhunter_29a a+c>c (p.14)
	Todhunter_29b a>c-b (p.14)
	Todhunter_30 a+b+c<2pi (p.15)
	Todhunter_32 A+B+C>2pi (p.15)
	Todhunter_33 controllare
	Todhunter_34 similitude criteria (p.17)
	Todhunter_35 isoscele spherical triangles A (p.17)
	Todhunter_36 isoscele spherical triangles B (p.17)
	Todhunter_47 controllare
		verticale da punto a vettore (ci dovrebbe essere)
		risolvere sistema 3 incognite (ci dovrebbe essere)
	controllare perchè Todhunter e dotProduct non coincidono esattamente
"""

# ------------- utilities -----------

def about(mw):
	"""
		Visualizza info sulla procedura
	"""
	QtGui.QMessageBox.about(
		mw,
		'About',
		"Spherical Trigonometry Educational plugin for qGis"
+ "\n----------------------------"
+ "\nversion:\t" + vers
+ "\nbuild_date:\t" + build_date
+ "\nauthor:\t\t" + author
+ "\nemail:\t\t"+ mail
+ "\ncopyright:\t" + copyright
+ "\nlicense:\t" + license
	)

def info(mw):
	"""
		First advice to users
	"""
	msg = """
Spherical Trigonometry Educational plugin for qGis.
-----------------------------------------------
The procedure allows the learning of Spherical
Trigonometry as esplained by I.Todhunter - J.G.Leathem,
McMillan, London 1914;

this procedure is EXPERIMENTAL; it might
contains many bugs, few duplications and some
mistakes, only in part known to the author;

please inform me about any encountered problems.
	"""

	QtGui.QMessageBox.about(
		mw,
		'Info',
		msg + mail
	)

# ----- trigonometric utilities ----------

def sessad2rad(vIn):
	"""
		converte gradi sessadecimali dalla notazione corrente (a_giro) in radianti
	"""
	return (vIn/180)*math.pi

def rad2sessad(angRad):
	"""
		Converte i radianti nella notazione corrente (a_giro)
	"""
	a = angRad * 180 / math.pi
#	if a < 0:
#		a = a + 360
	return a

def sessad2sessag(a):
	"""
		trasforma da sessadecimali in sessagesimali
	"""
	# calcola gradi
	b = int(a)
	a = 60.*(a-b)
	res = str(b)+'d'
	#print b,res,a
	# calcola minuti
	b = int(a)
	a = 60.*(a-b)	
	res = res+str(b)+'\''
	#print b,res,a
	# calcola secondi
	res = res+'%5.3f"' % (a)
	#print res
	return res

def sessages2sessadecim(a):
	"""
		trasforma da sessagesimali in sessadecimali
	"""
	# preleva gradi
	pos = a.index('d')
	res = float(a[0:pos])
	a = a[pos+1:]
	# preleva minuti primi
	pos = a.index('\'')
	res = res + float(a[0:pos])/60
	a = a[pos+1:]
	# preleva minuti secondi
	pos = a.index('"')
	res = res + float(a[0:pos])/3600
	return res

# ----- spherical trigonometry ----------

def todhunter_73_1(a,b):
	"""
		Todhunter/Leathem, art. 73 (1)
		
		NB: questo angolo si può calcolare con
		il dotProduct() fra i vettori OA e Oy)
	"""
	return math.acos(math.cos(a)*math.cos(b))


# ----- vector&matrix functions ----------------

def nullMatrix(nr,nc):
	"""
		resituisce la matrice nulla di nr x nc
	"""
	mat = []
	for i in range(nr):
		tmp = []
		for j in range(nc):
			tmp.append(0.0)
		mat.append(tmp)
	return mat

def interpolationPoint(a,b,k):
	"""
		restituisce le cds del punto P tale che AP = k*AB
	"""
	xa,ya,za = a
	xb,yb,zb = b
	x,y,z = (1-k)*xa+k*xb,(1-k)*ya+k*yb,(1-k)*za+k*zb
#	print x,y,z
	return x,y,z

def dotProduct(u,v):
	ux,uy,uz = u
	vx,vy,vz = v
	return ux*vx + uy*vy + uz*vz

def crossProduct(u,v):
	ux,uy,uz = u
	vx,vy,vz = v
	return uy*vz - uz*vy,-ux*vz + uz*vx,ux*vy - uy*vx

def normalFromPtoAB(p,a,b):
	"""
		genera la normale da P al vettore AB
	"""
	xp,yp,zp = p
	xa,ya,za = a
	xb,yb,zb = b
	ux,uy,uz = xp-xa,yp-ya,zp-za
	vx,vy,vz = xb-xa,yb-ya,zb-za
	print 'vettore PA',ux,uy,uz
	print 'vettore AB',vx,vy,vz
	mu = math.sqrt(dotProduct([ux,uy,uz],[ux,uy,uz]))
	mv = math.sqrt(dotProduct([vx,vy,vz],[vx,vy,vz]))
	# piede H della normale da A al vettore OB
	if mu > 0:
		if mv > 0:
			k = dotProduct([ux,uy,uz],[vx,vy,vz])/(mu*mv)
			print 'k=',k
			x,y,z = xa+k*(mu/mv)*vx,ya+k*(mu/mv)*vy,za+k*(mu/mv)*vz
			print "piede H",x,y,z
			return x,y,z
		else:
			print 'vettore %s-%s nullo' % (a,b)
	return 0.,0.,0.

def normalToTwoVectors(u,v):
	"""
		normale ai vettori u,v

		condizione di normalitÃ :
			nx*ux+ny*uy+nz*uz=0
			nx*vx+ny*vy+nz*vz=0
		che diventa:
			|ux uy uz||nx| |0|
			|vx vy vz||ny|=|0|
			          |nz|
		risolvendo si ha, a meno del divisore ux*vy-vx*uy:
			|nx| |uy*vz-vy*uz|
			|ny|=|vx*uz-ux*vz| = crossProduct(u,v)
			|nz| |ux*vy-vx*uy|
	"""
	a,b,c = crossProduct(u,v)
#	print a,b,c
	m = math.sqrt(dotProduct([a,b,c],[a,b,c]))
	if m:
		return a/m,b/m,c/m
	else:
		return 0,0,0

def normaleDaPuntoAPiano(p,a,b):
	"""
		calcola il piede H della normale da P al piano a,b 
		p,a,b sono terne [x,y,z]
	"""
	# calcolo elementi matrice
	aa = dotProduct(a,a)
	bb = dotProduct(b,b)
	ab = dotProduct(a,b)
	pa = dotProduct(p,a)
	pb = dotProduct(p,b)
	print 'matrice'
	print aa,ab,'k',pa
	print ab,bb,'j',pb
	# determinante
	d = aa*bb-ab*ab
	if d:
		a1 = pa*bb-ab*pb
		b1 = aa*pb-pa*ab
		# soluzione
		k = a1/d
		j = b1/d
	else:
		k,j = 0.,0.
#	print 'k=',k,'j=',j
	return k,j	

def intersection3Dvectors(U,V):
	"""
		compute the intersection, if any, between
		the 3D vectors U = A + ku and V = B + jv
	"""
	A,u = U
	B,v = V
	xa,ya,za = A
	print 'A:',xa,ya,za
	xb,yb,zb = B
	print 'B:',xb,yb,zb
	x1,y1,z1 = u
	print 'u:',x1,y1,z1
	x2,y2,z2 = v
	print 'v:',x2,y2,z2
	r = xa-xb
	s = ya-yb
	t = za-zb
	print 'r:',r,'s:',s,'t:',t
	# parametri
	a,b,c = crossProduct(u,v)
	print 'cross product',a,b,c
	# condizione di esistenza
	if r*a+s*b+t*c == 0:
		k,j = 'na','na'
		if a != 0:
			j = (s*z1-t*y1)/a
			k = (s*z2-t*y2)/a
		elif b != 0:
			j = (r*z1-t*x1)/b
			k = (r*z2-t*x2)/b
		elif c != 0:
			j = (r*y1-s*x1)/c
			k = (r*y2-s*x2)/c
		else:
			return -1
		print 'k:',k,'j:',j
		if k != 'na':
			s1x,s1y,s1z = xa+k*x1,ya+k*y1,za+k*z1
			print 'soluzione 1',s1x,s1y,s1z
		if j != 'na':
			s2x,s2y,s2z = xb+j*x2,yb+j*y2,zb+j*z2
			print 'soluzione 2',s2x,s2y,s2z
		if k != 'na':
			return s1x,s1y,s1z
		elif j != 'na':
			return s2x,s2y,s2z
		else:
			return -1
	else:
		return -1

def matRotation3D(axys,ang):
	"""
		rotazione dell'angolo a [radians] rispetto a axys
	"""
	c = math.cos(ang)
	s = math.sin(ang)
	if axys in ('X','x'):
		return [
			[1.,0.,0.,0.],
			[0., c,-s,0.],
			[0., s, c,0.],
			[0.,0.,0.,1.]
		]
	if axys in ('Y','y'):		# forse occorre cambiare il segno del seno
		return [
			[ c,0.,-s,0.],
			[0.,1.,0.,0.],
			[ s,0., c,0.],
			[0.,0.,0.,1.]
		]
	else:
		return [
			[ c,-s,0.,0.],
			[ s, c,0.,0.],
			[0.,0.,1.,0.],
			[0.,0.,0.,1.]
		]

def trasformaPunti3D(x,y,z,mat):
	"""
		trasforma i x,y,z secondo la matrice mat
		NB: così non tiene conto della coordinata omogenea!
	"""
	num = len(x)
	for i in range(num):
#		print "il nodo %s (%7.3f %7.3f%7.3f)" % (x[i],y[i],z[i]),
		cds = []
		for [a,b,c,d] in mat:
			cds.append(a*x[i]+b*y[i]+c*z[i]+d)
		x[i],y[i],z[i] = cds[0],cds[1],cds[2]
#		print 'diventa (%7.3f %7.3f%7.3f)' % (x[i],y[i],z[i])
	return x,y,z

def matrixMultiplication(mat1,mat2):
	"""
		Calcola il prodotto delle matrici mat1*mat2
	"""
#	print "matrice 1:"
#	printMatrix(mat1)
#	print "matrice 2:"
#	printMatrix(mat2)
	nr1 = len(mat1)
	nc1 = len(mat1[0])
	nr2 = len(mat2)
	nc2 = len(mat2[0])
	mat = []
	if nc1 == nr2:
		for i in range(nr1):
			tmp = []
			for j in range (nc2):
				val = 0
				for k in range (nc1):
			  		val += mat1[i][k] * mat2[k][j]
				tmp.append(val)
			mat.append(tmp)
		return mat
	else:
		print "matrici non congruenti per la moltiplicazione"
		return -1

def matrixTranspose(mat):
	"""
		esegue la trasposta di una matrice
	"""
	nr = len(mat)
	nc = len(mat[0])
	matT = nullMatrix(nc,nr)
	for r in range(nr):
		for c in range(nc):
			matT[c][r] = mat[r][c]
	return matT

# ----- model ---------------------

def pointSphericCds(rad,lat,lon):
	"""
		genera un punto in coordinate sferiche
	"""
	x = rad*math.cos(lat)*math.cos(lon)
	y = rad*math.cos(lat)*math.sin(lon)
	z = rad*math.sin(lat)
# 	print x,y,z
	return x,y,z

def pointRectCds(x,y,z):
	"""
		genera un punto in coordinate sferiche
		se il punto è sulla verticale Z lon è indefinita

		-90 <= lat <= 90
		0 <= lon <= 360
	"""
#	print x,y,z,'diventa',
	h = math.sqrt(x**2+y**2)
	rad = math.sqrt(h**2+z**2)
	# calcola e controlla la latitudine
	lat = math.atan2(z,h)
#	print 'inizio',lat,
	if lat >0 and lat <= math.pi:
		if lat > math.pi/2:
			lat = lat-math.pi/2
	else:
		if lat > 0:
			lat = lat-2*math.pi
		if lat < -math.pi/2:
			lat = lat+math.pi/2
#	print 'infine',lat
	# calcola e controlla la longitudine
	if h:
		lon = math.atan2(y,x)
		if lon < 0.:
			lon = 2*math.pi+lon
	else:
		lon = 'n.a.'
#	print rad,rad2sessad(lat),rad2sessad(lon)
	return rad,lat,lon

def isOnSphere(rad,x,y,z):
	"""
		controlla se la curva x,y,z giace sulla sfera
	"""
	for i in range(len(x)):
		d = math.sqrt(x[i]**2+y[i]**2+z[i]**2) - 100.
		if d > 0.5:
			print i,x[i],y[i],z[i],'esterno alla sfera di',d

def circle(rad,num):
	"""
		genera un cerchio di raggio rad e num spicchi
	"""
	x = []
	y = []
	da = 2*math.pi/num
	a = 0.
	for i in range (num):
		x.append(rad*math.cos(a))
		y.append(rad*math.sin(a))
		a += da
	# chiude sul primo
	x.append(x[1])	# in prima posizione c'è il centro
	y.append(y[1])
	return x,y

def arco2D(rad,ang,num):
	"""
		genera un arco orizzontale di raggio rad, ampiezza ang e num spicchi
	"""
	x = []
	y = []
	da = ang/num
	a = 0.
	for i in range (num+1):
		x.append(rad*math.cos(a))
		y.append(rad*math.sin(a))
		a += da
	return x,y

def arco3D(rad,a,b,num):
	"""
		disegna arco di ampiezza AOB;
		calcola i vari angoli per portare la normale
		N sulla verticale Z e il primo punto A sull'asse X;
		quindi calcola l'ampiezza AOB dell'arco
		infine genera l'arco 3D;
	"""
	xa,ya,za = a
	xb,yb,zb = b
#	print 'A:',xa,ya,za
#	print 'B:',xb,yb,zb
#	calcolo normale a AOB
	xn,yn,zn = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
#	print 'N:',xn,yn,zn
	# raggio orizzontale OS
	xs,ys,zs = -yn,xn,0.
#	print 'raggio orizzontale OS',xs,ys,zs
	# --- prepara matrice di trasformazione ----
	# rotazione su Z per portare S sull'asse X
	a = math.atan2(ys,xs)
	mat = matRotation3D('Z',-a)
#	print 'rotazione Z:',-a,rad2sessad(-a),mat
	# rotazione su X per portare N su Z
	xy = math.sqrt(xn**2+yn**2)
	a = math.atan2(xy,zn)
	tmp = matRotation3D('X',-a)
	mat = matrixMultiplication(tmp,mat)
#	print 'rot. X',-a,rad2sessad(-a),mat
	# trasformazione
	x,y,z = [xa,xb,xn,xs],[ya,yb,yn,ys],[za,zb,zn,zs]
	x,y,z = trasformaPunti3D(x,y,z,mat)
	xa2,xb2,xn2,xs2 = x
	ya2,yb2,yn2,ys2 = y
	za2,zb2,zn2,zs2 = z
#	print 'A2',xa2,ya2,za2
#	print 'B2',xb2,yb2,zb2
	# angolo aFirst = AOS
	aFirst = math.atan2(ya2,xa2)
#	print 'aFirst:',aFirst,rad2sessad(aFirst)
	# rotazione su Z per portare A su X pari a AOS
	tmp = matRotation3D('Z',-aFirst)
	mat = matrixMultiplication(tmp,mat)
	# terza trasformazione
	x,y,z = [xa,xb,xn,xs],[ya,yb,yn,ys],[za,zb,zn,zs]
	x,y,z = trasformaPunti3D(x,y,z,mat)
	xa3,xb3,xn3,xs3 = x
	ya3,yb3,yn3,ys3 = y
	za3,zb3,zn3,zs3 = z
	# angolo AOB
	aob = math.atan2(yb3,xb3)
#	print 'ampiezza arco',yb3,xb3,aob,rad2sessad(aob)
	# genera arco
	x,y = arco2D(rad,aob,num)
	z = [0.]*len(x)
	# matrice inversa
	mat = matrixTranspose(mat)	# matrice ortogonale
#	print 'matrice inversa',mat
	# trasforma
	x,y,z = trasformaPunti3D(x,y,z,mat)
	return x,y,z

def circle3D(rad,a,n,num):	
	"""
		draw a circle of radius rad by A with normal N
	"""
	xa,ya,za = a
	xn,yn,zn = n
	# --- prepara matrice di trasformazione ----
	xy = math.sqrt(xn**2+yn**2)
	a = math.atan2(xy,zn)
	mat = matRotation3D('Y',-a)
#	print 'rot Y:',-a,rad2sessad(a),mat
	a = math.atan2(yn,xn)
	tmp = matRotation3D('Z',a)
	mat = matrixMultiplication(tmp,mat)	
#	print 'rot Z:',a,rad2sessad(a),mat
	# genera e trasforma cerchio
	x,y,z = parallelo(0.,rad,num)
	x,y,z = trasformaPunti3D(x,y,z,mat)
	return x,y,z

def parallelo(z,rad,num):
	"""
		parallelo di raggio rad e quota z
	"""
	p_x,p_y = circle(rad,num)
	p_z = [z] * len(p_x)
	return p_x,p_y,p_z

def meridiano(lon,rad,num):
	"""
		meridiano lon e raggio rad
		
		valutare l'adozione della routine circle() e della rotazione
	"""
	p_x = []
	p_y = []
	p_z = []
	c = math.cos(lon)
	s = math.sin(lon)
# 	print lon,c,s
	for i in range (num):
		a = 2*math.pi*i/num
		x0 = rad*math.cos(a)
		y0 = 0.0
		z0 = rad*math.sin(a)
		# ruota dell'angolo lon
		x = x0*c - y0*s
		y = x0*s + y0*c
# 		print x,y
		p_x.append(x)
		p_y.append(y)
		p_z.append(z0)
	# chiudere sul primo
	p_x.append(p_x[0])
	p_y.append(p_y[0])
	p_z.append(p_z[0])
	return p_x,p_y,p_z

def sphere(rad,num):
	"""
		genera il modello di sfera costituito da linestring che rappresentano:
		- paralleli
		- meridiani
		parametri:
		- rad	raggio
		- num	numero di intervalli (num_parall -1 = num_merid -1)
	"""
	parall = []
	meridian = []
	# paralleli
	da = 2*math.pi/num
	lat = -math.pi/2
	for i in range(num/2):	# in questo modo hanno lo stesso intervalo dei meridiani
		z0 = rad*math.sin(lat)
#		print lat,z0
		parall.append(parallelo(z0,rad*math.cos(lat),num))
		lat += da
	# meridiani
	lon = 0.
	for i in range(num):
		meridian.append(meridiano(lon,rad,num))
		lon += da
	# thicks
	thicks = []
	a = 0
	da = 360/num
	for i in range(num):
		thicks.append(a)
		a += da
	return [parall,meridian,thicks]

def arcPolar(rad,a,b):
	xa,ya,za = a
	xb,yb,zb = b
	xn,yn,zn = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
	xN,yN,zN = rad*xn,rad*yn,rad*zn
	print 'normale a AOB',xN,yN,zN
	return xN,yN,zN

# ====== class genericDlg =============

class genericDlg(QtGui.QDialog):
	""" Dialogo per N parametri """

	def __init__(self,title,fields):
		"""
			Inizializza la maschera per l'editing di N parametri.
			I nomi dei parametri sono dati dalla lista fields.
		"""
		QtGui.QDialog.__init__(self)
		# impostazione interfaccia utente
		self.setWindowTitle(title)
		num = len(fields)
		self.resize(280, 40+30*num)

		self.entries = []
		for i,l in enumerate(fields):
			label = QtGui.QLabel(self)
			label.setGeometry(QtCore.QRect(10, 10+30*i, 100, 20))
			label.setText(l)
			entry = QtGui.QLineEdit(self)
			entry.setGeometry(QtCore.QRect(120, 10+30*i, 150, 20))
			self.entries.append(entry)

		buttonBox = QtGui.QDialogButtonBox(self)
		buttonBox.setGeometry(QtCore.QRect(120, 40+30*i, 150, 30))
		buttonBox.setOrientation(QtCore.Qt.Horizontal)
		buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
		buttonBox.setObjectName("buttonBox")
		QtCore.QObject.connect(buttonBox,QtCore.SIGNAL("accepted()"),self.accept)
		QtCore.QObject.connect(buttonBox,QtCore.SIGNAL("rejected()"),self.reject)
		QtCore.QMetaObject.connectSlotsByName(self)

	def setValues(self,params):
		""" Inizializza i parametri nella maschera di input  """
		for i,v in enumerate(self.entries):
			v.setText(str(params[i]))

	def clearValues(self):
		""" Azzera i parametri nella maschera di input  """
		for i,v in enumerate(self.entries):
			v.setText('')

	def getValues(self):
		""" Restituisce il parametro  """
		val = []
		for v in self.entries:
			val.append(v.text())
		return val


# ====== class Main Window =============

class MainWindow(QtGui.QMainWindow):

	myRad = 100.
	myNum = 36		# così il passo è di 10deg

	sessagesimalFlag = False
	thickFlag = False

	cGray   = '#aaaaaa'
	cBlack  = '#000000'
	cRed    = '#ff0000'
	cBlue   = '#0000ff'
	cGreen  = '#00ff00'
	cYellow = '#ffff00'
	cViolet = '#ff00ff'

	mpl.rcParams['legend.fontsize'] = 10

	def __init__(self):

#		----------- mainwindow ---------------
		QtGui.QMainWindow.__init__(self)
		self.resize(500,20)
		self.setWindowTitle('Spherical Trigonometry')

#		----------- menu ---------------
		mbar = self.menuBar()
		mbFile = mbar.addMenu('File')
		mOpts = mbar.addMenu('Options')

		tmp = QtGui.QAction(QtGui.QIcon(''),'Redraw',self)        
		tmp.triggered.connect(self.redraw)
		mbar.addAction(tmp)

		mDb = mbar.addMenu('Database')
		mInsert = mbar.addMenu('Insert')
		mInqui = mbar.addMenu('Inquiry')
		mTodh = mbar.addMenu('Todhunter')
		mHelp = mbar.addMenu('Help')

#		----------- file menu ---------------
		tmp = QtGui.QAction(QtGui.QIcon(prefix_icons+"actionExit.png"), "Quit", self)
		tmp.setShortcut("Ctrl+Q")
		tmp.setStatusTip("Quit application")
		self.connect(tmp, QtCore.SIGNAL('triggered()'), self.finish)	#QtCore.SLOT('close()'))
		mbFile.addAction(tmp)

		# -------- help menu -----
		tmp = QtGui.QAction(QtGui.QIcon(''),'About',self)        
		tmp.triggered.connect(self.about)
		mHelp.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Info',self)        
		tmp.triggered.connect(self.info)
		mHelp.addAction(tmp)

		# -------- options menu -----
		self.menuAngle = QtGui.QAction(QtGui.QIcon(''),'dd.mmm',self)        
		self.menuAngle.triggered.connect(self.angleNotation)
		mOpts.addAction(self.menuAngle)

		self.menuThicks = QtGui.QAction(QtGui.QIcon(''),'Thicks on',self)        
		self.menuThicks.triggered.connect(self.thiksOnOff)
		mOpts.addAction(self.menuThicks)

		# -------- database menu -----
		tmp = QtGui.QAction(QtGui.QIcon(''),'List all',self)        
		tmp.triggered.connect(self.dbListAll)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Redraw',self)        
		tmp.triggered.connect(self.redraw)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete vertex',self)        
		tmp.triggered.connect(self.vertexDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete vector',self)        
		tmp.triggered.connect(self.vectorDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete parallel',self)        
		tmp.triggered.connect(self.parallelDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete meridian',self)        
		tmp.triggered.connect(self.meridianDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete great circle',self)        
		tmp.triggered.connect(self.gCircleDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete arc',self)        
		tmp.triggered.connect(self.arcDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete all arcs',self)        
		tmp.triggered.connect(self.arcDelAll)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete triangle',self)        
		tmp.triggered.connect(self.triangleDel)
		mDb.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Delete all',self)        
		tmp.triggered.connect(self.dbClearAll)
		mDb.addAction(tmp)

		# -------- insert menu -----
		tmp = QtGui.QAction(QtGui.QIcon(''),'Point in lat/lon',self)        
		tmp.triggered.connect(self.sphericalCds)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Point in x,y,z',self)        
		tmp.triggered.connect(self.rectangularCds)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Antipodal point',self)        
		tmp.triggered.connect(self.antipodalPoint)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Interpolation point',self)        
		tmp.triggered.connect(self.interpolationPoint)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Parallel by a point',self)        
		tmp.triggered.connect(self.parallelByPoint)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Meridian by a point',self)        
		tmp.triggered.connect(self.meridianByPoint)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'New vector',self)        
		tmp.triggered.connect(self.newVector)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Two vectors intersection',self)        
		tmp.triggered.connect(self.twoVectorIntersection)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Normal from point to vector',self)        
		tmp.triggered.connect(self.normalFromPtoV)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Tangent from point to vector',self)        
		tmp.triggered.connect(self.tangentFromPtoV)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Great Circle by 2 points',self)        
		tmp.triggered.connect(self.greatCircle2points)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Intersection of two great circle',self)        
		tmp.triggered.connect(self.greatCircleIntersection)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Geodetic between 2 points',self)        
		tmp.triggered.connect(self.geodetic2points)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Geodetic polar',self)        
		tmp.triggered.connect(self.geodeticPolar)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Antipodal geodetic',self)        
		tmp.triggered.connect(self.antipodalGeodetic)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Spherical triangle',self)        
		tmp.triggered.connect(self.sphericalTriangle)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Polar triangle',self)        
		tmp.triggered.connect(self.polarTriangle)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Antipodal triangle',self)        
		tmp.triggered.connect(self.antipodalTriangle)
		mInsert.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'test',self)        
		tmp.triggered.connect(self.prova)
		mInsert.addAction(tmp)

		# -------- inquiry menu -----
		tmp = QtGui.QAction(QtGui.QIcon(''),'Cartesian cds of a vertex',self)        
		tmp.triggered.connect(self.recCdsVertexInquiry)
		mInqui.addAction(tmp)

# fare inquiry cds geografiche

		tmp = QtGui.QAction(QtGui.QIcon(''),'Vector',self)        
		tmp.triggered.connect(self.vectorInquiry)
		mInqui.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'Angle',self)        
		tmp.triggered.connect(self.angleInquiry)
		mInqui.addAction(tmp)

		# -------- Todhunter menu -----
		tmp = QtGui.QAction(QtGui.QIcon(''),'Advice',self)        
		tmp.triggered.connect(self.todhunter_000_advice)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'001 sphere',self)        
		tmp.triggered.connect(self.todhunter_001_sphere)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'002a small circle',self)        
		tmp.triggered.connect(self.todhunter_002a_smallCircle)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'002b parallels',self)        
		tmp.triggered.connect(self.todhunter_002b_parallels)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'003a great circle',self)        
		tmp.triggered.connect(self.todhunter_003a_greatCircle)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'003b meridians',self)        
		tmp.triggered.connect(self.todhunter_003b_meridians)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'003c base meridian',self)        
		tmp.triggered.connect(self.todhunter_003c_baseMeridian)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'003d equator',self)        
		tmp.triggered.connect(self.todhunter_003d_equator)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'003e spherical coordinates',self)        
		tmp.triggered.connect(self.todhunter_003e_sphericalCds)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'004 arc',self)        
		tmp.triggered.connect(self.todhunter_004_arc)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'005 polar',self)        
		tmp.triggered.connect(self.todhunter_005_polar)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'006 polar theorem',self)        
		tmp.triggered.connect(self.todhunter_006_polarTheorem)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'007 quadrant',self)        
		tmp.triggered.connect(self.todhunter_007_quadrant)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'008 great circle inclination theorem',self)        
		tmp.triggered.connect(self.todhunter_008_gcInclinationTheorem)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'009 great circle inclination definition',self)        
		tmp.triggered.connect(self.todhunter_009_gcInclination)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'010_greatCircleBisection',self)        
		tmp.triggered.connect(self.todhunter_010_greatCircleBisection)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'014 arc comparison',self)        
		tmp.triggered.connect(self.todhunter_014_arcComparison)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'015 solid angle',self)        
		tmp.triggered.connect(self.todhunter_015_solidAngle)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'016 spherical triangle',self)        
		tmp.triggered.connect(self.todhunter_016_sphericalTriangle)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'017 sides & angles',self)        
		tmp.triggered.connect(self.todhunter_017_sidesAndAngles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'018 angles',self)        
		tmp.triggered.connect(self.todhunter_018_angles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'024a lune',self)        
		tmp.triggered.connect(self.todhunter_024a_lune)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'024b colunar triangles',self)        
		tmp.triggered.connect(self.todhunter_024b_colunarTriangles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'024c_antipodalTriangles',self)        
		tmp.triggered.connect(self.todhunter_024c_antipodalTriangles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'025_polarTriangles',self)        
		tmp.triggered.connect(self.todhunter_025_polarTriangles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'026 simmetry of polarity',self)        
		tmp.triggered.connect(self.todhunter_026_simmetryOfPolarity)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'027a supplemental triangles',self)        
		tmp.triggered.connect(self.todhunter_027a_supplementalTriangles)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'027b supplemental triangles 2',self)        
		tmp.triggered.connect(self.todhunter_027b_supplementalTriangles_2)
		mTodh.addAction(tmp)

		tmp = QtGui.QAction(QtGui.QIcon(''),'solid angle',self)        
		tmp.triggered.connect(self.euclid_XI_20_solidAngle)
		mTodh.addAction(tmp)









		tmp = QtGui.QAction(QtGui.QIcon(''),'Art. 47',self)        
		tmp.triggered.connect(self.todhunter_047)
		mTodh.addAction(tmp)

		# initialization
		self.dbClearAll()

	def about(self):
		about(self)

	def info(self):
		info(self)

	def finish(self):
		res = QtGui.QMessageBox.question(
			self,
			'quit application',
			'vuoi chiudere la sessione?',
			QtGui.QMessageBox.Yes,
			QtGui.QMessageBox.No
		)
		if res == QtGui.QMessageBox.Yes:
			sys.exit()

#	---- utilità -------------

	def angleNotation(self):
		if self.sessagesimalFlag:
			self.menuAngle.setText('dd.mmm')
			self.sessagesimalFlag = False
		else:
			self.menuAngle.setText('dd mm ss.sss')
			self.sessagesimalFlag = True

	def thiksOnOff(self):
		if self.thickFlag:
			self.thickFlag = False
			self.menuThicks.setText('Thicks on')
		else:
			self.thickFlag = True
			self.menuThicks.setText('Thicks off')

#	============= struttura dati ================

#	----------- funzioni generali ---------------

	def dbClearAll(self):
		self.myTitle   = ''
		self.nextCod   = 65
		self.points    = {}
		self.points    = {}	# punti as
									# lab -> [type,par] depending on type:
									# geo: rad,lat,lon
									# rec: x,y,z
									# int: gc1,gc2
									# int2v: v1,v2
									# itp: p1,p2,k
									# pol:
									#		sc,a,n
									#		gc,a,b
		self.points['O']     = ['rec',[0.,0.,0.]]
		self.points['Nord']  = ['rec',[0.,0.,self.myRad]]
		self.points['Sud']   = ['rec',[0.,0.,-self.myRad]]
		self.vectors   = {}	# tipi:
									# 'vec'
									# 'tan'
		self.circles   = {}	# tipi:
									# 'gc'
									# 'sc'
									# 'par'
									# 'mer'
		self.arcs      = {}	# tipi:
									# 'arc': [a,b]
									# 'scarc': [a,b,c]

	def dbListAll(self):
		print '----------- point list ---------'
		print 'label    parameters'
		print '--------------------------------'
		for k in self.points.keys():
			print k,self.points[k]
		print '--------------------------------'
		print '----------- vector list ---------'
		print 'start    end'
		print '--------------------------------'
		for k in self.vectors.keys():
			print k,self.vectors[k]
		print '--------------------------------'
		print '------ circles list ---------'
		print 'label    parameters'
		print '--------------------------------'
		for k in self.circles.keys():
			print k,self.circles[k]
		print '--------------------------------'
		print '----------- arc list ---------'
		print 'label   start    end'
		print '--------------------------------'
		for k in self.arcs.keys():
			print k,self.arcs[k]
		print '--------------------------------'

#	-------- gestione del codice per le label -----

	def dbGetNextCod(self):
		"""
			gestisce il codice che genera la label del punto
		"""
		if self.nextCod == 79:
			self.nextCod += 1
		p = str(unichr(self.nextCod))
		# aggiorna il codice
		self.nextCod += 1
		return p

#	---------- punti --------------

	def dbAddPointLatLon(self,lat,lon):
		"""
			aggiunge al DB un punto in cds sferiche (radianti)
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['geo',[self.myRad,lat,lon]]
		return p

	def dbAddPointRect(self,x,y,z):
		"""
			aggiunge al DB un punto in cds rettangolari
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['rec',[x,y,z]]
		return p

	def dbAddAntipodalPoint(self,a):
		"""
			insert a point antipodal of the point A
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['ant',a]
		return p

	def dbAddInterpolation(self,p1,p2,k):
		"""
			add an interpolation point between p1,p2 at k from p1
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['itp',[p1,p2,k]]
		return p

	def dbAddPolar(self,cod,a,b):
		"""
			add a polar point to circle ref
			cod = 'sc' or 'gc'
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['pol',cod+a+b]
		return p

	def dbAddIntersection(self,gc1,gc2):
		"""
			add an intersection of two great circles
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['int',[gc1,gc2]]
		return p

	def dbAddIntersection2Vec(self,v1,v2):
		"""
			add an intersection of two vectors
		"""
		# preleva il codice disponibile
		p = self.dbGetNextCod()
		# salva i dati
		self.points[p] = ['int2v',[v1,v2]]
		return p

	def dbGetPointSpherCds(self,p):
		"""
		 restituisce le coordinate sferiche del punto p
		"""
		if p in self.points.keys():
#			print p,'esiste',
			type,par = self.points[p]
			if type == 'geo':
#				print 'in cds sferiche'
				return par
			elif type == 'rec':
#				print 'in cds cartesiane'
				return pointRectCds(par)
		else:
			print 'NB:',p,'NON esiste in archivio'
			return -1

	def dbGetPointRectCds(self,p):
		"""
			restituisce le coordinate rettangolari del punto p
			la funzione è ricorsiva (type == 'ant')
		"""
		if p in self.points.keys():
#			print p,'esiste',
			type,par = self.points[p]
			if type == 'rec':
#				print 'in cds cartesiane'
				return par
			elif type == 'geo':
#				print 'in cds geografiche'
				return pointSphericCds(*par)
			elif type == 'ant':
#				print 'punto antipodale di',par
				tmp = self.dbGetPointRectCds(par)	# è ricorsiva !
				if tmp != -1:
					x,y,z = tmp
					return -x,-y,-z
			elif type == 'itp':
#				print 'punto interpolato fra',par
				a,b,k = par
#				print a,b,k
				tmp = self.dbGetPointRectCds(a)	# è ricorsiva !
				if tmp != -1:
					xa,ya,za = tmp
					tmp = self.dbGetPointRectCds(b)	# è ricorsiva !
					if tmp != -1:
						xb,yb,zb = tmp
						x,y,z = (1-k)*xa+k*xb,(1-k)*ya+k*yb,(1-k)*za+k*zb
#						print x,y,z
						return x,y,z
			elif type == 'pol':
				type,a,b = self.dbGetPolar(p)
#				print 'ricevuto',type,a,b
				tmp = self.dbGetPointRectCds(a)
				if tmp != -1:
					xa,ya,za = tmp
					tmp = self.dbGetPointRectCds(b)
					if tmp != -1:
						xb,yb,zb = tmp
						if type == 'gc':
							# calcolo normale a AOB
							xn,yn,zn = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
#							print 'tipo _gc_: normale:',xn,yn,zn
						else:
							xn,yn,zn = xb,yb,zb
#							print 'tipo _sc_: normale:',xn,yn,zn
						# il polo è l'intersezione della normale con la sfera
						mn = math.sqrt(dotProduct([xn,yn,zn],[xn,yn,zn]))
						k = self.myRad/mn
						xn,yn,zn = k*xn,k*yn,k*zn
						return xn,yn,zn
					else:
						print 'polare: punto',b,'non in archivio'
				else:
					print 'polare: punto',a,'non in archivio'
			elif type == 'int':
#				print '====punto',p,type,par,'======'
				c1,c2 = par
				# prende il primo great circle
				tmp = self.circles[c1]
				if tmp[0] == 'gc':
					a,b = tmp[1],tmp[2]
					tmp = self.dbGetPointRectCds(a)
					if tmp != -1:
						xa,ya,za = tmp
						tmp = self.dbGetPointRectCds(b)
						if tmp != -1:
							xb,yb,zb = tmp
							# calcolo normale a AOB (non è aggiunta al Database)
							xn1,yn1,zn1 = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
#							print 'polo 1',xn1,yn1,zn1
							# prende il secondo great circle
							tmp = self.circles[c2]
							if tmp[0] == 'gc':
								a,b = tmp[1],tmp[2]
								tmp = self.dbGetPointRectCds(a)
								if tmp != -1:
									xa,ya,za = tmp
									tmp = self.dbGetPointRectCds(b)
									if tmp != -1:
										xb,yb,zb = tmp
										# calcolo normale a AOB (non è aggiunta al Database)
										xn2,yn2,zn2 = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
#										print 'polo 2',xn2,yn2,zn2
										# normale alle due normali
										xn,yn,zn = normalToTwoVectors([xn1,yn1,zn1],[xn2,yn2,zn2])
										mn = math.sqrt(dotProduct([xn,yn,zn],[xn,yn,zn]))
										k = self.myRad/mn
										xn,yn,zn = k*xn,k*yn,k*zn
#										print 'intersezione:',xn,yn,zn
										return xn,yn,zn
									else:
										print 'intersezione: punto',b,'non in archivio'
								else:
									print 'intersezione: punto',a,'non in archivio'
							else:
								print 'intersezione: secondo cerchio non valido'
						else:
							print 'polare: punto',b,'non in archivio'
					else:
						print 'polare: punto',a,'non in archivio'
				else:
					print 'intersezione: primo cerchio non valido'
			elif type == 'int2v':
#				print '====punto',p,type,par,'======'
				u,v = par
				# prende il primo vettore u = AB
				a,b = u[0],u[1]
				# prende il secondo vettore v = CD
				c,d = v[0],v[1]
				tmp = self.dbGetPointRectCds(a)
				if tmp != -1:
					xa,ya,za = tmp
					tmp = self.dbGetPointRectCds(b)
					if tmp != -1:
						xb,yb,zb = tmp
						tmp = self.dbGetPointRectCds(c)
						if tmp != -1:
							xc,yc,zc = tmp
							tmp = self.dbGetPointRectCds(d)
							if tmp != -1:
								xd,yd,zd = tmp
								# prepara i vettori
								u = [xb-xa,yb-ya,zb-za]
								v = [xd-xc,yd-yc,zd-zc]
								U = [[xa,ya,za],u]
								V = [[xc,yc,zc],v]
								# calcola intersezione
								tmp = intersection3Dvectors(U,V)
								if tmp != -1:
									x,y,z = tmp
									print 'intersezione:',x,y,z
									return x,y,z
			else:
				print 'tipo di punto',type,'non riconosciuto'
		else:
			print 'NB:',p,'NON esiste in archivio'
			return -1

	def dbGetPolar(self,p):
		"""
		 restituisce i parametri della polare
		"""
		if p in self.points.keys():
#			print p,'esiste',
			type,par = self.points[p]
			if type == 'pol':
				ctype,a,b = par[0:2],par[2],par[3]
#				print 'ok, polare del cerchio:',ctype,'per i punti',a,b
				return ctype,a,b
		else:
			print 'NB:',p,'NON esiste in archivio'
			return -1

#	------- vettori -------------

	def vectorExists(self,v):
		"""
			controlla l'esistenza del vettore v
		"""
		if v in self.vectors.keys():
			return 1
		else:
			return 0

	def dbAddVector(self,p1,p2):
		"""
			aggiunge un vettore
		"""
		self.vectors[p1+p2] = ['vec',[p1,p2]]

	def dbAddTangent(self,p,v1,v2):
		"""
			aggiunge una tangente ai vettori
		"""
		self.vectors[p+v1+v2] = ['tan',[p,v1,v2]]

	def dbGetVector(self,ref):
		"""
			fornisce le cds rettangolari del vettore ref
		"""
		if ref in self.vectors.keys():
			typ,par = self.vectors[ref]
			if typ == 'vec':
				a,b = par	
				tmp = self.dbGetPointRectCds(a)
				if tmp != -1:
					xa,ya,za = tmp
#					print 'punto',a,xa,ya,za
					tmp = self.dbGetPointRectCds(b)
					if tmp != -1:
						xb,yb,zb = tmp
#						print 'punto',b,xb,yb,zb
						return [xa,ya,za],[xb,yb,zb]
			elif typ == 'tan':
				p,a,b = par
				tmp = self.dbGetPointRectCds(p)
				if tmp != -1:
					xp,yp,zp = tmp
#					print 'punto',p,xp,yp,zp
					mp = math.sqrt(dotProduct([xp,yp,zp],[xp,yp,zp]))
					tmp = self.dbGetPointRectCds(a)
					if tmp != -1:
						xa,ya,za = tmp
#						print 'punto',a,xa,ya,za	
						tmp = self.dbGetPointRectCds(b)
						if tmp != -1:
							xb,yb,zb = tmp
#							print 'punto',b,xb,yb,zb
							xv,yv,zv = xb-xa,yb-ya,zb-za
							# piede della tangente in P sul vettore AB
							k = mp**2/dotProduct([xp,yp,zp],[xv,yv,zv])
#							print mp**2,dotProduct([xp,yp,zp],[xv,yv,zv]),k
							xh,yh,zh = xa+k*(xb-xa),ya+k*(yb-ya),za+k*(zb-za)
#							print 'H',xh,yh,zh
							return [xp,yp,zp],[xh,yh,zh]
			print 'vettore',ref,'inesistente'
			return -1

#	---------- archi --------------

	def dbAddArc(self,a,b):
		"""
			aggiunge un arco (great circle)
		"""
		self.arcs[str(a)+str(b)] = ['arc',[str(a),str(b)]]

	def dbAddSmallCircleArc(self,a,b,c):
		"""
			aggiunge un arco (small circle)
		"""
		self.arcs[str(a)+str(b)] = ['scarc',[str(a),str(b),str(c)]]

	def dbGetArc(self,ref):
		"""
			restituisce l'arco ref
		"""
		return self.arcs[ref]

#	------- cerchi ----------------

	def dbAddCircle(self,typ,a,b):
		"""
			aggiunge un great|small circle
		"""
		self.circles[typ+str(a)+str(b)] = [typ,a,b]

#	----- funzioni di disegno ----------

	def drawBase(self,ax):
		"""
			disegna la sfera di base
		"""
		# sfera di base
		parall,meridian,thicks = sphere(self.myRad,self.myNum)
		# paralleli
		for tmp in (parall):
			x,y,z = tmp
			ax.plot(x,y,z,color=self.cGray,linewidth=0.5)
		# equator
		x,y,z = parallelo(0.,self.myRad,self.myNum)
		ax.plot(x,y,z,color=self.cBlack,linewidth=0.5)
		# equator thicks
		if self.thickFlag:
			x,y,z = parall[int(self.myNum/4)]
			for i,t in enumerate(thicks):
				ax.text(x[i],y[i],z[i],t)
		# meridian
		for tmp in (meridian):
			x,y,z = tmp
			ax.plot(x,y,z,color=self.cGray,linewidth=0.5)
		# base meridian
		x,y,z = meridiano(0.,self.myRad,self.myNum)
		ax.plot(x,y,z,color=self.cBlack,linewidth=0.5)
		# meridian thicks
		if self.thickFlag:
			x,y,z = meridian[0]
			for i,t in enumerate(thicks):
				ax.text(x[i],y[i],z[i],t)
		# axys labels
		plt.xlabel('X')
		plt.ylabel('Y')

	def drawPoint(self,ax,a):
		"""
			draw the point a
		"""
		tmp = self.dbGetPointRectCds(a)
		if tmp != -1:
			x,y,z = tmp
			ax.plot([0.,x],[0.,y],[0.,z],'r-',color=self.cRed)
			ax.text(x,y,z,a)

	def drawVector(self,ax,ref):
		"""
			draw the vector ref
		"""
		tmp = self.dbGetVector(ref)
		if tmp != -1:
			[xa,ya,za],[xb,yb,zb] = tmp
			ax.plot([xa,xb],[ya,yb],[za,zb],'r-',color=self.cViolet,linewidth=1.5)

	def drawGreatCircle(self,ax,a,b):
		"""
			great circle
		"""
		tmp = self.dbGetPointRectCds(a)
		if tmp != -1:
			xa,ya,za = tmp
			tmp = self.dbGetPointRectCds(b)
			if tmp != -1:
				xb,yb,zb = tmp
				# calcolo normale a AOB (non è aggiunta al Database)
				xn,yn,zn = normalToTwoVectors([xa,ya,za],[xb,yb,zb])
#				print 'normale a AOB',xn,yn,zn
				x,y,z = circle3D(self.myRad,[xa,ya,za],[xn,yn,zn],self.myNum)
				# visualizzazione cerchio
				ax.plot(x,y,z,color=self.cBlue,linewidth=1.5)

	def drawSmallCircle(self,ax,a,n):
		"""
			small circle through A normal to N
		"""
		tmp = self.dbGetPointRectCds(a)
		if tmp != -1:
			xa,ya,za = tmp
			tmp = self.dbGetPointRectCds(n)
			if tmp != -1:
				xn,yn,zn = tmp
				# distanza del centro dall'origine
				nx = math.sqrt(dotProduct([xn,yn,zn],[xn,yn,zn]))
				k = (xa*xn+ya*yn+za*zn)/nx
#				print 'k=',k
				# centro
				xc,yc,zc = k*xn/nx,k*yn/nx,k*zn/nx
# non serve salvarlo				c = self.dbAddPointRect(xc,yc,zc)
				ax.plot([0.,xc],[0.,yc],[0.,zc],color=self.cRed,linewidth=1.0)
				ax.plot([xa,xc],[ya,yc],[za,zc],color=self.cRed,linewidth=1.0)
				ax.text(xc,yc,zc,'C')
#				print 'centro:',xc,yc,zc
				# raggio
				r = math.sqrt(self.myRad**2-k**2)	# math.sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
#				print 'raggio=',r
				x,y,z = circle3D(r,[xa,ya,za],[xn,yn,zn],self.myNum)
				# traslazione in C
				for i in range(len(x)):
					x[i] += xc
					y[i] += yc
					z[i] += zc
				# visualizzazione cerchio
				ax.plot(x,y,z,color=self.cBlue,linewidth=1.0)

	def drawArc(self,ax,a,b):
		"""
			disegna l'arco AB (great circle)
		"""
		tmp = self.dbGetPointRectCds(a)
#		print 'ricevo',tmp
		if tmp != -1:
			xa,ya,za = tmp
			tmp = self.dbGetPointRectCds(b)
#			print 'ricevo',tmp
			if tmp != -1:
				xb,yb,zb = tmp
				x,y,z = arco3D(self.myRad,[xa,ya,za],[xb,yb,zb],self.myNum)
				ax.plot(x,y,z,color=self.cGreen,linewidth=2.0)

	def drawSmallCircleArc(self,ax,a,b,c):
		"""
			disegna l'arco AB (small circle)
		"""
		tmp = self.dbGetPointRectCds(a)
#		print 'ricevo',tmp
		if tmp != -1:
			xa,ya,za = tmp
			tmp = self.dbGetPointRectCds(b)
#			print 'ricevo',tmp
			if tmp != -1:
				xb,yb,zb = tmp
				tmp = self.dbGetPointRectCds(c)
#				print 'ricevo',tmp
				if tmp != -1:
					xc,yc,zc = tmp
					# sposta il centro C nell'origine
					xa1,ya1,za1 = xa-xc,ya-yc,za-zc
					xb1,yb1,zb1 = xb-xc,yb-yc,zb-zc
					ax.text(xa1,ya1,za1,'A1')
					ax.text(xb1,yb1,zb1,'B1')
					rad = math.sqrt(dotProduct([xa1,ya1,za1],[xa1,ya1,za1]))
					print 'rad',rad
					x,y,z = arco3D(rad,[xa1,ya1,za1],[xb1,yb1,zb1],self.myNum)
					# sposta il centro C nella posizione iniziale
					for i in range(len(x)):
						x[i] += xc
						y[i] += yc
						z[i] += zc
					# disegna
					ax.plot(x,y,z,color=self.cGreen,linewidth=2.0)

	def redraw(self):
		"""
			ridisegna il modello
		"""
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		# titolo
		plt.title(self.myTitle)
		# sfera di base
		self.drawBase(ax)
		# disegna i punti in archivio
		for k in self.points.keys():
			self.drawPoint(ax,k)
		# disegna i vettori
		for k in self.vectors.keys():
			self.drawVector(ax,k)
		# disegna circle
		for k in self.circles.keys():
			tmp = self.circles[k]
			if tmp[0] == 'gc':
				self.drawGreatCircle(ax,tmp[1],tmp[2])
			elif tmp[0] == 'sc':
				self.drawSmallCircle(ax,tmp[1],tmp[2])
			else:
				tmp2 = self.dbGetPointSpherCds(tmp[1])
				if tmp2 != -1:
					rad,lat,lon = tmp2
					if tmp[0] == 'par':
						z0 = rad*math.sin(lat)
						r = rad*math.cos(lat)
						x,y,z = parallelo(z0,r,self.myNum)
					elif tmp[0] == 'mer':
						x,y,z = meridiano(lon,rad,self.myNum)
				# visualizzazione cerchio
				ax.plot(x,y,z,color=self.cBlue,linewidth=1.0)
		# disegna archi
		for k in self.arcs.keys():
			a,b = self.dbGetArc(k)
			if a == 'arc':
				a,b = b
				self.drawArc(ax,a,b)
			elif a == 'scarc':
				a,b,c = b
				self.drawSmallCircleArc(ax,a,b,c)
		# visualizzazione modello
		plt.show()

#	--------- funzioni di inserimento ------------

	def sphericalCds(self):
		"""
			dialogo per l'immissione e disegno di un punto in cds sferiche
		"""
		dlg = genericDlg('spherical coordinates',['lat','lon'])
		if self.sessagesimalFlag:
			dlg.setValues(['xxd mm\' ss.sss"','xxd mm\' ss.sss"'])
		else:
			dlg.setValues([0.,0.])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				lat,lon = tmp
				if self.sessagesimalFlag:
					lat,lon = str(lat),str(lon)
					lat = sessages2sessadecim(lat)
					lon = sessages2sessadecim(lon)
				else:
					lat = float(lat)
					lon = float(lon)
#				print lat,lon
				lat = sessad2rad(lat)
				lon = sessad2rad(lon)
#				print lat,lon
				l = self.dbAddPointLatLon(lat,lon)
			else:
				print "errore: numero di parametri inadeguato"

	def rectangularCds(self):
		"""
			dialogo per l'immissione e disegno di un punto in cds rettangolari
		"""
		dlg = genericDlg('rectangular coordinates',['x','y','z'])
		dlg.setValues([0.,0.,0.])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				x,y,z = tmp
				x = float(x)
				y = float(y)
				z = float(z)
#				print x,y,z
				l = self.dbAddPointRect(x,y,z)
			else:
				print "errore: numero di parametri inadeguato"

	def antipodalPoint(self):
		dlg = genericDlg('Point',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				p = str(tmp.pop())
#				print p,type(p)
				l = self.dbAddAntipodalPoint(p)
			else:
				print "errore: numero di parametri inadeguato"

	def interpolationPoint(self):
		"""
			inserisce un punto interpolato linearmente fra P1 e P2
		"""
		dlg = genericDlg('Interpolated point',['P1','P2','k'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,k = tmp
				a,b,k = str(a),str(b),float(k)
#				print a,b,k
				l = self.dbAddInterpolation(a,b,k)
			else:
				print "errore: numero di parametri inadeguato"

	def newVector(self):
		"""
			dialogo per l'inserimento di un vettore
		"""
		dlg = genericDlg('two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				a = str(a)
				b = str(b)
				if a in self.points.keys():
					if b in self.points.keys():
						self.dbAddVector(a,b)
					else:
						print 'vertice',b,'inesistente'
				else:
					print 'vertice',a,'inesistente'
			else:
				print "errore: numero di parametri inadeguato"

	def twoVectorIntersection(self):
		"""
			gestisce l'intersezione, se esiste, di 2 vettori 3D
		"""
		dlg = genericDlg('Two vectors',['v1','v2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				u,v = tmp
				u,v = str(u),str(v)
				# controlla i vettori
				if u != v:
					if self.vectorExists(u):
						if self.vectorExists(v):
							# salva l'intersezione
							self.dbAddIntersection2Vec(u,v)
						else:
							print 'il vettore',v,'non esiste in archivio'
					else:
						print 'il vettore',u,'non esiste in archivio'
				else:
					print 'i vettori non possono coincidere'
			else:
				print "errore: numero di parametri inadeguato"

	def parallelByPoint(self):
		dlg = genericDlg('Point',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				p = str(tmp.pop())
#				print p,type(p)
				self.circles['par'+str(p)] = ['par',str(p)]
#				self.redraw()
			else:
				print "errore: numero di parametri inadeguato"

	def meridianByPoint(self):
		dlg = genericDlg('Point',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				p = str(tmp.pop())
#				print p,type(p)
				self.circles['mer'+str(p)] = ['mer',str(p)]
#				self.redraw()
			else:
				print "errore: numero di parametri inadeguato"

	def normalFromPtoV(self):
		"""
			dialogo per la normale dal punto P al vettore P1-P2
		"""
		dlg = genericDlg('normal from P to v',['P','vect (P1,P2)'])
		dlg.setValues(['',','])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				p,v = tmp
				p = str(p)	# altrimenti diventa una QString
				a,b = v.split(',')
				a,b = str(a),str(b)	# altrimenti diventano QString
				if p in self.points.keys():
					if a+b in self.vectors.keys():
						# si omette la verifica sui vertici a,b perchè si presume esistenti
						xp,yp,zp = self.points[p][1]
						xa,ya,za = self.points[a][1]
						xb,yb,zb = self.points[b][1]
#						print 'P',xp,yp,zp
#						print 'A',xa,ya,za
#						print 'B',xb,yb,zb
						x,y,z = normalFromPtoAB([xp,yp,zp],[xa,ya,za],[xb,yb,zb])
#						print "normale",x,y,z
						l = self.dbAddPointRect(x,y,z)
					else:
						print 'vettore',a+b,'inesistente'
				else:
					print 'vertice',p,'inesistente'

	def tangentFromPtoV(self):
		dlg = genericDlg('normal from P to v',['P','vect (P1,P2)'])
		dlg.setValues(['',','])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				p,v = tmp
				p = str(p)	# altrimenti diventa una QString
				if p in self.points.keys():
					a,b = v.split(',')
					a,b = str(a),str(b)	# altrimenti diventano QString
					if a in self.points.keys():
						if b in self.points.keys():
							self.dbAddTangent(p,a,b)
						else:
							print 'errore: punto',b,'non in archivio'
					else:
						print 'errore: punto',a,'non in archivio'
				else:
					print 'errore: punto',p,'non in archivio'
			else:
				print "errore: numero di parametri inadeguato"

	def greatCircle2points(self):
		"""
			dialogo per il calcolo del great circle passante
			da due punti in cds sferiche
		"""
		dlg = genericDlg('two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
#				self.redraw()
			else:
				print "errore: numero di parametri inadeguato"

	def greatCircleIntersection(self):
		"""
			dialogo per il calcolo dell'intersezione
			fra due great circles
		"""
		dlg = genericDlg('two great circles',['AB','CD'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				c1 = str('gc'+a)
				c2 = str('gc'+b)
				l = self.dbAddIntersection(c1,c2)
				m = self.dbAddAntipodalPoint(l)
			else:
				print "errore: numero di parametri inadeguato"

	def geodetic2points(self):
		"""
			dialogo per il calcolo della geodetica fra
			due punti in cds sferiche
		"""
		dlg = genericDlg('Two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				self.dbAddArc(str(a),str(b))
			else:
				print "errore: numero di parametri inadeguato"

	def geodeticPolar(self):
		"""
			dialogo per il calcolo del polo della geodetica per 2 punti
		"""
		dlg = genericDlg('Two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				a = str(a)
				b = str(b)
				self.dbAddPolar('gc',a,b)
			else:
				print "errore: numero di parametri inadeguato"

	def antipodalGeodetic(self):
		"""
			dialogo per il calcolo della geodetica antipodale
			(questa non mi sembra definita da Todhunter, però
			risulta facilmente derivabile dal suo lavoro)
		"""
		dlg = genericDlg('Arc vertices',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				a,b = str(a),str(b)
				xa,ya,za = self.points[a][1]
				xb,yb,zb = self.points[b][1]
				a1 = self.dbAddPointRect(-xa,-ya,-za)
				b1 = self.dbAddPointRect(-xb,-yb,-zb)
				self.dbAddArc(str(a),str(b))
				self.redraw()
			else:
				print "errore: numero di parametri inadeguato"

	def sphericalTriangle(self):
		"""
			dialogo per il calcolo di un triangolo sferico
		"""
		dlg = genericDlg('Three points',['P1','P2','P3'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,c = tmp
#				print a,b,c
				if a in self.points.keys():
					if b in self.points.keys():
						if c in self.points.keys():
							self.dbAddArc(str(a),str(b))
							self.dbAddArc(str(b),str(c))
							self.dbAddArc(str(a),str(c))
						else:
							print 'il punto',c,'non esiste in archivio'
					else:
						print 'il punto',b,'non esiste in archivio'
				else:
					print 'il punto',a,'non esiste in archivio'
			else:
				print "errore: numero di parametri inadeguato"

	def polarTriangle(self):
		"""
			dialogo per l'inserimento del triangolo polare
		"""
		dlg = genericDlg('Polar triangle of',['P1','P2','P3'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,c = tmp
				a,b,c = str(a),str(b),str(c)
#				print a,b,c
				if a in self.points.keys():
					if b in self.points.keys():
						if c in self.points.keys():
							xa,ya,za = self.points[a][1]
							xb,yb,zb = self.points[b][1]
							xc,yc,zc = self.points[c][1]
							x1c,y1c,z1c = arcPolar(self.myRad,[xa,ya,za],[xb,yb,zb])
							x1a,y1a,z1a = arcPolar(self.myRad,[xb,yb,zb],[xc,yc,zc])
							x1b,y1b,z1b = arcPolar(self.myRad,[xc,yc,zc],[xa,ya,za])
							a1 = self.dbAddPointRect(x1a,y1a,z1a)
							b1 = self.dbAddPointRect(x1b,y1b,z1b)
							c1 = self.dbAddPointRect(x1c,y1c,z1c)
							self.dbAddArc(str(a1),str(b1))
							self.dbAddArc(str(b1),str(c1))
							self.dbAddArc(str(c1),str(a1))
						else:
							print 'il punto',c,'non esiste in archivio'
					else:
						print 'il punto',b,'non esiste in archivio'
				else:
					print 'il punto',a,'non esiste in archivio'
			else:
				print "errore: numero di parametri inadeguato"

	def antipodalTriangle(self):
		"""
			dialogo per l'inserimento del triangolo antipodale
		"""
		dlg = genericDlg('Polar triangle of',['P1','P2','P3'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,c = tmp
				a,b,c = str(a),str(b),str(c)
#				print a,b,c
				xa,ya,za = self.points[a][1]
				xb,yb,zb = self.points[b][1]
				xc,yc,zc = self.points[c][1]
				a1 = self.dbAddPointRect(-xa,-ya,-za)
				b1 = self.dbAddPointRect(-xb,-yb,-zb)
				c1 = self.dbAddPointRect(-xc,-yc,-zc)
				self.dbAddArc(str(a1),str(b1))
				self.dbAddArc(str(b1),str(c1))
				self.dbAddArc(str(c1),str(a1))
			else:
				print "errore: numero di parametri inadeguato"

#	-------- inquiry functions -----------

	def recCdsVertexInquiry(self):
		"""
			dialogo per l'inquiry delle coordinate cartesiane di un vertice
		"""
		dlg = genericDlg('A vertex',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				a = tmp.pop()
				a = str(a)
				"""
				if a in self.points.keys():
					lat,lon = self.points[a][0]
					x,y,z = self.points[a][1]
				"""
				tmp = self.dbGetPointRectCds(a)
				if tmp != -1:
					x,y,z = tmp
					print '------------vertex----------------'
					print 'vertex',a
					print 'rectangular coordinates %7.3f %7.3f %7.3f' % (x,y,z)
				else:
					print 'vertice',a,'inesistente'
			else:
				print 'numero di parametri inadeguato'

# inquiry cds geografiche
#					print 'spherical cds           %8.5f %8.5f' % (rad2sessad(lat),rad2sessad(lon))


	def vectorInquiry(self):
		"""
			dialogo per l'inquiry di un vettore
			(consente di interrogare tutti i vettori tra due estremi,
			non solo i vettori archiviati)
		"""
		dlg = genericDlg('Vector properties',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				a,b = str(a),str(b)
				if a in self.points.keys():
					if b in self.points.keys():
						xa,ya,za = self.points[a][1]
						xb,yb,zb = self.points[b][1]
						ux,uy,uz = xb-xa,yb-ya,zb-za
						mu = math.sqrt(dotProduct([ux,uy,uz],[ux,uy,uz]))
						print '------------vector----------------'
						print 'estremo A %s %7.3f %7.3f %7.3f' % (a,xa,ya,za)
						print 'estremo B %s %7.3f %7.3f %7.3f' % (b,xb,yb,zb)
						print 'vettore      %7.3f %7.3f %7.3f' % (ux,uy,uz)
						print 'modul        %7.3f' % (mu)
					else:
						print 'vertice',b,'inesistente'
				else:
					print 'vertice',a,'inesistente'
			else:
				print "errore: numero di parametri inadeguato"

	def angleInquiry(self):
		"""
			dialogo per l'inquiry di un vettore
			(consente di interrogare tutti i vettori tra due estremi,
			non solo i vettori archiviati)
		"""
		dlg = genericDlg('Three points angle',['P1','P2','P3'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,c = tmp
				a,b,c = str(a),str(b),str(c)
#				print a,b,c
				"""
				if a in self.points.keys():
					if b in self.points.keys():
						if c in self.points.keys():
							xa,ya,za = self.points[a][1]
							xb,yb,zb = self.points[b][1]
							xc,yc,zc = self.points[c][1]
				"""
				tmp = self.dbGetPointRectCds(a)
				if tmp != -1:
					xa,ya,za = tmp
					tmp = self.dbGetPointRectCds(b)
					if tmp != -1:
						xb,yb,zb = tmp
						tmp = self.dbGetPointRectCds(c)
						if tmp != -1:
							xc,yc,zc = tmp
							ux,uy,uz = xa-xb,ya-yb,za-zb
							vx,vy,vz = xc-xb,yc-yb,zc-zb
							mu = math.sqrt(dotProduct([ux,uy,uz],[ux,uy,uz]))
							if mu != 0:
								mv = math.sqrt(dotProduct([vx,vy,vz],[vx,vy,vz]))
								if mv != 0:
									ang = math.acos(dotProduct([ux,uy,uz],[vx,vy,vz])/(mu*mv))
									print '------------angle----------------'
									print 'angle          %s-%s-%s' % (a,b,c)
									print 'radianti       %8.5f' % (ang)
									print 'gradi decimali %8.5f' % (rad2sessad(ang))
								else:
									print 'vettore',c,b,'nullo'
							else:
								print 'vettore',a,b,'nullo'
						else:
							print 'vertice',c,'inesistente'
					else:
						print 'vertice',b,'inesistente'
				else:
					print 'vertice',a,'inesistente'
			else:
				print "errore: numero di parametri inadeguato"

#	-------- delete functions -----------

	def vertexDel(self):
		"""
			dialogo per la cancellazione di un vertice
		"""
		dlg = genericDlg('Vertex to delete',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				a = tmp.pop()
				a = str(a)
#				print 'cancello vertice',a
				# cerco nei vertici
				if a in self.points.keys():
					print 'trovato il punto',a
					del self.points[a]
					# il punto c'è e quindi lo cerco nei great circles
					for k in self.circles.keys():
						if a in self.circles[k]:
							print 'trovato great circle',k
							del self.circles[k]
					# .... negli archi
					for k in self.arcs.keys():
						if a in self.arcs[k]:
							print 'trovato arco',k
							del self.arcs[k]
					# ..... nei vettori
					for k in self.vectors.keys():
						if a in self.vectors[k]:
							print 'trovato vettore',k
							del self.vectors[k]
				else:
					print 'il punto',a,'non esiste in archivio'

	def vectorDel(self):
		dlg = genericDlg('two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				a,b = str(a),str(b)
				if a+b in self.vectors.keys():
					del self.vectors[a+b]
#					self.redraw()
				else:
					print 'Vettore',a+b,'inesistente'
			else:
				print "errore: numero di parametri inadeguato"

	def parallelDel(self):
		dlg = genericDlg('Point',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				p = str(tmp.pop())
#				print p,type(p)
				del self.circles['par'+str(p)]
#				self.redraw()

	def meridianDel(self):
		dlg = genericDlg('Point',['P'])
		dlg.setValues([''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 1:
				p = str(tmp.pop())
#				print p,type(p)
				del self.circles['mer'+str(p)]
#				self.redraw()

	def gCircleDel(self):
		dlg = genericDlg('two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				del self.circles['gc'+str(a)+str(b)]
#				self.redraw()

	def arcDel(self):
		dlg = genericDlg('two points',['P1','P2'])
		dlg.setValues(['',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 2:
				a,b = tmp
				del self.arcs[str(a)+str(b)]
#				self.redraw()

	def arcDelAll(self):
		self.arcs = {}
#		self.redraw()

	def triangleDel(self):
		"""
			dialogo per la cancellazione di un triangolo sferico
		"""
		dlg = genericDlg('Three points',['P1','P2','P3'])
		dlg.setValues(['','',''])
		dlg.show()
		result = dlg.exec_()
		if result:
			tmp = dlg.getValues()
			if len(tmp) == 3:
				a,b,c = tmp
#				print a,b,c
				del self.arcs[str(a)+str(b)]
				del self.arcs[str(b)+str(c)]
				del self.arcs[str(c)+str(a)]
#				self.redraw()
			else:
				print "errore: numero di parametri inadeguato"

#	---- algoritmi di Todhunter-Leathem -------

	def todhunter_000_advice(self):
		QtGui.QMessageBox.information(
			self,
			'000_Advice',
			'''
That is a free interpretation of the book
       "Spherical Trigonometry"
    by I.Todhunter and J.G.Leathem,
        McMillan, London, 1914;
all errors, mistakes, misunderstandings, etc.,
are of mine responsability'
			'''
		)

	def todhunter_001_sphere(self):
		QtGui.QMessageBox.information(
			self,
			'001_sphere',
			'''
Given a 3D euclidean space
with a norm defined by |X| = sqrt(x**2+y**2+z**2),
let choose a point C = (xc,yc,zc),
then the set of all P = (xp,yp,zp) such that
|CP|**2 = (xp-xc)**2+(yp-yc)**2+(zp-zc)**2 <= R**2 = cost
is said SPHERE of radius R;
the set of all P such that |CP|**2 = R**2
is said SURFACE of the sphere;
if C = O the condition simplify in:
xp**2+yp**2+zp**2 <= R**2.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 1 pag. 1'
		# visualizza
		self.redraw()

	def todhunter_002a_smallCircle(self):
		QtGui.QMessageBox.information(
			self,
			'002a_smallCircle',
			'''
Theorem:
The intersection of a sphere of radius r and a plane S
passing through a point P and normal to a vector N
is a circle;
Proof:
let C the intersection of the normal N from O and the
plane S, we have |OC| = k; for all points X belonging
to intersection we have
a) OX = r (X belongs to surface of the sphere)
b) CX is perpendicular to OC by construction,
thus, by Pytagora, |CX|**2 = r**2-|OC|**2 = cost.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 2 pag. 1'
		# punto P
		lat = sessad2rad(70.)
		lon = sessad2rad(10.)
		p = self.dbAddPointLatLon(lat,lon)
		# normale
		xn,yn,zn = self.myRad/math.sqrt(3),self.myRad/math.sqrt(3),self.myRad/math.sqrt(3)
		n = self.dbAddPointRect(xn,yn,zn)
		# cerchio
		self.circles['sc'+str(p)+str(n)] = ['sc',str(p),str(n)]
		# visualizza
		self.redraw()

	def todhunter_002b_parallels(self):
		QtGui.QMessageBox.information(
			self,
			'002b_parallels',
			'''
Among others, important small circles are these
parallel to horizontal plane XY by wich the name
PARALLELS
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914'
		# paralleli
		da = 2*math.pi/self.myNum
		lat = -math.pi/2
		lon = 0.
		for i in range(self.myNum/2):	# in questo modo hanno lo stesso intervalo dei meridiani
			p = self.dbAddPointLatLon(lat,lon)
			self.circles['par'+str(p)] = ['par',str(p)]
			lat += da
		# visualizza
		self.redraw()

	def todhunter_003a_greatCircle(self):
		QtGui.QMessageBox.information(
			self,
			'003a_greatCircle',
			'''
When the secant plane S pass through the center O
of the sphere we have a GREAT CIRCLE
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 3 pag. 1'
		# punto A
		lat = sessad2rad(20.)
		lon = sessad2rad(10.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(60.)
		lon = sessad2rad(40.)
		b = self.dbAddPointLatLon(lat,lon)
		# great circle
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# visualizza
		self.redraw()

	def todhunter_003b_meridians(self):
		QtGui.QMessageBox.information(
			self,
			'003b_meridians',
			'''
Among others, important great circles are these
vertical named MERIDIANS
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914'
		# meridiani
		da = 2*math.pi/self.myNum
		lat = 0.
		lon = 0.
		for i in range(self.myNum):
			p = self.dbAddPointLatLon(lat,lon)
			self.circles['mer'+str(p)] = ['mer',str(p)]
			lon += da
		# visualizza
		self.redraw()

	def todhunter_003c_baseMeridian(self):
		QtGui.QMessageBox.information(
			self,
			'003c_baseMeridian',
			'''
Among the great circles, one of the two most important is
the BASE MERIDIAN
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914'
		# meridiano base
		lat = 0.
		lon = 0.
		p = self.dbAddPointLatLon(lat,lon)
		self.circles['mer'+str(p)] = ['mer',str(p)]
		# visualizza
		self.redraw()

	def todhunter_003d_equator(self):
		QtGui.QMessageBox.information(
			self,
			'003d_equator',
			'''
Among the great circles, one of the two most important is
the EQUATOR
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914'
		# equatore
		lat = 0.
		lon = 0.
		p = self.dbAddPointLatLon(lat,lon)
		self.circles['par'+str(p)] = ['par',str(p)]
		# visualizza
		self.redraw()

	def todhunter_003e_sphericalCds(self):
		QtGui.QMessageBox.information(
			self,
			'003e_sphericalCds',
			'''
The set of parallels and meridians draws a grid
on the sphere surface by wich we can identifify
their points; given a point B and their projection
C on equator along a meridian and D on base meridian
along a parallel, the angle AOC is said LONGITUDE,
the angle AOD is said LATITUDE; the (lat,lon) are
known as SPHERICAL COORDINATES, no other point share
the same spherical coordinates
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914'
		# punto A
		lat = sessad2rad(0.)
		lon = sessad2rad(0.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(20.)
		b = self.dbAddPointLatLon(lat,lon)
		# projection C
		lat = sessad2rad(0.)
		lon = sessad2rad(20.)
		c = self.dbAddPointLatLon(lat,lon)
		# projection D
		lat = sessad2rad(30.)
		lon = sessad2rad(0.)
		d = self.dbAddPointLatLon(lat,lon)
		# visualizza
		self.redraw()

	def todhunter_004_arc(self):
		QtGui.QMessageBox.information(
			self,
			'004_arc',
			'''
The portion between two points A and B
along a great circle is said ARC
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 4 pag. 2'
		# punto A
		lat = sessad2rad(20.)
		lon = sessad2rad(10.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(60.)
		lon = sessad2rad(40.)
		b = self.dbAddPointLatLon(lat,lon)
		# great circle
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# arc
		self.dbAddArc(str(a),str(b))
		# visualizza
		self.redraw()

	def todhunter_005_polar(self):
		QtGui.QMessageBox.information(
			self,
			'005_polar',
			'''
Given a circle and its normal, the intercepts of the last
with the sphere are said POLARS;
in this context we name POLE the intercept nearer to the
center C of the circle
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 5, pag. 3'
		# genera il punto P
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		p = self.dbAddPointLatLon(lat,lon)
		# normale
		xn,yn,zn = 10/math.sqrt(3),10/math.sqrt(3),10/math.sqrt(3)
		n = self.dbAddPointRect(xn,yn,zn)
		# cerchio
		self.circles['sc'+str(p)+str(n)] = ['sc',str(p),str(n)]
		# add polar
		self.dbAddPolar('sc',p,n)
		# visualizza
		self.dbListAll()
		self.redraw()

	def todhunter_006_polarTheorem(self):
		QtGui.QMessageBox.information(
			self,
			'006_polarTheorem',
			'''
Theorem:
The pole is equidistant from the points of the circle;
Proof:
- the distance CP = CP-OC is constant;
- for all point A of the circle, CA is normal by definition to CP
- by Pytagoras PA^2 = CP^2+CA^2 = constant
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 6, pag. 3'
		# genera il punto P
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		p = self.dbAddPointLatLon(lat,lon)
		# normale
		xn,yn,zn = 10/math.sqrt(3),10/math.sqrt(3),10/math.sqrt(3)
		n = self.dbAddPointRect(xn,yn,zn)
		# cerchio
		self.circles['sc'+str(p)+str(n)] = ['sc',str(p),str(n)]
		# polo
		a = self.dbAddPolar('sc',p,n)
		# vector
		self.dbAddVector(p,a)
		# visualizza
		self.redraw()

	def todhunter_007_quadrant(self):
		QtGui.QMessageBox.information(
			self,
			'007_quadrant',
			'''
Theorem:
The arc from a point A of a great circle and its pole P is a quadrant.
Proof:
- for all point A of the circle, CA is normal by definition to CP
- CA = CP = r
- the arc AP = 2math.pi*rad/4 =  math.pi*rad/2
- the angle AP/rad = math.pi/2.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 7, pag. 4'
		# genera il punto A
		lat = sessad2rad(40.)
		lon = sessad2rad(10.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(20.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# cerchio
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# add polar
		p = self.dbAddPolar('gc',a,b)
		# add arc
		self.dbAddArc(str(a),str(p))
		# visualizza
		self.redraw()

	def todhunter_008_gcInclinationTheorem(self):
		QtGui.QMessageBox.information(
			self,
			'008_gcInclinationTheorem',
			'''
Theorem:
The angle subtended by the arc between the poles of two
great circles is equal to the angle between the planes
of the great circles.
Proof:
Let the great circle passing through the poles E, F and
define as G, H the intersections with the two great circles;
because the line OE, OF, OG, OH are pairwise normal, the
angle GOH is equal to EOF.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 8, pag. 4'
		# punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(20.)
		lon = sessad2rad(60.)
		b = self.dbAddPointLatLon(lat,lon)
		# great circle AB
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# punto C
		lat = sessad2rad(30.)
		lon = sessad2rad(20.)
		c = self.dbAddPointLatLon(lat,lon)
		# punto D
		lat = sessad2rad(50.)
		lon = sessad2rad(60.)
		d = self.dbAddPointLatLon(lat,lon)
		# great circle CD
		self.circles['gc'+str(c)+str(d)] = ['gc',str(c),str(d)]
		# polo E
		e = self.dbAddPolar('gc',a,b)
		# polo F
		f = self.dbAddPolar('gc',c,d)
		# great circle EF
		self.circles['gc'+str(e)+str(f)] = ['gc',str(e),str(f)]
		# arco EF
		self.dbAddArc(str(e),str(f))
		# fare intersezioni M,N
		c1 = 'gcEF'
		c2 = 'gcAB'
		m = self.dbAddIntersection(c1,c2)
		c1 = 'gcEF'
		c2 = 'gcCD'
		n = self.dbAddIntersection(c1,c2)
		# fare arco MN
		self.dbAddArc(str(m),str(n))
		# visualizza
		self.redraw()

	def todhunter_009_gcInclination(self):
		QtGui.QMessageBox.information(
			self,
			'009_gcInclination',
			'''
When two circles intersects, the angle between the tangents
at either of their points of intersection is called the angle
between the circles.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 9, pag. 5'
		# punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(20.)
		lon = sessad2rad(60.)
		b = self.dbAddPointLatLon(lat,lon)
		# great circle AB
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# punto C
		lat = sessad2rad(30.)
		lon = sessad2rad(20.)
		c = self.dbAddPointLatLon(lat,lon)
		# punto D
		lat = sessad2rad(50.)
		lon = sessad2rad(60.)
		d = self.dbAddPointLatLon(lat,lon)
		# great circle CD
		self.circles['gc'+str(c)+str(d)] = ['gc',str(c),str(d)]
		# intersection: valutare se completare il diamtero con l'antipodale
		c1 = 'gcAB'
		c2 = 'gcCD'
		l = self.dbAddIntersection(c1,c2)
		m = self.dbAddAntipodalPoint(l)
		# tangents
		self.dbAddTangent(l,'O',a)
		self.dbAddTangent(l,'O',c)
		# visualizza
		self.redraw()

	def todhunter_010_greatCircleBisection(self):
		QtGui.QMessageBox.information(
			self,
			'010_greatCircleBisection',
			'''
Theorem:
Two great circles bisect each other.
Proof:
(to be done)
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 10, pag. 5'
		# punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# punto B
		lat = sessad2rad(20.)
		lon = sessad2rad(60.)
		b = self.dbAddPointLatLon(lat,lon)
		# great circle AB
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		# punto C
		lat = sessad2rad(30.)
		lon = sessad2rad(20.)
		c = self.dbAddPointLatLon(lat,lon)
		# punto D
		lat = sessad2rad(50.)
		lon = sessad2rad(60.)
		d = self.dbAddPointLatLon(lat,lon)
		# great circle CD
		self.circles['gc'+str(c)+str(d)] = ['gc',str(c),str(d)]
		# intersection: valutare se completare il diametro con l'antipodale
		c1 = 'gcAB'
		c2 = 'gcCD'
		l = self.dbAddIntersection(c1,c2)
		m = self.dbAddAntipodalPoint(l)
		# visualizza
		self.redraw()

	def todhunter_014_arcComparison(self):
		QtGui.QMessageBox.information(
			self,
			'014_arcComparison',
			'''
The arc of small circle  as function of the arc of
great circles subtending the same angle:
CD = AB * cos(lat).
Caveat:
the arco CD is a small circle arc, (se ci si dimentica
si può inferire l'eguaglianza CD/rad = angleCD =
AB/rad * cos = lon * cos(lat) che è sbagliata)

NB: arco errato: deve essere un arco di small circle,
non di great circle!
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 14, pag. 7'
		# genera il punto A sull'equatore
		lat = sessad2rad(0.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B sull'equatore
		lat = sessad2rad(0.)
		lon = sessad2rad(60.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(50.)
		lon = sessad2rad(20.)
		c = self.dbAddPointLatLon(lat,lon)
		# genera il punto D
		lat = sessad2rad(50.)
		lon = sessad2rad(60.)
		d = self.dbAddPointLatLon(lat,lon)
		# genera il punto E
		x,y,z = 0.,0.,self.myRad*math.sin(sessad2rad(50.))
		e = self.dbAddPointRect(x,y,z)
		# genera i meridiani 
		self.circles['mer'+str(a)] = ['mer',str(a)]
		self.circles['par'+str(b)] = ['mer',str(b)]
		# genera gli archi
		self.dbAddSmallCircleArc(str(a),str(b),'O')
		self.dbAddSmallCircleArc(str(c),str(d),str(e))
		# genera i vettori
		self.dbAddVector(e,c)
		self.dbAddVector(e,d)
		# visualizza
		self.redraw()

	def todhunter_015_solidAngle(self):
		QtGui.QMessageBox.information(
			self,
			'015_solidAngle',
			'''
Planes intersecting at a common point define a
solid angle; if the point is a center of a sphere
every plane define a great circle; that defines
a spherical polygon A, B, C, C, E, etc...
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 15, pag. 8'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto c
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# genera il punto d
		lat = sessad2rad(30.)
		lon = sessad2rad(-20.)
		d = self.dbAddPointLatLon(lat,lon)
		# genera il punto e
		lat = sessad2rad(10.)
		lon = sessad2rad(-10.)
		e = self.dbAddPointLatLon(lat,lon)
		# greate circle
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		self.circles['gc'+str(b)+str(c)] = ['gc',str(b),str(c)]
		self.circles['gc'+str(c)+str(d)] = ['gc',str(c),str(d)]
		self.circles['gc'+str(d)+str(e)] = ['gc',str(d),str(e)]
		self.circles['gc'+str(e)+str(a)] = ['gc',str(e),str(a)]
		# arcs
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(c),str(d))
		self.dbAddArc(str(d),str(e))
		self.dbAddArc(str(e),str(a))
		# visualizza
		self.redraw()

	def todhunter_016_sphericalTriangle(self):
		QtGui.QMessageBox.information(
			self,
			'016_sphericalTriangle',
			'''
In the case of three planes we have a spherical
triangle A, B, C.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 15, pag. 8'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto c
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# greate circle
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		self.circles['gc'+str(b)+str(c)] = ['gc',str(b),str(c)]
		self.circles['gc'+str(c)+str(a)] = ['gc',str(c),str(a)]
		# arcs
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(c),str(a))
		# visualizza
		self.redraw()

	def todhunter_017_sidesAndAngles(self):
		QtGui.QMessageBox.information(
			self,
			'017_sidesAndAngles',
			'''
In a spherical triangle A, B, C we have:
- the sides defined as "a", "b" and "c"
- the angles defined by "A", "B" and "C";
the sides counter???? the angle, that is
"a" = "BC"/R.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 15, pag. 8'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto c
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# greate circle
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		self.circles['gc'+str(b)+str(c)] = ['gc',str(b),str(c)]
		self.circles['gc'+str(c)+str(a)] = ['gc',str(c),str(a)]
		# arcs
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(c),str(a))
		# visualizza
		self.redraw()

	def todhunter_018_angles(self):
		QtGui.QMessageBox.information(
			self,
			'018_angles',
			'''
Let A, B, C a spherical triangle, take the tangents
from A to OB and OC; the angle between the tangents
is the angle at A.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 15, pag. 8'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto c
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# arcs
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(a),str(c))
		# tangents
		self.dbAddTangent(a,'O',b)
		self.dbAddTangent(a,'O',c)
		# visualizza
		self.redraw()

	def todhunter_024a_lune(self):
		QtGui.QMessageBox.information(
			self,
			'024a_lune',
			'''
Let A, B and C, D two great circle; thei divide the
sphere surface in four section each said LUNE.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 24, pag. 12'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# genera il punto D
		lat = sessad2rad(80.)
		lon = sessad2rad(20.)
		d = self.dbAddPointLatLon(lat,lon)
		# great circles
		self.circles['gc'+str(a)+str(b)] = ['gc',str(a),str(b)]
		self.circles['gc'+str(c)+str(d)] = ['gc',str(c),str(d)]
		# visualizza
		self.redraw()

	def todhunter_024b_colunarTriangles(self):
		QtGui.QMessageBox.information(
			self,
			'024b_colunarTriangles',
			'''
Let A, B, C and D, B, C two spherical triangle where BC is in common
and D is antipodal of A; the triangles together form a lune, thus
they are said COLUNAR.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 24, pag. 12'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(30.)
		lon = sessad2rad(30.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto c
		lat = sessad2rad(50.)
		lon = sessad2rad(0.)
		c = self.dbAddPointLatLon(lat,lon)
		# punto D antipodo di A
		d = self.dbAddAntipodalPoint(a)
		# arcs
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		self.dbAddArc(str(d),str(b))
		self.dbAddArc(str(d),str(c))
		# visualizza
		self.redraw()

	def todhunter_024c_antipodalTriangles(self):
		QtGui.QMessageBox.information(
			self,
			'024c_antipodalTriangles',
			'''
Let A, B, C three point and D, E, F their antipodes;
the triangles A, B, C and D, E, F are said  ANTIPODAL.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 24, pag. 12'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# antipolare
		a1 = self.dbAddAntipodalPoint(a)
		b1 = self.dbAddAntipodalPoint(b)
		c1 = self.dbAddAntipodalPoint(c)
		# genera il triangolo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# antipolare
		self.dbAddArc(str(a1),str(b1))
		self.dbAddArc(str(b1),str(c1))
		self.dbAddArc(str(c1),str(a1))
		self.redraw()

	def todhunter_025_polarTriangles(self):
		QtGui.QMessageBox.information(
			self,
			'025_polarTriangles',
			'''
Let A, B, C three point and D, E, F poles respectively of
BC, AC and AB, the triangles ABC and DEF are POLARS.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 25, pag. 12'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# polari
		a1 = self.dbAddPolar('gc',b,c)
		b1 = self.dbAddPolar('gc',a,c)
		c1 = self.dbAddPolar('gc',a,b)
		# genera il triangolo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# polare
		self.dbAddArc(str(a1),str(b1))
		self.dbAddArc(str(b1),str(c1))
		self.dbAddArc(str(c1),str(a1))
		self.redraw()

	def todhunter_026_simmetryOfPolarity(self):
		QtGui.QMessageBox.information(
			self,
			'026_simmetryOfPolarity',
			'''
Theorem:
Let ABC a spherical triangle and DEF it's polar
triangle, thus ABC is the polar triangle of DEF.
Proof:
D polar of BC imply DB and DC are quadrants,
E polar of AC imply EA and EC are quadrants,
F polar of AB imply FA and FB are quadrants,
thus A is polar of EF, B of DF and C of DE.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 26, pag. 13'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# polari
		a1 = self.dbAddPolar('gc',b,c)
		b1 = self.dbAddPolar('gc',a,c)
		c1 = self.dbAddPolar('gc',a,b)
		# genera il triangolo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# polare
		self.dbAddArc(str(a1),str(b1))
		self.dbAddArc(str(b1),str(c1))
		self.dbAddArc(str(c1),str(a1))
		self.redraw()

	def todhunter_027a_supplementalTriangles(self):
		QtGui.QMessageBox.information(
			self,
			'027a_supplementalTriangles',
			'''
Theorem:
Let ABC and DEF two polar triangle, the sides og
one are supplemental of the angles of the other;
the angles the supplemental of the sides.
Proof:
Let G and H the intersection of great circles BA
and BC with great circle DF; DH and FG are by
definition quadrant, thus DH + FG = DF + GH = 2pi;
by the simmetry of polarity D is a pole of DF, thus
the ar GH is equals to the angle in B, from wich
DF = 2pi - GH = 2pi - B.
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 26, pag. 13'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# polari
		a1 = self.dbAddPolar('gc',b,c)
		b1 = self.dbAddPolar('gc',a,c)
		c1 = self.dbAddPolar('gc',a,b)
		# genera il triangolo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# polare
		self.dbAddArc(str(a1),str(b1))
		self.dbAddArc(str(b1),str(c1))
		self.dbAddArc(str(a1),str(c1))
		# great circles
		typ = 'gc'
		self.circles[typ+str(b)+str(a)] = [typ,b,a]
		self.circles[typ+str(b)+str(c)] = [typ,b,c]
		self.circles[typ+str(a1)+str(c1)] = [typ,a1,c1]
		# intersezioni
		c1 = 'gc'+a1+c1
		c2 = 'gc'+b+a
		l = self.dbAddIntersection(c1,c2)
		c2 = 'gc'+b+c
		m = self.dbAddIntersection(c1,c2)
		# visualizza
		self.redraw()

	def todhunter_027b_supplementalTriangles_2(self):
		QtGui.QMessageBox.information(
			self,
			'027b_supplementalTriangles_2',
			'''
Proposition:
The sides (angles) of the polar triangle are equal
or supplemental of the angles (sides) of the
primitive triangle;
Proof:

			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 26, pag. 13'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# polari
		a1 = self.dbAddPolar('gc',b,c)
		b1 = self.dbAddPolar('gc',a,c)
		c1 = self.dbAddPolar('gc',a,b)
		# genera il triangolo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# polare
		self.dbAddArc(str(a1),str(b1))
		self.dbAddArc(str(b1),str(c1))
		self.dbAddArc(str(a1),str(c1))
		# great circles
		typ = 'gc'
		self.circles[typ+str(a)+str(b)] = [typ,a,b]
		self.circles[typ+str(a)+str(c)] = [typ,a,c]
		self.circles[typ+str(b1)+str(c1)] = [typ,b1,c1]
		# intersezioni
		c1 = 'gc'+b1+c1
		c2 = 'gc'+a+b
		l = self.dbAddIntersection(c1,c2)
		c2 = 'gc'+a+c
		m = self.dbAddIntersection(c1,c2)
		# visualizza
		self.redraw()

	def euclid_XI_20_solidAngle(self):
		QtGui.QMessageBox.information(
			self,
			'XI_20_solidAngle',
			'''
Proposition:
In a solid angle defined by the three planes
AOB, BOC, COA, the sum of any two angles
AOB + BOC is greater than the third angle AOC ;
Proof:
If the angle AOC be less than or equal to either
of the other angles the proposition is evident.
If not, suppose it greater: take any point D in
AC such that AOD = AOB, thus AB = AD.
In plane triangle ABC we have AB+BC > AC = AD+DC,
thus subtracting the same quantity AB = AD, we
find BC > DC; from ??? we have BOC > DOC, thus
AOB+BOC > AOB+DOC = AOD+DOC = AOC. 
			'''
		)
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'Euclid\n Elements, book XI, proposition 20'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		xa,ya,za = pointSphericCds(self.myRad,lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(50.)
		b = self.dbAddPointLatLon(lat,lon)
		xb,yb,zb = pointSphericCds(self.myRad,lat,lon)
		# genera il punto C
		lat = sessad2rad(10.)
		lon = sessad2rad(70.)
		c = self.dbAddPointLatLon(lat,lon)
		xc,yc,zc = pointSphericCds(self.myRad,lat,lon)
		# genera il triangolo piano
		self.dbAddVector(a,b)
		self.dbAddVector(b,c)
		self.dbAddVector(c,a)
		# genera il punto D tale che AB=AD
		mab = math.sqrt(dotProduct([xa-xb,ya-yb,za-zb],[xa-xb,ya-yb,za-zb]))
		print 'AB:',mab
		mac = math.sqrt(dotProduct([xc-xa,yc-ya,zc-za],[xc-xa,yc-ya,zc-za]))
		print 'AC:',mac
		l = self.dbAddInterpolation(a,c,mab/mac)
		self.dbAddVector('O',l)
		self.dbAddVector(b,l)
		# visualizza
		self.redraw()





	def todhunter_047(self):
		"""
			spherical triangle
		"""
		# pulisce tutto
		self.dbClearAll()
		# titolo
		self.myTitle = 'I.Todhunter J.G.Leathem\n Spherical Trigonometry\n McMillan, London 1914, art. 47, pag. 25'
		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(20.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(40.)
		lon = sessad2rad(70.)
		b = self.dbAddPointLatLon(lat,lon)
		# genera il punto C
		lat = sessad2rad(60.)
		lon = sessad2rad(30.)
		c = self.dbAddPointLatLon(lat,lon)
		# triangolo primitivo
		self.dbAddArc(str(a),str(b))
		self.dbAddArc(str(b),str(c))
		self.dbAddArc(str(a),str(c))
		# legge le coordinate
		xa,ya,za = self.points[a][1]
		xb,yb,zb = self.points[b][1]
		xc,yc,zc = self.points[c][1]
		# punto P a metà OC
		xp,yp,zp = interpolationPoint([0.,0.,0.],[xc,yc,zc],0.5)
		p = self.dbAddPointRect(xp,yp,zp)
		k,j = normaleDaPuntoAPiano([xp,yp,zp],[xa,ya,za],[xb,yb,zb])
		xh,yh,zh = k*xa+j*xb,k*ya+j*yb,k*za+j*zb
		print 'punto H',xh,yh,zh
		h = self.dbAddPointRect(xh,yh,zh)
		# normale E da H a OA
		xe,ye,ze = normalFromPtoAB([xp,yp,zp],[0.,0.,0.],[xa,ya,za])
		e = self.dbAddPointRect(xe,ye,ze)
		# normale F da H a OB
		xf,yf,zf = normalFromPtoAB([xp,yp,zp],[0.,0.,0.],[xb,yb,zb])
		f = self.dbAddPointRect(xf,yf,zf)
		self.redraw()

#	--- routine di prova-------

	def prova(self):

		# genera il punto A
		lat = sessad2rad(10.)
		lon = sessad2rad(0.)
		a = self.dbAddPointLatLon(lat,lon)
		# genera il punto B
		lat = sessad2rad(70.)
		lon = sessad2rad(0.)
		b = self.dbAddPointLatLon(lat,lon)
		# origine
		c = 'O'
		# genera il punto C
		lat = sessad2rad(30.)
		lon = sessad2rad(0.)
		d = self.dbAddPointLatLon(lat,lon)

		self.dbAddVector(a,b)
		self.dbAddVector(c,d)







# ======================= runtime =================================

myApp = QtGui.QApplication(sys.argv)
myWin = MainWindow()
myWin.show()
sys.exit(myApp.exec_())




