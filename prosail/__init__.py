import numpy as np
from prosail_fortran import run_sail, run_prosail, prospect_5b
from prosail_fortran import mod_dataspec_p5b

def trans_prosail ( N, cab, car, cbrown, cw, cm, lai, lidfa, lidfb, rsoil, psoil, \
        hspot, tts, tto, psi, typelidf):
    """A version of PROSAIL that uses transformed parameters to quasi-linearise
    the   model. See http://dx.doi.org/10.1016/j.rse.2011.12.027"""
    # Define the constants
    slai = -2.0
    skab = -100.0
    skar = -100.0
    skw =  -1./50.
    skm =  -1./100.
    # Transform the parameters to real units
    xlai = slai * np.log ( lai )
    xkab = skab * np.log ( cab )
    xkar = skar * np.log ( car )
    xkw = skw * np.log ( cw )
    xdm = skm * np.log ( cm )
    # Run the PROSAIL model
    retval = run_prosail ( N, xkab, xkar, cbrown, xkw, xdm, xlai, \
            lidfa, lidfb, rsoil, psoil, hspot, tts, tto, psi, typelidf )
    return retval

# nested stuff
'''
tau average
'''
def tav_abs(theta,refr):
    thetarad=np.pi*theta/180.
    if theta == 0:
        res=4.*refr/(refr+1.)**2
        return res

    refr2=refr*refr
    ax=(refr+1.)**2/2.
    bx=-(refr2-1.)**2/4.

    if thetarad == np.pi/2.:
        b1=0.
    else:
        b1=np.sqrt((np.sin(thetarad)**2-(refr2+1.)/2.)**2+bx)

    b2=np.sin(thetarad)**2-(refr2+1.)/2.
    b0=b1-b2
    ts=(bx**2/(6.*b0**3)+bx/b0-b0/2.)-(bx**2/(6.*ax**3)+bx/ax-ax/2.)
    tp1=-2.*refr2*(b0-ax)/(refr2+1.)**2
    tp2=-2.*refr2*(refr2+1.)*np.log(b0/ax)/(refr2-1.)**2
    tp3=refr2*(1./b0-1./ax)/2.
    tp4=16.*refr2**2*(refr2**2+1.)*np.log((2.*(refr2+1.)*b0-(refr2-1.)**2)/ \
            (2.*(refr2+1.)*ax-(refr2-1.)**2))/((refr2+1.)**3*(refr2-1.)**2)
    tp5=16.*refr2**3*(1./(2.*(refr2+1.)*b0-((refr2-1.)**2))-1./(2.*(refr2+1.) \
            *ax-(refr2-1.)**2))/(refr2+1.)**3
    tp=tp1+tp2+tp3+tp4+tp5
    res=(ts+tp)/(2.*np.sin(thetarad)**2)

    return res
    
'''
Define soil
'''

charsoil = '''
350 0.04,360 0.04,370 0.04,380 0.04,390 0.04,400 0.04,410 0.04,420 0.04,430 0.04,440 0.04,450 0.04,460 0.04,470 0.04,480 0.04,490 0.04,500 0.04,510 0.04,520 0.04,530 0.04,540 0.04,550 0.04,560 0.04,570 0.04,580 0.04,590 0.04,600 0.04,610 0.04,620 0.04,630 0.04,640 0.04,650 0.04,660 0.04,670 0.04,680 0.04,690 0.04,700 0.04,710 0.04,720 0.04,730 0.04,740 0.04,750 0.04,760 0.04,770 0.04,780 0.04,790 0.04,800 0.04,810 0.04,820 0.04,830 0.04,840 0.04,850 0.04,860 0.04,870 0.04,880 0.04,890 0.04,900 0.04,910 0.04,920 0.04,930 0.04,940 0.04,950 0.04,960 0.04,970 0.04,980 0.04,990 0.04,1000 0.04,1010 0.04,1020 0.04,1030 0.04,1040 0.04,1060 0.04,1060 0.04,1070 0.04,1080 0.04,1090 0.04,1100 0.04,1110 0.04,1120 0.04,1130 0.04,1140 0.04,1150 0.04,1160 0.04,1170 0.04,1180 0.04,1190 0.04,1200 0.04,1210 0.04,1220 0.05,1230 0.05,1240 0.05,1250 0.05,1260 0.05,1270 0.05,1280 0.05,1290 0.05,1300 0.05,1310 0.05,1320 0.05,1330 0.05,1340 0.05,1350 0.05,1360 0.05,1370 0.05,1380 0.05,1390 0.05,1400 0.05,1410 0.05,1420 0.05,1430 0.05,1440 0.05,1450 0.05,1460 0.05,1470 0.05,1480 0.05,1490 0.05,1500 0.05,1510 0.05,1520 0.05,1530 0.05,1540 0.05,1550 0.05,1560 0.05,1570 0.05,1580 0.05,1590 0.05,1600 0.05,1610 0.05,1620 0.05,1630 0.05,1640 0.05,1650 0.05,1660 0.05,1670 0.05,1680 0.05,1690 0.05,1700 0.05,1710 0.05,1720 0.05,1730 0.05,1740 0.05,1750 0.05,1760 0.05,1770 0.05,1780 0.05,1790 0.06,1800 0.06,1810 0.06,1820 0.06,1830 0.06,1840 0.06,1850 0.06,1860 0.06,1870 0.06,1880 0.06,1890 0.06,1900 0.06,1910 0.06,1920 0.06,1930 0.06,1940 0.06,1950 0.06,1960 0.06,1970 0.06,1980 0.06,1990 0.06,2000 0.06,2010 0.06,2020 0.06,2030 0.06,2040 0.06,2050 0.06,2060 0.06,2070 0.06,2080 0.06,2090 0.06,2100 0.06,2110 0.06,2120 0.06,2130 0.06,2140 0.06,2150 0.06,2160 0.06,2170 0.06,2180 0.06,2190 0.06,2200 0.06,2210 0.06,2220 0.06,2230 0.06,2240 0.06,2250 0.06,2260 0.06,2270 0.06,2280 0.06,2290 0.06,2300 0.07,2310 0.07,2320 0.07,2330 0.07,2340 0.07,2350 0.07,2360 0.07,2370 0.07,2380 0.07,2390 0.07,2400 0.07,2410 0.07,2420 0.07,2430 0.07,2440 0.07,2450 0.07,2460 0.07,2470 0.07,2480 0.07,2490 0.07,2500 0.06'''

def soil(x,scale=None,trans=False):
   '''
    x is a dictionary
    
    ************************************
    x['spectra'] should contain spectra
    ************************************
    
    e.g.
    
    x['spectra']['dry']
    x['spectra']['char']
    
    etc.
    
    If this doesnt exist, it is loaded from the fortran data files
    and provides 'cab', 'car', 'cbrown', 'cw', 'cm'
    
    ************************************
    x['params'] should contain concentration values
    ************************************
    
     e.g.
    
    x['params']['cbrown']
    x['params']['cab']
    
    etc.
    
    '''
    # wavelength
    if not 'lamdba' in x:
        #try:
        #    x['lambda'] = eval('mod_dataspec_p5b.lambda')
        #except:
        x['lambda'] = np.arange(400.,2501.)

    x['nw'] = len(x['lambda'])

    if not 'params' in x:
        # return zero arrays of correct size if no parameters specified
        return np.zeros(x['nw']),np.zeros(x['nw']),x

    if not 'spectra' in x:
        x['spectra'] = {}

    # load defaults
    result = np.zeros(x['nw'])
    for p in x['params']:
        if not p in x['spectra']:
            if p is 'dry':
              x['spectra'][p] = mod_dataspec_p5b.rsoil1
            if p is 'wet':
              x['spectra'][p] = mod_dataspec_p5b.rsoil2
            if p is 'char':
              cchar = np.array([i.split() for i in np.split(charsoil,',')]).T
              x['spectra'][p] = \
                scipy.interpolate.interp1d(cchar[0],cchar[1])(x['lambda'])
        try:
                result += x['params'][p] * x['spectra'][p]
        



'''
Define leaf absorbing constituents
'''
def leaf(x,theta2=40.,scale=None,trans=False):
    '''
    x is a dictionary
    
    ************************************
    x['spectra'] should contain spectra
    ************************************
    
    e.g.
    
    x['spectra']['cbrown']
    x['spectra']['cab']
    
    etc.
    
    If this doesnt exist, it is loaded from the fortran data files
    and provides 'cab', 'car', 'cbrown', 'cw', 'cm'
    
    ************************************
    x['params'] should contain concentration values
    ************************************
    
     e.g.
    
    x['params']['cbrown']
    x['params']['cab']
    
    etc.
    
    '''
    # number of leaf layers
    if not 'n' in x:
        x['n'] = 1.0

    # wavelength
    if not 'lamdba' in x:
        #try:
        #    x['lambda'] = eval('mod_dataspec_p5b.lambda')
        #except:
        x['lambda'] = np.arange(400.,2501.)

    x['nw'] = len(x['lambda'])

    if not 'params' in x:
        # return zero arrays of correct size if no parameters specified
        return np.zeros(x['nw']),np.zeros(x['nw']),x

    if not 'spectra' in x:
        x['spectra'] = {}

    # get refractive index
    if not 'refractive' in x:
        x['refractive'] = mod_dataspec_p5b.refractive

    # load defaults
    for p in x['params']:
        if not p in x['spectra']:
            try:
                kterm = 'k_c' + p[1:]
                # potential for security issues using eval? find a
                # better way to do this
                x['spectra'][p] = eval('mod_dataspec_p5b'+'.'+kterm)
            except:
                pass

    x['k'] = np.zeros(x['nw'])
    x['tau'] = np.ones(x['nw'])

    # check scaling
    scale = scale or {'cab':100.,'car':100.,'cw':1/50.,'cm':1/100.}
    if not 'scale' in x:
      x['scale'] = {}

    for p in scale.keys():
      if p in x['params']:
        x['scale'][p] = scale[p]


    for p in x['params']:
        if p in x['spectra']:
            v = x['params'][p]
            if trans:
                # try a transformation
                if p in x['scale']:
                    v =  - x['scale'][p] * np.log ( v )
            x['k'] += v * x['spectra'][p]
    x['k'] /= x['n']

    # upper limit
    ww = np.where(x['k'] >= 85)[0]
    x['tau'][ww] = 0.

    # lower limit
    ww = np.where(x['k'] <= 4)[0]
    if len(ww):
        xx=0.5*x['k'][ww]-1.0
        yy=(((((((((((((((-3.60311230482612224e-13 \
            *xx+3.46348526554087424e-12)*xx-2.99627399604128973e-11) \
            *xx+2.57747807106988589e-10)*xx-2.09330568435488303e-9) \
            *xx+1.59501329936987818e-8)*xx-1.13717900285428895e-7) \
            *xx+7.55292885309152956e-7)*xx-4.64980751480619431e-6) \
            *xx+2.63830365675408129e-5)*xx-1.37089870978830576e-4) \
            *xx+6.47686503728103400e-4)*xx-2.76060141343627983e-3) \
            *xx+1.05306034687449505e-2)*xx-3.57191348753631956e-2) \
            *xx+1.07774527938978692e-1)*xx-2.96997075145080963e-1
        yy=(yy*xx+8.64664716763387311e-1)*xx+7.42047691268006429e-1
        yy=yy-np.log(x['k'][ww])
        x['tau'][ww] = (1.0-x['k'][ww])*np.exp(-x['k'][ww])+x['k'][ww]**2*yy

    ww = np.where((x['k'] > 4) * (x['k'] <= 85))[0]
    if len(ww):
        xx=14.5/(x['k'][ww]+3.25)-1.0
        yy=(((((((((((((((-1.62806570868460749e-12 \
                *xx-8.95400579318284288e-13)*xx-4.08352702838151578e-12) \
                *xx-1.45132988248537498e-11)*xx-8.35086918940757852e-11) \
                *xx-2.13638678953766289e-10)*xx-1.10302431467069770e-9) \
                *xx-3.67128915633455484e-9)*xx-1.66980544304104726e-8) \
                *xx-6.11774386401295125e-8)*xx-2.70306163610271497e-7) \
                *xx-1.05565006992891261e-6)*xx-4.72090467203711484e-6) \
                *xx-1.95076375089955937e-5)*xx-9.16450482931221453e-5) \
                *xx-4.05892130452128677e-4)*xx-2.14213055000334718e-3
        yy=((yy*xx-1.06374875116569657e-2)*xx-8.50699154984571871e-2)*xx+\
                9.23755307807784058e-1
        yy=np.exp(-x['k'][ww])*yy/x['k'][ww]
        x['tau'][ww]=(1.0-x['k'][ww])*np.exp(-x['k'][ww])+x['k'][ww]**2*yy

    # transmissivity of the layer
    tau = x['tau']

    theta1=90.
    t1 = tav_abs(theta1,x['refractive'])
    t2 = tav_abs(theta2,x['refractive'])
    x['t1'] = t1
    x['t2'] = t2
    x1=1-t1
    x2=t1**2*tau**2*(x['refractive']**2-t1)
    x3=t1**2*tau*x['refractive']**2
    x4=x['refractive']**4-tau**2*(x['refractive']**2-t1)**2
    x5=t2/t1
    x6=x5*(t1-1)+1-t2
    r=x1+x2/x4
    t=x3/x4
    ra=x5*r+x6
    ta=x5*t
    return ra,ta,x
