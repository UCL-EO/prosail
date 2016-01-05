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
    res=(ts+tp)/(2.*sin(thetarad)**2)
    
    return res

'''
Define leaf absorbing constituents
'''
def leaf(x,theta2=40.):
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
    if not 'N' in x:
        x['N'] = 1.0
        
    # wavelength
    if not 'lamdba' in x:
        x['lambda'] = eval('mod_dataspec_p5b.lambda')
        # yuk = bacuase of reserved name 'lambda'
        
    x['nw'] = len(x['lambda'])
    
    if not 'params' in x:
        # return zero arrays of correct size if no parameters specified
        return np.zeros(x['nw']),np.zeros(x['nw'])
    
    if not 'spectra' in x: 
        x['spectra'] = {}
        
    # get refractive index
    if not 'n' in x:
        x['n'] = mod_dataspec_p5b.refractive
        
    # load defaults
    for p in x['params']:
        if not p in x['spectra']:
            try:
                kterm = 'k_' + p[1:]
                # potential for security issues using eval? find a
                # better way to do this
                x['spectra'][p] = eval('mod_dataspec_p5b'+'.'+kterm)
            except:
                pass

    x['k'] = np.zeros(mod_dataspec_p5b.nw)
    x['k'][x['k']<0] = 0.
    
    x['tau'] = np.ones(mod_dataspec_p5b.nw)
    for p in x['params']:
        if p in x['spectra']:
            x['k'] += p * x['spectra']
    x['k'] /= x['N']
    
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
        x['tau']=(1.0-x['k'][ww])*np.exp(-x['k'][ww])+x['k'][ww]**2*yy
    
    # transmissivity of the layer
    tau = x['tau']
    
    theta1=90.
    t1 = tav_abs(theta1,x['n'])
    t2 = tav_abs(theta2,x['n'])
    x1=1-t1
    x2=t1**2*tau**2*(x['n']**2-t1)
    x3=t1**2*tau*x['n']**2
    x4=x['n']**4-tau**2*(x['n']**2-t1)**2
    x5=t2/t1
    x6=x5*(t1-1)+1-t2
    r=x1+x2/x4
    t=x3/x4
    ra=x5*r+x6
    ta=x5*t
    return ra,ta

    
