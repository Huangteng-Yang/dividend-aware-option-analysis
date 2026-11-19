import numpy as np
import pandas as pd

from scipy.stats import norm 


def bs_calculator(S0,K,T,r,q,sigma,option_type):
    '''
    Calculate option price with continuous dividend yield based on Black Scholes formula
    ================================================
    r: float
        risk free rate
    q: float 
        continuous dividend yield
    option_type: str 'call'/'put'
        confirm option type whether call or put
    ================================================
    return: float
        option price based on Black Scholes formula
    '''
    if (S0<=0) or (K<=0) or (T<=0):
        raise ValueError('Invalid market parameters')
    d1=(np.log(S0/K)+(r-q+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    d2=d1-sigma*np.sqrt(T)

    # Calculate
    if option_type=='call':
        c=S0*np.exp(-q*T)*norm.cdf(d1)-K*np.exp(-r*T)*norm.cdf(d2)
        return c
    elif option_type=='put':
        p=K*np.exp(-r*T)*norm.cdf(-d2)-S0*np.exp(-q*T)*norm.cdf(-d1)
        return p 
    else:
        raise ValueError('Invalid option type')

