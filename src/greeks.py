import numpy as np
from scipy.stats import norm

def bs_greeks(S0,K,T,r,q,sigma,option_type='put'):
    '''
    Calculate Greeks letter with continuous dividend yield based on Black Scholes formula
    ================================================
    r: float
        risk free rate
    q: float 
        continuous dividend yield
    option_type: str 'call'/'put'
        confirm option type whether call or put
    ================================================
    return: float,float
        Greek letters Delta, Vega respectively
    '''

    if S0 <= 0 or K <= 0 or T <= 0 or sigma <= 0:
            raise ValueError("S0, K, T and sigma must be positive")
    
    d1 = (
        np.log(S0/K) + (r-q+0.5*sigma**2)*T
    ) / (sigma*np.sqrt(T))

    if option_type == 'call':
        delta = np.exp(-q*T)*norm.cdf(d1)
    elif option_type == 'put':
        delta = np.exp(-q*T)*(norm.cdf(d1)-1)
    else:
        raise ValueError("Invalid option type")

    vega = (
        S0*np.exp(-q*T)
        * norm.pdf(d1)
        * np.sqrt(T)
    )

    return delta, vega
