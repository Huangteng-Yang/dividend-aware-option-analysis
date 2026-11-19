import black_scholes as bs
import greeks as gk

import numpy as np
import math

from scipy.optimize import brentq


def bisection(S0,K,T,r,q,mkt_price,option_type='put',sigma_tol=1e-8,price_tol=1e-8,max_sigma=8,max_iteration=500):
    '''
    mkt_price: market price
    tol: tolerance of error
    max_iteration: acceptable maximum times of iteration
    =========================
    return: float,int,[float]
        implied volatility from bisection method, number of iterations, iteration history of sigma
    '''

    # bounds of market price
    discount_S=S0*np.exp(-q*T)
    discount_K=K*np.exp(-r*T)

    if option_type=='call':
        low_bound=max(discount_S-discount_K,0)
        high_bound=discount_S
    elif option_type=='put':
        low_bound=max(discount_K-discount_S,0)
        high_bound=discount_K
    else:
        raise ValueError('Invalid option type')

    # Check market price in price interval 
    if not (low_bound<mkt_price<high_bound):
        raise ValueError('market price outside valid BS put bounds')

    # set interval of sigma
    low=1e-8
    high=1
        
    f_low=bs.bs_calculator(S0,K,T,r,q,low,option_type=option_type)-mkt_price
    if f_low>=0:
        raise ValueError('Lower bound does not bracket the root')

    while bs.bs_calculator(S0,K,T,r,q,high,option_type=option_type)-mkt_price<0:
        high*=2
        if high>max_sigma:
            raise ValueError('Invalid sigma upper bound')
    
    history=[]

    for i in range(max_iteration):
        mid=(low+high)/2
        history.append(mid)
        error=bs.bs_calculator(S0,K,T,r,q,mid,option_type=option_type)-mkt_price
        if (high-low)/2<sigma_tol or abs(error)<price_tol:
            return mid,i+1,history
        if error>0:
            high=mid 
        else:
            low=mid 
    
    raise RuntimeError('Bisection did not converge')


def Newton(S0,K,T,r,q,mkt_price,option_type='put',initial_sigma=0.3,tol=1e-8,max_iteration=50):
    '''
    mkt_price: market price
    initial_sigma: start point of sigma
    tol: acceptable error range
    max_iteration: acceptable maximum times of iteration
    =================
    return: float,int,[float]
        implied volatility from Newton's method, number of iterations, iteration history of sigma
    '''

    if initial_sigma<=0:
        raise ValueError('Invalid negative initial sigma')
    
    sigma=initial_sigma
    history=[]

    for i in range(max_iteration):
        error=bs.bs_calculator(S0,K,T,r,q,sigma,option_type=option_type)-mkt_price
        history.append(sigma)
        if abs(error)<tol:
            return sigma,i+1,history

        *_, vega=gk.bs_greeks(S0,K,T,r,q,sigma,option_type=option_type)
        if not np.isfinite(vega) or vega<=1e-12:
            raise RuntimeError('vega be too small or invalid')

        new_sigma=sigma-error/vega
        if not(math.isfinite(new_sigma) and new_sigma>0):
            raise RuntimeError('Invalid Newton iteration')

        sigma=new_sigma

    raise RuntimeError("Newton's method did not converge")


def Newton_tracking(S0,K,T,r,q,mkt_price,option_type='put',initial_sigma=0.3,tol=1e-8,max_iteration=50):
    '''
    mkt_price: market price
    initial_sigma: start point of sigma
    tol: acceptable error range
    max_iteration: acceptable maximum times of iteration
    =================
    return: float,int,str
        implied volatility from Newton's method, number of iterations, status of outcomes
    '''

    if initial_sigma<=0:
        raise ValueError('Invalid negative initial sigma')
    
    sigma=initial_sigma
    
    for i in range(max_iteration):
        error=bs.bs_calculator(S0,K,T,r,q,sigma,option_type=option_type)-mkt_price
        
        if abs(error)<tol:
            return sigma,i+1,'price_tolerance'

        *_, vega=gk.bs_greeks(S0,K,T,r,q,sigma,option_type=option_type)
        if not np.isfinite(vega) or vega<=1e-12:
            return np.nan,i+1,'vega_too_small'

        new_sigma=sigma-error/vega
        if not(math.isfinite(new_sigma) and new_sigma>0):
            return np.nan,i+1,'invalid_step'

        sigma=new_sigma

    return np.nan,max_iteration,'max_iter'
    

def iv_calculator(S0,K,T,r,q,price,option_type):
    '''
    Calculate implied volatiltiy based on continuous dividend yield option
    ================================================
    r: float
        risk free rate
    q: float 
        continuous dividend yield
    option_type: str 'call'/'put'
        confirm option type whether call or put
    ================================================
    return: float
        implied volatility
    '''
    if (S0<=0) or (K<=0) or (T<=0) or (price<=0) or not np.isfinite(price):
        return np.nan

    # Check price interval
    discount_S=S0*np.exp(-q*T)
    discount_K=K*np.exp(-r*T)
    if option_type=='call':
        lower=max(discount_S-discount_K,0)
        upper=discount_S
    elif option_type=='put':
        lower=max(discount_K-discount_S,0)
        upper=discount_K
    else:
        return np.nan

    if not (lower<price<upper):
        return np.nan 

    def objective(sigma):
        return bs.bs_calculator(S0,K,T,r,q,sigma,option_type)-price

    sigma_min=1e-6
    sigma_max=5

    if not objective(sigma_min)*objective(sigma_max)<0:
        return np.nan 

    try:
        return brentq(objective,sigma_min,sigma_max,xtol=1e-10)
    
    except ValueError:
        return np.nan 