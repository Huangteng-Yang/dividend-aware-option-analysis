import numpy as np

def mc_put(S0,K,T,r,q,sigma,N):
    '''
    Calculate option price based on Monte Carlo simulation
    ================================================
    r: float
        risk free rate
    q: float 
        continuous dividend yield
    N: int
        number of random sample path
    ================================================
    return: list of floats and interval
        European put price via monte carlo simulation, std of sample path, se of estimator, confidence interval
    '''

    Z=np.random.randn(N)
    ST=S0*np.exp((r-q-1/2*sigma**2)*T)*np.exp((sigma*np.sqrt(T))*Z)
    payoff=np.maximum(K-ST,0)

    # Monte Carlo price
    price=np.exp(-r*T)*np.mean(payoff)

    #Sample standard deviation
    std=np.std(payoff,ddof=1)

    # Standard error of the estimator
    se=np.exp(-r*T)*std/np.sqrt(N)

    # 95% confidence interval
    CI=(price-1.96*se,price+1.96*se) 

    return price,std,se,CI 
    