import implied_volatility as iv

import numpy as np
import pandas as pd

def calculate_iv_surface(options,max_spread=0.3):
    df=options.copy()

    # Filter overfull spread price
    if max_spread is not None:
        df=df.loc[df['relative_spread']<=max_spread]

    df['iv_bs']=df.apply(
        lambda row:
        iv.iv_calculator(
            row['S0'],
            row['strike'],
            row['T'],
            row['r'],
            row['q'],
            row['mid_price'],
            row['option_type']
        ),axis=1
    )

    # Check numerical results
    total=len(df)
    success=df['iv_bs'].notna().sum()

    print('Contracts after spread filtering:',total)
    print('Sucessful IV calculations:',success)
    print('Success rate:',f'{success/total:.2%}' if total else 'N/A')

    return df 

