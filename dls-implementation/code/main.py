import numpy as np
from scipy.optimize import minimize, root_scalar
import pandas as pd


def load_data(path):
    df = pd.read_csv(path)
    return df


def preprocess_data(data):
    
    data = data[(data['Runs.Remaining'] >= 0) & (data['Wickets.in.Hand'] > 0)]
    data = data.copy()

    data['Overs.Remaining'] = 50 - data['Over']
    data['Wickets.Down'] = 10 - data['Wickets.in.Hand']
    filtered_columns = ['Innings','Overs.Remaining', 'Wickets.Down','Runs.Remaining']
    data = data[filtered_columns]
    data.dropna(inplace=True, axis=0)

    return data


def Z_std_u_w(params, u, w):
    G_50, b, a1, a2, a3 = params

    Z_0 = G_50/(1 - np.exp(-50*b))
    F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)

    Z_std = Z_0 * F_w * (1 - np.exp((-u*b)/F_w))

    return Z_std 


def Z_pro_u_w(params, u, w, lmda):
    G_50, b, a1, a2, a3 = params

    n0 = 5

    Z_0 = G_50/(1 - np.exp(-50*b))
    F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)

    n_w = n0 * F_w

    return Z_0 * F_w * (lmda**(n_w+1)) * (1 - np.exp((-u*b)/(F_w * (lmda**n_w)))) 


def Z_dls(params, args, u, w, lmbda):
    c1, c2, c3, c4 = params
    G_50, b, a1, a2, a3 = args
    
    Z_0 = G_50/(1-np.exp(-50*b))
    F_w = 1 + a1*w + a2*(w**2) + a3*(w**3)

    n_w = 5 * F_w

    alpha = -1 / (1+c1*(lmbda-1)*np.exp(-c2*(lmbda-1)))
    beta = -c3*(lmbda-1)*np.exp(-c4*(lmbda-1))

    g_u_lmbda = 0

    if u != 0:
        g_u_lmbda = np.power(u/50, -(1+alpha+beta*u))

    # if u != 50:
    #     print(f"{Z_0 * F_w * (lmbda**(n_w+1)) = }")
    #     print(f"{(1-np.exp((-u*b*g_u_lmbda)/(F_w * (lmbda**n_w)))) = }")

    return Z_0 * F_w * (lmbda**(n_w+1)) * (1-np.exp((-u*b*g_u_lmbda)/(F_w * (lmbda**n_w))))


def calculate_loss(params, data):
    total_loss = 0
    count = 0
    for over in range(1, 50):
        for wicket in range(10):
            temp_data = data[(data['Overs.Remaining'] == over) & (data['Wickets.Down'] == wicket)]
            z_val = Z_std_u_w(params, over, wicket)

            if temp_data.shape[0] > 0:
                count += 1
                loss = (temp_data['Runs.Remaining'].values - z_val)**2
                total_loss += np.mean(loss)

    return total_loss / count


def lmda_equation(x, params, team1_score):

    G_50, b, a1, a2, a3 = params

    u, w, n0 = 50, 0, 5

    Z_0 = G_50/(1 - np.exp(-50*b))
    F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)
    n_w = n0 * F_w

    return Z_0 * F_w * (x**(n_w+1)) * (1 - np.exp((-u*b)/(F_w * (x**n_w)))) - team1_score


def calculate_lambda(params, team1_score):

    solution = root_scalar(lambda x: lmda_equation(x, params, team1_score), bracket=[1,100])

    return solution.root


if __name__ == '__main__':

    input_path = '../data/04_cricket_1999to2011.csv'

    df_data = load_data(input_path)

    cleaned_data = preprocess_data(df_data)

    print(f"{cleaned_data.shape = }")

    #bounds = [(200, None), (1e-3, None),(1e-3, None),(1e-3, None),(1e-3, None)]
    initial_guess = [200,0.1,0.1,0.1,0.1]

    result = minimize(calculate_loss, initial_guess, args = cleaned_data, method='L-BFGS-B')
    print(f"{result = }")     

