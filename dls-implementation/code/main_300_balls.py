import numpy as np
from scipy.optimize import minimize, root_scalar
import pandas as pd
import matplotlib.pyplot as plt


def load_data(path):
    df = pd.read_csv(path)
    return df


def preprocess_data(data):
    data = data.copy()
    data = data[(data['season'] >= 2000) & (data['season'] <= 2006)]
    filtered_columns = ['innings','balls_remaining', 'wickets_down','runs_remaining']
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



def calculate_loss(params, data):
    total_loss = 0
    count = 0
    for ball in range(1, 301):
        for wicket in range(10):
            temp_data = data[(data['balls_remaining'] == ball) & (data['wickets_down'] == wicket)]
            z_val = Z_std_u_w(params, ball/6, wicket)

            if temp_data.shape[0] > 0:
                count += 1
                loss = (temp_data['runs_remaining'].values - z_val)**2
                total_loss += np.mean(loss)

    return total_loss / count


def g_equation(params, lmda):
    c1, c2, c3, c4 = params
    u = 50

    alpha_lmda = -1/(1 + c1 * (lmda-1)* np.exp(-c2*(lmda - 1)))
    beta_lmda = -c3 * (lmda - 1) * np.exp(-c4 * (lmda - 1))
    exponent = -(1 + alpha_lmda + beta_lmda * u)

    g = (u/50) ** exponent
    return g - 1


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

    input_path = '../data/modified_ODI_data.csv'

    df_data = load_data(input_path)
    cleaned_data = preprocess_data(df_data)
    print(f"{cleaned_data.shape = }")

    #bounds = [(200, None), (1e-3, None),(1e-3, None),(1e-3, None),(1e-3, None)]
    initial_guess = [200,.1,.1,.1,.1]

    result = minimize(calculate_loss, initial_guess, args = cleaned_data, method='Powell')
    print(f"{result = }")     

