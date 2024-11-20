from main import Z_dls, calculate_lambda
import numpy as np
import scipy as sp
from scipy.optimize import minimize, root_scalar
import pandas as pd
import matplotlib.pyplot as plt


def load_data(path):
    df = pd.read_csv(path)
    return df


def preprocess_data(data, params):
    # data = data[(data['target_score'] >290)]
    data = data[(data['season'] >= 2010) & (data['season'] <= 2015)]
    data = data[(data['total_runs'] > 300) | (data['target_score'] >300) ]
    data = data.copy()
    runs = np.where(data['innings'] == 1,data['total_runs'], data['target_score'] - 1)  
    lambdas = []
    for i in runs:
        lambdas.append(calculate_lambda(params, i))
    data['match_factor'] =lambdas
    # print(data)
    # print("max lambda", data["Match.Factor"].max())
    filtered_columns = ['innings','balls_remaining', 'wickets_down','runs_remaining', 'total_runs', 'match_factor']
    data = data[filtered_columns]
    print(data.shape)
    data.dropna(inplace=True, axis=0)
    return data


# def Z_std_u_w(params, u, w):
#     G_50, b, a1, a2, a3 = params

#     Z_0 = G_50/(1 - np.exp(-50*b))
#     F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)

#     Z_std = Z_0 * F_w * (1 - np.exp((-u*b)/F_w))

#     return Z_std 


# def Z_pro_u_w(params, u, w, lmda):
#     G_50, b, a1, a2, a3 = params

#     n0 = 5

#     Z_0 = G_50/(1 - np.exp(-50*b))
#     F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)

#     n_w = n0 * F_w

#     return Z_0 * F_w * (lmda**(n_w+1)) * (1 - np.exp((-u*b)/(F_w * (lmda**n_w)))) 

# def Z_dls(params, args, u, w, lmbda):
#     c1, c2, c3, c4 = params
#     G_50, b, a1, a2, a3 = args
    
#     Z_0 = G_50/(1-np.exp(-50*b))
#     F_w = 1 + a1*w + a2*(w**2) + a3*(w**3)
#     n_w = 5 * F_w

#     alpha = -1 / (1+c1*(lmbda-1)*np.exp(-c2*(lmbda-1)))
#     beta = -c3*(lmbda-1)*np.exp(-c4*(lmbda-1))
#     g_u_lmbda = np.power(u/50, -(1+alpha+beta*u))
#     # print("In Z_dls")
#     # print(lmbda.shape)
#     # print(f"In Z_Dls {params=}, {args=}, {Z_0=}, {u=}, {w=}, {lmbda=},{g_u_lmbda=}, {n_w=}, {(-u*b*g_u_lmbda)/(F_w * (lmbda**n_w))}")
#     # print("pred", Z_0 * F_w * (lmbda**(n_w+1)) * (1-np.exp((-u*b*g_u_lmbda)/F_w * (lmbda**n_w))))

#     return Z_0 * F_w * (lmbda**(n_w+1)) * (1-np.exp((-u*b*g_u_lmbda)/(F_w * (lmbda**n_w))))


def calculate_loss1(params, args, data):
    total_loss = 0
    count = 0
    for ball in range(1, 301):
        for wicket in range(10):
            # Filter the data based on overs remaining and wickets down
            temp_data = data[(data['balls_remaining'] == ball) & (data['wickets_down'] == wicket)]
            # print(temp_data)

            if temp_data.shape[0] > 0:
                count += 1

                lambda_val = temp_data['match_factor'].to_numpy()
                predicted_value = Z_dls(params, args, ball/6, wicket, lambda_val)
                
                # print(f"{params=}, {args=}, {over=}, {wicket=}, {lambda_val=},{predicted_value=}")
                for i in predicted_value:
                    assert(i!=float('inf'))
                
                loss = (temp_data['runs_remaining'].values - predicted_value)**2

                total_loss += np.mean(loss)
        
    return total_loss/count

            # if temp_data.shape[0] > 0:
            #     for index, row in temp_data.iterrows():
        
            #         overs_remaining = row['Overs.Remaining']
            #         wickets_down = row['Wickets.Down']
            #         runs_remaining = row['Runs.Remaining']

            #         runs = row['Innings.Total.Runs'] if row['Innings'] == 1 else row['Target.Score'] - 1

            #         lmbda = calculate_lambda(args, runs)
            #         # print(lmbda)

            #         # Calculate the loss for each row using Z_dls function
            #         predicted_value = Z_dls(params, args, overs_remaining, wickets_down, lmbda)
            #         loss += (runs_remaining - predicted_value) ** 2
            #         count += 1

    #return loss / count 




# def equation(x, params, team1_score):

#     G_50, b, a1, a2, a3 = params

#     u, w, n0 = 50, 0, 5

#     Z_0 = G_50/(1 - np.exp(-50*b))
#     F_w = 1 + a1*w + a2 * (w**2) + a3 * (w**3)

#     n_w = n0 * F_w

#     return Z_0 * F_w * (x**(n_w+1)) * (1 - np.exp((-u*b)/(F_w * (x**n_w)))) - team1_score


# def calculate_lambda(params, team1_score):

#     solution = root_scalar(lambda x: equation(x, params, team1_score), bracket=[1,100])

#     return solution.root

    


if __name__ == '__main__':

    input_path = '../data/modified_ODI_data.csv'
    G_50,b, a1, a2, a3 =  2.609e+02,  1.204e-02, -1.252e-01,  2.624e-03,  1.758e-04
    args = (G_50,b, a1, a2, a3)
    df_data = load_data(input_path)

    cleaned_data = preprocess_data(df_data, args)
    
    # print("Clean")
    # print(cleaned_data)

    # print(f"{cleaned_data.shape = }")

    #bounds = [(200, None), (1e-3, None),(1e-3, None),(1e-3, None),(1e-3, None)]

    initial_guess = [-0.1, 1, -0.1, 1]
    # args = 2.808e+02, -3.020e-02,  4.439e-01, -4.038e-01,  1.082e-01
    # initial_guess = [200,0.1,0.1,0.1,0.1]
    # bounds = [(0.01, 10), (0.01, 10), (0.01, 10), (0.01, 10)]

    result = minimize(calculate_loss1, initial_guess, args = (args,cleaned_data), method='Nelder-Mead')
    print(f"{result = }")     













