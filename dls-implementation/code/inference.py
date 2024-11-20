import numpy as np

from main import Z_dls, lmda_equation, calculate_lambda, Z_std_u_w, Z_pro_u_w


def inference_ODI(Z_dls_params, Z_std_params,  team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method = 'DLS'):
    
    G_50, _, _, _, _ = Z_std_params

    u_1 = 50 - team2_overs_completed
    u_2 = 50 - (team2_overs_completed + lost_overs_due_to_rain)

    lambda_val = None

    if team1_score <= G_50:
        lambda_val = 1
    else:
        lambda_val = calculate_lambda(Z_std_params, team1_score)
        # print(f"Lambda value: {lambda_val}")
        # print(f"Equation solver value: {lmda_equation(lambda_val, Z_std_params, team1_score)}")

    P_u1_w, P_u2_w = None, None
    if method == 'DL-STD':
        P_u1_w = Z_std_u_w(Z_std_params, u_1, team2_wicktes_down)/Z_std_u_w(Z_std_params, 50, 0)
        P_u2_w = Z_std_u_w(Z_std_params, u_2, team2_wicktes_down)/Z_std_u_w(Z_std_params, 50, 0)
    
    elif method == 'DL-PRO':
        P_u1_w = Z_pro_u_w(Z_std_params, u_1, team2_wicktes_down, lambda_val)/Z_pro_u_w(Z_std_params, 50, 0, 1)
        P_u2_w = Z_pro_u_w(Z_std_params, u_2, team2_wicktes_down, lambda_val)/Z_pro_u_w(Z_std_params, 50, 0, 1)
    
    else:
        P_u1_w = Z_dls(Z_dls_params,Z_std_params, u_1, team2_wicktes_down, lambda_val)/Z_dls(Z_dls_params,Z_std_params, 50, 0, lambda_val)
        P_u2_w = Z_dls(Z_dls_params,Z_std_params, u_2, team2_wicktes_down, lambda_val)/Z_dls(Z_dls_params,Z_std_params, 50, 0, lambda_val)

    R_1 = 1 * 100
    R_2 = (1 - P_u1_w + P_u2_w)*100

    team2_target = None

    if R_2 <= R_1:
        team2_target = team1_score * R_2/R_1
    else:
        team2_target = team1_score + G_50*(R_2 - R_1)

    return team2_target, int(np.ceil(team2_target))


def inference_T20I(Z_dls_params, Z_std_params,  team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain,method = 'DLS'):

    initial_run = Z_std_u_w(Z_std_params, 30, 0)
    new_team1_score = initial_run + team1_score
    new_team2_overs_completed = 30 + team2_overs_completed
    predicted_par_ODI, _ = inference_ODI(Z_dls_params, Z_std_params, new_team1_score, team2_wicktes_down, new_team2_overs_completed, lost_overs_due_to_rain, method)

    return (predicted_par_ODI - initial_run), int(np.ceil(predicted_par_ODI - initial_run))


if __name__ == '__main__':

    Z_dls_params = [1.101e+00,  1.540e+00, -1.404e-01 , 1.083e+01]
    Z_std_params = [ 2.609e+02,  1.204e-02, -1.252e-01,  2.624e-03,  1.758e-04]

    # ODI Scenario
    print("\nODI Scenario: ")
    team1_score = 350
    team2_wicktes_down = 0
    team2_overs_completed = 10
    lost_overs_due_to_rain = 5

    print("\nDL-STD: ",inference_ODI(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DL-STD'))
    print("\nDL-PRO: ",inference_ODI(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DL-PRO'))
    print("\nDLS: ",inference_ODI(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DLS'))
    #print(inference_ODI_final(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DL'))


    #T20I scenario
    print("\nT20I Scenario: ")
    team1_score = 219
    team2_wicktes_down = 3
    team2_overs_completed = 10
    lost_overs_due_to_rain = 5
    print("\nDL-STD: ",inference_T20I(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DL-STD'))
    print("\nDL-PRO: ",inference_T20I(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DL-PRO'))
    print("\nDLS: ",inference_T20I(Z_dls_params, Z_std_params, team1_score, team2_wicktes_down, team2_overs_completed, lost_overs_due_to_rain, method='DLS'))
    

