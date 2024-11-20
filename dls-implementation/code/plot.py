from inference import *
from main import calculate_lambda,Z_dls,Z_pro_u_w,Z_std_u_w
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    path='../data/modified_ODI_data.csv'
    data=pd.read_csv(path)
    data=data[['over','ball','team_score','total_runs','wickets_down','target_score','total_balls','balls_remaining','season']]
    
    # dls_data=data[(data['total_runs']>300 ) & ( data['target_score']!=-1) & (data['wickets_down']==0) ]
    # dl_data=data[ (data['total_runs']<300) & (data['total_runs']>=200) & (data['target_score']!=-1) & (data['wickets_down']==0) ]

    dls_data=data[(data['total_runs']>300 ) & ( data['target_score']!=-1) & (data['wickets_down']==0) &(data['season']>2010) &(data['season']<2016)  ]
    dl_data=data[ (data['total_runs']<300) & (data['total_runs']>=200) & (data['target_score']!=-1) & (data['wickets_down']==0)&(data['season']>2010) &(data['season']<2016) ]
    print(dls_data.shape)
    print(dl_data.shape)
    print(dls_data['total_runs'].mean())
    print(dl_data['total_runs'].mean())
    dls_act=[]
    dl_act=[]
    for i in range(0,301):
        dls_i=dls_data[ (dls_data['total_balls']-dls_data['balls_remaining'])==i]
        dls_runstogo= ((dls_i['total_runs']-dls_i['team_score']) / dls_i['total_runs'])
        dls_act.append( (dls_runstogo.mean(),(300-i)/6))
        
        dl_i=dl_data[(dl_data['total_balls']-dl_data['balls_remaining'])==i]
        dl_runstogo=((  dl_i['total_runs']-dl_i['team_score']) / dl_i['total_runs'] )
        dl_act.append((dl_runstogo.mean(),(300-i)/6))

    # Z_dls_params = [1.101e+00,  1.540e+00, -1.404e-01 , 1.083e+01]
    # Z_std_params = [ 2.609e+02,  1.204e-02, -1.252e-01,  2.624e-03,  1.758e-04]
    Z_dls_params=[ 1.150e+00,  2.113e+00 , 9.041e-02, -4.502e-01]
    Z_std_params=[ 2.609e+02,  1.204e-02, -1.252e-01,  2.624e-03,  1.758e-04]
    lamda= calculate_lambda(Z_std_params,331)
    
    arr_dls=[]
    arr_dl=[]
    
    for i in np.linspace(0,50,200):
        arr_dls.append((Z_dls(Z_dls_params,Z_std_params,50-i,0,lamda)/Z_dls(Z_dls_params,Z_std_params,50,0,lamda),50-i))
        arr_dl.append((Z_std_u_w(Z_std_params,50-i,0)/Z_std_u_w(Z_std_params,50,0),50-i))

    y1, x1 = zip(*arr_dls)
    y2, x2 = zip(*arr_dl)
    y3, x3 = zip(*dls_act)
    y4, x4 = zip(*dl_act)
    plt.figure(figsize=(10, 6))

    plt.plot(x1, y1, label='dls_pred  ')
    plt.plot(x2, y2, label='dl_pred  ')
    plt.scatter(x3, y3, label='Avg_325 ', marker='o')
    plt.scatter(x4, y4, label='Avg_250', marker='o')

    # Adding labels and legend
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plot of Four Lists of (Y, X) Pairs')
    plt.legend()
    plt.savefig('plot.png', format='png', dpi=300)
    # Show the plot
    
        
