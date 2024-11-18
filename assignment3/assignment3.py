import pandas as pd
import numpy as np
from numpy.linalg import matrix_rank, pinv
from scipy.stats import f
import matplotlib.pyplot as plt

df=pd.read_csv('../data/Raw Data_GeneSpring.txt',sep='\t')
# print(df.head)
# total_columns = len(df.columns)
# print(f"Total number of columns: {total_columns}")

X =df.iloc[:,1:-3]

X = X.to_numpy()
X=np.exp2(X)
# print(X.shape)
# total_columns = len(df_cleaned.columns)
# print(f"Total number of columns: {total_columns}")

m = [1, 0, 1, 0]
fe = [1, 0, 0, 1]
s = [0, 1, 1, 0]
ns =[0, 1, 0, 1]

m1 = np.tile(m, (12, 1))
f1 = np.tile(fe, (12, 1))
s1 = np.tile(s, (12, 1))
ns1 = np.tile(ns, (12, 1))

N = np.vstack([m1, f1, s1, ns1])

# print(n_final)

ms =  [1, 0, 0, 0]
mns = [0, 1, 0, 0]
fs =  [0, 0, 1, 0]
fns = [0, 0, 0, 1]

ms1 =  np.tile(ms, (12, 1))
mns1 = np.tile(mns, (12, 1))
fs1 =  np.tile(fs, (12, 1))
fns1 = np.tile(fns, (12, 1))

D=np.vstack([ms1, mns1, fs1 , fns1])

n=48
rank_n=matrix_rank(N)
rank_d=matrix_rank(D)
df1 = rank_d - rank_n  
df2 = n - rank_d  
# print(rank_d)
# print(rank_n)

const= (1 / (rank_d - rank_n)) / (1 / (48 - rank_d))
# print(const)

I = np.eye(48)

N_proj = I - (np.matmul(np.matmul(N, pinv(np.matmul(N.T, N))), N.T))
D_proj = I - (np.matmul(np.matmul(D, pinv(np.matmul(D.T, D))), D.T))
# print(N_proj.shape)
# print(D_proj.shape)
f_scores = np.zeros(X.shape[0])
# print(X.shape)
for i in range(X.shape[0]):
    row =np.array( X[i, :])
    # print(row.shape)
    numerator = np.matmul(np.matmul(row, N_proj), row.T)  
    denominator = np.matmul(np.matmul(row, D_proj), row.T)
    
    f_scores[i]=(const*((numerator/denominator)-1))
p_values = 1 - f.cdf(f_scores, df1, df2)
# print(p_values.shape)
plt.figure(figsize=(10, 6))
plt.hist(p_values, bins=50, color='skyblue', edgecolor='black')
plt.title('Histogram of p-values')
plt.xlabel('p-value')
plt.ylabel('Frequency')
plt.grid(True)


plt.savefig('p_values_histogram.png')


plt.show()
print("Histogram plotted of p_values")
