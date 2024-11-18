import numpy as np
import multiprocessing as mp
from datetime import datetime
import math
import matplotlib.pyplot as plt

def get_times(data):
    times=[]
    for row in data:
        date0=datetime(year=row[0], month=row[1], day=row[2], hour=row[3], minute=row[4]) 
        date1=datetime(1580,11,18,1,31)
        diff=date0-date1
        days=diff.days
        sec=diff.seconds
        times.append(days+sec/(24*60*60))
    return np.array(times)
    
# Year,Month,Day,Hour,Minute,ZodiacIndex5,Degree,Minute,Second
def get_oppositions(data):
    oppositions=[]
    for row in data:
        deg=(row[5]*30) +row[6]+ (row[7]/60)+(row[8]/3600)
        oppositions.append(deg)
    return np.array(oppositions)


def polar_to_cartesian(r, theta):
    return r * math.cos(theta), r * math.sin(theta)


def calculate_slope(x, y):
    angle_degrees = math.degrees(math.atan2(y, x))
    if angle_degrees < 0:
        angle_degrees += 360
    return angle_degrees


def MarsEquantModel(c,r,e1,e2,z,s,times,oppositions):
    h,k=polar_to_cartesian(1,math.radians(c))
    x_1,y_1=polar_to_cartesian(e1,math.radians(e2))
    
    error=[]
    for i in range(12):
        slope=(z+s*times[i])%360
        m=math.tan(math.radians(slope))
        c=y_1-m*x_1
        
        A = 1 + m**2
        B = -2 * h + 2 * m * (c - k)
        C = h**2 + (c - k)**2 - r**2
        
        discriminant = B**2 - 4*A*C
        
        if discriminant < 0:
            print("No real intersection points")
            return []  
        elif discriminant == 0:
            X = -B / (2 * A)
            Y = m * X + c
                        
        else:
            sqrt_discriminant = np.sqrt(discriminant)
            x1 = (-B + sqrt_discriminant) / (2 * A)
            x2 = (-B - sqrt_discriminant) / (2 * A)
            y1 = m * x1 + c
            y2 = m * x2 + c
            slope1=calculate_slope(x1,y1)
            slope2=calculate_slope(x2,y2)
            # print(f"inside mars equant, {slope1 = }, {slope2 = }")
            if(abs(slope1-oppositions[i]) < abs(slope2-oppositions[i])):
                X=x1
                Y=y1
            else:
                X=x2
                Y=y2

        error.append(abs(calculate_slope(X,Y)-oppositions[i]))
    return error,max(error)
# def MarsEquantModel2(c,r,e1,e2,z,s,times,oppositions):
#     #print(f"Mars Equant model called")
#     h,k=polar_to_cartesian(1,math.radians(c))
#     x_1,y_1=polar_to_cartesian(e1,math.radians(e2))
    
#     predicted_x=[]
#     predicted_y=[]
#     for i in range(12):
#         slope=(z+s*times[i])%360
#         m=math.tan(math.radians(slope))
#         c=y_1-m*x_1
        
#         A = 1 + m**2
#         B = -2 * h + 2 * m * (c - k)
#         C = h**2 + (c - k)**2 - r**2
        
#         discriminant = B**2 - 4*A*C
        
#         if discriminant < 0:
#             print("No real intersection points")
#             return []  
#         elif discriminant == 0:
#             X = -B / (2 * A)
#             Y = m * X + c
                        
#         else:
#             sqrt_discriminant = np.sqrt(discriminant)
#             x1 = (-B + sqrt_discriminant) / (2 * A)
#             x2 = (-B - sqrt_discriminant) / (2 * A)
#             y1 = m * x1 + c
#             y2 = m * x2 + c
#             slope1=calculate_slope(x1,y1)
#             slope2=calculate_slope(x2,y2)
#             # print(f"inside mars equant, {slope1 = }, {slope2 = }")
#             if(abs(slope1-oppositions[i]) < abs(slope2-oppositions[i])):
#                 X=x1
#                 Y=y1
#             else:
#                 X=x2
#                 Y=y2
#         predicted_x.append(X)
#         predicted_y.append(Y)
        
#     return predicted_x,predicted_y

def bestOrbitInnerParams_only_e2_z(params):
    c, r, e1, s, times, oppositions = params
    best_e2, best_z = None, None
    min_max_error = float('inf')
    best_errors = []
    e2_values = np.linspace(140, 160, 60)  
    z_values = np.linspace(40, 60, 60)  
    for e2 in e2_values:
        for z in z_values:
            errors, max_error = MarsEquantModel(c, r, e1, e2, z, s, times, oppositions)
            if max_error < min_max_error:
                min_max_error = max_error
                best_e2, best_z = e2, z
                best_errors = errors

    return c, e1, best_e2, best_z, best_errors, min_max_error 

def bestOrbitInnerParams(r, s, times, oppositions):
    print(f"best orbit inner params called")
    best_c, best_e1, best_e2, best_z = None, None, None, None
    min_max_error = float('inf')

    c_values = np.linspace(140, 160, 60)  
    e1_values = np.linspace(1, 2, 50)  
    
    # Generate a list of all (e2, z) combinations to be evaluated
    params_list = [(c, r, e1, s, times, oppositions) for e1 in e1_values for c in c_values]

    # Use multiprocessing to evaluate these combinations in parallel
    with mp.Pool(processes=180) as pool:
        results = pool.map(bestOrbitInnerParams_only_e2_z, params_list)
    
    for c, e1, e2, z, errors, max_error in results:
        if max_error < min_max_error:
            min_max_error = max_error
            best_c, best_e1, best_e2, best_z = c, e1, e2, z
            best_errors = errors

    return best_c, best_e1, best_e2, best_z, best_errors, min_max_error

def bestS(r, times, oppositions):
    print(f"bestS")
    best_s = None
    min_max_error = float('inf')
    best_errors = None
    
    s_values = np.linspace((360-0.04)/687, (360+0.04)/687, 50)
    
    for s in s_values:
        c, e1, e2, z, errors, max_error = bestOrbitInnerParams(r, s, times, oppositions)
        if (max(errors) != max_error):
            print(f"my error : {errors = }, {max_error = }")
        if max_error < min_max_error:
            min_max_error = max_error
            best_s = s
            best_errors = errors

    return best_s, best_errors, min_max_error

def bestR(s, times, oppositions):
    print(f"best R called")
    best_r = None
    min_max_error = float('inf')
    best_errors = None
    r_values = np.linspace(8.5, 9.5, 50)  
    

    for r in r_values:
        c, e1, e2, z, errors, max_error = bestOrbitInnerParams(r, s, times, oppositions)
        if max_error < min_max_error:
            min_max_error = max_error
            best_r = r
            best_errors = errors

    return best_r, best_errors, min_max_error

def bestMarsOrbitParams(times, oppositions):
    print(f"best mars orbit params called")
    r_guess=10
    s, errors, max_error = bestS(r_guess, times, oppositions)
    print(s, errors, max_error)
    
    r, errors, max_error = bestR(s, times, oppositions)
    c, e1, e2, z, errors, max_error = bestOrbitInnerParams(r, s, times, oppositions)
    return r, s, c, e1, e2, z, errors, max_error



def plot_mars_orbit(r, s, z, e1, e2, c, times, oppositions):
    
    fig, ax = plt.subplots()
    h, k = polar_to_cartesian(1, math.radians(c))

    circle = plt.Circle((h, k), r, color='blue', fill=False, linestyle='dashed')
    ax.add_artist(circle)

    original_x, original_y = [], []
    
    for i in range(len(oppositions)):
        
        m=math.tan(math.radians(oppositions[i]%360))
        c2=0
        
        A = 1 + m**2
        B = -2 * h + 2 * m * (c2 - k)
        C = h**2 + (c2 - k)**2 - r**2
        discriminant = B**2 - 4*A*C
        if discriminant < 0:
            print("No real intersection points")
            return []  
        elif discriminant == 0:
            X=-B / (2 * A)
            original_x.append(X )
            original_y.append(m * X + c2)
        else:
            sqrt_discriminant = np.sqrt(discriminant)
            x1 = (-B + sqrt_discriminant) / (2 * A)
            x2 = (-B - sqrt_discriminant) / (2 * A)
            y1 = m * x1 + c2
            y2 = m * x2 + c2
            slope1=calculate_slope(x1,y1)
            slope2=calculate_slope(x2,y2)
            # print(f"actual, {slope1 = }, {slope2 = }")
            if(abs(slope1-oppositions[i]) <= abs(slope2-oppositions[i])):
                original_x.append(x1)
                original_y.append(y1)
            else:
                original_x.append(x2)
                original_y.append(y2)
    # Predicted points using the Mars Equant model
    predicted_x, predicted_y = [], []
    x_1,y_1=polar_to_cartesian(e1,math.radians(e2))
    for i in range(12):
        slope=(z+s*times[i])%360
        m=math.tan(math.radians(slope))
        c1=y_1-m*x_1
        
        A = 1 + m**2
        B = -2 * h + 2 * m * (c1 - k)
        C = h**2 + (c1 - k)**2 - r**2
        discriminant = B**2 - 4*A*C
        
        if discriminant < 0:
            print("No real intersection points")
            return []  
        elif discriminant == 0:
            X = -B / (2 * A)
            Y = m * X + c1
                        
        else:
            sqrt_discriminant = np.sqrt(discriminant)
            x1 = (-B + sqrt_discriminant) / (2 * A)
            x2 = (-B - sqrt_discriminant) / (2 * A)
            y1 = m * x1 + c1
            y2 = m * x2 + c1
            slope1=calculate_slope(x1,y1)
            slope2=calculate_slope(x2,y2)
            # print(f"{slope1 = }, {slope2 = }")
            # print(abs(slope1-oppositions[i]))
            # print(abs(slope2-oppositions[i]))
            if(abs(slope1-oppositions[i]) < abs(slope2-oppositions[i])):
                X=x1
                Y=y1
            else:
                X=x2
                Y=y2
        predicted_x.append(X)
        predicted_y.append(Y)
        
    # for i in range(0,12):
    #     print(predicted_x[i],predicted_y[i])
    #     print(original_x[i],original_y[i])

    plt.scatter(original_x, original_y, color='red', label='Original Oppositions', marker = 'o')
    
    # Plot predicted points
    plt.scatter(predicted_x, predicted_y, color='blue', label='Predicted Points', marker='x')
    plt.scatter(h,k, color='blue', label='Center', s=100, zorder=5)
    plt.scatter(0,0, color='orange', label='Sun', s=100, zorder=5)
    plt.scatter(x_1,y_1, color='green', label='equant', s=100, zorder=5)

    # Set axis limits and labels
    ax.set_xlim(-r - 1, r + 1)
    ax.set_ylim(-r - 1, r + 1)
    ax.set_aspect('equal', 'box')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Mars Orbit: Original vs Predicted Points')
    plt.legend(loc='best',fontsize=7)
    plt.grid(True)
    plt.savefig("plot1")


