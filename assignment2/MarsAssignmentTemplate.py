import numpy as np
from functions import *

if __name__ == "__main__":

    # Import oppositions data from the CSV file provided
    data = np.genfromtxt(
        "../data/01_data_mars_opposition_updated.csv",
        delimiter=",",
        skip_header=True,
        dtype="int",
    )
    
    # Extract times from the data in terms of number of days.
    # "times" is a numpy array of length 12. The first time is the reference
    # time and is taken to be "zero". That is times[0] = 0.0
    
    times = get_times(data)
    assert len(times) == 12, "times array is not of length 12"


    oppositions = get_oppositions(data)
    assert len(oppositions) == 12, "oppositions array is not of length 12"

    # s=0.5240799553026627
    # r= 8.908163265306122
    # c = 148.27586206896552
    # e1 = 1.6530612244897958
    # e2 = 148.79310344827587
    # z = 55.86206896551724
    # print(MarsEquantModel(c,r,e1,e2,z,s,times,oppositions))
    # errors = array([0.0385409 , 0.01037564, 0.01638428, 0.01167598, 0.03651423,
    #    0.03044545, 0.04764352, 0.0149071 , 0.00210762, 0.00658089,
    #    0.047372  , 0.01033101]), maxError = 0.047643517415735914
    r, s, c, e1, e2, z, errors, maxError = bestMarsOrbitParams(times, oppositions)
    
    assert max(list(map(abs, errors))) == maxError, "maxError is not computed properly!"
    print("Fit parameters: r = {:.4f}, s = {:.4f}, c = {:.4f}, e1 = {:.4f}, e2 = {:.4f}, z = {:.4f}".format(r, s, c, e1, e2, z))
    print("The maximum angular error = {:2.4f}".format(maxError))

   
    # print(MarsEquantModel(c,r,e1,e2,z,s,times,oppositions))
    plot_mars_orbit(r, s, z, e1, e2, c, times, oppositions)