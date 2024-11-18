import os
import pickle
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Union
import scipy
import scipy.optimize 


if not os.path.exists('../models'):
    os.makedirs('../models')
if not os.path.exists('../plots'):
    os.makedirs('../plots')


class DLModel:
    """
        Model Class to approximate the Z function as defined in the assignment.
    """

    def __init__(self):
        """Initialize the model."""
        self.Z0 = [None] * 10
        self.L = None
        self.L_list=[None]*10
        self.pre_z=[None] * 10
    
    def get_predictions(self, X, Z_0=None, w=10, L=None) -> np.ndarray:
        """Get the predictions for the given X values.

        Args:
            X (np.array): Array of overs remaining values.
            Z_0 (float, optional): Z_0 as defined in the assignment.
                                   Defaults to None.
            w (int, optional): Wickets in hand.
                               Defaults to 10.
            L (float, optional): L as defined in the assignment.
                                 Defaults to None.

        Returns:
            np.array: Predicted score possible
        """
        a=-L*X
        b=a/Z_0

        Z=Z_0*(1-np.exp(b))
        return Z

    def calculate_loss(self, Params, X, Y, w=10) -> float:
        """ Calculate the loss for the given parameters and datapoints.
        Args:
            Params (list): List of parameters to be optimized.
            X (np.array): Array of overs remaining values.
            Y (np.array): Array of actual average score values.
            w (int, optional): Wickets in hand.
                               Defaults to 10.

        Returns:
            float: Mean Squared Error Loss for the model parameters 
                   over the given datapoints.
        """
        
        Z_0,L= Params
        # print("in Loss",Z_0,L)
        y_pred= self.get_predictions(X,Z_0,w=w,L=L)
# (y'+1) log((y'+1)/(y+1)) - y' + y,
        loss=(y_pred+1)*np.log((y_pred+1)/(Y+1)  ) -y_pred + Y
        return np.mean(loss)
#     def calculate_loss_l(self, L, X, Y, w=10,) -> float:
#         y_pred= self.get_predictions(X,Z_0,w=w,L=L)
# # (y'+1) log((y'+1)/(y+1)) - y' + y,
#         loss=(y_pred+1)*np.log((y_pred+1)/(Y+1)  ) -y_pred + Y
#         return np.mean(loss)
    
    def save(self, path):
        """Save the model to the given path.

        Args:
            path (str): Location to save the model.
        """
        with open(path, 'wb') as f:
            pickle.dump((self.L, self.Z0), f)
    
    def load(self, path):
        """Load the model from the given path.

        Args:
            path (str): Location to load the model.
        """
        with open(path, 'rb') as f:
            (self.L, self.Z0) = pickle.load(f)


def get_data(data_path) -> Union[pd.DataFrame, np.ndarray]:
    """
    Loads the data from the given path and returns a pandas dataframe.

    Args:
        path (str): Path to the data file.

    Returns:
        pd.DataFrame, np.ndarray: Data Structure containing the loaded data
    """
    return pd.read_csv(data_path)
    


def preprocess_data(data: Union[pd.DataFrame, np.ndarray]) -> Union[pd.DataFrame, np.ndarray]:
    """Preprocesses the dataframe by
    (i)   removing the unnecessary columns,
    (ii)  loading date in proper format DD-MM-YYYY,
    (iii) removing the rows with missing values,
    (iv)  anything else you feel is required for training your model.

    Args:
        data (pd.DataFrame, nd.ndarray): Pandas dataframe containing the loaded data

    Returns:
        pd.DataFrame, np.ndarray: Datastructure containing the cleaned data.
    """
    data=data[['Over','Runs.Remaining','Wickets.in.Hand','Innings']]
    # print(data[:5])
    
    data = data[
        (data['Runs.Remaining'] > 0) & 
        (data['Wickets.in.Hand'] > 0) & 
        (data['Innings'] == 1)]
    
    
    data['Overs.REMAINING'] = 50 - data['Over']
    data = data[['Overs.REMAINING', 'Wickets.in.Hand', 'Runs.Remaining']]

    # print(data)
    data.dropna()
    return data


def train_model(data: Union[pd.DataFrame, np.ndarray], model: DLModel) -> DLModel:
    """Trains the model

    Args:
        data (pd.DataFrame, np.ndarray): Datastructure containg the cleaned data
        model (DLModel): Model to be trained
    """
    Z0_list=[]
    L_list=[]
    l_sum=0
    for i in range(1,11):
        data1=data[data['Wickets.in.Hand'] == i]
        x=data1['Overs.REMAINING'].to_numpy()
        y=data1['Runs.Remaining'].to_numpy()
        # print("shape",x.shape, y.shape)
        ans=scipy.optimize.minimize(model.calculate_loss,[3.0,1.0],args=(x,y,i),bounds=[(0.001, float("inf")),(0.001, float("inf"))])
        Z0_list.append(ans.x[0])
        l_sum+=(ans.x[1])*len(x)
        L_list.append(ans.x[1])
        # print(ans.x)
        # break
    # print(data.shape)
    l_sum=l_sum/data.shape[0]
    # print(l_sum)
    model.pre_z=Z0_list
    model.L=l_sum
    model.L_list=L_list

    Z0_list=[]
    L_list=[]
    
    for i in range(1,11):
        data1=data[data['Wickets.in.Hand'] == i]
        x=data1['Overs.REMAINING'].to_numpy()
        y=data1['Runs.Remaining'].to_numpy()
        # print("shape",x.shape, y.shape)
        # ans=scipy.optimize.minimize(model.calculate_loss_l,[3.0,1.0],args=(x,y,i),bounds=[(0.001, float("inf")),(0.001, float("inf"))])
        ans=scipy.optimize.minimize(lambda z0: model.calculate_loss([z0, l_sum], x, x, i), 1, bounds=[(0.01, None)])
        Z0_list.append(ans.x[0])
        
        # print(ans.x)
        # break
    # print(data.shape)
    
    # print(l_sum)
    model.Z0=Z0_list
   
    # print("zpre",model.pre_z)
    # print("lpre",model.L_list)
    # print("zo",model.Z0)
    # print("lfinal",model.L)


    # print(Z0_list)
    # print(l_sum)
    # print(L_list)
    return model


def plot(model: DLModel, plot_path: str) -> None:
    """ Plots the model predictions against the number of overs
        remaining according to wickets in hand.

    Args:
        model (DLModel): Trained model
        plot_path (str): Path to save the plot
    """
    overs_array = np.arange(0, 50.25, 0.25)
    plot_score=[]
    for i in range(1,11):
        score=model.get_predictions(overs_array,model.Z0[i-1],w=i,L=model.L)
        plot_score.append(score)

    
    _, axs = plt.subplots(5, 2, figsize=(12, 15))  
    
    axs = axs.flatten()  
    
   
    for i, (score, ax) in enumerate(zip(plot_score, axs), start=1):
        ax.plot(overs_array, score, label=f'{i} wickets in hand')
        ax.set_xlabel('Overs remaining')
        ax.set_ylabel('Predicted Score')
        ax.set_title(f'Wickets in hand: {i}')
        ax.legend()

    
    plt.tight_layout()
    
    
    plt.savefig(plot_path+'1.png')
    plt.show()

    overs_array = np.arange(0, 50.25, 0.25)
    plot_score=[]
    for i in range(1,11):
        score=model.get_predictions(overs_array,model.Z0[i-1],w=i,L=model.L_list[i-1])
        plot_score.append(score)

    _, axs = plt.subplots(5, 2, figsize=(12, 15))  
    
    axs = axs.flatten()  
    
   
    for i, (score, ax) in enumerate(zip(plot_score, axs), start=1):
        ax.plot(overs_array, score, label=f'{i} wickets in hand')
        ax.set_xlabel('Overs remaining')
        ax.set_ylabel('Predicted Score')
        ax.set_title(f'Wickets in hand: {i}')
        ax.legend()

    
    plt.tight_layout()
    
    
    plt.savefig(plot_path+'2.png')
    plt.show()
def print_model_params(model: DLModel) -> List[float]:
    '''
    Prints the 11 (Z_0(1), ..., Z_0(10), L) model parameters

    Args:
        model (DLModel): Trained model
    
    Returns:
        array: 11 model parameters (Z_0(1), ..., Z_0(10), L)

    '''
    print(model.Z0+[model.L])
    return model.Z0+[model.L]
    
    
    


def calculate_loss(model: DLModel, data: Union[pd.DataFrame, np.ndarray]) -> float:
    '''
    Calculates the normalised squared error loss for the given model and data

    Args:
        model (DLModel): Trained model
        data (pd.DataFrame or np.ndarray): Data to calculate the loss on
    
    Returns:
        float: Normalised squared error loss for the given model and data
    '''
    loss_array=[]
    
    for i in range(1,11):
        data1=data[data['Wickets.in.Hand'] == i]
        x=data1['Overs.REMAINING'].to_numpy()
        y=data1['Runs.Remaining'].to_numpy()
        temp=(model.Z0[i-1],model.L)
        loss=model.calculate_loss(temp,x,y,i)
        loss_array.append(loss)
        # print("shape",x.shape, y.shape)
        # ans=scipy.optimize.minimize(model.calculate_loss,[3.0,1.],args=(x,y,i),bounds=[(0.001, float("inf")),(0.001, float("inf"))])
    # print(np.mean(loss_array))
    return np.mean(loss_array)


def main(args):
    """Main Function"""

    data = get_data(args['data_path'])  # Loading the data
    print("Data loaded.")
    
    # Preprocess the data
    data = preprocess_data(data)
    print("Data preprocessed.")
    
    model = DLModel()  # Initializing the model

    model = train_model(data, model)  # Training the model
    model.save(args['model_path'])  # Saving the model
    
    plot(model, args['plot_path'])  # Plotting the model
    
    # Printing the model parameters
    print_model_params(model)

    # Calculate the normalised squared error
    calculate_loss(model, data)


if __name__ == '__main__':
    args = {
        "data_path": "../data/04_cricket_1999to2011.csv",
        "model_path": "../models/model.pkl",  # ensure that the path exists
        "plot_path": "../plots/plot.png",  # ensure that the path exists
    }
    main(args)
