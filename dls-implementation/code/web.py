import streamlit as st
from inference import inference


    
st.write('''
         # Duckworth Lewis(Standard Edition) Method
         ''')

st.title("Cricket Target Score Prediction")

# Input fields for parameters and values
st.header("Input Parameters")
# G_50 = st.number_input("G_50", value=200.0)
# a = st.number_input("Parameter a", value=0.0)
# b = st.number_input("Parameter b", value=0.0)
# c = st.number_input("Parameter c", value=0.0)
# d = st.number_input("Parameter d", value=0.0)

team1_score = st.number_input("Team 1 Score", value=250)
team2_wicktes_down = st.number_input("Team 2 Wickets Down", value=2)
team2_overs_down = st.number_input("Team 2 Overs Down", value=20)
lost_overs_due_to_rain = st.number_input("Lost Overs Due to Rain", value=5)

# Compute the target score
# params = [G_50, a, b, c, d]
params = [ 2.808e+02, -3.020e-02, 4.439e-01, -4.038e-01, 1.082e-01]

if st.button("Calculate Target"):
    target_score = inference(params, team1_score, team2_wicktes_down, team2_overs_down, lost_overs_due_to_rain)
    st.success(f"Predicted Target Score for Team 2: {target_score}")