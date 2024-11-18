# Data_analytics
---
# 🌌 **Assignment 2 -Analyzing Kepler's Data on Mars**  


This repository provides materials, datasets, and instructions for Analyzing Kepler's Data on Mars.  

---

## 📚 **Contents**
1. [Overview](#overview)
2. [Datasets](#datasets)
3. [Assignments](#assignments)
4. [Resources](#resources)
5. [Contact](#contact)  

---

## 🌍 **Overview**  
This module focuses on analyzing Mars's orbital data using trigonometric and regression techniques, including:  
- Exploring Mars opposition data to identify orbital parameters.  
- Computing orbital projections using triangulation data.  
- Applying regression and nonlinear optimization.  

---

## 📊 **Datasets**  

### **Mars Opposition Data**  
- **File:** `mars_opposition_data.csv`  
- **Description:** Contains Mars's heliocentric longitude, geocentric latitude, and Keplerian approximations during its opposition with the Sun.  
- **Key Details:**  
  - **Longitude Formula:** `s * 30 + Degree + Minute/60 + Second/3600` (degrees)  
  - Includes columns for date and Mars's mean longitude (Kepler's equant approximation).  

### **Triangulation Data**  
- **File:** `triangulation_data.csv`  
- **Description:** Provides longitudinal Earth-Mars observational data for triangulation analysis.  
- **Key Details:**  
  - Heliocentric longitude of Earth and geocentric longitude of Mars.  
  - Includes observation indices and corresponding dates.  

---

## 📝 **Assignments**  

### **1. Mars Opposition Analysis**  
- Derive Mars's projected position and identify the best-fit parameters using a custom loss function.  

### **2. Triangulation Data Analysis**  
- Use paired observations to triangulate Mars's ecliptic projections and determine its best-fit circular orbit around the Sun.  

### **3. Orbital Plane Identification**  
- Apply iterative regression to determine Mars's orbital plane inclination using heliocentric latitudes.  

### **4. Orbital Ellipse Calculation**  
- Compute Mars's 3D positions, fit a circle and ellipse on the orbital plane, and visualize the orbital elements.  

---

## 🔗 **Resources**  
- [Mars Opposition Data (`mars_opposition_data.csv`)](https://example.com/mars_opposition_data.csv)  
- [Triangulation Data (`triangulation_data.csv`)](https://example.com/triangulation_data.csv)  

---
