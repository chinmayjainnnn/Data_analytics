# 🌌 E0 259 - Data Analytics: Module 2 (August 2019)

**Analyzing Kepler's data on Mars**  
This repository contains all the materials, datasets, and instructions for Module 2 of the E0 259 Data Analytics course.

---

## 📚 **Contents**
1. [Overview](#overview)
2. [Lectures](#lectures)
3. [Datasets](#datasets)
   - [Mars Opposition Data](#mars-opposition-data)
   - [Triangulation Data](#triangulation-data)
4. [Assignment 2](#assignment-2)
   - [Questions](#assignment-questions)
5. [Submission Guidelines](#submission-guidelines)
6. [Resources and Links](#resources-and-links)
7. [Contact](#contact)

---

## 🌍 **Overview**
Module 2 of E0 259 focuses on analyzing Mars's orbital data using trigonometric and regression techniques. This involves exploring:
- Mars's opposition data to identify orbital parameters.
- Triangulation data to compute orbital projections.
- Application of regression and nonlinear optimization techniques.

---

## 🧑‍🏫 **Lectures**
1. **Lecture 1:** Introduction to Mars orbit analysis (Ramesh Hariharan)
2. **Lectures 2 & 3:** Techniques for regression and nonlinear optimization (Ramesh Hariharan)

Lecture slides and additional materials can be found on the course Moodle platform.

---

## 📊 **Datasets**

### **Mars Opposition Data**
- **File:** `mars_opposition_data.csv`
- **Description:** Contains Mars's longitude and latitude during its opposition with the Sun.
- **Key Details:**
  - **Columns A/B/C:** Day, Month, Year
  - **Columns D/E/F/G:** ZodiacIndex, Degree, Minute, Second (heliocentric longitude)
    - **Formula:** Longitude = `s * 30 + Degree + Minute/60 + Second/3600` (degrees)
  - **Columns H/I:** Degree, Minute (geocentric latitude)
  - **Columns J/K/L/M:** Mars's mean longitude (Kepler's equant approximation)

### **Triangulation Data**
- **File:** `triangulation_data.csv`
- **Description:** Longitudinal data for Earth-Mars observations.
- **Key Details:**
  - **Column A:** Observation pair index
  - **Columns B/C/D:** Day, Month, Year
  - **Columns E/F:** Heliocentric longitude of Earth
  - **Columns G/H:** Geocentric longitude of Mars

---

## 📝 **Assignment 2**
### **Deadlines**
- **Part 1 (Question 1):** Monday, 09 September 2019 (11:55 PM)
- **Part 2 (Questions 2-4):** Monday, 23 September 2019 (11:55 PM)

### **Assignment Questions**

#### **1. Mars Opposition Analysis**
- Use the Mars opposition data to:
  - Derive an expression for Mars's projected position.
  - Identify the best-fit parameters using a custom loss function.

#### **2. Triangulation Data Analysis**
- Use paired observations to:
  - Triangulate Mars's ecliptic projections.
  - Identify the best-fit circular orbit centered at the Sun.

#### **3. Orbital Plane Identification**
- Use heliocentric latitudes to:
  - Develop an iterative regression algorithm.
  - Determine the inclination of Mars's orbital plane.

#### **4. Orbital Ellipse Calculation**
- Identify:
  - Mars's 3D positions on the orbital plane.
  - The best-fit circle and ellipse on the plane.
- Visualize results with plots of orbital elements.

---

## 📥 **Submission Guidelines**
- **Platform:** Course Moodle
- **Instructions:**
  1. Submit your code as a single package.
  2. Accept the submission statement.
  3. Click the **Submit** button.
- **Important Notes:**
  - Late or email submissions will **not** be accepted.
  - Collaborative discussion is allowed, but your code must be original.

---

## 🔗 **Resources and Links**
- [Mars Opposition Data (`mars_opposition_data.csv`)](https://example.com/mars_opposition_data.csv)
- [Triangulation Data (`triangulation_data.csv`)](https://example.com/triangulation_data.csv)
- [Course Moodle](https://example.com/moodle-login) *(Login required)*

---

## 💡 **Getting Started**
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/e0-259-mars-data-analytics.git
   ```
2. Access datasets in the `data/` directory.
3. Review lecture slides in the `lectures/` directory.
4. Follow the assignment instructions in `assignments/`.

---

## 📧 **Contact**
For queries, reach out to:  
**Course Coordinator:** Ramesh Hariharan  
**Email:** [course-email@example.com](mailto:course-email@example.com)

---

Start early, plan your approach, and enjoy the journey of unraveling Mars's mysteries! 🚀