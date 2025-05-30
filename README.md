# 📡 Log-Periodic Dipole Antenna Design & Analysis (1200–2520 MHz)

This repository contains MATLAB code for designing and analyzing a **Log-Periodic Dipole Antenna (LPDA)**. It was developed as part of an assignment for the *Antenna and Microwave Engineering* course at VIT Vellore.

📄 **Assignment PDF**: [View on Google Drive](https://drive.google.com/file/d/1tACaPMnFt6kfZKexkk-ji3y1cQ20gHHQ/view?usp=sharing)

---

## 📌 Objectives

* Design a Log-Periodic Dipole Antenna covering **1200 MHz to 2520 MHz**
* Compute key antenna parameters like:

  * Apex angle
  * Number of elements
  * Length of the antenna
  * Spacing and impedance values
* Visualize:

  * Antenna structure
  * Radiation patterns (2D and 3D)
  * Return loss, gain, VSWR, and impedance vs. frequency

---

## 🧮 Parameters Used

| Parameter                     | Value         |
| ----------------------------- | ------------- |
| Frequency Range               | 1200–2520 MHz |
| Directivity                   | 8.7 dB        |
| Tau (τ)                       | 0.925         |
| Sigma (σ)                     | 0.1738        |
| Z₀ (Characteristic Impedance) | 50 Ω          |
| Speed of Light (c)            | 3 × 10⁸ m/s   |

---

## 📊 Outputs and Visualizations

The MATLAB script computes and visualizes:

* 📈 **Return Loss vs Frequency**
* 🛰 **Antenna Structure**
* 🌀 **2D & 3D Radiation Patterns**
* ⚡ **Gain vs Frequency**
* 📉 **VSWR vs Frequency**
* 🧮 **Real & Imaginary Impedance vs Frequency**
* 🎯 **Directivity vs Frequency**

---

## 📂 File Structure

```
📁 Antenna DA
├── LogPeriodicAntenna.m  # MATLAB source code
└── README.md
```

---

## ▶️ How to Run

1. Open MATLAB.
2. Run the script: `LogPeriodicAntenna.m`
3. View outputs in the command window and plotted figures.

---

## 🧠 Key Concepts Covered

* Log-Periodic Antenna Theory
* Impedance Matching and Return Loss
* Directivity and Gain Calculations
* Radiation Pattern Modeling
* VSWR and Input Impedance Analysis

---

## 📝 Author

**Prarthna Puhan**
3rd-year B.Tech. in Electronics and Communication Engineering
VIT Vellore

---

## 📘 References

* Assignment Brief [PDF](https://drive.google.com/file/d/1tACaPMnFt6kfZKexkk-ji3y1cQ20gHHQ/view?usp=sharing)
* Balanis, C. A., *Antenna Theory: Analysis and Design*