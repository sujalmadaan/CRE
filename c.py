# -*- coding: utf-8 -*-
"""
Created on Fri May  9 12:21:59 2025
@author: sujal
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title("Non-Adiabatic CSTR Volume Calculator")

with st.form("cstr_form"):
    st.header("Input Parameters")

    X_user = st.number_input("Desired Conversion (X)", min_value=0.001, max_value=0.999, value=0.8, step=0.01)

    F_A0 = st.number_input("Molar feed rate of A (mol/s)", min_value=1e-10, value=1.0, step=0.1)
    F_I0 = st.number_input("Molar feed rate of Inert (mol/s)", min_value=0.0, value=1.0, step=0.1)
    F_B0 = st.number_input("Molar feed rate of B (mol/s)", min_value=0.0, value=1.0, step=0.1)

    C_A0 = st.number_input("Feed concentration of A (mol/L)", min_value=1e-10, value=1.0, step=0.1)
    C_I0 = st.number_input("Feed concentration of Inert (mol/L)", min_value=0.0, value=1.0, step=0.1)

    T_0 = st.number_input("Feed temperature T0 (K)", min_value=1e-5, value=300.0)
    T_r = st.number_input("Reference temperature Tr (K)", min_value=1e-5, value=300.0)
    T_a = st.number_input("Ambient temperature Ta (K)", min_value=1e-5, value=300.0)

    delta_Hr0 = st.number_input("Heat of reaction ΔHr0 (J/mol)", value=-50000.0)
    Cp_a = st.number_input("Heat capacity of A Cp_a (J/mol·K)", min_value=1e-5, value=100.0)
    Cp_b = st.number_input("Heat capacity of B Cp_b (J/mol·K)", min_value=1e-5, value=100.0)
    Cp_i = st.number_input("Heat capacity of inert Cp_i (J/mol·K)", min_value=1e-5, value=100.0)

    k_0 = st.number_input("Rate constant at reference temperature k₀ (1/s)", min_value=1e-10, value=1e6)
    E = st.number_input("Activation energy E (J/mol)", min_value=0.0, value=80000.0)

    U = st.number_input("Overall heat transfer coefficient U (J/s·m²·K)", min_value=0.0, value=500.0)
    A_heat = st.number_input("Heat exchange area A (m²)", min_value=0.0, value=1.0)

    submitted = st.form_submit_button("Calculate")

if submitted:
    theta_a = 1
    theta_b = F_B0 / F_A0
    theta_i = F_I0 / F_A0
    R = 8.314  # J/mol·K

    # Overall heat capacity (J/mol·K)
    Cp_0 = theta_a * Cp_a + theta_b * Cp_b + theta_i * Cp_i

    # Kappa: heat transfer parameter
    kappa = (U * A_heat) / (F_A0 * Cp_0)

    # Coolant temperature
    T_c = (kappa * T_a + T_0) / (1 + kappa)

   

  

    # Final results
    final_T = T_c + (-delta_Hr0 * X_user) / (Cp_0 * (1 + kappa))
    final_k = k_0 * np.exp(-E / R * (1 / final_T - 1 / T_r))
    final_rA = final_k * C_A0 * (1 - X_user)
    final_V = (F_A0 * X_user) / final_rA

    st.subheader("Results at Desired Conversion")
    st.write(f"Reactor Temperature T: **{final_T:.2f} K**")
    st.write(f"Rate Constant k: **{final_k:.4e} 1/s**")
    st.write(f"Reaction Rate -rA: **{final_rA:.4e} mol/L·s**")
    st.write(f"Required CSTR Volume V: **{final_V:.4f} L**")
    st.success("Calculation and plotting complete!")
