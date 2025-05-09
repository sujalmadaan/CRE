# -*- coding: utf-8 -*-
"""
Created on Fri May  9 23:13:08 2025

@author: sujal
"""





import streamlit as st
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.optimize import brentq

st.title("Non-Adiabatic CSTR Sizing")
st.sidebar.link_button(" Home", "https://cll-122-project.vercel.app")
st.sidebar.link_button(" Adiabatic Reactor", "https://cll-122-project.vercel.app/adiabatic")
st.sidebar.link_button(" Team", "https://cll-122-project.vercel.app/team")

# Sidebar selection
option = st.sidebar.radio(
    "Choose Calculation Mode",
    ("Specify X", "Specify T", "Specify V")
)

if option == "Specify X":
    with st.form("cstr_form_x"):
        st.header("Input Parameters for Desired Conversion (X)")

        X_user = st.number_input("Desired Conversion (X)", min_value=0.001, max_value=0.999, value=0.8, step=0.01)
        F_A0 = st.number_input("Molar feed rate of A (mol/s)", min_value=1e-10, value=1.0, step=0.1)
        F_I0 = st.number_input("Molar feed rate of Inert (mol/s)", min_value=0.0, value=1.0, step=0.1)
        F_B0 = st.number_input("Molar feed rate of B (mol/s)", min_value=0.0, value=1.0, step=0.1)
       
        T_0 = st.number_input("Feed temperature T₀ (K)", min_value=1e-5, value=300.0)
        T_r = st.number_input("Reference temperature Tr (K)", min_value=1e-5, value=300.0)
        T_a = st.number_input("Ambient temperature Ta (K)", min_value=1e-5, value=300.0)
        delta_Hr0 = st.number_input("Heat of reaction ΔHr₀ (J/mol)", value=-50000.0)
        Cp_a = st.number_input("Heat capacity of A Cₚₐ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_b = st.number_input("Heat capacity of B Cₚᵦ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_i = st.number_input("Heat capacity of inert Cₚᵢ (J/mol·K)", min_value=1e-5, value=100.0)
        k_0 = st.number_input("Rate constant at reference temperature k₀ (1/s)", min_value=1e-10, value=0.001)
        E = st.number_input("Activation energy E (J/mol)", min_value=0.0, value=80000.0)
        U = st.number_input("Overall heat transfer coefficient U (J/s·m²·K)", min_value=0.0, value=500.0)
        A_heat = st.number_input("Heat exchange area A (m²)", min_value=0.0, value=1.0)

        submitted = st.form_submit_button("Calculate")

    if submitted:
        theta_a = 1
        theta_b = F_B0 / F_A0
        theta_i = F_I0 / F_A0
        R = 8.314

        Cp_0 = theta_a * Cp_a + theta_b * Cp_b + theta_i * Cp_i
        C_A0=F_A0/(F_A0+F_B0+F_I0)
        kappa = (U * A_heat) / (F_A0 * Cp_0)
        T_c = (kappa * T_a + T_0) / (1 + kappa)
        final_T = T_c + (-delta_Hr0 * X_user) / (Cp_0 * (1 + kappa))

        final_k = k_0 * np.exp(-E / R * (1 / final_T - 1 / T_r))
        final_rA = final_k * C_A0 * (1 - X_user)
        final_V = (F_A0 * X_user) / final_rA

        st.subheader("Results at Desired Conversion")
        st.write(f"Reactor Temperature T: **{final_T:.2f} K**")
        st.write(f"Rate Constant k: **{final_k:.4e} 1/s**")
        st.write(f"Reaction Rate -rA: **{final_rA:.4e} mol/L·s**")
        st.write(f"Required CSTR Volume V: **{final_V:.4f} L**")
        st.success("Calculation complete!")

elif option == "Specify T":
    with st.form("cstr_form_t"):
        st.header("Input Parameters for Specified Reactor Temperature (T)")

        T = st.number_input("Reactor temperature T (K)", min_value=200.0, max_value=1000.0, value=300.0, step=5.0)
        F_A0 = st.number_input("Molar feed rate of A (mol/s)", min_value=1e-10, value=1.0, step=0.1)
        F_I0 = st.number_input("Molar feed rate of Inert (mol/s)", min_value=0.0, value=1.0, step=0.1)
        F_B0 = st.number_input("Molar feed rate of B (mol/s)", min_value=0.0, value=1.0, step=0.1)
        
        T_0 = st.number_input("Feed temperature T0 (K)", min_value=1e-5, value=300.0)
        T_r = st.number_input("Reference temperature Tr (K)", min_value=1e-5, value=300.0)
        T_a = st.number_input("Ambient temperature Ta (K)", min_value=1e-5, value=300.0)
        delta_Hr0 = st.number_input("Heat of reaction ΔHr₀ (J/mol)", value=-50000.0)
        Cp_a = st.number_input("Heat capacity of A Cₚₐ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_b = st.number_input("Heat capacity of B Cₚᵦ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_i = st.number_input("Heat capacity of inert Cₚᵢ (J/mol·K)", min_value=1e-5, value=100.0)
        k_0 = st.number_input("Rate constant at reference temperature k₀ (1/s)", min_value=1e-10, value=1e6)
        E = st.number_input("Activation energy E (J/mol)", min_value=0.0, value=80000.0)
        U = st.number_input("Overall heat transfer coefficient U (J/s·m²·K)", min_value=0.0, value=500.0)
        A_heat = st.number_input("Heat exchange area A (m²)", min_value=0.0, value=1.0)

        submitted = st.form_submit_button("Calculate")

    if submitted:
        theta_a = 1
        theta_b = F_B0 / F_A0
        theta_i = F_I0 / F_A0
        R = 8.314
        C_A0=F_A0/(F_A0+F_B0+F_I0)

        Cp_0 = theta_a * Cp_a + theta_b * Cp_b + theta_i * Cp_i
        kappa = (U * A_heat) / (F_A0 * Cp_0)
        T_c = (kappa * T_a + T_0) / (1 + kappa)

        X = -Cp_0 * (1 + kappa) * (T - T_c) / delta_Hr0
        k = k_0 * np.exp(-E / R * (1 / T - 1 / T_r))
        rA = k * C_A0 * (1 - X)
        V = (F_A0 * X) / rA

        st.subheader("Results at Specified Temperature")
        st.write(f"Conversion X: **{X:.4f}**")
        st.write(f"Rate Constant k: **{k:.4e} 1/s**")
        st.write(f"Reaction Rate -rA: **{rA:.4e} mol/L·s**")
        st.write(f"Required CSTR Volume V: **{V:.4f} L**")
        st.success("Calculation complete!")

elif option == "Specify V":
    with st.form("volume_form"):
        st.header("Input Parameters for Specified Reactor Volume (V)")

        V = st.number_input("Reactor Volume V (L)", min_value=1e-6, value=10.0)
        F_A0 = st.number_input("Molar feed rate of A (mol/s)", min_value=1e-10, value=1.0)
        F_I0 = st.number_input("Molar feed rate of Inert (mol/s)", min_value=0.0, value=1.0)
        F_B0 = st.number_input("Molar feed rate of B (mol/s)", min_value=0.0, value=1.0)

      
        T_0 = st.number_input("Feed temperature T₀ (K)", min_value=1e-5, value=300.0)
        T_r = st.number_input("Reference temperature Tr (K)", min_value=1e-5, value=300.0)
        T_a = st.number_input("Ambient temperature Ta (K)", min_value=1e-5, value=300.0)

        delta_Hr0 = st.number_input("Heat of reaction ΔHr₀ (J/mol)", value=-50000.0)
        Cp_a = st.number_input("Heat capacity of A Cₚₐ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_b = st.number_input("Heat capacity of B Cₚᵦ (J/mol·K)", min_value=1e-5, value=100.0)
        Cp_i = st.number_input("Heat capacity of inert Cₚᵢ (J/mol·K)", min_value=1e-5, value=100.0)

        k_0 = st.number_input("Rate constant at reference temperature k₀ (1/s)", min_value=1e-10, value=1e6)
        E = st.number_input("Activation energy E (J/mol)", min_value=0.0, value=80000.0)

        U = st.number_input("Overall heat transfer coefficient U (J/s·m²·K)", min_value=0.0, value=500.0)
        A_heat = st.number_input("Heat exchange area A (m²)", min_value=0.0, value=1.0)

        submitted = st.form_submit_button("Calculate")
    if submitted:
        R = 8.314
        theta_a = 1
        theta_b = F_B0 / F_A0
        theta_i = F_I0 / F_A0
        C_A0=F_A0/(F_A0+F_B0+F_I0)

        Cp_0 = theta_a * Cp_a + theta_b * Cp_b + theta_i * Cp_i
        kappa = (U * A_heat) / (F_A0 * Cp_0)
        T_c = (kappa * T_a + T_0) / (1 + kappa)
        tau = V / (F_A0 + F_B0 + F_I0)  # residence time

        def X_mb(T):
            k = k_0 * np.exp(-E / R * (1 / T - 1 / T_r))
            return (tau * k) / (1 + tau * k)

        def X_eb(T):
            return -Cp_0 * (1 + kappa) * (T - T_c) / delta_Hr0

        def func_to_solve(T):
            return X_eb(T) - X_mb(T)

        # Temperature range to search for roots
        T_min, T_max = 250, 800
        T_vals = np.linspace(T_min, T_max, 5000)  # finer grid for better root detection
        f_vals = func_to_solve(T_vals)

        # Find intervals where f changes sign
        sign_changes = np.where(np.diff(np.sign(f_vals)))[0]

        roots = []
        for idx in sign_changes:
            T_left = T_vals[idx]
            T_right = T_vals[idx + 1]
            try:
                root = brentq(func_to_solve, T_left, T_right)
                # Avoid duplicates (due to numerical noise)
                if not any(np.isclose(root, r, atol=1e-4) for r in roots):
                    roots.append(root)
            except ValueError:
                # No root in this interval, skip
                pass

        # Calculate conversions at roots
        roots = sorted(roots)
        X_roots = [X_mb(T) for T in roots]

        st.subheader("Steady-State Solutions")
        if roots:
            for i, (T_root, X_root) in enumerate(zip(roots, X_roots), start=1):
                st.write(f"Solution {i}: Temperature = **{T_root:.2f} K**, Conversion = **{X_root:.4f}**")
        else:
            st.write("No steady-state solutions found in the specified temperature range.")

        # Plot curves and all intersection points
        X_mb_vals = [np.clip(X_mb(T), 0, 1) for T in T_vals]
        X_eb_vals = [np.clip(X_eb(T), 0, 1) for T in T_vals]

        fig, ax = plt.subplots()
        ax.plot(T_vals, X_mb_vals, label="X_mb (Mole Balance)", color='blue')
        ax.plot(T_vals, X_eb_vals, label="X_eb (Energy Balance)", color='red')

        for T_root, X_root in zip(roots, X_roots):
            ax.axvline(T_root, linestyle='--', color='green', alpha=0.7)
            ax.axhline(X_root, linestyle=':', color='grey', alpha=0.7)
            ax.plot(T_root, X_root, 'ko')  # mark intersection points

        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Conversion X")
        ax.set_title("Steady-State Conversion from X_mb and X_eb")
        ax.set_ylim(0, 1)
        ax.legend()
        ax.grid(True)

        st.pyplot(fig)
