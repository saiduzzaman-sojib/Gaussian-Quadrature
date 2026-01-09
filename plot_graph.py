import matplotlib.pyplot as plt

# --- 1. Data Setup ---
evals_gauss = [2, 3, 4]
errors_gauss = [1.82e-02, 5.34e-05, 2.10e-07]

evals_nc = [5]
error_trap = [1.19e-01]
error_simp = [2.15e-03]

# --- 2. Plotting ---
plt.figure(figsize=(10, 6))

plt.semilogy(evals_gauss, errors_gauss, 'o-', color='blue', linewidth=2, markersize=8, label='Gaussian Quadrature')
plt.semilogy(evals_nc, error_simp, 's', color='green', markersize=10, label="Simpson's 1/3 (n=4)")
plt.semilogy(evals_nc, error_trap, '^', color='red', markersize=10, label="Trapezoidal (n=4)")

# --- 3. Formatting ---
plt.title('Convergence Analysis: Gaussian vs. Newton-Cotes', fontsize=14)
plt.xlabel('Number of Function Evaluations (Computational Cost)', fontsize=12)
plt.ylabel('Absolute Error (Log Scale)', fontsize=12)
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend(fontsize=11)
plt.xticks([2, 3, 4, 5])

# --- 4. Save Output ---
plt.savefig('error_plot.png', dpi=300, bbox_inches='tight')
print("✅ Success! 'error_plot.png' has been generated.")
