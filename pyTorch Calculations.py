import torch

# -----------------------------------------------------------------------------
#  Defining and Organising S-Parameters
# -----------------------------------------------------------------------------

# Assigning the provided values
frequency = torch.tensor(5.000)  # GHz

# S11 Parameters
S11_mag = torch.tensor(0.6555)
S11_ang_deg = torch.tensor(126.9)

# S21 Parameters
S21_mag = torch.tensor(2.011)
S21_ang_deg = torch.tensor(22.5)

# S12 Parameters
S12_mag = torch.tensor(0.1142)
S12_ang_deg = torch.tensor(7.7)

# S22 Parameters
S22_mag = torch.tensor(0.2001)
S22_ang_deg = torch.tensor(-138.1)

# -----------------------------------------------------------------------------
# Stage 2: Converting S-Parameters from Polar to Rectangular Form
# -----------------------------------------------------------------------------

def polar_to_rect(magnitude, angle_deg):
    angle_rad = torch.deg2rad(angle_deg)  # Convert degrees to radians
    real = magnitude * torch.cos(angle_rad)
    imag = magnitude * torch.sin(angle_rad)
    return torch.complex(real, imag)

# Convert each S-parameter to rectangular form
S11 = polar_to_rect(S11_mag, S11_ang_deg)
S21 = polar_to_rect(S21_mag, S21_ang_deg)
S12 = polar_to_rect(S12_mag, S12_ang_deg)
S22 = polar_to_rect(S22_mag, S22_ang_deg)

print("----- S-Parameters (5GHz) in Rectangular Form -----")
print(f"S11: {S11:.3f}")
print(f"S21: {S21:.3f}")
print(f"S12: {S12:.3f}")
print(f"S22: {S22:.3f}\n")

# -----------------------------------------------------------------------------
# Calculating the Determinant (Δ)
# -----------------------------------------------------------------------------

# Calculate Δ = S11 * S22 - S12 * S21
Delta = S11 * S22 - S12 * S21
Delta_magnitude = torch.abs(Delta)

print("----- Determinant (Δ) -----")
print(f"Δ = S11 * S22 - S12 * S21 = {Delta:.3f}")
print(f"|Δ| = {Delta_magnitude.item():.3f}\n")

# -----------------------------------------------------------------------------
# Calculating the Stability Factor (K)
# -----------------------------------------------------------------------------

# Calculate numerator:
S11_mag_2 = torch.abs(S11)**2
S22_mag_2 = torch.abs(S22)**2
Delta_magnitude_2 = Delta_magnitude**2
numerator = 1 - S11_mag_2 - S22_mag_2 + Delta_magnitude_2

# Calculate denominator:
S12_S21 = S12 * S21
S12_S21_magnitude = torch.abs(S12_S21)
denominator = 2 * S12_S21_magnitude

# Calculate Stability Factor (K)
K = numerator / denominator

print("----- Stability Factor (K) -----")

print(f"Stability Factor K = {K.item():.3f}\n")

# -----------------------------------------------------------------------------
# Assessing Transistor Stability
# -----------------------------------------------------------------------------

print("----- Stability Check -----")
if (K > 1) and (Delta_magnitude < 1):
    print(f"K = {K.item():.3f} > 1 and |Δ| = {Delta_magnitude.item():.3f} < 1")
    print("The transistor is unconditionally stable at 5 GHz.\n")
    unconditionally_stable = True
else:
    print(f"K = {K.item():.3f} <= 1 or |Δ| = {Delta_magnitude.item():.3f} >= 1")
    print("The transistor is not unconditionally stable at 5 GHz.\n")
    unconditionally_stable = False

# -----------------------------------------------------------------------------
# Calculating Maximum Available Gain (G_Tmax)
# -----------------------------------------------------------------------------

if unconditionally_stable:
    # Calculate G_Tmax using the formula:
    # G_Tmax = (|S21| / |S12|) * (K - sqrt(K^2 - 1))
    
    # Calculate |S21| / |S12|
    S21_over_S12 = S21_mag / S12_mag
    
    # Calculate sqrt(K^2 - 1)
    sqrt_term = torch.sqrt(K**2 - 1)
    
    # Calculate (K - sqrt(K^2 - 1))
    K_minus_sqrt = K - sqrt_term
    
    # Calculate G_Tmax
    G_Tmax = S21_over_S12 * K_minus_sqrt
    
    # Convert G_Tmax to decibels (dB)
    G_Tmax_dB = 10 * torch.log10(G_Tmax)
    
    print("----- Maximum Available Transducer Gain (G_Tmax) -----")
    print(f"G_Tmax (linear) = {G_Tmax.item():.3f}")
    print(f"G_Tmax (dB) = {G_Tmax_dB.item():.3f} dB\n")
else:
    print("Since the transistor is not unconditionally stable, G_Tmax cannot be reliably calculated.\n")

# -----------------------------------------------------------------------------
# Summary of Results
# -----------------------------------------------------------------------------

if unconditionally_stable:
    print("===== Summary =====")
    print(f"Transistor Stability: Unconditionally Stable")
    print(f"Maximum Transducer Gain (G_Tmax): {G_Tmax.item():.3f} (linear)")
    print(f"Maximum Transducer Gain (G_Tmax): {G_Tmax_dB.item():.3f} dB")
else:
    print("===== Summary =====")
    print(f"Transistor Stability: Not Unconditionally Stable")
    print(f"Maximum Transducer Gain (G_Tmax): Cannot be calculated reliably.")
