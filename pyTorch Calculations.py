import torch

# -----------------------------------------------------------------------------
#  Defining S-Parameters of BFP540 @ 5GHz
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
# Converting S-Parameters from Polar to Cartesian Form
# -----------------------------------------------------------------------------

def polar_to_cartesian(magnitude, angle_deg):
    angle_rad = torch.deg2rad(angle_deg)  # Convert degrees to radians
    real = magnitude * torch.cos(angle_rad)
    imag = magnitude * torch.sin(angle_rad)
    return torch.complex(real, imag)

# Convert each S-parameter to cartesian form
S11 = polar_to_cartesian(S11_mag, S11_ang_deg)
S21 = polar_to_cartesian(S21_mag, S21_ang_deg)
S12 = polar_to_cartesian(S12_mag, S12_ang_deg)
S22 = polar_to_cartesian(S22_mag, S22_ang_deg)

# -----------------------------------------------------------------------------
# Calculating the Determinant (Δ)
# -----------------------------------------------------------------------------

Delta = S11 * S22 - S12 * S21
Delta_magnitude = torch.abs(Delta)

#------------------------------------------------------------------------------
# Calculating B1, B2, C1, C2
#------------------------------------------------------------------------------

S11_mag_2 = torch.abs(S11)**2
S22_mag_2 = torch.abs(S22)**2
Delta_magnitude_2 = Delta_magnitude**2

B1 = 1 + S11_mag_2 - S22_mag_2 - Delta_magnitude_2 
B2 = 1 + S22_mag_2 - S11_mag_2 - Delta_magnitude_2
C1 = S11 - (Delta * (torch.conj(S22)))
C2 = S22 - (Delta * (torch.conj(S11)))

#------------------------------------------------------------------------------\
# Reflection Coefficients 
#------------------------------------------------------------------------------

# Source Reflection Coeffcient (Positive and Negative)
Gamma_s_pos = (B1 + torch.sqrt((B1 **2) - (4 * torch.abs(C1)**2))) / (2 * C1)
Gamma_s_neg = (B1 - torch.sqrt((B1 **2) - (4 * torch.abs(C1)**2))) / (2 * C1)

# Load Reflection Coefficients (Positive and Negative)
Gamma_l_pos = (B2 + torch.sqrt((B2 **2) - (4 * torch.abs(C2)**2))) / (2 * C2)
Gamma_l_neg = (B2 - torch.sqrt((B2 **2) - (4 * torch.abs(C2)**2))) / (2 * C2)

# Checking which Absolute Values of Reflection Coefficients are less than 1

if torch.abs(Gamma_s_pos) < torch.abs(Gamma_s_neg):
    Gamma_in = S11 + ((S12 * S21 * Gamma_l_pos)/(1 - (S11 * Gamma_l_pos)))
else:
     Gamma_in = S11 + ((S12 * S21 * Gamma_l_neg)/(1 - (S11 * Gamma_l_neg)))

if torch.abs(Gamma_l_pos) < torch.abs(Gamma_l_neg):
    Gamma_out = S22 + ((S12 * S21 * Gamma_s_pos)/(1 - (S11 * Gamma_s_pos)))
else:
    Gamma_out = S22 + ((S12 * S21 * Gamma_s_neg)/(1 - (S11 * Gamma_s_neg)))
   
#------------------------------------------------------------------------------
# Calculating Gs, Go, GL and G_T
#------------------------------------------------------------------------------

gain_source = 1 / (1 - (torch.abs(Gamma_s_neg)** 2))
gain_internal = torch.abs(S21)**2
gain_load = (1-torch.abs(Gamma_l_neg)**2) / (torch.abs(1- S22*Gamma_l_neg)**2)

gain_max = gain_source * gain_internal * gain_load
gain_max_dB = 10 * torch.log10(gain_max)

# ------------------------------------------------------------------------------
# Outputting Results
# ------------------------------------------------------------------------------

print(f"""
===== S-Parameters (5GHz) in Rectangular Form =====
S11: {S11:.3f}
S21: {S21:.3f}
S12: {S12:.3f}
S22: {S22:.3f}

===== Determinant (Δ) =====
Δ: {Delta:.3f}
|Δ|: {Delta_magnitude.item():.3f} 

===== B1, B2, C1, C2 =====
B1: {B1:.3f}
B2: {B2:.3f}
C1: {C1:.3f}
C2: {C2:.3f} 

===== Reflection Coefficients ===== 
Positive Source Reflection Coefficient: {Gamma_s_pos:.3f}
Negative Source Reflection Coefficient: {Gamma_s_neg:.3f}
Positive Load Reflection Coefficient: {Gamma_l_pos:.3f}
Negative Load Reflection Coefficient: {Gamma_l_neg:.3f}
Input Reflection Coefficient: {Gamma_in:.4f}
Output Reflection Coefficient: {Gamma_out:.3f}

===== Amplifier Gains =====
Gain Source (Gs): {gain_source.item():.3f}
Gain Internal (Go): {gain_internal.item():.3f}
Gain Load (Gl): {gain_load.item():.3f}

========= Final Results ==========
Max Gain (Linear): {gain_max.item():.3f}
Max Gain (dB): {gain_max_dB.item():.3f} dB
""")