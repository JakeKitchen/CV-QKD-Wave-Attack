import pennylane as qml
import numpy as np

# Constants
N0 = 1.0  # Shot noise unit
I_LO = 1e4  # Local oscillator amplitude (scaled for realism)
num_trials = 1000

wavelengths = {
    "normal": {"signal": 0.5, "local_oscillator": 0.5},
    "wave_attack": {"signal": 0.5144, "local_oscillator": 0.5155}
}

def calculate_noise(eta_ch, noise_real, r1, r2):
    """
    Calculates total noise under channel and attenuation effects.
    """
    return (2 + noise_real) * N0 / eta_ch + (r1 + r2 - 2) * N0 / eta_ch

dev = qml.device("default.gaussian", wires=2)

@qml.qnode(dev)
def homodyne_measure(alpha, trans_signal, trans_LO, quadrature='x'):
    """
    Simulates a homodyne measurement using a 50:50 beamsplitter.
    
    Args:
        alpha (float): Signal displacement.
        trans_signal (float): Signal transmittance.
        trans_LO (float): LO transmittance.
        quadrature (str): 'x' or 'p' for the quadrature type.

    Returns:
        float: Expectation value of chosen quadrature.
    """
    # Prepare signal and LO
    qml.Displacement(alpha * np.sqrt(trans_signal), 0.0, wires=0)
    qml.Displacement(I_LO * np.sqrt(trans_LO), 0.0, wires=1)
    
    # Mix them
    qml.Beamsplitter(np.pi / 4, 0, wires=[0, 1])
    
    # Measure quadrature on wire 0
    if quadrature == 'x':
        return qml.expval(qml.ops.cv.X(0))
    elif quadrature == 'p':
        return qml.expval(qml.ops.cv.P(0))
    else:
        raise ValueError("Quadrature must be 'x' or 'p'")

def simulate_scenario(alpha_vals, scenario_params, eta_ch=0.9):
    """
    Runs simulated homodyne measurements for given alpha values and scenario.
    
    Args:
        alpha_vals (np.ndarray): Array of signal amplitudes.
        scenario_params (dict): Contains transmittance values.
        eta_ch (float): Channel transmission efficiency.

    Returns:
        list: Simulated noise values.
    """
    noises = []
    for alpha in alpha_vals:
        r1, r2 = np.random.choice([0.1, 0.5, 1.0], 2, replace=False)
        
        X = homodyne_measure(alpha * np.sqrt(r1),
                             scenario_params["signal"],
                             scenario_params["local_oscillator"],
                             quadrature='x')
        P = homodyne_measure(alpha * np.sqrt(r2),
                             scenario_params["signal"],
                             scenario_params["local_oscillator"],
                             quadrature='p')
        
        # Simulated shot noise
        shot_noise = N0 * np.random.uniform(0.9, 1.1)
        scenario_noise = calculate_noise(eta_ch, shot_noise, r1, r2)
        noises.append(scenario_noise)
    return noises

alpha_vals = np.random.normal(scale=0.5, size=num_trials)

normal_noises = simulate_scenario(alpha_vals, wavelengths["normal"])
attack_noises = simulate_scenario(alpha_vals, wavelengths["wave_attack"])

print(f'--- Normal Scenario ---')
print(f'Mean Noise Level: {np.mean(normal_noises):.4f}')
print(f'Standard Deviation of Noise: {np.std(normal_noises):.4f}')

print(f'\n--- Wave Attack Scenario ---')
print(f'Mean Noise Level: {np.mean(attack_noises):.4f}')
print(f'Standard Deviation of Noise: {np.std(attack_noises):.4f}')