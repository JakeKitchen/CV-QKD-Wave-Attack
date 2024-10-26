import pennylane as qml
import numpy as np

N0 = 1.0
I_LO = 1e8
num_trials = 1000

wavelengths = {
    "normal": {"signal": 0.5, "local_oscillator": 0.5},
    "wave_attack": {"signal": 0.5144, "local_oscillator": 0.5155}
}

def calculate_noise(eta_ch, noise_real, r1, r2):
    """
    Calculates noise for the given scenario with attenuation.

    Args:
        eta_ch (float): Channel efficiency.
        noise_real (float): Shot noise level.
        r1, r2 (float): Attenuation ratios.

    Returns:
        float: Calculated noise.
    """
    return (2 + noise_real) * N0 / eta_ch + (r1 + r2 - 2) * N0 / eta_ch

dev = qml.device("default.gaussian", wires=2)

@qml.qnode(dev)
def measure_x_quadrature(alpha, params):
    """
    Measures the X quadrature for the given scenario.

    Args:
        alpha (float): Displacement amplitude.
        params (dict): Contains transmittance values.

    Returns:
        float: Expectation value of the X quadrature.
    """
    transmittance = params["signal"]
    qml.Displacement(alpha * np.sqrt(transmittance), 0, wires=0)
    qml.Displacement(I_LO * np.sqrt(transmittance), 0, wires=1)
    return qml.expval(qml.QuadX(0))

@qml.qnode(dev)
def measure_p_quadrature(alpha, params):
    """
    Measures the P quadrature for the given scenario.

    Args:
        alpha (float): Displacement amplitude.
        params (dict): Contains transmittance values.

    Returns:
        float: Expectation value of the P quadrature.
    """
    transmittance = params["local_oscillator"]
    qml.Displacement(alpha * np.sqrt(transmittance), 0, wires=0)
    qml.Displacement(I_LO * np.sqrt(transmittance), 0, wires=1)
    return qml.expval(qml.QuadP(1))

def simulate_scenario(alpha_vals, scenario_params):
    """
    Simulates measurements for either normal or attack scenario.

    Args:
        alpha_vals (array): Array of initial alpha values.
        scenario_params (dict): Parameters for either the normal or wave attack scenario.

    Returns:
        list: Noise levels for the given scenario.
    """
    noises = []
    for alpha in alpha_vals:
        r1, r2 = np.random.choice([0.1, 0.5, 1.0], 2, replace=False)
        X = measure_x_quadrature(alpha * np.sqrt(r1), scenario_params)
        P = measure_p_quadrature(alpha * np.sqrt(r2), scenario_params)
        shot_noise = N0 * np.random.uniform(0.9, 1.1)
        scenario_noise = calculate_noise(0.9, shot_noise, r1, r2)
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
