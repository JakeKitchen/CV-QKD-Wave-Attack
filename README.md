# Overview

This Continuous-Variable Quantum Key Distribution (CV-QKD) simulation demonstrates the effectiveness of a wave attack on noise levels during quantum key distribution. The project allows a comparison between normal and attack scenarios, illustrating how the wave attack maintains low noise levels, thus showing its effectiveness in going undetected.

## Features

- **Wave Attack Simulation**: Simulates a wave attack to assess its effectiveness in maintaining noise levels similar to a normal scenario.
- **Normal Scenario Comparison**: Provides a baseline noise comparison without attack.
- **Noise Calculation**: Calculates noise levels under randomized attenuation conditions, giving realistic noise distributions.
- **Statistical Output**: Outputs mean noise levels and standard deviation for both normal and attack scenarios.

## Requirements

Install the required packages via pip:
```bash
pip install pennylane numpy
```

## Usage

1. **Run the Simulation**: Use the provided code to simulate noise in both normal and attack scenarios.

    ```python
    # Run simulations for normal and wave attack scenarios
    normal_noises = simulate_scenario(alpha_vals, wavelengths["normal"])
    attack_noises = simulate_scenario(alpha_vals, wavelengths["wave_attack"])
    ```

2. **Analyze the Results**: The code prints the mean and standard deviation of noise for each scenario, allowing you to observe how the wave attack maintains similar noise levels to the normal scenario.

    ```python
    print(f'--- Normal Scenario ---')
    print(f'Mean Noise Level: {np.mean(normal_noises):.4f}')
    print(f'Standard Deviation of Noise: {np.std(normal_noises):.4f}')

    print(f'\n--- Wave Attack Scenario ---')
    print(f'Mean Noise Level: {np.mean(attack_noises):.4f}')
    print(f'Standard Deviation of Noise: {np.std(attack_noises):.4f}')
    ```

## Functions

- **`calculate_noise`**: Computes noise based on channel efficiency, shot noise, and attenuation ratios.
- **`measure_x_quadrature` & `measure_p_quadrature`**: Separate QNodes for measuring quadrature values under specified transmittance conditions.
- **`simulate_scenario`**: Conducts the main simulation, calculating noise levels for either normal or wave attack scenarios.

## Example Output

Sample output from the simulation:

```
--- Normal Scenario ---
Mean Noise Level: 0.5123
Standard Deviation of Noise: 0.0234

--- Wave Attack Scenario ---
Mean Noise Level: 0.5121
Standard Deviation of Noise: 0.0232
```
