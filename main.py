import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BIN_WIDTH_M = 0.32
BIN_LENGTH_M = 0.26
BIN_HEIGHT_M = 0.08

TUBE_DIAMETER_M = 0.00794
TUBE_LENGTHS_M = [0.2, 0.3, 0.4, 0.6]

END_HEIGHT_M = 0.02
TOTAL_HEIGHT_CHANGE_M = 0.08
START_HEIGHT_M = END_HEIGHT_M + TOTAL_HEIGHT_CHANGE_M

WATER_MU = 8.90e-4
WATER_RHO = 1.0
GRAVITATIONAL_CONSTANT = 9.8

SIN_THETA_RADS = 1 / 150

TIMESTEP_S = 1e-1

EXPERIMENTAL_RESULTS = []

AX_INDICES = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)}


def get_bin_volume_m3():
    return BIN_WIDTH_M * BIN_LENGTH_M * (BIN_HEIGHT_M + END_HEIGHT_M)


def get_pipe_area():
    return np.pi * TUBE_DIAMETER_M**2 / 4


def get_tube_volume_m3(L):
    return get_pipe_area() * L


def get_start_relative_height_m(L):
    return SIN_THETA_RADS * L + START_HEIGHT_M


def get_velocity_v2_m_per_s(v1, h):
    return np.sqrt(v1**2 + 2 * GRAVITATIONAL_CONSTANT * h)


def get_reynolds_number(rho, u, L, mu):
    return rho * u * L / mu


def run_simulation(L):
    bin_volume_m3 = get_bin_volume_m3()
    tube_volume_m3 = get_tube_volume_m3(L)
    
    total_volume_m3 = bin_volume_m3 + tube_volume_m3

    start_relative_height_m = get_start_relative_height_m(L)

    elapsed_time_s = 0
    v1 = 0
    current_height_m = start_relative_height_m
    previous_height_m = current_height_m
    total_delta_h_m = 0
    PIPE_AREA_M2 = get_pipe_area()
    dt = TIMESTEP_S
    previous_volume_m3 = total_volume_m3

    columns = ["t", "v1", "v2", "h", "V", "dh", "dV"]
    # data = [[elapsed_time_s, v1, 0, current_height_m, total_volume_m3, 0, 0]]
    data = []

    while total_delta_h_m < TOTAL_HEIGHT_CHANGE_M and elapsed_time_s < 500:
        v2 = get_velocity_v2_m_per_s(v1, current_height_m)
        dV = v2 * PIPE_AREA_M2 * dt
        current_volume_m3 = previous_volume_m3 - dV
        dh = dV / (BIN_WIDTH_M * BIN_LENGTH_M)
        # print(dh, current_volume_m3, tube_volume_m3, BIN_WIDTH_M, BIN_LENGTH_M, BIN_LENGTH_M * BIN_WIDTH_M)
        # exit()
        current_height_m = previous_height_m - dh
        v1 = dh / dt

        data.append([elapsed_time_s, v1, v2, current_height_m, current_volume_m3, dh, dV])

        previous_volume_m3 = current_volume_m3
        previous_height_m = current_height_m
        elapsed_time_s += dt
        total_delta_h_m += dh

    df = pd.DataFrame(data=data, columns=columns)
    return df


def main():
    fig, ax = plt.subplots(2, 2)

    for i, L in enumerate(TUBE_LENGTHS_M):
        df = run_simulation(L)

        x = df["t"]
        y = df["h"]
        title = f"{round(L, 1)}m"
        ax_indices = AX_INDICES[i]
        ax[ax_indices].plot(x, y)
        ax[ax_indices].set_title(title)
    
    plt.show()



        


if __name__ == "__main__":
    main()