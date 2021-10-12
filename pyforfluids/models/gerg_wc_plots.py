import matplotlib.pyplot as plt
import numpy as np


def conditions_plot(pressure, temperature):

    # Invalid area
    x_invalid = np.arange(0.0, 1000.0, 0.01)

    # Extended range area
    x_extended = np.arange(60.0, 700.0, 0.01)

    # Normal range area
    x_normal = np.arange(90.0, 450.0, 0.01)

    # Fill plot
    fig, ax = plt.subplots(sharex=True, figsize=(9, 6))
    ax.fill_between(
        x_invalid, 100.0, color="r", alpha=0.6, label="Invalid area"
    )
    ax.fill_between(
        x_extended,
        70.0,
        color="orange",
        alpha=0.5,
        label="Extended range of validity",
    )
    ax.fill_between(
        x_normal,
        35.0,
        color="green",
        alpha=0.5,
        label="Normal range of validity",
    )

    # Ticks
    xticks = [60.0, 90.0, 450.0, 700.0, 1000.0, temperature]
    yticks = [35.0, 70.0, 100.0, pressure]
    plt.xticks(xticks)
    plt.yticks(yticks)

    plt.axhline(pressure, linestyle="dashed", color="black", linewidth=0.8)
    plt.axvline(temperature, linestyle="dashed", color="black", linewidth=0.8)

    # Plot
    plt.scatter(temperature, pressure)
    plt.annotate("You are here", (temperature + 4, pressure + 2))
    plt.legend(loc="upper right")
    plt.title("Pressure - Temeperature condition plot")
    plt.xlabel(r"$Temperature[K]$", fontsize=12)
    plt.ylabel(r"$Pressure [MPa]$", fontsize=12)
    plt.xlim(0.0, 1000.0)
    plt.ylim(0.0, 100.0)

    plt.show()

    return
