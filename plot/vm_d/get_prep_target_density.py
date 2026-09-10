#!/usr/bin/env python3
"""Plot EPICS sensor readings from a comma-delimited CSV file.

Expected columns:
time, unix time, top temperature, bottom temperature, exhaust pressure, fill pressure

The time column must use the format YYYY-MM-DD HH:MM:SS.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import numpy as np


TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
MIN_TEMPERATURE = 20.7
MAX_PRESSURE = 25.0
MISSING_PERIOD_START = datetime(2021, 11, 23, 12, 0, 0)
MISSING_PERIOD_END = datetime(2021, 12, 15, 0, 0, 0)
MISSING_PERIOD_DURATION = MISSING_PERIOD_END - MISSING_PERIOD_START


@dataclass
class EpicsReadings:
    timestamps: List[datetime]
    unix_time: List[float]
    top_temperature: List[float]
    bottom_temperature: List[float]
    fill_pressure: List[float]
    exhaust_pressure: List[float]
    average_temperature: List[float]
    average_pressure: List[float]
    density: List[float]
    density_nearest: List[float]
    density_low: List[float]
    density_high: List[float]
    density_uncertainty_low_pct: List[float]
    density_uncertainty_high_pct: List[float]


@dataclass
class DensityTable:
    temperatures: np.ndarray
    pressures: np.ndarray
    densities: np.ndarray


def load_density_table(table_path: Path) -> DensityTable:
    raw_table = np.loadtxt(table_path, usecols=(0, 1, 2))
    if raw_table.ndim != 2 or raw_table.shape[1] != 3:
        raise ValueError(f"Unexpected EOS table shape in {table_path}")

    temperatures = np.unique(raw_table[:, 0])
    pressures = np.unique(raw_table[:, 1])
    density_grid = np.full((temperatures.size, pressures.size), np.nan, dtype=float)

    temperature_index = {value: index for index, value in enumerate(temperatures)}
    pressure_index = {value: index for index, value in enumerate(pressures)}

    for temperature, pressure, density in raw_table:
        density_grid[temperature_index[temperature], pressure_index[pressure]] = density

    if np.isnan(density_grid).any():
        raise ValueError(f"EOS grid is incomplete in {table_path}")

    return DensityTable(temperatures=temperatures, pressures=pressures, densities=density_grid)


def interpolate_density(table: DensityTable, temperature: float, pressure: float) -> float:
    temperature = float(np.clip(temperature, table.temperatures[0], table.temperatures[-1]))
    pressure = float(np.clip(pressure, table.pressures[0], table.pressures[-1]))

    temp_index = np.searchsorted(table.temperatures, temperature, side="right") - 1
    press_index = np.searchsorted(table.pressures, pressure, side="right") - 1

    temp_index = int(np.clip(temp_index, 0, table.temperatures.size - 2))
    press_index = int(np.clip(press_index, 0, table.pressures.size - 2))

    temp_low = table.temperatures[temp_index]
    temp_high = table.temperatures[temp_index + 1]
    press_low = table.pressures[press_index]
    press_high = table.pressures[press_index + 1]

    temp_fraction = 0.0 if temp_high == temp_low else (temperature - temp_low) / (temp_high - temp_low)
    press_fraction = 0.0 if press_high == press_low else (pressure - press_low) / (press_high - press_low)

    density_ll = table.densities[temp_index, press_index]
    density_lh = table.densities[temp_index, press_index + 1]
    density_hl = table.densities[temp_index + 1, press_index]
    density_hh = table.densities[temp_index + 1, press_index + 1]

    return (
        density_ll * (1.0 - temp_fraction) * (1.0 - press_fraction)
        + density_hl * temp_fraction * (1.0 - press_fraction)
        + density_lh * (1.0 - temp_fraction) * press_fraction
        + density_hh * temp_fraction * press_fraction
    )


def nearest_density(table: DensityTable, temperature: float, pressure: float) -> float:
    temperature_index = int(np.abs(table.temperatures - temperature).argmin())
    pressure_index = int(np.abs(table.pressures - pressure).argmin())
    return float(table.densities[temperature_index, pressure_index])


def read_readings(csv_path: Path, density_table_path: Path) -> EpicsReadings:
    density_table = load_density_table(density_table_path)

    timestamps: List[datetime] = []
    unix_time: List[float] = []
    top_temperature: List[float] = []
    bottom_temperature: List[float] = []
    fill_pressure: List[float] = []
    exhaust_pressure: List[float] = []
    average_temperature: List[float] = []
    average_pressure: List[float] = []
    density: List[float] = []
    density_nearest: List[float] = []
    density_low: List[float] = []
    density_high: List[float] = []
    density_uncertainty_low_pct: List[float] = []
    density_uncertainty_high_pct: List[float] = []
    skipped_rows = 0

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) != 6:
                raise ValueError(
                    f"Expected 6 comma-separated fields on line {row_number}, got {len(row)}"
                )

            timestamp = datetime.strptime(row[0].strip(), TIME_FORMAT)
            if MISSING_PERIOD_START <= timestamp < MISSING_PERIOD_END:
                skipped_rows += 1
                continue

            if timestamp >= MISSING_PERIOD_END:
                timestamp -= MISSING_PERIOD_DURATION

            top_value = float(row[2])
            bottom_value = float(row[3])
            exhaust_value = float(row[4])
            fill_value = float(row[5])

            if (
                top_value < MIN_TEMPERATURE
                or bottom_value < MIN_TEMPERATURE
                or fill_value > MAX_PRESSURE
                or exhaust_value > MAX_PRESSURE
            ):
                skipped_rows += 1
                continue

            timestamps.append(timestamp)
            unix_time.append(float(row[1]))
            top_temperature.append(top_value)
            bottom_temperature.append(bottom_value)
            fill_pressure.append(fill_value)
            exhaust_pressure.append(exhaust_value)

            avg_temperature = 0.5 * (top_value + bottom_value)
            avg_pressure = 0.5 * (fill_value + exhaust_value)
            average_temperature.append(avg_temperature)
            average_pressure.append(avg_pressure)
            density.append(interpolate_density(density_table, avg_temperature, avg_pressure))
            density_nearest.append(nearest_density(density_table, avg_temperature, avg_pressure))

            density_samples = [
                interpolate_density(density_table, top_value, fill_value),
                interpolate_density(density_table, top_value, exhaust_value),
                interpolate_density(density_table, bottom_value, fill_value),
                interpolate_density(density_table, bottom_value, exhaust_value),
            ]
            density_min = min(density_samples)
            density_max = max(density_samples)
            density_low.append(density_min)
            density_high.append(density_max)

            density_value = density[-1]
            density_uncertainty_low_pct.append(100.0 * max(0.0, density_value - density_min) / density_value)
            density_uncertainty_high_pct.append(100.0 * max(0.0, density_max - density_value) / density_value)

    if not timestamps:
        raise ValueError(f"No data rows found in {csv_path}")

    if skipped_rows:
        print(f"Skipped {skipped_rows} outlier row(s) while reading {csv_path.name}")

    return EpicsReadings(
        timestamps=timestamps,
        unix_time=unix_time,
        top_temperature=top_temperature,
        bottom_temperature=bottom_temperature,
        fill_pressure=fill_pressure,
        exhaust_pressure=exhaust_pressure,
        average_temperature=average_temperature,
        average_pressure=average_pressure,
        density=density,
        density_nearest=density_nearest,
        density_low=density_low,
        density_high=density_high,
        density_uncertainty_low_pct=density_uncertainty_low_pct,
        density_uncertainty_high_pct=density_uncertainty_high_pct,
    )


def plot_density(ax: plt.Axes, readings: EpicsReadings) -> None:
    ax.plot(readings.timestamps, readings.density, color="tab:purple", linewidth=1.5)
    ax.set_ylabel("Density (g/mL)")
    ax.grid(True, alpha=0.3)


def plot_density_uncertainty(ax: plt.Axes, readings: EpicsReadings) -> None:
    ax.fill_between(
        readings.timestamps,
        readings.density_uncertainty_low_pct,
        readings.density_uncertainty_high_pct,
        color="tab:purple",
        alpha=0.25,
    )
    ax.plot(
        readings.timestamps,
        readings.density_uncertainty_low_pct,
        color="tab:purple",
        linewidth=1.0,
        linestyle="--",
    )
    ax.plot(
        readings.timestamps,
        readings.density_uncertainty_high_pct,
        color="tab:purple",
        linewidth=1.0,
        linestyle="--",
    )
    ax.set_ylabel("Uncertainty (%)")
    ax.grid(True, alpha=0.3)


def get_running_periods(readings: EpicsReadings) -> List[np.ndarray]:
    if not readings.top_temperature:
        return []

    period_indices: List[np.ndarray] = []
    start_index = 0
    top_temperature = np.asarray(readings.top_temperature, dtype=float)

    for index in range(1, len(top_temperature)):
        if abs(top_temperature[index] - top_temperature[index - 1]) > 1.0:
            period_indices.append(np.arange(start_index, index))
            start_index = index

    period_indices.append(np.arange(start_index, len(top_temperature)))
    return period_indices


def plot_density_variation(ax: plt.Axes, readings: EpicsReadings) -> None:
    density_values = np.asarray(readings.density, dtype=float)
    periods = get_running_periods(readings)

    for period_number, indices in enumerate(periods, start=1):
        if indices.size == 0:
            continue

        period_density = density_values[indices]
        period_average_density = float(np.mean(period_density))
        if period_average_density == 0.0:
            continue

        variation_pct = 100.0 * (period_density - period_average_density) / period_average_density
        ax.plot(
            np.asarray(readings.timestamps, dtype=object)[indices],
            variation_pct,
            linewidth=1.2,
            color="tab:purple",
        )

    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
    ax.set_ylabel("Variation (%)")
    ax.set_ylim(-0.4, 0.1)
    ax.grid(True, alpha=0.3)


def plot_density_vs_temperature(density_table: DensityTable, title: str | None = None) -> plt.Figure:
    fig, ax = plt.subplots(
        1,
        1,
        figsize=(10, 7),
        dpi=300,
        constrained_layout=True,
    )

    temperature_mask = (density_table.temperatures >= 20.5) & (density_table.temperatures <= 23.0)
    selected_temperatures = density_table.temperatures[temperature_mask]
    selected_densities = density_table.densities[temperature_mask, :]

    for pressure_index, pressure_value in enumerate(density_table.pressures):
        ax.plot(
            selected_temperatures,
            selected_densities[:, pressure_index],
            linewidth=1.2,
            label=f"P = {pressure_value:.2f}",
        )

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Density (g/mL)")
    ax.set_title(title or "Density vs Temperature at Constant Pressure")
    ax.set_xlim(20.5, 23.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)
    return fig


def plot_readings(readings: EpicsReadings, title: str | None = None) -> plt.Figure:
    fig, (ax_temp, ax_pressure, ax_density, ax_uncertainty, ax_variation) = plt.subplots(
        5,
        1,
        figsize=(14, 14),
        dpi=300,
        sharex=True,
        constrained_layout=False,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0, 0.5, 0.5]},
    )

    fig.suptitle(title or "EPICS Sensor Readings", fontsize=14)
    fig.subplots_adjust(hspace=0.0)

    ax_temp.plot(readings.timestamps, readings.top_temperature, label="Top temperature", color="tab:red", linewidth=1.5)
    ax_temp.plot(
        readings.timestamps,
        readings.bottom_temperature,
        label="Bottom temperature",
        color="tab:orange",
        linewidth=1.5,
    )
    ax_temp.set_ylabel("Temperature (K)")
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend(loc="best")

    ax_pressure.plot(readings.timestamps, readings.fill_pressure, label="Fill pressure", color="tab:blue", linewidth=1.5)
    ax_pressure.plot(
        readings.timestamps,
        readings.exhaust_pressure,
        label="Exhaust pressure",
        color="tab:green",
        linewidth=1.5,
    )
    ax_pressure.set_ylabel("Pressure (psia)")
    ax_pressure.set_xlabel("Time")
    ax_pressure.grid(True, alpha=0.3)
    ax_pressure.legend(loc="best")

    plot_density(ax_density, readings)
    ax_density.tick_params(labelbottom=False)

    plot_density_uncertainty(ax_uncertainty, readings)
    ax_uncertainty.set_xlabel("Time")

    plot_density_variation(ax_variation, readings)
    ax_variation.set_xlabel("Time")

    def format_actual_time(x_value: float, _position: int) -> str:
        tick_time = mdates.num2date(x_value).replace(tzinfo=None)
        if tick_time >= MISSING_PERIOD_END - MISSING_PERIOD_DURATION:
            tick_time += MISSING_PERIOD_DURATION
        return tick_time.strftime("%Y-%m-%d\n%H:%M:%S")

    ax_pressure.xaxis.set_major_formatter(FuncFormatter(format_actual_time))
    fig.autofmt_xdate()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot EPICS sensor readings from a CSV file.")
    parser.add_argument(
        "csv_file",
        nargs="?",
        default="input/target_density_epics_readings.csv",
        help="Input CSV file (default: input/target_density_epics_readings.csv in the current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output/plots_phi_d_target_density.pdf",
        help="Output PDF file (default: output/plots_phi_d_target_density.pdf)",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title",
    )
    parser.add_argument(
        "--density-table",
        default=str(Path(__file__).resolve().parent / "input" / "target_density_deuterium_eos.csv"),
        help="EOS table for density lookup (default: input/target_density_deuterium_eos.csv next to this script)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    csv_path = Path(args.csv_file).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    density_table_path = Path(args.density_table).expanduser().resolve()

    density_table = load_density_table(density_table_path)
    readings = read_readings(csv_path, density_table_path)
    density_values = np.asarray(readings.density, dtype=float)
    density_mean = float(np.mean(density_values))
    density_stddev = float(np.std(density_values))
    running_periods = get_running_periods(readings)

    print(f"Density over time: mean = {density_mean:.6f} g/mL, std dev = {density_stddev:.6f} g/mL")
    for period_number, indices in enumerate(running_periods, start=1):
        if indices.size == 0:
            continue
        period_mean_density = float(np.mean(np.asarray(readings.density, dtype=float)[indices]))
        print(f"Period {period_number} average density = {period_mean_density:.6f} g/mL")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_path) as pdf:
        pdf.savefig(plot_readings(readings, title=args.title), bbox_inches="tight")
        pdf.savefig(plot_density_vs_temperature(density_table, title=args.title), bbox_inches="tight")
        plt.close("all")

    print(f"Saved plots to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())