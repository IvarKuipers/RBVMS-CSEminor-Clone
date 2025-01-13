# Re-importing the necessary components after environment reset

import os
import re
import csv
import matplotlib.pyplot as plt
import pandas as pd

# Define the functions again

def parse_log_to_csv(log_file, csv_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()

    # Regex patterns for relevant information
    step_pattern = re.compile(r'\s*step = (\d+) :')
    setup_time_pattern = re.compile(r'Setup time for new preconditioner:\s+([\d\.e\-]+)')
    iteration_time_pattern = re.compile(r'Time spent:\s+([\d\.]+) seconds')
    total_time_pattern = re.compile(r'Time taken for one Time step:\s+([\d\.]+) seconds')
    newton_cycle_pattern = re.compile(r'Newton iteration\s+(\d+)')
    newton_time_pattern = re.compile(r'\|\|r\|\|.*Time spent:\s+([\d\.]+) seconds')

    # Storage for extracted data
    data = []

    current_step = None
    setup_time = 0
    iteration_times = []
    total_time = 0
    newton_cycles = 0
    newton_times = []

    for line in lines:
        # Match step
        step_match = step_pattern.match(line)
        if step_match:
            # Save previous step data if exists
            if current_step is not None:
                solve_time = sum(iteration_times)
                time_per_iter = solve_time / len(iteration_times) if iteration_times else 0
                data.append([current_step, setup_time, solve_time, time_per_iter, total_time, newton_cycles, newton_times])

            # Start new step
            current_step = int(step_match.group(1))
            setup_time = 0
            iteration_times = []
            total_time = 0
            newton_cycles = 0
            newton_times = []

        # Match setup time
        setup_time_match = setup_time_pattern.search(line)
        if setup_time_match:
            setup_time += float(setup_time_match.group(1))

        # Match iteration time
        iteration_time_match = iteration_time_pattern.search(line)
        if iteration_time_match:
            iteration_times.append(float(iteration_time_match.group(1)))

        # Match total time
        total_time_match = total_time_pattern.search(line)
        if total_time_match:
            total_time = float(total_time_match.group(1))

        # Match newton cycle
        newton_cycle_match = newton_cycle_pattern.match(line)
        if newton_cycle_match:
            newton_cycles += 1

        # Match newton time
        newton_time_match = newton_time_pattern.search(line)
        if newton_time_match:
            newton_times.append(float(newton_time_match.group(1)))

    # Save last step data
    if current_step is not None:
        solve_time = sum(iteration_times)
        time_per_iter = solve_time / len(iteration_times) if iteration_times else 0
        data.append([current_step, setup_time, solve_time, time_per_iter, total_time, newton_cycles, newton_times])

    # Write to CSV
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Step", "Setup Time", "Solve Time", "Time Per Iteration", "Total Time", "Newton Cycles", "Newton Times"])
        for row in data:
            # Flatten Newton times into a string for CSV
            row[-1] = "|".join(map(str, row[-1]))
            writer.writerow(row)

log_file_path = '../cases/von-karman-re100/log0'  
csv_output_path = 'timings/run1.csv'    # Save the CSV in the graphs folder

# Parse the log file and generate graphs
parse_log_to_csv(log_file_path, csv_output_path)


