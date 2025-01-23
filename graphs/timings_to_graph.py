import os
import re
import csv
import matplotlib.pyplot as plt
import pandas as pd

def generate_combined_graphs(csv_files, output_dir):
    """
    Generate combined graphs from multiple CSV files.
    
    Parameters:
        csv_files (list): List of CSV file paths.
        output_dir (str): Directory to save the output graphs.
    """
    #os.makedirs(output_dir, exist_ok=True)
    
    # Initialize a dictionary to store data from all files
    combined_data = {}
    averages = []  # To store averages for the bar chart
    
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df['Execution Time'] = df['Newton Time'] - df['Jacobi Time']
        # Calculate averages for this run
        avg_setup_time = df['Setup Time'].mean()
        avg_total_time = df['Execution Time'].mean()
        run_label = os.path.basename(csv_file).split('.')[0]
        Total_time = df['Total Time'].sum()
        print('Total time spent on this run: ', Total_time, 'Seconds')
        averages.append({
            'Run': run_label,
            'Average Setup Time': avg_setup_time,
            'Average Execution Time': avg_total_time
        })
        
        label = os.path.basename(csv_file).split('.')[0]  # Use file name as label
        combined_data[label] = df
        
        
    
    # Define a helper function to plot combined data
    def plot_combined(x, y, xlabel, ylabel, title, output_file):
        plt.figure()
        for label, df in combined_data.items():
            plt.plot(df[x], df[y], marker='o', label=label, linestyle='None')
        plt.xlabel('Newton iteration')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(output_dir, output_file))
        plt.close()
    
    # Generate combined graphs
    plot_combined('Step', 'Jacobi Time','Newton Iterations', 'Jacobi Setup Time (s)', 'Total Time per Step', 'Jacobi_time_vs_step.png')
    plot_combined('Step', 'Setup Time','Newton Iterations', 'Setup Time (s)', 'Setup Time per Step', 'setup_time_vs_step.png')
    plot_combined('Step', 'Newton Time','Newton Iterations', 'Newton Iteration Time (s)', 'Newton Iteration Time per Step', 'Newton_time_vs_step.png')
    plot_combined('Step', 'FGMRES Iterations','Newton Iterations', 'Krylov Iterations', 'Krylov Iterations per Step', 'Krylov_vs_step.png')
    plot_combined('Step', 'Execution Time', 'Newton Iterations', 'Exectution Time Per Iteration (s)', 'Average Time Per Iteration per Step', 'Execution_vs_step.png')


    # Create a bar chart for average setup and execution times
    avg_df = pd.DataFrame(averages)
    ax = avg_df.plot(kind='bar', x='Run', y=['Average Setup Time', 'Average Execution Time'], 
                title='Average Setup and Execution Time per Run', 
                ylabel='Time (s)', figsize=(8, 6))

    for container in ax.containers:
        ax.bar_label(container, fmt='%.4f', label_type='edge')  # Format with 4 decimal places

    # Customize legend
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, title='Metrics')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'average_setup_execution_time.png'))
    plt.close()
    
# Example usage
csv_files = [
    #'timings/run1.csv',  # Path to the first CSV file
    'timings/HypreParasails_SetFilter(0,055).csv',  # Path to the second CSV file
]
output_dir = 'generated_graphs'

# Generate combined graphs
generate_combined_graphs(csv_files, output_dir)