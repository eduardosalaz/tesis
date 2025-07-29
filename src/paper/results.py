#!/usr/bin/env python3
"""
Optimization Results Parser
Processes optimization run results from text files and generates summary statistics
"""

import os
import re
import csv
from collections import defaultdict
from pathlib import Path
import argparse
from tabulate import tabulate

class OptimizationResult:
    """Class to store optimization result data"""
    def __init__(self, instance_id, p_value, bus, centers, instance_set, 
                 objective, bound, gap, runtime):
        self.instance_id = instance_id
        self.p_value = p_value
        self.bus = bus
        self.centers = centers
        self.instance_set = instance_set
        self.objective = objective
        self.bound = bound
        self.gap = gap * 100  # Convert to percentage
        self.runtime = runtime

def parse_filename(filename):
    """Extract information from filename"""
    # Pattern: results_BUs_centers_p_first_instance_set_size_runtime_s_id
    pattern = r'results_(\d+)_(\d+)_(\d+)_first_instance_set_(\w+)_[\d.]+s_(\d+)'
    match = re.match(pattern, filename)
    
    if match:
        bus = int(match.group(1))
        centers = int(match.group(2))
        p_value = int(match.group(3))
        instance_set = match.group(4)
        instance_id = int(match.group(5))
        
        return {
            'bus': bus,
            'centers': centers,
            'p_value': p_value,
            'instance_set': instance_set,
            'instance_id': instance_id
        }
    return None

def read_result_file(filepath):
    """Read the 4-line result file"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        if len(lines) >= 4:
            objective = float(lines[0].strip())
            bound = float(lines[1].strip())
            gap = float(lines[2].strip())
            runtime = float(lines[3].strip())
            
            return {
                'objective': objective,
                'bound': bound,
                'gap': gap,
                'runtime': runtime
            }
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
    
    return None

def process_results_folder(folder_path):
    """Process all result files in a folder"""
    results = []
    
    for filename in os.listdir(folder_path):
        if not filename.startswith('results_'):
            continue
            
        filepath = os.path.join(folder_path, filename)
        if not os.path.isfile(filepath):
            continue
            
        # Parse filename
        file_info = parse_filename(filename)
        if not file_info:
            print(f"Could not parse filename: {filename}")
            continue
            
        # Read file content
        file_data = read_result_file(filepath)
        if not file_data:
            print(f"Could not read file: {filename}")
            continue
            
        # Create result object
        result = OptimizationResult(
            instance_id=file_info['instance_id'],
            p_value=file_info['p_value'],
            bus=file_info['bus'],
            centers=file_info['centers'],
            instance_set=file_info['instance_set'],
            objective=file_data['objective'],
            bound=file_data['bound'],
            gap=file_data['gap'],
            runtime=file_data['runtime']
        )
        
        results.append(result)
    
    return results

def organize_results(results):
    """Organize results by instance ID"""
    organized = defaultdict(dict)
    
    for result in results:
        organized[result.instance_id][result.p_value] = result
    
    return organized

def calculate_statistics(results, by_p_value=False):
    """Calculate statistics for each metric, optionally grouped by p value"""
    if by_p_value:
        # Group results by p value
        p_groups = defaultdict(list)
        for result in results:
            p_groups[result.p_value].append(result)
        
        # Calculate stats for each p value
        stats_by_p = {}
        for p_value, p_results in p_groups.items():
            stats_by_p[p_value] = calculate_single_stats(p_results)
        
        return stats_by_p
    else:
        return calculate_single_stats(results)

def calculate_single_stats(results):
    """Calculate statistics for a single group of results"""
    stats = {
        'objective': {'values': [], 'avg': 0, 'min': float('inf'), 'max': float('-inf')},
        'bound': {'values': [], 'avg': 0, 'min': float('inf'), 'max': float('-inf')},
        'gap': {'values': [], 'avg': 0, 'min': float('inf'), 'max': float('-inf')},
        'runtime': {'values': [], 'avg': 0, 'min': float('inf'), 'max': float('-inf')}
    }
    
    for result in results:
        stats['objective']['values'].append(result.objective)
        stats['bound']['values'].append(result.bound)
        stats['gap']['values'].append(result.gap)
        stats['runtime']['values'].append(result.runtime)
    
    for metric in stats:
        if stats[metric]['values']:
            stats[metric]['avg'] = sum(stats[metric]['values']) / len(stats[metric]['values'])
            stats[metric]['min'] = min(stats[metric]['values'])
            stats[metric]['max'] = max(stats[metric]['values'])
    
    return stats

def write_csv(organized_results, output_file, instance_set_info):
    """Write results to CSV file"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        p_values = sorted(set(p for instance in organized_results.values() 
                            for p in instance.keys()))
        header = ['Instance_ID']
        for p in p_values:
            header.extend([f'p{p}_Objective', f'p{p}_Bound', f'p{p}_Gap(%)', f'p{p}_Runtime(s)'])
        writer.writerow(header)
        
        # Data rows
        all_results = []
        for instance_id in sorted(organized_results.keys()):
            row = [instance_id]
            for p in p_values:
                if p in organized_results[instance_id]:
                    result = organized_results[instance_id][p]
                    row.extend([
                        f"{result.objective:.2f}",
                        f"{result.bound:.2f}",
                        f"{result.gap:.2f}",
                        f"{result.runtime:.2f}"
                    ])
                    all_results.append(result)
                else:
                    row.extend(['N/A', 'N/A', 'N/A', 'N/A'])
            writer.writerow(row)
        
        # Statistics
        if all_results:
            # Statistics by p value
            stats_by_p = calculate_statistics(all_results, by_p_value=True)
            writer.writerow([])
            writer.writerow(['STATISTICS BY P VALUE'])
            
            for p_value in sorted(stats_by_p.keys()):
                writer.writerow([])
                writer.writerow([f'p = {p_value}'])
                writer.writerow(['Metric', 'Average', 'Min', 'Max'])
                stats = stats_by_p[p_value]
                writer.writerow(['Objective', f"{stats['objective']['avg']:.2f}", 
                               f"{stats['objective']['min']:.2f}", f"{stats['objective']['max']:.2f}"])
                writer.writerow(['Bound', f"{stats['bound']['avg']:.2f}", 
                               f"{stats['bound']['min']:.2f}", f"{stats['bound']['max']:.2f}"])
                writer.writerow(['Gap (%)', f"{stats['gap']['avg']:.2f}", 
                               f"{stats['gap']['min']:.2f}", f"{stats['gap']['max']:.2f}"])
                writer.writerow(['Runtime (s)', f"{stats['runtime']['avg']:.2f}", 
                               f"{stats['runtime']['min']:.2f}", f"{stats['runtime']['max']:.2f}"])
            
            # Total statistics
            total_stats = calculate_statistics(all_results)
            writer.writerow([])
            writer.writerow(['TOTAL STATISTICS (ALL P VALUES)'])
            writer.writerow(['Metric', 'Average', 'Min', 'Max'])
            writer.writerow(['Objective', f"{total_stats['objective']['avg']:.2f}", 
                           f"{total_stats['objective']['min']:.2f}", f"{total_stats['objective']['max']:.2f}"])
            writer.writerow(['Bound', f"{total_stats['bound']['avg']:.2f}", 
                           f"{total_stats['bound']['min']:.2f}", f"{total_stats['bound']['max']:.2f}"])
            writer.writerow(['Gap (%)', f"{total_stats['gap']['avg']:.2f}", 
                           f"{total_stats['gap']['min']:.2f}", f"{total_stats['gap']['max']:.2f}"])
            writer.writerow(['Runtime (s)', f"{total_stats['runtime']['avg']:.2f}", 
                           f"{total_stats['runtime']['min']:.2f}", f"{total_stats['runtime']['max']:.2f}"])
            
            # Instance set info
            if instance_set_info:
                writer.writerow([])
                writer.writerow(['Instance Set Info'])
                writer.writerow(['Set', 'BUs', 'Centers'])
                writer.writerow([instance_set_info['set'], instance_set_info['bus'], 
                               instance_set_info['centers']])

def print_table(organized_results):
    """Print results as formatted table"""
    p_values = sorted(set(p for instance in organized_results.values() 
                        for p in instance.keys()))
    
    # Prepare table data
    table_data = []
    all_results = []
    
    for instance_id in sorted(organized_results.keys()):
        row = [instance_id]
        for p in p_values:
            if p in organized_results[instance_id]:
                result = organized_results[instance_id][p]
                row.append(f"Obj: {result.objective:.2e}\n"
                          f"Bnd: {result.bound:.2e}\n"
                          f"Gap: {result.gap:.2f}%\n"
                          f"Time: {result.runtime:.0f}s")
                all_results.append(result)
            else:
                row.append("N/A")
        table_data.append(row)
    
    # Create headers
    headers = ['Instance'] + [f'p = {p}' for p in p_values]
    
    # Print main table
    print("\n" + "="*80)
    print("OPTIMIZATION RESULTS")
    print("="*80)
    print(tabulate(table_data, headers=headers, tablefmt='grid'))
    
    # Print statistics
    if all_results:
        # Statistics by p value
        stats_by_p = calculate_statistics(all_results, by_p_value=True)
        
        print("\n" + "="*80)
        print("STATISTICS BY P VALUE")
        print("="*80)
        
        for p_value in sorted(stats_by_p.keys()):
            print(f"\np = {p_value}:")
            stats = stats_by_p[p_value]
            
            stats_table = [
                ['Objective', f"{stats['objective']['avg']:.2e}", 
                 f"{stats['objective']['min']:.2e}", f"{stats['objective']['max']:.2e}"],
                ['Bound', f"{stats['bound']['avg']:.2e}", 
                 f"{stats['bound']['min']:.2e}", f"{stats['bound']['max']:.2e}"],
                ['Gap (%)', f"{stats['gap']['avg']:.2f}", 
                 f"{stats['gap']['min']:.2f}", f"{stats['gap']['max']:.2f}"],
                ['Runtime (s)', f"{stats['runtime']['avg']:.2f}", 
                 f"{stats['runtime']['min']:.2f}", f"{stats['runtime']['max']:.2f}"]
            ]
            
            print(tabulate(stats_table, headers=['Metric', 'Average', 'Min', 'Max'], 
                          tablefmt='grid'))
        
        # Total statistics
        total_stats = calculate_statistics(all_results)
        
        print("\n" + "="*80)
        print("TOTAL STATISTICS (ALL P VALUES)")
        print("="*80)
        
        total_stats_table = [
            ['Objective', f"{total_stats['objective']['avg']:.2e}", 
             f"{total_stats['objective']['min']:.2e}", f"{total_stats['objective']['max']:.2e}"],
            ['Bound', f"{total_stats['bound']['avg']:.2e}", 
             f"{total_stats['bound']['min']:.2e}", f"{total_stats['bound']['max']:.2e}"],
            ['Gap (%)', f"{total_stats['gap']['avg']:.2f}", 
             f"{total_stats['gap']['min']:.2f}", f"{total_stats['gap']['max']:.2f}"],
            ['Runtime (s)', f"{total_stats['runtime']['avg']:.2f}", 
             f"{total_stats['runtime']['min']:.2f}", f"{total_stats['runtime']['max']:.2f}"]
        ]
        
        print(tabulate(total_stats_table, headers=['Metric', 'Average', 'Min', 'Max'], 
                      tablefmt='grid'))

def main():
    parser = argparse.ArgumentParser(description='Process optimization results')
    parser.add_argument('--folder', type=str, default='.', 
                       help='Folder containing result files (default: current directory)')
    parser.add_argument('--csv', type=str, help='Output CSV filename')
    parser.add_argument('--no-table', action='store_true', 
                       help='Do not print table output')
    
    args = parser.parse_args()
    
    # Process results
    print(f"Processing files in: {args.folder}")
    results = process_results_folder(args.folder)
    
    if not results:
        print("No valid result files found!")
        return
    
    print(f"Found {len(results)} result files")
    
    # Organize results
    organized = organize_results(results)
    
    # Get instance set info from first result
    instance_set_info = None
    if results:
        instance_set_info = {
            'set': results[0].instance_set,
            'bus': results[0].bus,
            'centers': results[0].centers
        }
    
    # Output results
    if args.csv:
        write_csv(organized, args.csv, instance_set_info)
        print(f"\nCSV file written to: {args.csv}")
    
    if not args.no_table:
        print_table(organized)
    
    # Print instance set info
    if instance_set_info:
        print(f"\nInstance Set: {instance_set_info['set'].upper()}")
        print(f"BUs: {instance_set_info['bus']}, Centers: {instance_set_info['centers']}")

if __name__ == "__main__":
    main()