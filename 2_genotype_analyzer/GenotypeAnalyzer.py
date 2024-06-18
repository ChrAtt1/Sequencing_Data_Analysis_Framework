# Import necessary libraries
import tkinter as tk
import pandas as pd
from tkinter import filedialog, ttk, messagebox
from pathlib import Path
import time
from datetime import timedelta
# Import necessary modules for genotype analysis
import step_1_build_genotype
import step_2_merge_files
import step_3_analyse_file
import step_4_color_count

# Load the Excel file containing the settings
setting_file = pd.ExcelFile('Setting_Genotype-Analyzer.xlsx')
# Parse 'Sequencing Method Setting' sheet into a DataFrame
df1 = setting_file.parse('Sequencing Method Setting')
# Convert the 'Sequencing Method' column to a list of methods
methods = df1['Sequencing Method'].tolist()
# Parse 'Amplicon Setting' sheet into another DataFrame
df2 = setting_file.parse('Amplicon Setting')


# Function to open a file dialog and get the selected file path
# Inserts the file path into the provided entry field
def browse_file(entry, filetypes=(("all files", "*.*"),)):
    file_path = filedialog.askopenfilename(filetypes=filetypes)
    entry.delete(0, tk.END)
    entry.insert(0, file_path)


# Function to open a directory dialog and get the selected directory path
# Inserts the directory path into the provided entry field
def browse_folder(entry):
    folder_path = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, folder_path)


# Dictionary to store the StringVar references
limit_vars = {}


# Function to create two entry fields in the GUI for a given feature (name)
# These entry fields are used to specify the lower and upper limits for the feature
def create_limit_entries(name, row, root_data):
    tk.Label(root_data, text=f"{name}:").grid(row=row, column=0, sticky="e")
    lower_limit_var = tk.StringVar(value="0")
    upper_limit_var = tk.StringVar(value="15000")
    tk.Entry(root_data, textvariable=lower_limit_var).grid(row=row, column=1)
    tk.Entry(root_data, textvariable=upper_limit_var).grid(row=row, column=2)
    limit_vars[name] = {"lower": lower_limit_var, "upper": upper_limit_var}


# Function to check whether the input paths provided in the GUI are valid
# If the paths are valid, it calls the function to start the genotype building process
def check_paths_and_build_genotypes():
    ref_file_path = Path(ref_path_var.get())
    output_folder = output_folder_var.get()

    # Validate the paths and build genotypes based on the selected option
    # For option "A", VCF files are processed
    # For option "B", Excel files are processed
    if not ref_file_path.is_file() or not ref_file_path.name.lower().endswith(".fa"):
        messagebox.showwarning("Warning", "Please enter a valid reference sequence (.fa) file.")
        return
    if output_folder == '':
        messagebox.showwarning("Warning", "The provided output folder path is invalid. Please check it.")
        return
    if option.get() == "A":
        vcf_method1 = Path(vcf_path_entry_method1.get())
        vcf_method2 = Path(vcf_path_entry_method2.get())
        if not vcf_method1.is_dir() or not vcf_method2.is_dir():
            messagebox.showwarning("Warning", "The provided paths are invalid. Please check them.")
            return
        vcf_files_method1 = list(vcf_method1.glob("*.vcf"))
        vcf_files_method2 = list(vcf_method2.glob("*.vcf"))
        if not vcf_files_method1 or not vcf_files_method2:
            messagebox.showwarning("Warning", "Please enter valid vcf folder files.")
            return
    elif option.get() == "B":
        excel_path_method1 = Path(excel_path_var_method1.get())
        excel_path_method2 = Path(excel_path_var_method2.get())
        if not all(path.is_file() and path.name.lower().endswith(".xlsx") for path in
                   [excel_path_method1, excel_path_method2]):
            messagebox.showwarning("Warning", "Please enter valid excel (.xlsx) files.")
            return
    build_genotypes()


# Function to orchestrate the genotype building process
def build_genotypes():
    # Dictionary to store amplicon settings data
    amplicon_data = {}
    # Populate amplicon_data with rows from df2
    for i, row in df2.iterrows():
        amplicon_data[row['Amplicon Name']] = {
            'Lower Limit': row['Lower Limit'],
            'Upper Limit': row['Upper Limit'],
            'Sequencing Direction': row['Gene Direction (forward/backward)'],
            'Offset NCBI': row['Offset NCBI'],
            'Offset CDS': row['Offset CDS']
        }

    # Helper function to measure execution time of a function
    def measure_execution_time(func, *args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        minutes, seconds = divmod(execution_time, 60)
        print(f"{func.__name__} execution time: {int(minutes)} minutes {seconds} seconds")
        return result, execution_time

    path_output = output_folder_var.get()
    total_execution_time = 0.0

    # If the selected option is 'A', run the genotype building for methods
    if option.get() == "A":
        excel_path_method1, exec_time1 = measure_execution_time(step_1_build_genotype.setting_phasing_tool,
                                                                vcf_path_var_method1.get(), path_output,
                                                                methods[0], amplicon_data)
        excel_path_method2, exec_time2 = measure_execution_time(step_1_build_genotype.setting_phasing_tool,
                                                                vcf_path_var_method2.get(),
                                                                path_output, methods[1], amplicon_data)
        total_execution_time += exec_time1 + exec_time2
    else:
        # Otherwise, just retrieve the path for the Excel files
        excel_path_method1 = excel_path_var_method1.get()
        excel_path_method2 = excel_path_var_method2.get()

    # Update progress after genotype building
    time.sleep(5)  # pause for 5 seconds
    update_progress()
    time.sleep(5)  # pause for 5 seconds


    # Merge files and update progress
    file_name_merged, exec_time3 = measure_execution_time(step_2_merge_files.merging_tool_setting, amplicon_data,
                                                          methods, excel_path_method1, excel_path_method2, path_output)
    total_execution_time += exec_time3
    time.sleep(5)  # pause for 5 seconds
    update_progress()
    time.sleep(5)  # pause for 5 seconds

    # Analyze the merged file and update progress
    file_name_analysed, exec_time4 = measure_execution_time(step_3_analyse_file.analysis_tool, amplicon_data, methods,
                                                            file_name_merged, path_output, ref_path_var.get(),
                                                            qual_filter_var_method1.get(), dp_filter_var_method1.get(),
                                                            qual_filter_var_method2.get(), dp_filter_var_method2.get())

    total_execution_time += exec_time4
    time.sleep(5)  # pause for 5 seconds
    update_progress()
    time.sleep(5)  # pause for 5 seconds
    # Count colors in the analyzed file and update progress
    exec_time5 = measure_execution_time(step_4_color_count.fill_row_and_color_count_setting, amplicon_data, methods,
                                        file_name_analysed, path_output)
    total_execution_time += exec_time5[1]
    time.sleep(5)  # pause for 5 seconds
    update_progress()
    time.sleep(5)  # pause for 5 seconds

    # Print total analysis completion and total execution time
    print('Analysing finished!')
    print(f"Total execution time: " + str(timedelta(seconds=total_execution_time)) + " seconds")


# Function to increase the progress bar value by 25 units
# Ensures that the GUI updates to reflect the new progress bar value
def update_progress():
    # Add 25 units to the progress bar value
    progress['value'] += 25
    root.update_idletasks()


# Create a GUI window
root = tk.Tk()
root.title("Genotype Analyser")

# Declare variables for path entry fields
vcf_path_var = tk.StringVar()
ref_path_var = tk.StringVar()
excel_path_var_method1 = tk.StringVar()
excel_path_var_method2 = tk.StringVar()

# Option A: Phase VCF file
option = tk.StringVar(value="A")
tk.Radiobutton(root, text="Option A: Phase VCF files", variable=option, value="A").grid(
    row=0, column=0, columnspan=2, sticky="w")

# VCF file path entries for method1 and method2, with Browse button to select file path
tk.Label(root, text="Path " + methods[0] + ' (.vcf-file)' + ":").grid(row=1, column=0, sticky="e")
vcf_path_var_method1 = tk.StringVar()
vcf_path_entry_method1 = tk.Entry(root, textvariable=vcf_path_var_method1)
vcf_path_entry_method1.grid(row=1, column=1)
tk.Button(root, text="Browse", command=lambda: browse_folder(vcf_path_entry_method1)).grid(row=1, column=2)

# VCF file path method2
tk.Label(root, text="Path " + methods[1] + ' (.vcf-file)' + ":").grid(row=2, column=0, sticky="e")
vcf_path_var_method2 = tk.StringVar()
vcf_path_entry_method2 = tk.Entry(root, textvariable=vcf_path_var_method2)
vcf_path_entry_method2.grid(row=2, column=1)
tk.Button(root, text="Browse", command=lambda: browse_folder(vcf_path_entry_method2)).grid(row=2, column=2)

# Option B
tk.Radiobutton(root, text="Option B: Use previous Genotype Analyzer output excel files", variable=option,
               value="B").grid(row=4, column=0,
                               columnspan=2,
                               sticky="w")

# Excel file path method1
tk.Label(root, text="Path " + methods[0] + ' (.xlsx-file)' + ":").grid(row=5, column=0, sticky="e")
excel_entry_method1 = tk.Entry(root, textvariable=excel_path_var_method1)
excel_entry_method1.grid(row=5, column=1)
tk.Button(root, text="Browse", command=lambda: browse_file(excel_entry_method1)).grid(row=5, column=2)

# Excel file path method2
tk.Label(root, text="Path " + methods[1] + ' (.xlsx-file)' + ":").grid(row=6, column=0, sticky="e")
excel_entry_method2 = tk.Entry(root, textvariable=excel_path_var_method2)
excel_entry_method2.grid(row=6, column=1)
tk.Button(root, text="Browse", command=lambda: browse_file(excel_entry_method2)).grid(row=6, column=2)

# Filter settings
tk.Label(root, text="Filter settings:").grid(row=7, column=0, sticky="w")

# Reference sequence file path
tk.Label(root, text="Reference sequence file path (.fa-file):").grid(row=9, column=0, sticky="e")
ref_path_entry = tk.Entry(root, textvariable=ref_path_var)
ref_path_entry.grid(row=9, column=1)
tk.Button(root, text="Browse", command=lambda: browse_file(ref_path_entry)).grid(row=9, column=2)

# Define a bold font
bold_font = ("Helvetica", 10, "bold")

# Qual and DP entries
tk.Label(root, text=methods[0] + ":", font=bold_font).grid(row=10, column=1, sticky="w")
tk.Label(root, text=methods[1] + ":", font=bold_font).grid(row=10, column=2, sticky="w")

# Qual and DP entries
qual_filter_var_method1 = tk.StringVar(value="0")
dp_filter_var_method1 = tk.StringVar(value="0")
tk.Entry(root, textvariable=qual_filter_var_method1).grid(row=11, column=1)
tk.Entry(root, textvariable=dp_filter_var_method1).grid(row=12, column=1)

# Qual and DP entries
tk.Label(root, text="Qual:").grid(row=11, column=0, sticky="e")
tk.Label(root, text="DP:").grid(row=12, column=0, sticky="e")
qual_filter_var_method2 = tk.StringVar(value="0")
dp_filter_var_method2 = tk.StringVar(value="0")
tk.Entry(root, textvariable=qual_filter_var_method2).grid(row=11, column=2)
tk.Entry(root, textvariable=dp_filter_var_method2).grid(row=12, column=2)

# Merge settings
tk.Label(root, text="Haplotype merge settings:").grid(row=17, column=0, sticky="w")

# Qual and DP entries
tk.Label(root, text="Qual: ").grid(row=17, column=0, sticky="e")
tk.Label(root, text="DP: ").grid(row=18, column=0, sticky="e")
qual_merge_var_method1 = tk.StringVar(value="0")
dp_merge_var_method1 = tk.StringVar(value="0")
tk.Entry(root, textvariable=qual_merge_var_method1).grid(row=17, column=1)
tk.Entry(root, textvariable=dp_merge_var_method1).grid(row=18, column=1)

# Qual and DP entries
tk.Label(root, text="Qual: ").grid(row=17, column=0, sticky="e")
tk.Label(root, text="DP: ").grid(row=18, column=0, sticky="e")
qual_merge_var_method2 = tk.StringVar(value="0")
dp_merge_var_method2 = tk.StringVar(value="0")
tk.Entry(root, textvariable=qual_merge_var_method2).grid(row=17, column=2)
tk.Entry(root, textvariable=dp_merge_var_method2).grid(row=18, column=2)

# Output folder
tk.Label(root, text="Output folder: ").grid(row=19, column=0, sticky="w")
output_folder_var = tk.StringVar()
output_folder_entry = tk.Entry(root, textvariable=output_folder_var)
output_folder_entry.grid(row=20, column=1)
tk.Button(root, text="Browse", command=lambda: browse_folder(output_folder_entry)).grid(row=20, column=2)

# Progress bar
progress = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
progress.grid(row=21, column=0, columnspan=3, pady=10)

# Button to start analysis
tk.Button(root, text="Analyse Genotypes", command=check_paths_and_build_genotypes).grid(row=22, column=1, sticky="e")

# Start GUI loop
root.mainloop()
