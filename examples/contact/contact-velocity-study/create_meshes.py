import shutil
import subprocess
import os

cubit_path = '/Applications/Cubit-16.18/Cubit.app/Contents/MacOS/cubit'

def run_cubit_with_replacement(jou_file, nel_value, output_folder):
    # Define the new jou file name
    new_jou_file = 'modified_' + os.path.basename(jou_file)

    # Copy the original .jou file to a new file
    shutil.copy(jou_file, new_jou_file)

    # Read the contents of the new .jou file and replace nel values
    with open(new_jou_file, 'r') as file:
        content = file.read()

    # Replace nel_fine and nel_coarse with the input value
    content = content.replace('${nel_fine = 1}', f'${{nel_fine = {nel_value}}}')
    content = content.replace('${nel_coarse = 1}', f'${{nel_coarse = {nel_value}}}')

    # Write the modified content back to the new .jou file
    with open(new_jou_file, 'w') as file:
        file.write(content)

    # Run the Cubit command
    try:
        subprocess.run([cubit_path, '-batch', '-nojournal', new_jou_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Cubit: {e}")
        return

    # Move the output files to the specified output folder
    output_files = ['bar-1.g', 'bar-2.g']
    for output_file in output_files:
        if os.path.exists(output_file):
            shutil.move(output_file, os.path.join(output_folder, output_file))
            print(f"Moved {output_file} to {output_folder}")
        else:
            print(f"{output_file} does not exist.")

# Example usage
if __name__ == "__main__":
    jou_file_path = 'bars.jou'  # Path to your .jou file

    nel_value = [ 1, 2, 4, 8 ]
    output_folder = [ f'mesh_{n}' for n in range(4) ]

    for nv, fol in zip(nel_value, output_folder):
        # Create the output folder if it doesn't exist
        os.makedirs(fol, exist_ok=True)

        # Run the function
        run_cubit_with_replacement(jou_file_path, nv, fol)