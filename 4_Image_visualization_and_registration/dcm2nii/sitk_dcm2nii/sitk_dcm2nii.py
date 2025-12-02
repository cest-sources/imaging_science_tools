from pathlib import Path
from tkinter import filedialog
import os
import SimpleITK as sitk

def dcm2nii(folder_path = None):
    """
    Converts DICOM series in a specified folder to NIfTI (.nii) format using SimpleITK.

    Parameters
    ----------
    folder_path : str or None, optional
        Path to the parent folder containing a 'dcm' subfolder with DICOM files. If None, a folder selection dialog will open.

    Behavior
    --------
    - Searches for DICOM series in the 'dcm' subfolder of the provided folder_path.
    - Converts each DICOM series found to a NIfTI file and saves it in a 'nii' subfolder within folder_path.
    - Prints the type of the converted image and the output filename for each successful conversion.
    - Prints an error message if a series cannot be converted.
    """
    
    if folder_path == None:
        import tkinter as tk
        root = tk.Tk()
        folder_path = filedialog.askdirectory()
        root.destroy()

    dicom_directory = folder_path + '/dcm'
    output_folder = folder_path + '/nii'

    os.makedirs(output_folder, exist_ok=True)
        # Walk through the directory and find DICOM series
    for root, _, _ in os.walk(dicom_directory):
        reader = sitk.ImageSeriesReader()
        dicom_series = reader.GetGDCMSeriesIDs(root)

        if not dicom_series:
            continue  # Skip if no DICOM series found in this folder

        for series_id in dicom_series:
            dicom_files = reader.GetGDCMSeriesFileNames(root, series_id)
            reader.SetFileNames(dicom_files)
                
            try:
                image = reader.Execute()
                output_filename = os.path.join(output_folder, f"{os.path.basename(root)}.nii")
                sitk.WriteImage(image, output_filename)
                print(type(image))
                print(f"Converted: {output_filename}")              
                
            except Exception as e:
                print(f"Error converting {series_id} in {root}: {e}")
    
    return
