import pydicom
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.widgets as widgets

# --- Configuration Constants ---
# Define the directory containing your DICOM series.
# Use a raw string literal (r'') to handle backslashes correctly on Windows paths.
SERIES_PATH = r'path/to/your/dicom/series'

# Default Window/Level settings for displaying CT images (common for soft tissue).
DEF_WC = 40.0
DEF_WW = 400.0

# Matplotlib plot configuration.
FIG_TITLE = "DICOM Series 3D Viewer"
FIG_SIZE = (13, 7) # Inches (width, height)

# Slider layout parameters.
SL_W = 0.17
SL_H = 0.03
SL_ROW1_BOT = 0.18
SL_ROW2_BOT = 0.08
SL_HSPACE = 0.15
SL_LEFT_MARGIN = 0.12
SL_COL = 'lightgoldenrodyellow'

# --- DICOM Loading and Processing Functions ---

def load_dicom_files(directory_path):
    """
    Loads valid DICOM image files from the specified directory.

    Args:
        directory_path (str): The path to the directory containing DICOM files.

    Returns:
        list: A list of pydicom.Dataset objects, each representing a valid DICOM slice.
    """
    if not os.path.isdir(directory_path):
        print(f"Error: DICOM series directory not found at {directory_path}")
        print("Please ensure the path is correct and the directory exists.")
        return None

    print(f"Loading DICOM files from: {directory_path}")
    dicom_files = []
    for filename in os.listdir(directory_path):
        filepath = os.path.join(directory_path, filename)

        if os.path.isfile(filepath) and (filename.endswith('.dcm') or '.' not in filename):
            try:
                ds = pydicom.dcmread(filepath)
                # Ensure it's an image slice with spatial positioning information.
                if 'PixelData' in ds and hasattr(ds, 'ImagePositionPatient'):
                    dicom_files.append(ds)
            except Exception:
                # Silently skip files that are not valid DICOMs or are corrupted.
                pass

    if not dicom_files:
        print("No valid DICOM image files found in the specified directory.")
        return None

    print(f"Found {len(dicom_files)} valid DICOM image files.")
    return dicom_files

def sort_dicom_slices(dicom_files):
    """
    Sorts DICOM slices primarily by ImagePositionPatient[2] (Z-axis),
    falling back to SliceLocation or SOPInstanceUID if necessary.

    Args:
        dicom_files (list): A list of pydicom.Dataset objects.
    """
    try:
        dicom_files.sort(key=lambda x: float(x.ImagePositionPatient[2]))
        print("Slices sorted by ImagePositionPatient[2].")
    except AttributeError:
        try:
            dicom_files.sort(key=lambda x: float(x.SliceLocation))
            print("Slices sorted by SliceLocation.")
        except AttributeError:
            print("Warning: Neither ImagePositionPatient nor SliceLocation found for sorting. Slices may be misordered.")
            dicom_files.sort(key=lambda x: x.SOPInstanceUID) # Fallback for consistent (but not spatial) order

def build_3d_volume(dicom_files):
    """
    Extracts pixel data from sorted DICOM files and reconstructs a 3D volume.
    Applies Rescale Slope and Intercept for Hounsfield Unit conversion.

    Args:
        dicom_files (list): Sorted list of pydicom.Dataset objects.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: The 3D reconstructed volume (Z, Y, X) in Hounsfield Units.
            - dict: A dictionary of aspect ratios for axial, sagittal, and coronal views.
    """
    first_slice = dicom_files[0]
    image_height = first_slice.Rows
    image_width = first_slice.Columns
    num_slices = len(dicom_files)

    # Determine voxel spacing.
    pixel_spacing = [float(s) for s in first_slice.PixelSpacing]
    try:
        slice_thickness = float(first_slice.SliceThickness)
    except AttributeError:
        print("Warning: SliceThickness not found. Defaulting to 1.0 for Z-spacing.")
        slice_thickness = 1.0

    # Calculate aspect ratios for correct display proportions in imshow.
    aspect_ratios = {
        'axial': pixel_spacing[1] / pixel_spacing[0],  # Y-spacing / X-spacing
        'sagittal': pixel_spacing[0] / slice_thickness, # X-spacing / Z-spacing
        'coronal': pixel_spacing[1] / slice_thickness  # Y-spacing / Z-spacing
    }

    volume = np.zeros((num_slices, image_height, image_width), dtype=np.float32)

    for i, ds in enumerate(dicom_files):
        slice_data = ds.pixel_array.astype(np.float32)

        if 'RescaleIntercept' in ds and 'RescaleSlope' in ds:
            rescale_intercept = ds.RescaleIntercept
            rescale_slope = ds.RescaleSlope
            slice_data = slice_data * rescale_slope + rescale_intercept

        volume[i, :, :] = slice_data

    print(f"3D volume created with dimensions: {volume.shape} (slices, height, width)")
    print(f"Voxel spacing: X={pixel_spacing[0]:.2f}, Y={pixel_spacing[1]:.2f}, Z={slice_thickness:.2f}")

    return volume, aspect_ratios

def apply_window_level(image_data, window_center, window_width):
    """
    Applies window and level transformation to image data for display.

    Args:
        image_data (np.ndarray): The input image data (e.g., Hounsfield Units).
        window_center (float): The center of the display window.
        window_width (float): The width of the display window.

    Returns:
        np.ndarray: The image data scaled to the [0, 1] range for display.
    """
    min_val = window_center - window_width / 2
    max_val = window_center + window_width / 2

    display_image = np.clip(image_data, min_val, max_val)

    if window_width == 0:
        return np.zeros_like(display_image)

    return (display_image - min_val) / (max_val - min_val)

# --- Main Application Logic ---

def main():
    """
    Main function to load DICOM data, build a 3D volume, and display it
    with interactive slicing and window/level adjustment using Matplotlib.
    """
    dicom_files = load_dicom_files(SERIES_PATH)
    if dicom_files is None:
        return # Exit if no DICOM files are found or directory is invalid

    sort_dicom_slices(dicom_files)

    volume, aspect_ratios = build_3d_volume(dicom_files)

    # Calculate initial slice indices for display (middle of each dimension).
    num_slices, image_height, image_width = volume.shape
    initial_axial_idx = num_slices // 2
    initial_sagittal_idx = image_width // 2
    initial_coronal_idx = image_height // 2

    # --- Setup Matplotlib Plot ---
    fig, axes = plt.subplots(1, 3, figsize=FIG_SIZE, num=FIG_TITLE, )
    fig.subplots_adjust(bottom=0.25, wspace=0.3)

    ax_axial, ax_sagittal, ax_coronal = axes

    # Display initial slices.
    # Note: im_ objects store the Image object, allowing for data updates later.
    im_axial = ax_axial.imshow(
        apply_window_level(volume[initial_axial_idx, :, :], DEF_WC, DEF_WW),
        cmap='gray',
        aspect=aspect_ratios['axial']
    )
    ax_axial.set_title(f"Axial View (Slice {initial_axial_idx+1}/{num_slices})")
    ax_axial.axis('off')

    im_sagittal = ax_sagittal.imshow(
        apply_window_level(volume[:, initial_sagittal_idx, :], DEF_WC, DEF_WW),
        cmap='gray',
        aspect=aspect_ratios['sagittal']
    )
    ax_sagittal.set_title(f"Sagittal View (Slice {initial_sagittal_idx+1}/{image_width})")
    ax_sagittal.axis('off')

    im_coronal = ax_coronal.imshow(
        apply_window_level(volume[:, :, initial_coronal_idx], DEF_WC, DEF_WW),
        cmap='gray',
        aspect=aspect_ratios['coronal']
    )
    ax_coronal.set_title(f"Coronal View (Slice {initial_coronal_idx+1}/{image_height})")
    ax_coronal.axis('off')

    # --- Create Sliders ---
    # Axial Slice Slider
    axial_slider_ax = plt.axes([SL_LEFT_MARGIN, SL_ROW1_BOT, SL_W, SL_H], facecolor=SL_COL)
    axial_slider = widgets.Slider(
        ax=axial_slider_ax, label='Axial Slice',
        valmin=0, valmax=num_slices - 1, valinit=initial_axial_idx, valstep=1
    )

    # Sagittal Slice Slider
    sagittal_slider_ax = plt.axes([SL_LEFT_MARGIN + SL_W + SL_HSPACE, SL_ROW1_BOT, SL_W, SL_H], facecolor=SL_COL)
    sagittal_slider = widgets.Slider(
        ax=sagittal_slider_ax, label='Sagittal Slice',
        valmin=0, valmax=image_width - 1, valinit=initial_sagittal_idx, valstep=1
    )

    # Coronal Slice Slider
    coronal_slider_ax = plt.axes([SL_LEFT_MARGIN + 2 * (SL_W + SL_HSPACE), SL_ROW1_BOT, SL_W, SL_H], facecolor=SL_COL)
    coronal_slider = widgets.Slider(
        ax=coronal_slider_ax, label='Coronal Slice',
        valmin=0, valmax=image_height - 1, valinit=initial_coronal_idx, valstep=1
    )

    # Window Center Slider
    wc_slider_ax = plt.axes([SL_LEFT_MARGIN, SL_ROW2_BOT, SL_W, SL_H], facecolor=SL_COL)
    wc_slider = widgets.Slider(
        ax=wc_slider_ax, label='Window Center',
        valmin=volume.min(), valmax=volume.max(), valinit=DEF_WC
    )

    # Window Width Slider
    ww_slider_ax = plt.axes([SL_LEFT_MARGIN + SL_W + SL_HSPACE, SL_ROW2_BOT, SL_W, SL_H], facecolor=SL_COL)
    ww_slider = widgets.Slider(
        ax=ww_slider_ax, label='Window Width',
        valmin=1, valmax=(volume.max() - volume.min()) * 2, valinit=DEF_WW
    )

    # --- Update Function for Sliders ---
    def update(val):
        """
        Callback function to update all displayed images when any slider value changes.
        """
        current_wc = wc_slider.val
        current_ww = ww_slider.val

        # Update Axial View
        axial_idx = int(axial_slider.val)
        im_axial.set_data(apply_window_level(volume[axial_idx, :, :], current_wc, current_ww))
        ax_axial.set_title(f"Axial View (Slice {axial_idx+1}/{num_slices})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")

        # Update Sagittal View
        sagittal_idx = int(sagittal_slider.val)
        im_sagittal.set_data(apply_window_level(volume[:, sagittal_idx, :], current_wc, current_ww))
        ax_sagittal.set_title(f"Sagittal View (Slice {sagittal_idx+1}/{image_width})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")

        # Update Coronal View
        coronal_idx = int(coronal_slider.val)
        im_coronal.set_data(apply_window_level(volume[:, :, coronal_idx], current_wc, current_ww))
        ax_coronal.set_title(f"Coronal View (Slice {coronal_idx+1}/{image_height})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")

        fig.canvas.draw_idle()

    # Connect the update function to each slider's 'on_changed' event.
    axial_slider.on_changed(update)
    sagittal_slider.on_changed(update)
    coronal_slider.on_changed(update)
    wc_slider.on_changed(update)
    ww_slider.on_changed(update)

    # Initialize the display with current slider values (important for initial title consistency).
    update(None) # Pass None as a dummy argument, as 'val' is not used directly.

    plt.show()

if __name__ == "__main__":
    main()