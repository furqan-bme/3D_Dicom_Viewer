import pydicom                  # Import the pydicom library for DICOM file handling.
import matplotlib.pyplot as plt # Import the pyplot module for plotting and visualization.
import numpy as np              # Import numpy for numerical operations, especially with arrays.
import os                       # Import the os module for interacting with the operating system (e.g., file paths).
import matplotlib.widgets as widgets # Import widgets for interactive elements like sliders.

# --- 1. Define the directory containing your DICOM series ---
# This line defines the file path to the directory holding the DICOM slices.
# Syntax: Raw string literal (r'') is used to prevent backslashes (\) from being interpreted
# as escape sequences (e.g., \n for newline, \t for tab), which is a common issue on Windows paths.
dicom_series_path = r'D:\VSCode\Projects\DATA\DICOMSamples\series-00000'


# Checks if the specified directory exists.
# os.path.isdir(path): Returns True if 'path' is an existing directory.
if not os.path.isdir(dicom_series_path):
    # If the directory does not exist, print an error message.
    print(f"Error: DICOM series directory not found at {dicom_series_path}")
    print("Please create this directory and place your 30 DICOM slices inside it,")
    print("or update the 'dicom_series_path' variable to your actual directory.")
    # Exit the script if the directory is not found, as further execution would fail.
    exit()

# --- 2. Load all DICOM files from the directory ---
print(f"Loading DICOM files from: {dicom_series_path}") # Informative message.
dicom_files = [] # Initialize an empty list to store loaded pydicom Dataset objects.

# Iterate over each filename found in the specified directory.
# os.listdir(path): Returns a list containing the names of the entries in the directory given by path.
for filename in os.listdir(dicom_series_path):
    # Construct the full file path by joining the directory path and the filename.
    # os.path.join(path, *paths): Concatenates path components intelligently.
    filepath = os.path.join(dicom_series_path, filename)
    
    # Check if the current item is a file (not a subdirectory) and looks like a DICOM file.
    # os.path.isfile(path): Returns True if 'path' is an existing regular file.
    # filename.endswith('.dcm') or '.' not in filename: Simple heuristic for DICOM files,
    # as some have '.dcm' extension and others might have no extension.
    if os.path.isfile(filepath) and (filename.endswith('.dcm') or '.' not in filename):
        try:
            # Attempt to read the DICOM file using pydicom.dcmread().
            # ds: A pydicom.Dataset object containing all DICOM tags and pixel data.
            ds = pydicom.dcmread(filepath)
            
            # Filter for actual image slices: Check if 'PixelData' tag exists (it's an image)
            # and 'ImagePositionPatient' exists (crucial for spatial sorting in 3D reconstruction).
            # hasattr(object, name): Returns True if the object has the given named attribute.
            if 'PixelData' in ds and hasattr(ds, 'ImagePositionPatient'):
                dicom_files.append(ds) # Add the loaded DICOM dataset to the list.
        except Exception as e:
            # Catch any exceptions during file reading (e.g., corrupted files, non-DICOM files).
            # print(f"Skipping non-DICOM or corrupted file: {filename} - {e}") # Debugging line.
            pass # Silently skip non-DICOM or corrupted files for cleaner output.

# After iterating through all files, check if any valid DICOM image files were found.
if not dicom_files:
    print("No valid DICOM image files found in the specified directory.")
    exit() # Exit if no images are found.

print(f"Found {len(dicom_files)} valid DICOM image files.") # Report the number of files found.

# --- 3. Sort slices by their Z-position (SliceLocation or ImagePositionPatient[2]) ---
# Physical Principle: To reconstruct a 3D volume, the 2D slices must be ordered correctly along the Z-axis.
# ImagePositionPatient[2] (the Z-component) or SliceLocation are used for this.
try:
    # Attempt to sort the list of DICOM datasets using the Z-component of ImagePositionPatient.
    # lambda x: float(x.ImagePositionPatient[2]): An anonymous function used as a key for sorting.
    # It extracts the ImagePositionPatient tag (a list/tuple of floats), takes its third element (index 2, Z-coordinate),
    # and converts it to a float for numerical sorting.
    dicom_files.sort(key=lambda x: float(x.ImagePositionPatient[2]))
    print("Slices sorted by ImagePositionPatient[2].")
except AttributeError:
    # If ImagePositionPatient is not found or is problematic, try SliceLocation as a fallback.
    try:
        dicom_files.sort(key=lambda x: float(x.SliceLocation))
        print("Slices sorted by SliceLocation.")
    except AttributeError:
        # If neither primary spatial tag is found, issue a warning.
        print("Warning: Neither ImagePositionPatient nor SliceLocation found for sorting. Slices may be misordered.")
        # As a last resort, sort by SOPInstanceUID, which provides a consistent but not necessarily spatial order.
        dicom_files.sort(key=lambda x: x.SOPInstanceUID)


# --- 4. Extract pixel data and build 3D volume ---
# Get dimensions from the first (now sorted) slice to define the 3D volume's shape.
first_slice = dicom_files[0]
image_height = first_slice.Rows    # DICOM tag (0028,0010) - Number of rows (Y-dimension).
image_width = first_slice.Columns  # DICOM tag (0028,0011) - Number of columns (X-dimension).
num_slices = len(dicom_files)      # Number of slices determines the Z-dimension.

# Determine voxel spacing for correct aspect ratio display in Matplotlib.
# Physical Principle: The physical size of each pixel/voxel determines the true proportions.
# Pixel Spacing (0028,0030): Represents the physical distance between pixel centers in mm, typically [column_spacing, row_spacing].
pixel_spacing = [float(s) for s in first_slice.PixelSpacing] # Convert list of strings to floats.

# Slice Thickness (0018,0050): Represents the thickness of the acquired slice in mm (Z-dimension spacing).
try:
    slice_thickness = float(first_slice.SliceThickness)
except AttributeError:
    # If SliceThickness is not available, default to 1.0mm.
    print("Warning: SliceThickness not found. Defaulting to 1.0 for Z-spacing.")
    slice_thickness = 1.0

# Calculate aspect ratios for imshow to ensure correct display proportions.
# Matplotlib imshow's 'aspect' parameter for (M, N) data: (width_pixel_size / height_pixel_size)
# Axial View (Y vs X): Displaying a 2D slice where rows are Y and columns are X.
aspect_axial = pixel_spacing[1] / pixel_spacing[0] # Y-spacing / X-spacing
# Sagittal View (Z vs X): Displaying a 2D slice where rows are Z and columns are X.
aspect_sagittal = pixel_spacing[0] / slice_thickness # X-spacing / Z-spacing
# Coronal View (Z vs Y): Displaying a 2D slice where rows are Z and columns are Y.
aspect_coronal = pixel_spacing[1] / slice_thickness # Y-spacing / Z-spacing


# Create an empty 3D NumPy array to store the reconstructed volume.
# The shape is (num_slices, image_height, image_width), representing (Z, Y, X) dimensions.
# dtype=np.float32: Important because Hounsfield Units (HU) are floating-point numbers,
# even if raw DICOM pixel data is integer. Using float32 saves memory over float64 while maintaining precision.
volume = np.zeros((num_slices, image_height, image_width), dtype=np.float32)

# Populate the 3D volume by iterating through the sorted DICOM files.
for i, ds in enumerate(dicom_files):
    # Get the raw pixel data for the current slice.
    # .astype(np.float32): Convert the raw pixel array to float32 *before* rescaling to prevent
    # potential integer overflow if RescaleSlope is large or to ensure floating-point arithmetic.
    slice_data = ds.pixel_array.astype(np.float32)

    # Mathematical Principle: Apply Rescale Slope and Intercept for HU conversion.
    # DICOM tags (0028,1053) RescaleIntercept and (0028,1052) RescaleSlope.
    if 'RescaleIntercept' in ds and 'RescaleSlope' in ds:
        rescale_intercept = ds.RescaleIntercept
        rescale_slope = ds.RescaleSlope
        # HU = Pixel Value * RescaleSlope + RescaleIntercept
        slice_data = slice_data * rescale_slope + rescale_intercept
    
    # Place the processed 2D slice data into the corresponding Z-position in the 3D volume array.
    volume[i, :, :] = slice_data

print(f"3D volume created with dimensions: {volume.shape} (slices, height, width)")
print(f"Voxel spacing: X={pixel_spacing[0]:.2f}, Y={pixel_spacing[1]:.2f}, Z={slice_thickness:.2f}")


# --- Window/Level Function ---
# This function applies the window and level transformation to image data for display.
# Mathematical Principle: Linear scaling of a specified intensity range to the display's full grayscale range.
def apply_window_level(image_data, window_center, window_width):
    # Calculate the minimum and maximum HU values to be displayed within the window.
    min_val = window_center - window_width / 2
    max_val = window_center + window_width / 2
    
    # Clip (clamp) the image data values to this min_val and max_val range.
    # np.clip(array, min, max): Values less than min become min; values greater than max become max.
    # Values within the range remain unchanged.
    display_image = np.clip(image_data, min_val, max_val)
    
    # Handle the edge case where window_width is zero to prevent division by zero.
    # If width is zero, all values are treated as black (or another default).
    if window_width == 0:
        return np.zeros_like(display_image) # Returns an array of zeros with same shape and type.
    
    # Normalize the clipped values to the range [0, 1].
    # This maps the windowed intensity range [min_val, max_val] to the display range [0, 1].
    # Matplotlib's imshow expects values in [0, 1] for floating point images for grayscale.
    display_image = (display_image - min_val) / (max_val - min_val)
    
    return display_image

# --- Default Window/Level Settings ---
# These are clinically common default settings for viewing specific tissue types in CT.
# Adjust these based on the type of scan data you are viewing (e.g., Lung, Bone, Brain).
default_window_center = 40.0  # Center of the window (e.g., for soft tissue).
default_window_width = 400.0  # Width of the window (e.g., for soft tissue contrast).
print(f"Applying default Window Center: {default_window_center}, Window Width: {default_window_width}")


# --- 5. Setup Matplotlib Plot for 3D Slicing ---
# Create a figure and a grid of subplots (1 row, 3 columns).
# figsize=(18, 7): Sets the figure size in inches for a wider layout.
figure_name = "DICOM_Series_to_3D"
fig, axes = plt.subplots(1, 3, figsize=(13, 7), num= figure_name)

# Adjust the subplot parameters for layout.
# bottom=0.25: Increases the space at the bottom of the figure to accommodate sliders.
# wspace=0.3: Adds horizontal space between subplots to prevent titles/images from overlapping.
fig.subplots_adjust(bottom=0.25, wspace=0.3)
fig.tight_layout()

# Assign each subplot (Axes object) to a descriptive variable.
ax_axial = axes[0]     # Left subplot for Axial view.
ax_sagittal = axes[1]  # Middle subplot for Sagittal view.
ax_coronal = axes[2]   # Right subplot for Coronal view.

# Calculate initial slice indices to display, typically the middle of the volume.
initial_axial_idx = num_slices // 2    # Integer division for middle slice in Z.
initial_sagittal_idx = image_width // 2  # Middle slice in X (for Sagittal).
initial_coronal_idx = image_height // 2 # Middle slice in Y (for Coronal).

# Display initial slices for each view.
# im_axial: Stores the Image object returned by imshow, allowing for later updates.
im_axial = ax_axial.imshow(
    apply_window_level(volume[initial_axial_idx, :, :], default_window_center, default_window_width),
    cmap='gray',          # Use a grayscale colormap.
    aspect=aspect_axial   # Apply calculated aspect ratio to maintain physical proportions.
)
ax_axial.set_title(f"Axial View (Slice {initial_axial_idx+1}/{num_slices})") # Set subplot title.
ax_axial.axis('off') # Turn off axis ticks and labels for cleaner image display.

im_sagittal = ax_sagittal.imshow(
    apply_window_level(volume[:, initial_sagittal_idx, :], default_window_center, default_window_width),
    cmap='gray',
    aspect=aspect_sagittal
)
ax_sagittal.set_title(f"Sagittal View (Slice {initial_sagittal_idx+1}/{image_width})")
ax_sagittal.axis('off')

im_coronal = ax_coronal.imshow(
    apply_window_level(volume[:, :, initial_coronal_idx], default_window_center, default_window_width),
    cmap='gray',
    aspect=aspect_coronal
)
ax_coronal.set_title(f"Coronal View (Slice {initial_coronal_idx+1}/{image_height})")
ax_coronal.axis('off')

# Create sliders for each view and for window/level adjustment.
axcolor = 'lightgoldenrodyellow' # Background color for slider axes.

# --- Slider Positioning Variables (Crucial for Layout) ---
# These variables define common dimensions and spacing for slider positioning.
slider_width = 0.2         # Width of each slider as a fraction of figure width.
slider_height = 0.02        # Height of each slider as a fraction of figure height.
row1_bottom = 0.15          # Vertical position (bottom) for the top row of sliders.
row2_bottom = 0.08          # Vertical position (bottom) for the second row of sliders.
horizontal_spacing = 0.13   # Horizontal space between adjacent sliders in a row.
left_margin = 0.1          # Starting left margin for the first slider in a row.

# Top row of sliders (Slice selection for Axial, Sagittal, Coronal)
# Each plt.axes([left, bottom, width, height]) defines a new axis for the slider.
# Axial Slice Slider
axial_slider_ax = plt.axes([left_margin, row1_bottom, slider_width, slider_height], facecolor=axcolor)
axial_slider = widgets.Slider(
    ax=axial_slider_ax,             # The Axes object where the slider is drawn.
    label='Axial Slice',            # Label displayed next to the slider.
    valmin=0,                       # Minimum value of the slider.
    valmax=num_slices - 1,          # Maximum value (0 to N-1 for array indices).
    valinit=initial_axial_idx,      # Initial position of the slider knob.
    valstep=1                       # Step size for the slider (ensures integer slice selection).
)

# Sagittal Slice Slider (positioned to the right of the Axial slider)
sagittal_slider_ax = plt.axes([left_margin + slider_width + horizontal_spacing, row1_bottom, slider_width, slider_height], facecolor=axcolor)
sagittal_slider = widgets.Slider(
    ax=sagittal_slider_ax,
    label='Sagittal Slice',
    valmin=0,
    valmax=image_width - 1,
    valinit=initial_sagittal_idx,
    valstep=1
)

# Coronal Slice Slider (positioned to the right of the Sagittal slider)
coronal_slider_ax = plt.axes([left_margin + 2 * (slider_width + horizontal_spacing), row1_bottom, slider_width, slider_height], facecolor=axcolor)
coronal_slider = widgets.Slider(
    ax=coronal_slider_ax,
    label='Coronal Slice',
    valmin=0,
    valmax=image_height - 1,
    valinit=initial_coronal_idx,
    valstep=1
)

# Bottom row of sliders (Global Window/Level adjustment)
# Window Center Slider
wc_slider_ax = plt.axes([left_margin, row2_bottom, slider_width, slider_height], facecolor=axcolor)
wc_slider = widgets.Slider(
    ax=wc_slider_ax,
    label='Window Center',
    valmin=volume.min(),              # Min slider value is the min HU in the volume.
    valmax=volume.max(),              # Max slider value is the max HU in the volume.
    valinit=default_window_center     # Initial value from default settings.
)

# Window Width Slider (positioned to the right of the Window Center slider)
ww_slider_ax = plt.axes([left_margin + slider_width + horizontal_spacing, row2_bottom, slider_width, slider_height], facecolor=axcolor)
ww_slider = widgets.Slider(
    ax=ww_slider_ax,
    label='Window Width',
    valmin=1,                         # Minimum window width is 1 (cannot be zero or negative).
    valmax=(volume.max() - volume.min()) * 2, # Max width is twice the total range of HU values.
    valinit=default_window_width      # Initial value from default settings.
)


# --- 6. Update Function for Sliders ---
# This function is called every time ANY of the sliders are moved.
# The 'val' argument is the specific new value of the slider that triggered the call,
# but we read all slider values directly for robustness.
def update(val):
    # Always read the current values from the Window Center and Window Width sliders.
    # This ensures that all image views are updated with the latest W/L settings,
    # regardless of which slice slider or W/L slider was moved.
    current_wc = wc_slider.val
    current_ww = ww_slider.val

    # Update Axial View:
    axial_idx = int(axial_slider.val) # Get the current integer slice index from the axial slider.
    # Apply window/level to the selected axial slice and update the image data.
    im_axial.set_data(apply_window_level(volume[axial_idx, :, :], current_wc, current_ww))
    # Update the title to reflect the current slice number and W/L settings.
    ax_axial.set_title(f"Axial View (Slice {axial_idx+1}/{num_slices})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")

    # Update Sagittal View:
    sagittal_idx = int(sagittal_slider.val) # Get current slice index from sagittal slider.
    # Select a sagittal slice (all Z, specific Y, all X) from the 3D volume.
    im_sagittal.set_data(apply_window_level(volume[:, sagittal_idx, :], current_wc, current_ww))
    ax_sagittal.set_title(f"Sagittal View (Slice {sagittal_idx+1}/{image_width})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")

    # Update Coronal View:
    coronal_idx = int(coronal_slider.val) # Get current slice index from coronal slider.
    # Select a coronal slice (all Z, all Y, specific X) from the 3D volume.
    im_coronal.set_data(apply_window_level(volume[:, :, coronal_idx], current_wc, current_ww))
    ax_coronal.set_title(f"Coronal View (Slice {coronal_idx+1}/{image_height})\nWC: {current_wc:.0f}, WW: {current_ww:.0f}")
    
    # Redraw the figure canvas to display the updated images.
    # fig.canvas.draw_idle(): Efficiently redraws the figure only when the UI thread is idle,
    # preventing excessive redraws during rapid slider movements.
    fig.canvas.draw_idle()

# Connect the update function to each slider's 'on_changed' event.
# This registers 'update' as the callback function to be executed when any slider's value changes.
axial_slider.on_changed(update)
sagittal_slider.on_changed(update)
coronal_slider.on_changed(update)
wc_slider.on_changed(update)
ww_slider.on_changed(update)

# Display the Matplotlib figure. This call blocks program execution until the figure window is closed.
plt.show()