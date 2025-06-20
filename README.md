# My DICOM Series 3D Viewer

This is a little project that I put together: a simple, interactive 3D viewer for **DICOM medical image series**. I built it using `pydicom` for handling the data, `NumPy` for the heavy lifting with numbers, and `Matplotlib` for the interactive visuals.

I originally crafted this script to make exploring DICOM data a bit more intuitive, letting you slice through images and tweak their appearance on the fly.

![](3dViewer.gif)

## What Can It Do?

I've packed a few key features into this viewer:

  * **Easy DICOM Loading:** Just point it to a directory, and it'll load your DICOM files, understanding how they fit together spatially.
  * **3D Reconstruction:** It takes all those individual slices and stitches them into a full 3D volume. Plus, it converts them to **Hounsfield Units (HU)**, which is super helpful for consistent viewing across different scans.
  * **Interactive Slicing:** You can effortlessly slide through **axial, sagittal, and coronal** views.
  * **Window/Level Adjustment:** This is a big one for medical images. You can dynamically adjust the **Window Center (WC)** and **Window Width (WW)**. This lets you highlight different tissues—from bone to soft tissue—by changing the contrast and brightness.
  * **Correct Proportions:** I've made sure it automatically corrects for aspect ratios, so your the images always look geometrically accurate, no stretching or squishing.

-----

## Getting Started Is Simple

Here’s how you can get it running on your machine.

### 1\. Grab the Essentials

First things first, you'll need a few Python libraries. If you don't have them, you can quickly install them using pip:

```bash
pip install pydicom matplotlib numpy
```

### 2\. Point to Your DICOMs

This is the most crucial step\! You'll need to tell the script where your DICOM series are located. Open up the script and find this line:

```python
# Use a raw string literal (r'') to handle backslashes correctly on Windows paths.
SERIES_PATH = r'path/to/your/dicom/series' # <--- IMPORTANT: Update this path!
```

Just **change the path** inside the quotes `r'...'` to wherever your DICOM files live. For Windows, using an `r` before the path (a raw string) is a good habit to avoid issues with backslashes.

### 3\. Run It\!

Once the path is set, simply open your terminal or command prompt, navigate to where you saved the script, and run:

```bash
python your_script_name.py
```

(Replace `your_script_name.py` with whatever you named the file, like `dicom_viewer.py`).

-----

## How to Use It

Once it starts, you'll see a Matplotlib window pop up with three different views of your medical scan:

  * **Axial View:** The typical "slice through the body" view.
  * **Sagittal View:** A side-to-side cross-section.
  * **Coronal View:** A front-to-back cross-section.

Below these images, you'll find a set of sliders. These are the controls:

  * The **"Axial Slice," "Sagittal Slice," and "Coronal Slice"** sliders let you move through the different depths of your 3D scan.
  * The **"Window Center"** slider adjusts the overall brightness, helping you bring different details into focus.
  * The **"Window Width"** slider controls the contrast. Play with this to expand or compress the range of visible intensities, making specific tissues stand out.

I encourage you to play around with these\! It’s really the best way to get a feel for your data and see how different settings affect the image.

-----

## A Few Tips & Troubleshooting

I've tried to make it robust, but sometimes things happen\! Here are a few common isuues and what to do:

  * **"Error: DICOM series directory not found"**: This almost always means the `DICOM_SERIES_PATH` isn't quite right. Double-check your path, and remember that `r'...'` for Windows can save you a headache.
  * **"No valid DICOM image files found"**: Make sure the folder you're pointing to actually contains DICOM image files. Sometimes a directory might have other files mixed in, or the DICOMs might be in a subfolder.
  * **Images look completely black or white**: This is usually just the **Window Center** and **Window Width** needing a tweak. Different scans and tissues (like bone vs. soft tissue) require very different settings. Just move those sliders until you see something\!
  * **Slices seem out of order**: The script tries its best to sort slices based on standard DICOM tags like `ImagePositionPatient` or `SliceLocation`. If those aren't present, it falls back to a consistent but not necessarily spatial order. It's rare, but sometimes old or unusual DICOMs might behave this way.

-----

## Want to Contribute?

I'm always open to ideas and improvements\! If you have suggestions, spot a bug, or want to add a feature, feel free to fork the repository and send a pull request.