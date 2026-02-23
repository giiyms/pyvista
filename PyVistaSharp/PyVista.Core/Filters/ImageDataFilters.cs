using PyVista.Core;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.ImageDataFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="ImageData"/> and provide image processing,
/// contouring, resampling, morphological operations, and frequency-domain filters.
/// </para>
/// </summary>
public static class ImageDataFilters
{
    // ---------------------------------------------------------------
    //  Smoothing
    // ---------------------------------------------------------------

    /// <summary>
    /// Performs Gaussian smoothing on the image data.
    /// </summary>
    /// <param name="self">The image data to smooth.</param>
    /// <param name="radiusFactor">
    /// Number of standard deviations used to compute the Gaussian kernel extent.
    /// </param>
    /// <param name="standardDeviation">
    /// Standard deviation of the Gaussian kernel. Can be a single value (applied to
    /// all dimensions) or a 3-element array for per-axis control.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array to smooth. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A smoothed <see cref="ImageData"/>.</returns>
    public static ImageData GaussianSmooth(this ImageData self, double radiusFactor = 1.5, double standardDeviation = 2.0, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("GaussianSmooth requires VTK vtkImageGaussianSmooth.");
    }

    /// <summary>
    /// Performs Gaussian smoothing with per-axis standard deviations.
    /// </summary>
    /// <param name="self">The image data to smooth.</param>
    /// <param name="radiusFactor">
    /// Number of standard deviations used to compute the Gaussian kernel extent.
    /// </param>
    /// <param name="standardDeviations">
    /// Per-axis standard deviations as a 3-element array <c>[sx, sy, sz]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array to smooth. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A smoothed <see cref="ImageData"/>.</returns>
    public static ImageData GaussianSmooth(this ImageData self, double radiusFactor, double[] standardDeviations, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(standardDeviations);
        if (standardDeviations.Length != 3)
        {
            throw new ArgumentException("Standard deviations must have exactly 3 elements.", nameof(standardDeviations));
        }

        throw new NotImplementedException("GaussianSmooth requires VTK vtkImageGaussianSmooth.");
    }

    /// <summary>
    /// Performs median smoothing on the image data.
    /// </summary>
    /// <param name="self">The image data to smooth.</param>
    /// <param name="kernelSize">
    /// Size of the median filter kernel as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array to smooth. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A smoothed <see cref="ImageData"/>.</returns>
    public static ImageData MedianSmooth(this ImageData self, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("MedianSmooth requires VTK vtkImageMedian3D.");
    }

    // ---------------------------------------------------------------
    //  Threshold
    // ---------------------------------------------------------------

    /// <summary>
    /// Thresholds the image data by replacing values outside the given range.
    /// </summary>
    /// <param name="self">The image data to threshold.</param>
    /// <param name="lowerThreshold">Lower threshold value. When <c>null</c>, uses the data minimum.</param>
    /// <param name="upperThreshold">Upper threshold value. When <c>null</c>, uses the data maximum.</param>
    /// <param name="inValue">Replacement value for voxels inside the threshold range.</param>
    /// <param name="outValue">Replacement value for voxels outside the threshold range.</param>
    /// <param name="scalars">
    /// Name of the scalar array to threshold. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A thresholded <see cref="ImageData"/>.</returns>
    public static ImageData ImageThreshold(this ImageData self, double? lowerThreshold = null, double? upperThreshold = null, double inValue = 1.0, double outValue = 0.0, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ImageThreshold requires VTK vtkImageThreshold.");
    }

    // ---------------------------------------------------------------
    //  Contouring
    // ---------------------------------------------------------------

    /// <summary>
    /// Generates labeled contour surfaces from a labeled (segmented) image.
    /// </summary>
    /// <param name="self">The labeled <see cref="ImageData"/>.</param>
    /// <param name="nLabels">
    /// Number of label values. When <c>null</c>, all unique labels are used.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <param name="smoothingIterations">
    /// Number of smoothing iterations applied to the output surface.
    /// </param>
    /// <returns>A <see cref="PolyData"/> of contour surfaces, one per label.</returns>
    public static PolyData ContourLabeled(this ImageData self, int? nLabels = null, string? scalars = null, int smoothingIterations = 0)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ContourLabeled requires VTK vtkDiscreteFlyingEdges3D.");
    }

    // ---------------------------------------------------------------
    //  Slicing
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts a 2-D slice of the image data along one axis at a given index.
    /// </summary>
    /// <param name="self">The image data to slice.</param>
    /// <param name="axis">Axis perpendicular to the slice: 0 (x), 1 (y), or 2 (z).</param>
    /// <param name="index">
    /// Index along the axis. When <c>null</c>, the midpoint index is used.
    /// </param>
    /// <returns>An <see cref="ImageData"/> containing the 2-D slice.</returns>
    public static ImageData SliceIndex(this ImageData self, int axis = 2, int? index = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (axis < 0 || axis > 2)
        {
            throw new ArgumentOutOfRangeException(nameof(axis), "Axis must be 0 (x), 1 (y), or 2 (z).");
        }

        throw new NotImplementedException("SliceIndex requires VTK vtkExtractVOI.");
    }

    // ---------------------------------------------------------------
    //  Subset extraction
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts a sub-volume from the image data using VOI (volume of interest) indices.
    /// </summary>
    /// <param name="self">The image data to extract from.</param>
    /// <param name="voi">
    /// Volume of interest as a 6-element array
    /// <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c> in index space.
    /// </param>
    /// <param name="sampleRate">
    /// Sampling rate per axis as a 3-element array <c>[rx, ry, rz]</c>.
    /// Defaults to <c>[1, 1, 1]</c>.
    /// </param>
    /// <returns>The extracted <see cref="ImageData"/> sub-volume.</returns>
    public static ImageData ExtractSubset(this ImageData self, int[] voi, int[]? sampleRate = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(voi);
        if (voi.Length != 6)
        {
            throw new ArgumentException("VOI must have exactly 6 elements [xMin, xMax, yMin, yMax, zMin, zMax].", nameof(voi));
        }

        throw new NotImplementedException("ExtractSubset requires VTK vtkExtractVOI.");
    }

    /// <summary>
    /// Crops the image to the specified extent.
    /// </summary>
    /// <param name="self">The image data to crop.</param>
    /// <param name="bounds">
    /// Cropping bounds as a 6-element array <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </param>
    /// <returns>The cropped <see cref="ImageData"/>.</returns>
    public static ImageData Crop(this ImageData self, double[] bounds)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(bounds);
        if (bounds.Length != 6)
        {
            throw new ArgumentException("Bounds must have exactly 6 elements.", nameof(bounds));
        }

        throw new NotImplementedException("Crop requires VTK vtkExtractVOI.");
    }

    // ---------------------------------------------------------------
    //  Morphological operations
    // ---------------------------------------------------------------

    /// <summary>
    /// Performs binary dilation and/or erosion on the image data.
    /// </summary>
    /// <param name="self">The image data to process.</param>
    /// <param name="dilateValue">Value to use for dilation.</param>
    /// <param name="erodeValue">Value to use for erosion.</param>
    /// <param name="kernelSize">
    /// Kernel size as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A processed <see cref="ImageData"/>.</returns>
    public static ImageData ImageDilateErode(this ImageData self, double dilateValue = 1.0, double erodeValue = 0.0, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ImageDilateErode requires VTK vtkImageDilateErode3D.");
    }

    /// <summary>
    /// Performs binary dilation on the image data.
    /// </summary>
    /// <param name="self">The image data to dilate.</param>
    /// <param name="kernelSize">
    /// Kernel size as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A dilated <see cref="ImageData"/>.</returns>
    public static ImageData Dilate(this ImageData self, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Dilate requires VTK vtkImageDilateErode3D.");
    }

    /// <summary>
    /// Performs binary erosion on the image data.
    /// </summary>
    /// <param name="self">The image data to erode.</param>
    /// <param name="kernelSize">
    /// Kernel size as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>An eroded <see cref="ImageData"/>.</returns>
    public static ImageData Erode(this ImageData self, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Erode requires VTK vtkImageDilateErode3D.");
    }

    /// <summary>
    /// Performs a morphological open operation (erosion followed by dilation).
    /// </summary>
    /// <param name="self">The image data to process.</param>
    /// <param name="kernelSize">
    /// Kernel size as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A morphologically opened <see cref="ImageData"/>.</returns>
    public static ImageData MorphologicalOpen(this ImageData self, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("MorphologicalOpen requires VTK vtkImageDilateErode3D (chained).");
    }

    /// <summary>
    /// Performs a morphological close operation (dilation followed by erosion).
    /// </summary>
    /// <param name="self">The image data to process.</param>
    /// <param name="kernelSize">
    /// Kernel size as a 3-element array <c>[kx, ky, kz]</c>.
    /// Defaults to <c>[3, 3, 3]</c>.
    /// </param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A morphologically closed <see cref="ImageData"/>.</returns>
    public static ImageData MorphologicalClose(this ImageData self, int[]? kernelSize = null, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("MorphologicalClose requires VTK vtkImageDilateErode3D (chained).");
    }

    // ---------------------------------------------------------------
    //  Frequency-domain filters
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the Fast Fourier Transform (FFT) of the image data.
    /// </summary>
    /// <param name="self">The image data to transform.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>An <see cref="ImageData"/> with the FFT result (complex values).</returns>
    public static ImageData Fft(this ImageData self, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Fft requires VTK vtkImageFFT.");
    }

    /// <summary>
    /// Computes the reverse (inverse) Fast Fourier Transform of the image data.
    /// </summary>
    /// <param name="self">The FFT image data to invert.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>An <see cref="ImageData"/> with the inverse FFT result.</returns>
    public static ImageData Rfft(this ImageData self, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Rfft requires VTK vtkImageRFFT.");
    }

    /// <summary>
    /// Applies a low-pass filter in the frequency domain.
    /// </summary>
    /// <param name="self">The image data (in frequency domain) to filter.</param>
    /// <param name="xCutoff">Normalized cutoff frequency for the x-axis.</param>
    /// <param name="yCutoff">Normalized cutoff frequency for the y-axis.</param>
    /// <param name="zCutoff">Normalized cutoff frequency for the z-axis.</param>
    /// <param name="order">Order of the Butterworth filter.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A low-pass filtered <see cref="ImageData"/>.</returns>
    public static ImageData LowPass(this ImageData self, double xCutoff = 0.5, double yCutoff = 0.5, double zCutoff = 0.5, int order = 1, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("LowPass requires VTK vtkImageButterworthLowPass.");
    }

    /// <summary>
    /// Applies a high-pass filter in the frequency domain.
    /// </summary>
    /// <param name="self">The image data (in frequency domain) to filter.</param>
    /// <param name="xCutoff">Normalized cutoff frequency for the x-axis.</param>
    /// <param name="yCutoff">Normalized cutoff frequency for the y-axis.</param>
    /// <param name="zCutoff">Normalized cutoff frequency for the z-axis.</param>
    /// <param name="order">Order of the Butterworth filter.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A high-pass filtered <see cref="ImageData"/>.</returns>
    public static ImageData HighPass(this ImageData self, double xCutoff = 0.5, double yCutoff = 0.5, double zCutoff = 0.5, int order = 1, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("HighPass requires VTK vtkImageButterworthHighPass.");
    }

    // ---------------------------------------------------------------
    //  Resampling
    // ---------------------------------------------------------------

    /// <summary>
    /// Resamples the image data to a new set of dimensions.
    /// </summary>
    /// <param name="self">The image data to resample.</param>
    /// <param name="dimensions">
    /// Target dimensions as a 3-element array <c>[nx, ny, nz]</c>.
    /// </param>
    /// <returns>A resampled <see cref="ImageData"/>.</returns>
    public static ImageData Resample(this ImageData self, int[] dimensions)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(dimensions);
        if (dimensions.Length != 3)
        {
            throw new ArgumentException("Dimensions must have exactly 3 elements.", nameof(dimensions));
        }

        throw new NotImplementedException("Resample requires VTK vtkImageResample or vtkResampleToImage.");
    }

    // ---------------------------------------------------------------
    //  Data association conversion
    // ---------------------------------------------------------------

    /// <summary>
    /// Converts point-associated data to cell-associated data.
    /// </summary>
    /// <param name="self">The image data to convert.</param>
    /// <returns>An <see cref="ImageData"/> with cell-associated data.</returns>
    public static ImageData PointsToCells(this ImageData self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("PointsToCells requires VTK vtkPointDataToCellData.");
    }

    /// <summary>
    /// Converts cell-associated data to point-associated data.
    /// </summary>
    /// <param name="self">The image data to convert.</param>
    /// <returns>An <see cref="ImageData"/> with point-associated data.</returns>
    public static ImageData CellsToPoints(this ImageData self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellsToPoints requires VTK vtkCellDataToPointData.");
    }

    // ---------------------------------------------------------------
    //  Padding
    // ---------------------------------------------------------------

    /// <summary>
    /// Pads the image data with a constant value.
    /// </summary>
    /// <param name="self">The image data to pad.</param>
    /// <param name="padSize">
    /// Number of voxels to pad on each side, as a 3-element array <c>[px, py, pz]</c>
    /// or a single value for uniform padding.
    /// </param>
    /// <param name="padValue">The constant value to fill padded voxels with.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>A padded <see cref="ImageData"/>.</returns>
    public static ImageData PadImage(this ImageData self, int[] padSize, double padValue = 0.0, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(padSize);
        throw new NotImplementedException("PadImage requires VTK vtkImageConstantPad.");
    }

    // ---------------------------------------------------------------
    //  Connectivity labeling
    // ---------------------------------------------------------------

    /// <summary>
    /// Labels connected regions in a binary image.
    /// </summary>
    /// <param name="self">The binary image data.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>An <see cref="ImageData"/> with a label array identifying connected regions.</returns>
    public static ImageData LabelConnectivity(this ImageData self, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("LabelConnectivity requires VTK vtkImageConnectivityFilter.");
    }

    // ---------------------------------------------------------------
    //  Value selection
    // ---------------------------------------------------------------

    /// <summary>
    /// Selects specific scalar values and masks all others.
    /// </summary>
    /// <param name="self">The image data to filter.</param>
    /// <param name="values">
    /// Array of scalar values to keep. All other values are replaced with
    /// <paramref name="maskValue"/>.
    /// </param>
    /// <param name="maskValue">Replacement value for non-selected voxels.</param>
    /// <param name="scalars">
    /// Name of the scalar array. When <c>null</c>, active scalars are used.
    /// </param>
    /// <returns>An <see cref="ImageData"/> with only the selected values.</returns>
    public static ImageData SelectValues(this ImageData self, double[] values, double maskValue = 0.0, string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(values);
        throw new NotImplementedException("SelectValues requires custom implementation.");
    }

    // ---------------------------------------------------------------
    //  Concatenation
    // ---------------------------------------------------------------

    /// <summary>
    /// Concatenates multiple <see cref="ImageData"/> objects along a specified axis.
    /// </summary>
    /// <param name="self">The first image data.</param>
    /// <param name="other">The image data to concatenate.</param>
    /// <param name="axis">Axis along which to concatenate: 0 (x), 1 (y), or 2 (z).</param>
    /// <returns>A concatenated <see cref="ImageData"/>.</returns>
    public static ImageData Concatenate(this ImageData self, ImageData other, int axis = 0)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        if (axis < 0 || axis > 2)
        {
            throw new ArgumentOutOfRangeException(nameof(axis), "Axis must be 0 (x), 1 (y), or 2 (z).");
        }

        throw new NotImplementedException("Concatenate requires VTK vtkImageAppend.");
    }
}
