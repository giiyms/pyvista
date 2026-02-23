using PyVista.Core;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.DataObjectFilters</c> mixin.
/// <para>
/// These filters operate on any <see cref="DataObject"/> and provide geometric
/// transformations, clipping, slicing, and other common operations.
/// </para>
/// </summary>
public static class DataObjectFilters
{
    // ---------------------------------------------------------------
    //  Transform / Rotate / Translate / Scale
    // ---------------------------------------------------------------

    /// <summary>
    /// Transforms the mesh by a 4×4 transformation matrix.
    /// </summary>
    /// <param name="self">The data object to transform.</param>
    /// <param name="matrix">
    /// A 4×4 transformation matrix stored as a 16-element array in row-major order.
    /// </param>
    /// <param name="transformAllInputVectors">
    /// When <c>true</c>, all input vectors are transformed. Otherwise only points are transformed.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The transformed data object.</returns>
    public static T Transform<T>(this T self, double[] matrix, bool transformAllInputVectors = false, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(matrix);
        if (matrix.Length != 16)
        {
            throw new ArgumentException("Transformation matrix must have exactly 16 elements (4x4).", nameof(matrix));
        }

        // Stub: requires VTK vtkTransformFilter.
        throw new NotImplementedException("Transform requires VTK bindings.");
    }

    /// <summary>
    /// Reflects the dataset across a plane defined by a normal and a point.
    /// </summary>
    /// <param name="self">The data object to reflect.</param>
    /// <param name="normal">The normal of the plane to reflect across, e.g. <c>"x"</c>, <c>"y"</c>, or <c>"z"</c>.</param>
    /// <param name="point">A point on the reflection plane. Defaults to the origin.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <param name="transformAllInputVectors">When <c>true</c>, transforms all input vectors.</param>
    /// <returns>The reflected data object.</returns>
    public static T Reflect<T>(this T self, string normal, double[]? point = null, bool inplace = false, bool transformAllInputVectors = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(normal);

        throw new NotImplementedException("Reflect requires VTK bindings.");
    }

    /// <summary>
    /// Rotates the mesh about the x-axis.
    /// </summary>
    /// <param name="self">The data object to rotate.</param>
    /// <param name="angle">Angle in degrees to rotate about the x-axis.</param>
    /// <param name="point">Point to rotate about. Defaults to the origin.</param>
    /// <param name="transformAllInputVectors">When <c>true</c>, transforms all input vectors.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The rotated data object.</returns>
    public static T RotateX<T>(this T self, double angle, double[]? point = null, bool transformAllInputVectors = false, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("RotateX requires VTK bindings.");
    }

    /// <summary>
    /// Rotates the mesh about the y-axis.
    /// </summary>
    /// <param name="self">The data object to rotate.</param>
    /// <param name="angle">Angle in degrees to rotate about the y-axis.</param>
    /// <param name="point">Point to rotate about. Defaults to the origin.</param>
    /// <param name="transformAllInputVectors">When <c>true</c>, transforms all input vectors.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The rotated data object.</returns>
    public static T RotateY<T>(this T self, double angle, double[]? point = null, bool transformAllInputVectors = false, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("RotateY requires VTK bindings.");
    }

    /// <summary>
    /// Rotates the mesh about the z-axis.
    /// </summary>
    /// <param name="self">The data object to rotate.</param>
    /// <param name="angle">Angle in degrees to rotate about the z-axis.</param>
    /// <param name="point">Point to rotate about. Defaults to the origin.</param>
    /// <param name="transformAllInputVectors">When <c>true</c>, transforms all input vectors.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The rotated data object.</returns>
    public static T RotateZ<T>(this T self, double angle, double[]? point = null, bool transformAllInputVectors = false, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("RotateZ requires VTK bindings.");
    }

    /// <summary>
    /// Rotates the mesh about an arbitrary vector.
    /// </summary>
    /// <param name="self">The data object to rotate.</param>
    /// <param name="vector">The rotation axis as a 3-element array.</param>
    /// <param name="angle">Angle in degrees.</param>
    /// <param name="point">Point to rotate about. Defaults to the origin.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <param name="transformAllInputVectors">When <c>true</c>, transforms all input vectors.</param>
    /// <returns>The rotated data object.</returns>
    public static T RotateVector<T>(this T self, double[] vector, double angle, double[]? point = null, bool inplace = false, bool transformAllInputVectors = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(vector);
        if (vector.Length != 3)
        {
            throw new ArgumentException("Rotation vector must have exactly 3 elements.", nameof(vector));
        }

        throw new NotImplementedException("RotateVector requires VTK bindings.");
    }

    /// <summary>
    /// Translates the mesh by the given offset vector.
    /// </summary>
    /// <param name="self">The data object to translate.</param>
    /// <param name="xyz">Translation offset as a 3-element array <c>[dx, dy, dz]</c>.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The translated data object.</returns>
    public static T Translate<T>(this T self, double[] xyz, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(xyz);
        if (xyz.Length != 3)
        {
            throw new ArgumentException("Translation vector must have exactly 3 elements.", nameof(xyz));
        }

        throw new NotImplementedException("Translate requires VTK bindings.");
    }

    /// <summary>
    /// Scales the mesh by the given factors along each axis.
    /// </summary>
    /// <param name="self">The data object to scale.</param>
    /// <param name="xyzFactor">
    /// Scale factors. Can be a single-element array for uniform scaling or a
    /// 3-element array <c>[sx, sy, sz]</c>.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The scaled data object.</returns>
    public static T Scale<T>(this T self, double[] xyzFactor, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(xyzFactor);

        throw new NotImplementedException("Scale requires VTK bindings.");
    }

    /// <summary>
    /// Flips the mesh across the YZ-plane (reflects in the x direction).
    /// </summary>
    /// <param name="self">The data object to flip.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The flipped data object.</returns>
    public static T FlipX<T>(this T self, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("FlipX requires VTK bindings.");
    }

    /// <summary>
    /// Flips the mesh across the XZ-plane (reflects in the y direction).
    /// </summary>
    /// <param name="self">The data object to flip.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The flipped data object.</returns>
    public static T FlipY<T>(this T self, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("FlipY requires VTK bindings.");
    }

    /// <summary>
    /// Flips the mesh across the XY-plane (reflects in the z direction).
    /// </summary>
    /// <param name="self">The data object to flip.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The flipped data object.</returns>
    public static T FlipZ<T>(this T self, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("FlipZ requires VTK bindings.");
    }

    // ---------------------------------------------------------------
    //  Clipping and Slicing
    // ---------------------------------------------------------------

    /// <summary>
    /// Clips the dataset by a plane defined by a normal and origin.
    /// </summary>
    /// <param name="self">The data object to clip.</param>
    /// <param name="normal">
    /// Plane normal direction. Can be a 3-element array or a string
    /// (<c>"x"</c>, <c>"y"</c>, <c>"z"</c>, or their negations).
    /// </param>
    /// <param name="origin">
    /// Point on the clipping plane. When <c>null</c>, the center of the dataset is used.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <param name="value">
    /// Scalar value to use when clipping with an implicit function.
    /// </param>
    /// <returns>The clipped data object.</returns>
    public static DataObject Clip(this DataObject self, double[] normal, double[]? origin = null, bool inplace = false, double value = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(normal);
        if (normal.Length != 3)
        {
            throw new ArgumentException("Normal must have exactly 3 elements.", nameof(normal));
        }

        throw new NotImplementedException("Clip requires VTK vtkClipDataSet.");
    }

    /// <summary>
    /// Clips the dataset by a bounding box.
    /// </summary>
    /// <param name="self">The data object to clip.</param>
    /// <param name="bounds">
    /// Bounding box as a 6-element array <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// When <c>null</c>, the dataset bounds are used.
    /// </param>
    /// <param name="invert">When <c>true</c>, clips the region outside the box.</param>
    /// <returns>The clipped data object.</returns>
    public static DataObject ClipBox(this DataObject self, double[]? bounds = null, bool invert = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ClipBox requires VTK vtkClipDataSet.");
    }

    /// <summary>
    /// Slices the dataset with a plane defined by a normal and an origin.
    /// </summary>
    /// <param name="self">The data object to slice.</param>
    /// <param name="normal">Plane normal as a 3-element array. Defaults to <c>[1, 0, 0]</c>.</param>
    /// <param name="origin">Point on the slice plane. When <c>null</c>, the center is used.</param>
    /// <param name="generateTriangles">When <c>true</c>, the output contains triangles.</param>
    /// <param name="contour">When <c>true</c>, contours are generated instead of slices.</param>
    /// <returns>A <see cref="PolyData"/> representing the slice.</returns>
    public static PolyData Slice(this DataObject self, double[]? normal = null, double[]? origin = null, bool generateTriangles = false, bool contour = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Slice requires VTK vtkCutter.");
    }

    /// <summary>
    /// Creates three orthogonal slices through the dataset.
    /// </summary>
    /// <param name="self">The data object to slice.</param>
    /// <param name="x">X-coordinate of the YZ slice plane. When <c>null</c>, the center is used.</param>
    /// <param name="y">Y-coordinate of the XZ slice plane. When <c>null</c>, the center is used.</param>
    /// <param name="z">Z-coordinate of the XY slice plane. When <c>null</c>, the center is used.</param>
    /// <returns>A <see cref="MultiBlock"/> containing the three slices.</returns>
    public static MultiBlock SliceOrthogonal(this DataObject self, double? x = null, double? y = null, double? z = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("SliceOrthogonal requires VTK vtkCutter.");
    }

    /// <summary>
    /// Creates multiple slices of the dataset along a specified axis.
    /// </summary>
    /// <param name="self">The data object to slice.</param>
    /// <param name="n">Number of evenly-spaced slices.</param>
    /// <param name="axis">
    /// Axis along which to create slices: 0 for x, 1 for y, 2 for z.
    /// Also accepts <c>"x"</c>, <c>"y"</c>, or <c>"z"</c> via overload.
    /// </param>
    /// <param name="generateTriangles">When <c>true</c>, the output contains triangles.</param>
    /// <param name="contour">When <c>true</c>, contours are generated instead of slices.</param>
    /// <returns>A <see cref="MultiBlock"/> of <see cref="PolyData"/> slices.</returns>
    public static MultiBlock SliceAlongAxis(this DataObject self, int n = 5, int axis = 0, bool generateTriangles = false, bool contour = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (axis < 0 || axis > 2)
        {
            throw new ArgumentOutOfRangeException(nameof(axis), "Axis must be 0 (x), 1 (y), or 2 (z).");
        }

        throw new NotImplementedException("SliceAlongAxis requires VTK vtkCutter.");
    }

    /// <summary>
    /// Slices a dataset along a line defined by two endpoints.
    /// </summary>
    /// <param name="self">The data object to slice.</param>
    /// <param name="line">
    /// A <see cref="PolyData"/> representing the line to slice along.
    /// </param>
    /// <param name="generateTriangles">When <c>true</c>, the output contains triangles.</param>
    /// <param name="contour">When <c>true</c>, contours are generated instead of slices.</param>
    /// <returns>A <see cref="PolyData"/> representing the slice.</returns>
    public static PolyData SliceAlongLine(this DataObject self, PolyData line, bool generateTriangles = false, bool contour = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(line);
        throw new NotImplementedException("SliceAlongLine requires VTK vtkCutter.");
    }

    // ---------------------------------------------------------------
    //  Data Conversion
    // ---------------------------------------------------------------

    /// <summary>
    /// Averages cell data to point data.
    /// </summary>
    /// <param name="self">The data object to convert.</param>
    /// <param name="passThrough">When <c>true</c>, the original cell data is preserved.</param>
    /// <returns>A new data object with point data derived from cell data.</returns>
    public static DataObject CellDataToPointData(this DataObject self, bool passThrough = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellDataToPointData requires VTK vtkCellDataToPointData.");
    }

    /// <summary>
    /// Averages point data to cell data.
    /// </summary>
    /// <param name="self">The data object to convert.</param>
    /// <param name="passThrough">When <c>true</c>, the original point data is preserved.</param>
    /// <returns>A new data object with cell data derived from point data.</returns>
    public static DataObject PointDataToCellData(this DataObject self, bool passThrough = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("PointDataToCellData requires VTK vtkPointDataToCellData.");
    }

    // ---------------------------------------------------------------
    //  Geometry / Topology
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts the surface geometry from the dataset.
    /// </summary>
    /// <param name="self">The data object to extract the surface from.</param>
    /// <param name="passPointIds">When <c>true</c>, passes original point IDs to the output.</param>
    /// <param name="passCellIds">When <c>true</c>, passes original cell IDs to the output.</param>
    /// <returns>A <see cref="PolyData"/> representing the external surface.</returns>
    public static PolyData ExtractSurface(this DataObject self, bool passPointIds = true, bool passCellIds = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractSurface requires VTK vtkDataSetSurfaceFilter.");
    }

    /// <summary>
    /// Extracts all edges of the dataset as a <see cref="PolyData"/> of lines.
    /// </summary>
    /// <param name="self">The data object whose edges are extracted.</param>
    /// <returns>A <see cref="PolyData"/> containing the edges as line cells.</returns>
    public static PolyData ExtractAllEdges(this DataObject self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractAllEdges requires VTK vtkExtractEdges.");
    }

    /// <summary>
    /// Triangulates all cells in the dataset.
    /// </summary>
    /// <param name="self">The data object to triangulate.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The triangulated data object.</returns>
    public static DataObject Triangulate(this DataObject self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Triangulate requires VTK vtkDataSetTriangleFilter.");
    }

    // ---------------------------------------------------------------
    //  Elevation / Cell sizes
    // ---------------------------------------------------------------

    /// <summary>
    /// Generates scalar values based on elevation (height along a direction).
    /// </summary>
    /// <param name="self">The data object to compute elevation for.</param>
    /// <param name="lowPoint">The low point of the elevation range, as a 3-element array.</param>
    /// <param name="highPoint">The high point of the elevation range, as a 3-element array.</param>
    /// <param name="scalarRange">
    /// Scalar range to map the elevation to. Defaults to <c>[0, 1]</c>.
    /// </param>
    /// <returns>A data object with an <c>"Elevation"</c> point scalar array.</returns>
    public static DataObject Elevation(this DataObject self, double[]? lowPoint = null, double[]? highPoint = null, double[]? scalarRange = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Elevation requires VTK vtkElevationFilter.");
    }

    /// <summary>
    /// Computes the size of each cell (length for 1-D, area for 2-D, volume for 3-D).
    /// </summary>
    /// <param name="self">The data object to compute cell sizes for.</param>
    /// <param name="length">When <c>true</c>, computes the length of 1-D cells.</param>
    /// <param name="area">When <c>true</c>, computes the area of 2-D cells.</param>
    /// <param name="volume">When <c>true</c>, computes the volume of 3-D cells.</param>
    /// <returns>A data object with cell size arrays.</returns>
    public static DataObject ComputeCellSizes(this DataObject self, bool length = true, bool area = true, bool volume = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ComputeCellSizes requires VTK vtkCellSizeFilter.");
    }

    /// <summary>
    /// Returns a <see cref="PolyData"/> of points at the centers of each cell.
    /// </summary>
    /// <param name="self">The data object whose cell centers are computed.</param>
    /// <returns>A <see cref="PolyData"/> with one vertex per cell center.</returns>
    public static PolyData CellCenters(this DataObject self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellCenters requires VTK vtkCellCenters.");
    }

    // ---------------------------------------------------------------
    //  Sampling / Quality
    // ---------------------------------------------------------------

    /// <summary>
    /// Samples the values of a dataset at the point locations of another dataset.
    /// </summary>
    /// <param name="self">The data object to resample onto.</param>
    /// <param name="source">The source dataset whose values will be sampled.</param>
    /// <param name="tolerance">Tolerance for the sampling operation.</param>
    /// <param name="passThrough">When <c>true</c>, passes the original data arrays through.</param>
    /// <returns>A new data object with sampled values.</returns>
    public static DataObject Sample(this DataObject self, DataObject source, double? tolerance = null, bool passThrough = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(source);
        throw new NotImplementedException("Sample requires VTK vtkResampleWithDataSet.");
    }

    /// <summary>
    /// Computes per-cell quality metrics for the dataset.
    /// </summary>
    /// <param name="self">The data object to evaluate.</param>
    /// <param name="qualityMeasure">
    /// The quality measure to compute, e.g., <c>"area"</c>, <c>"aspect_ratio"</c>,
    /// <c>"condition"</c>, <c>"jacobian"</c>, <c>"scaled_jacobian"</c>.
    /// When <c>null</c>, a default measure is used for each cell type.
    /// </param>
    /// <returns>A data object with a <c>"CellQuality"</c> cell data array.</returns>
    public static DataObject CellQuality(this DataObject self, string? qualityMeasure = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellQuality requires VTK vtkMeshQuality.");
    }

    /// <summary>
    /// Resizes the mesh to fit within the specified bounds.
    /// </summary>
    /// <param name="self">The data object to resize.</param>
    /// <param name="newSize">
    /// Target bounding size as a 3-element array <c>[xSize, ySize, zSize]</c>.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The resized data object.</returns>
    public static T Resize<T>(this T self, double[] newSize, bool inplace = false)
        where T : DataObject
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(newSize);
        if (newSize.Length != 3)
        {
            throw new ArgumentException("New size must have exactly 3 elements.", nameof(newSize));
        }

        throw new NotImplementedException("Resize requires VTK bindings.");
    }
}
