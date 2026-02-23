using PyVista.Core;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.DataSetFilters</c> mixin.
/// <para>
/// These filters operate on any <see cref="DataSet"/> and provide the broadest set of
/// VTK-based algorithms including clipping, contouring, glyphing, warping, streamlines,
/// and many more.
/// </para>
/// </summary>
public static class DataSetFilters
{
    // ---------------------------------------------------------------
    //  Clipping
    // ---------------------------------------------------------------

    /// <summary>
    /// Clips a <see cref="DataSet"/> by a scalar value.
    /// </summary>
    /// <param name="self">The dataset to clip.</param>
    /// <param name="scalars">
    /// Name of the scalar array to clip with. When <c>null</c>, the active scalars are used.
    /// </param>
    /// <param name="value">
    /// Scalar value to clip at. When <c>null</c>, the midpoint of the scalar range is used.
    /// </param>
    /// <param name="invert">When <c>true</c>, clips the region above the value.</param>
    /// <returns>The clipped dataset.</returns>
    public static DataSet ClipScalar(this DataSet self, string? scalars = null, double? value = null, bool invert = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ClipScalar requires VTK vtkClipDataSet.");
    }

    /// <summary>
    /// Clips a <see cref="DataSet"/> using an implicit surface.
    /// </summary>
    /// <param name="self">The dataset to clip.</param>
    /// <param name="surface">
    /// A <see cref="PolyData"/> surface to clip against.
    /// </param>
    /// <param name="invert">When <c>true</c>, clips the region outside the surface.</param>
    /// <param name="value">Scalar value for the implicit function clip level.</param>
    /// <returns>The clipped dataset.</returns>
    public static DataSet ClipSurface(this DataSet self, PolyData surface, bool invert = true, double value = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(surface);
        throw new NotImplementedException("ClipSurface requires VTK vtkClipDataSet.");
    }

    // ---------------------------------------------------------------
    //  Threshold
    // ---------------------------------------------------------------

    /// <summary>
    /// Thresholds the dataset by a scalar value or range.
    /// </summary>
    /// <param name="self">The dataset to threshold.</param>
    /// <param name="value">
    /// Single value or two-element array <c>[lower, upper]</c> defining the threshold range.
    /// When <c>null</c>, the non-NaN data range is used.
    /// </param>
    /// <param name="scalars">Name of the scalar array. When <c>null</c>, active scalars are used.</param>
    /// <param name="invert">When <c>true</c>, inverts the threshold selection.</param>
    /// <param name="method">
    /// Threshold method: <c>"upper"</c>, <c>"lower"</c>, or <c>"between"</c>.
    /// </param>
    /// <returns>The thresholded dataset.</returns>
    public static DataSet Threshold(this DataSet self, double[]? value = null, string? scalars = null, bool invert = false, string method = "upper")
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Threshold requires VTK vtkThreshold.");
    }

    /// <summary>
    /// Thresholds the dataset by a percentage of the scalar range.
    /// </summary>
    /// <param name="self">The dataset to threshold.</param>
    /// <param name="percent">
    /// Percentage in the range <c>[0, 1]</c>. Values below this percentile are removed.
    /// </param>
    /// <param name="scalars">Name of the scalar array. When <c>null</c>, active scalars are used.</param>
    /// <param name="invert">When <c>true</c>, inverts the threshold selection.</param>
    /// <param name="method">
    /// Threshold method: <c>"upper"</c>, <c>"lower"</c>, or <c>"between"</c>.
    /// </param>
    /// <returns>The thresholded dataset.</returns>
    public static DataSet ThresholdPercent(this DataSet self, double percent = 0.5, string? scalars = null, bool invert = false, string method = "upper")
    {
        ArgumentNullException.ThrowIfNull(self);
        if (percent < 0.0 || percent > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(percent), "Percent must be between 0 and 1.");
        }

        throw new NotImplementedException("ThresholdPercent requires VTK vtkThreshold.");
    }

    // ---------------------------------------------------------------
    //  Outline
    // ---------------------------------------------------------------

    /// <summary>
    /// Produces an outline (bounding box wireframe) of the dataset.
    /// </summary>
    /// <param name="self">The dataset to outline.</param>
    /// <param name="generateFaces">
    /// When <c>true</c>, generates polygonal faces instead of wireframe edges.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the outline.</returns>
    public static PolyData Outline(this DataSet self, bool generateFaces = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Outline requires VTK vtkOutlineFilter.");
    }

    /// <summary>
    /// Produces an outline showing only the corners of the bounding box.
    /// </summary>
    /// <param name="self">The dataset to outline.</param>
    /// <param name="factor">
    /// Controls the relative size of the corners. Defaults to <c>0.2</c>.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the corner outline.</returns>
    public static PolyData OutlineCorners(this DataSet self, double factor = 0.2)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("OutlineCorners requires VTK vtkOutlineCornerFilter.");
    }

    // ---------------------------------------------------------------
    //  Contouring
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes isolines or isosurfaces from the dataset using a scalar array.
    /// </summary>
    /// <param name="self">The dataset to contour.</param>
    /// <param name="isosurfaces">
    /// Number of evenly-spaced isosurfaces to generate, or a list of explicit scalar values.
    /// </param>
    /// <param name="scalars">Name of the scalar array. When <c>null</c>, active scalars are used.</param>
    /// <param name="computeNormals">When <c>true</c>, computes normals for the output surface.</param>
    /// <param name="computeGradients">When <c>true</c>, computes gradients for the output surface.</param>
    /// <param name="computeScalars">When <c>true</c>, computes scalars for the output surface.</param>
    /// <param name="method">Contouring method: <c>"contour"</c> or <c>"marching_cubes"</c>.</param>
    /// <returns>A <see cref="PolyData"/> of the isosurface(s).</returns>
    public static PolyData Contour(this DataSet self, int isosurfaces = 1, string? scalars = null, bool computeNormals = false, bool computeGradients = false, bool computeScalars = true, string method = "contour")
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Contour requires VTK vtkContourFilter.");
    }

    /// <summary>
    /// Computes isolines or isosurfaces at explicit scalar values.
    /// </summary>
    /// <param name="self">The dataset to contour.</param>
    /// <param name="values">An array of explicit iso-values.</param>
    /// <param name="scalars">Name of the scalar array. When <c>null</c>, active scalars are used.</param>
    /// <param name="computeNormals">When <c>true</c>, computes normals for the output surface.</param>
    /// <param name="computeGradients">When <c>true</c>, computes gradients for the output surface.</param>
    /// <param name="computeScalars">When <c>true</c>, computes scalars for the output surface.</param>
    /// <param name="method">Contouring method: <c>"contour"</c> or <c>"marching_cubes"</c>.</param>
    /// <returns>A <see cref="PolyData"/> of the isosurface(s).</returns>
    public static PolyData ContourValues(this DataSet self, double[] values, string? scalars = null, bool computeNormals = false, bool computeGradients = false, bool computeScalars = true, string method = "contour")
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(values);
        throw new NotImplementedException("ContourValues requires VTK vtkContourFilter.");
    }

    // ---------------------------------------------------------------
    //  Geometry Extraction
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts the geometry (points, cells) from the dataset as an
    /// <see cref="UnstructuredGrid"/> or <see cref="PolyData"/>.
    /// </summary>
    /// <param name="self">The dataset to extract geometry from.</param>
    /// <returns>The extracted geometry.</returns>
    public static DataSet ExtractGeometry(this DataSet self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractGeometry requires VTK vtkGeometryFilter.");
    }

    /// <summary>
    /// Extracts the largest connected region from the dataset.
    /// </summary>
    /// <param name="self">The dataset to process.</param>
    /// <returns>The largest connected region.</returns>
    public static DataSet ExtractLargest(this DataSet self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractLargest requires VTK vtkConnectivityFilter.");
    }

    /// <summary>
    /// Splits the dataset into separate bodies (connected regions).
    /// </summary>
    /// <param name="self">The dataset to split.</param>
    /// <returns>A <see cref="MultiBlock"/> where each block is a connected region.</returns>
    public static MultiBlock SplitBodies(this DataSet self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("SplitBodies requires VTK vtkConnectivityFilter.");
    }

    /// <summary>
    /// Finds and labels connected regions in the dataset.
    /// </summary>
    /// <param name="self">The dataset to label.</param>
    /// <param name="extractionMode">
    /// Mode of extraction: <c>"all"</c>, <c>"largest"</c>, <c>"specified"</c>,
    /// <c>"closest"</c>, <c>"cell_seed"</c>, or <c>"point_seed"</c>.
    /// </param>
    /// <param name="scalars">Name of the scalar array for labeling.</param>
    /// <returns>The dataset with a <c>"RegionId"</c> cell data array.</returns>
    public static DataSet Connectivity(this DataSet self, string extractionMode = "all", string? scalars = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Connectivity requires VTK vtkConnectivityFilter.");
    }

    // ---------------------------------------------------------------
    //  Glyphs
    // ---------------------------------------------------------------

    /// <summary>
    /// Copies a geometric representation (glyph) to every point in the dataset.
    /// </summary>
    /// <param name="self">The dataset to glyph.</param>
    /// <param name="orient">Name of the vector array for glyph orientation.</param>
    /// <param name="scale">Name of the scalar array for glyph scaling.</param>
    /// <param name="factor">Global scale factor applied to each glyph.</param>
    /// <param name="geom">
    /// The glyph geometry as a <see cref="PolyData"/>. When <c>null</c>, an arrow is used.
    /// </param>
    /// <param name="tolerance">
    /// Tolerance for merging coincident points. When <c>null</c>, no merging.
    /// </param>
    /// <param name="absoluteScaling">
    /// When <c>true</c>, uses absolute scaling (ignores dataset dimensions).
    /// </param>
    /// <returns>A <see cref="PolyData"/> containing the glyphed geometry.</returns>
    public static PolyData Glyph(this DataSet self, string? orient = null, string? scale = null, double factor = 1.0, PolyData? geom = null, double? tolerance = null, bool absoluteScaling = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Glyph requires VTK vtkGlyph3D.");
    }

    // ---------------------------------------------------------------
    //  Warping
    // ---------------------------------------------------------------

    /// <summary>
    /// Warps the dataset geometry by a scalar value.
    /// </summary>
    /// <param name="self">The dataset to warp.</param>
    /// <param name="scalars">Name of the scalar array. When <c>null</c>, active scalars are used.</param>
    /// <param name="factor">Scale factor for the warping displacement.</param>
    /// <param name="normal">
    /// Direction of warping as a 3-element array. When <c>null</c>, the point normals are used.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The warped dataset.</returns>
    public static DataSet WarpByScalar(this DataSet self, string? scalars = null, double factor = 1.0, double[]? normal = null, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("WarpByScalar requires VTK vtkWarpScalar.");
    }

    /// <summary>
    /// Warps the dataset geometry by a vector field.
    /// </summary>
    /// <param name="self">The dataset to warp.</param>
    /// <param name="vectors">Name of the vector array. When <c>null</c>, active vectors are used.</param>
    /// <param name="factor">Scale factor for the displacement.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The warped dataset.</returns>
    public static DataSet WarpByVector(this DataSet self, string? vectors = null, double factor = 1.0, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("WarpByVector requires VTK vtkWarpVector.");
    }

    // ---------------------------------------------------------------
    //  Streamlines
    // ---------------------------------------------------------------

    /// <summary>
    /// Generates streamlines from a vector field.
    /// </summary>
    /// <param name="self">The dataset containing the vector field.</param>
    /// <param name="vectors">Name of the vector array. When <c>null</c>, active vectors are used.</param>
    /// <param name="sourceCenter">Center of the seeding sphere as a 3-element array.</param>
    /// <param name="sourceRadius">Radius of the seeding sphere.</param>
    /// <param name="nPoints">Number of seed points.</param>
    /// <param name="startPosition">
    /// Starting position as a 3-element array. When specified, a single streamline is generated.
    /// </param>
    /// <param name="maxTime">Maximum integration time.</param>
    /// <param name="maxSteps">Maximum number of integration steps.</param>
    /// <param name="terminalSpeed">Speed below which integration is terminated.</param>
    /// <param name="integrationDirection">
    /// Direction of integration: <c>"forward"</c>, <c>"backward"</c>, or <c>"both"</c>.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the streamlines.</returns>
    public static PolyData Streamlines(
        this DataSet self,
        string? vectors = null,
        double[]? sourceCenter = null,
        double sourceRadius = 1.0,
        int nPoints = 100,
        double[]? startPosition = null,
        double maxTime = 200.0,
        int maxSteps = 2000,
        double terminalSpeed = 1e-12,
        string integrationDirection = "both")
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Streamlines requires VTK vtkStreamTracer.");
    }

    /// <summary>
    /// Generates streamlines from an explicit seed source.
    /// </summary>
    /// <param name="self">The dataset containing the vector field.</param>
    /// <param name="source">A <see cref="DataSet"/> defining the seed points.</param>
    /// <param name="vectors">Name of the vector array. When <c>null</c>, active vectors are used.</param>
    /// <param name="maxTime">Maximum integration time.</param>
    /// <param name="maxSteps">Maximum number of integration steps.</param>
    /// <param name="terminalSpeed">Speed below which integration is terminated.</param>
    /// <param name="integrationDirection">
    /// Direction of integration: <c>"forward"</c>, <c>"backward"</c>, or <c>"both"</c>.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the streamlines.</returns>
    public static PolyData StreamlinesFromSource(
        this DataSet self,
        DataSet source,
        string? vectors = null,
        double maxTime = 200.0,
        int maxSteps = 2000,
        double terminalSpeed = 1e-12,
        string integrationDirection = "both")
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(source);
        throw new NotImplementedException("StreamlinesFromSource requires VTK vtkStreamTracer.");
    }

    // ---------------------------------------------------------------
    //  Sampling
    // ---------------------------------------------------------------

    /// <summary>
    /// Samples a dataset's values over a line between two points.
    /// </summary>
    /// <param name="self">The dataset to sample.</param>
    /// <param name="pointa">Start point of the sampling line as a 3-element array.</param>
    /// <param name="pointb">End point of the sampling line as a 3-element array.</param>
    /// <param name="resolution">Number of sampling points along the line.</param>
    /// <param name="tolerance">Tolerance for the probe filter.</param>
    /// <returns>A <see cref="PolyData"/> of sampled values.</returns>
    public static PolyData SampleOverLine(this DataSet self, double[] pointa, double[] pointb, int resolution = 1000, double? tolerance = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(pointa);
        ArgumentNullException.ThrowIfNull(pointb);
        if (pointa.Length != 3 || pointb.Length != 3)
        {
            throw new ArgumentException("Point arrays must have exactly 3 elements.");
        }

        throw new NotImplementedException("SampleOverLine requires VTK vtkProbeFilter.");
    }

    /// <summary>
    /// Samples a dataset's values over multiple line segments.
    /// </summary>
    /// <param name="self">The dataset to sample.</param>
    /// <param name="points">
    /// An array of 3-D points defining the polyline vertices. Must contain at least two points.
    /// </param>
    /// <param name="resolution">Number of sampling points per segment.</param>
    /// <param name="tolerance">Tolerance for the probe filter.</param>
    /// <returns>A <see cref="PolyData"/> of sampled values.</returns>
    public static PolyData SampleOverMultipleLines(this DataSet self, double[,] points, int resolution = 1000, double? tolerance = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(points);
        if (points.GetLength(1) != 3 || points.GetLength(0) < 2)
        {
            throw new ArgumentException("Points must be an (N, 3) array with at least 2 rows.", nameof(points));
        }

        throw new NotImplementedException("SampleOverMultipleLines requires VTK vtkProbeFilter.");
    }

    /// <summary>
    /// Probes the dataset at the point locations of another dataset.
    /// </summary>
    /// <param name="self">The source dataset to probe.</param>
    /// <param name="points">
    /// A <see cref="DataSet"/> whose point locations define where to sample.
    /// </param>
    /// <param name="tolerance">Tolerance for the probe operation.</param>
    /// <returns>A dataset with sampled values at the probe locations.</returns>
    public static DataSet ProbeFilter(this DataSet self, DataSet points, double? tolerance = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(points);
        throw new NotImplementedException("ProbeFilter requires VTK vtkProbeFilter.");
    }

    // ---------------------------------------------------------------
    //  Interpolation
    // ---------------------------------------------------------------

    /// <summary>
    /// Interpolates point data from a source dataset onto this dataset.
    /// </summary>
    /// <param name="self">The target dataset.</param>
    /// <param name="source">The source dataset whose scalars are interpolated.</param>
    /// <param name="sharpness">Sharpness parameter for the interpolation kernel.</param>
    /// <param name="radius">Radius of the interpolation kernel.</param>
    /// <param name="strategy">
    /// Strategy: <c>"null_value"</c>, <c>"mask_points"</c>, or <c>"closest_point"</c>.
    /// </param>
    /// <returns>A dataset with interpolated point data.</returns>
    public static DataSet Interpolate(this DataSet self, DataSet source, double sharpness = 2.0, double radius = 1.0, string strategy = "null_value")
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(source);
        throw new NotImplementedException("Interpolate requires VTK vtkPointInterpolator.");
    }

    // ---------------------------------------------------------------
    //  Decimation
    // ---------------------------------------------------------------

    /// <summary>
    /// Reduces the number of triangles in a triangulated surface.
    /// </summary>
    /// <param name="self">The dataset to decimate.</param>
    /// <param name="targetReduction">
    /// Fraction of triangles to remove, in the range <c>[0, 1)</c>.
    /// </param>
    /// <returns>The decimated dataset.</returns>
    public static DataSet Decimate(this DataSet self, double targetReduction = 0.5)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (targetReduction < 0.0 || targetReduction >= 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(targetReduction), "Target reduction must be in [0, 1).");
        }

        throw new NotImplementedException("Decimate requires VTK vtkDecimatePro.");
    }

    /// <summary>
    /// Decimates the boundary of a dataset.
    /// </summary>
    /// <param name="self">The dataset whose boundary is decimated.</param>
    /// <param name="targetReduction">
    /// Fraction of boundary triangles to remove, in the range <c>[0, 1)</c>.
    /// </param>
    /// <returns>The decimated dataset.</returns>
    public static DataSet DecimateBoundary(this DataSet self, double targetReduction = 0.5)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("DecimateBoundary requires VTK vtkDecimatePro.");
    }

    // ---------------------------------------------------------------
    //  Delaunay
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the 3-D Delaunay triangulation of the dataset.
    /// </summary>
    /// <param name="self">The dataset to triangulate.</param>
    /// <param name="alpha">Alpha value for the Delaunay triangulation.</param>
    /// <param name="tolerance">Tolerance for coincident point merging.</param>
    /// <param name="offset">Multiplier to control bounding triangulation size.</param>
    /// <returns>An <see cref="UnstructuredGrid"/> of tetrahedra.</returns>
    public static UnstructuredGrid Delaunay3D(this DataSet self, double alpha = 0.0, double tolerance = 0.001, double offset = 2.5)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Delaunay3D requires VTK vtkDelaunay3D.");
    }

    // ---------------------------------------------------------------
    //  Selection
    // ---------------------------------------------------------------

    /// <summary>
    /// Marks points as inside or outside an enclosed surface.
    /// </summary>
    /// <param name="self">The dataset containing the points to test.</param>
    /// <param name="surface">A closed <see cref="PolyData"/> surface.</param>
    /// <param name="tolerance">Tolerance for the enclosure test.</param>
    /// <param name="insideOut">When <c>true</c>, inverts the selection.</param>
    /// <returns>A dataset with a <c>"SelectedPoints"</c> point data array.</returns>
    public static DataSet SelectEnclosedPoints(this DataSet self, PolyData surface, double tolerance = 0.001, bool insideOut = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(surface);
        throw new NotImplementedException("SelectEnclosedPoints requires VTK vtkSelectEnclosedPoints.");
    }

    // ---------------------------------------------------------------
    //  Distance
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the implicit distance from the dataset to a surface.
    /// </summary>
    /// <param name="self">The dataset whose points are evaluated.</param>
    /// <param name="surface">A <see cref="PolyData"/> surface to compute distance to.</param>
    /// <returns>A dataset with an <c>"implicit_distance"</c> point data array.</returns>
    public static DataSet ComputeImplicitDistance(this DataSet self, PolyData surface)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(surface);
        throw new NotImplementedException("ComputeImplicitDistance requires VTK vtkImplicitPolyDataDistance.");
    }

    // ---------------------------------------------------------------
    //  Texture mapping
    // ---------------------------------------------------------------

    /// <summary>
    /// Maps texture coordinates to a plane.
    /// </summary>
    /// <param name="self">The dataset to map texture coordinates on.</param>
    /// <param name="origin">Origin of the texture plane as a 3-element array.</param>
    /// <param name="point1">
    /// First axis-defining point of the texture plane as a 3-element array.
    /// </param>
    /// <param name="point2">
    /// Second axis-defining point of the texture plane as a 3-element array.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>A dataset with <c>"Texture Coordinates"</c> point data.</returns>
    public static DataSet TextureMapToPlane(this DataSet self, double[]? origin = null, double[]? point1 = null, double[]? point2 = null, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("TextureMapToPlane requires VTK vtkTextureMapToPlane.");
    }

    /// <summary>
    /// Maps texture coordinates to a sphere.
    /// </summary>
    /// <param name="self">The dataset to map texture coordinates on.</param>
    /// <param name="center">
    /// Center of the mapping sphere. When <c>null</c>, the center of the dataset is used.
    /// </param>
    /// <param name="preventSeam">When <c>true</c>, prevents seam artifacts in the mapping.</param>
    /// <returns>A dataset with <c>"Texture Coordinates"</c> point data.</returns>
    public static DataSet TextureMapToSphere(this DataSet self, double[]? center = null, bool preventSeam = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("TextureMapToSphere requires VTK vtkTextureMapToSphere.");
    }

    // ---------------------------------------------------------------
    //  Alignment
    // ---------------------------------------------------------------

    /// <summary>
    /// Aligns the dataset to its principal axes via an iterative closest-point algorithm.
    /// </summary>
    /// <param name="self">The dataset to align.</param>
    /// <param name="target">The target dataset to align to.</param>
    /// <param name="maxIterations">Maximum number of ICP iterations.</param>
    /// <param name="tolerance">Convergence tolerance.</param>
    /// <returns>The aligned dataset.</returns>
    public static DataSet Align(this DataSet self, DataSet target, int maxIterations = 100, double tolerance = 1e-6)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(target);
        throw new NotImplementedException("Align requires VTK vtkIterativeClosestPointTransform.");
    }

    // ---------------------------------------------------------------
    //  Gaussian splatting
    // ---------------------------------------------------------------

    /// <summary>
    /// Performs Gaussian splatting to create an image representation of point data.
    /// </summary>
    /// <param name="self">The dataset to splat.</param>
    /// <param name="radius">Radius of the Gaussian kernel.</param>
    /// <param name="dimensions">
    /// Output image dimensions as a 3-element array. Defaults to <c>[100, 100, 100]</c>.
    /// </param>
    /// <returns>An <see cref="ImageData"/> of the splatted result.</returns>
    public static ImageData GaussianSplatting(this DataSet self, double radius = 1.0, int[]? dimensions = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("GaussianSplatting requires VTK vtkGaussianSplatter.");
    }

    // ---------------------------------------------------------------
    //  Elevation (DataSet-specific overload)
    // ---------------------------------------------------------------

    /// <summary>
    /// Generates scalar values based on elevation (height along a direction).
    /// </summary>
    /// <param name="self">The dataset to compute elevation for.</param>
    /// <param name="lowPoint">The low point of the elevation range, as a 3-element array.</param>
    /// <param name="highPoint">The high point of the elevation range, as a 3-element array.</param>
    /// <param name="scalarRange">
    /// Scalar range to map the elevation to. Defaults to <c>[0, 1]</c>.
    /// </param>
    /// <returns>A dataset with an <c>"Elevation"</c> point scalar array.</returns>
    public static DataSet Elevation(this DataSet self, double[]? lowPoint = null, double[]? highPoint = null, double[]? scalarRange = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Elevation requires VTK vtkElevationFilter.");
    }

    // ---------------------------------------------------------------
    //  Cell operations
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a <see cref="PolyData"/> of points at the centers of each cell.
    /// </summary>
    /// <param name="self">The dataset whose cell centers are computed.</param>
    /// <param name="vertexCells">
    /// When <c>true</c>, the output includes vertex cells (one per center point).
    /// </param>
    /// <returns>A <see cref="PolyData"/> with one vertex per cell center.</returns>
    public static PolyData CellCenters(this DataSet self, bool vertexCells = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellCenters requires VTK vtkCellCenters.");
    }

    /// <summary>
    /// Averages cell data to point data.
    /// </summary>
    /// <param name="self">The dataset to convert.</param>
    /// <param name="passThrough">When <c>true</c>, the original cell data is preserved.</param>
    /// <returns>A dataset with point data derived from cell data.</returns>
    public static DataSet CellDataToPointData(this DataSet self, bool passThrough = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellDataToPointData requires VTK vtkCellDataToPointData.");
    }

    /// <summary>
    /// Averages point data to cell data.
    /// </summary>
    /// <param name="self">The dataset to convert.</param>
    /// <param name="passThrough">When <c>true</c>, the original point data is preserved.</param>
    /// <returns>A dataset with cell data derived from point data.</returns>
    public static DataSet PointDataToCellData(this DataSet self, bool passThrough = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("PointDataToCellData requires VTK vtkPointDataToCellData.");
    }

    /// <summary>
    /// Computes per-cell quality metrics.
    /// </summary>
    /// <param name="self">The dataset to evaluate.</param>
    /// <param name="qualityMeasure">
    /// The quality measure to compute, e.g., <c>"area"</c>, <c>"aspect_ratio"</c>,
    /// <c>"condition"</c>, <c>"jacobian"</c>, <c>"scaled_jacobian"</c>.
    /// When <c>null</c>, a default measure is used for each cell type.
    /// </param>
    /// <returns>A dataset with a <c>"CellQuality"</c> cell data array.</returns>
    public static DataSet CellQuality(this DataSet self, string? qualityMeasure = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CellQuality requires VTK vtkMeshQuality.");
    }

    /// <summary>
    /// Computes the size of each cell (length, area, or volume).
    /// </summary>
    /// <param name="self">The dataset to compute cell sizes for.</param>
    /// <param name="length">When <c>true</c>, computes the length of 1-D cells.</param>
    /// <param name="area">When <c>true</c>, computes the area of 2-D cells.</param>
    /// <param name="volume">When <c>true</c>, computes the volume of 3-D cells.</param>
    /// <returns>A dataset with cell size arrays.</returns>
    public static DataSet ComputeCellSizes(this DataSet self, bool length = true, bool area = true, bool volume = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ComputeCellSizes requires VTK vtkCellSizeFilter.");
    }

    // ---------------------------------------------------------------
    //  Surface extraction
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts the external surface of the dataset.
    /// </summary>
    /// <param name="self">The dataset to extract the surface from.</param>
    /// <param name="passPointIds">When <c>true</c>, passes original point IDs.</param>
    /// <param name="passCellIds">When <c>true</c>, passes original cell IDs.</param>
    /// <returns>A <see cref="PolyData"/> representing the external surface.</returns>
    public static PolyData ExtractSurface(this DataSet self, bool passPointIds = true, bool passCellIds = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractSurface requires VTK vtkDataSetSurfaceFilter.");
    }

    /// <summary>
    /// Triangulates all cells in the dataset.
    /// </summary>
    /// <param name="self">The dataset to triangulate.</param>
    /// <param name="inplace">When <c>true</c>, modifies the dataset in place.</param>
    /// <returns>The triangulated dataset.</returns>
    public static DataSet Triangulate(this DataSet self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Triangulate requires VTK vtkDataSetTriangleFilter.");
    }
}
