using PyVista.Core;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.PolyDataFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="PolyData"/> and provide surface mesh operations
/// such as Boolean operations, smoothing, decimation, normal computation, and more.
/// </para>
/// </summary>
public static class PolyDataFilters
{
    // ---------------------------------------------------------------
    //  Boolean operations
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the Boolean union of two <see cref="PolyData"/> meshes.
    /// </summary>
    /// <param name="self">The first mesh.</param>
    /// <param name="other">The second mesh to union with.</param>
    /// <param name="tolerance">Tolerance for the Boolean operation.</param>
    /// <returns>A <see cref="PolyData"/> representing the union.</returns>
    public static PolyData BooleanUnion(this PolyData self, PolyData other, double tolerance = 1e-6)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        throw new NotImplementedException("BooleanUnion requires VTK vtkBooleanOperationPolyDataFilter.");
    }

    /// <summary>
    /// Computes the Boolean intersection of two <see cref="PolyData"/> meshes.
    /// </summary>
    /// <param name="self">The first mesh.</param>
    /// <param name="other">The second mesh to intersect with.</param>
    /// <param name="tolerance">Tolerance for the Boolean operation.</param>
    /// <returns>A <see cref="PolyData"/> representing the intersection.</returns>
    public static PolyData BooleanIntersection(this PolyData self, PolyData other, double tolerance = 1e-6)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        throw new NotImplementedException("BooleanIntersection requires VTK vtkBooleanOperationPolyDataFilter.");
    }

    /// <summary>
    /// Computes the Boolean difference (subtraction) of two <see cref="PolyData"/> meshes.
    /// </summary>
    /// <param name="self">The first mesh.</param>
    /// <param name="other">The mesh to subtract from <paramref name="self"/>.</param>
    /// <param name="tolerance">Tolerance for the Boolean operation.</param>
    /// <returns>A <see cref="PolyData"/> representing the difference.</returns>
    public static PolyData BooleanDifference(this PolyData self, PolyData other, double tolerance = 1e-6)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        throw new NotImplementedException("BooleanDifference requires VTK vtkBooleanOperationPolyDataFilter.");
    }

    // ---------------------------------------------------------------
    //  Intersection
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the intersection lines between two <see cref="PolyData"/> meshes.
    /// </summary>
    /// <param name="self">The first mesh.</param>
    /// <param name="other">The second mesh.</param>
    /// <param name="splitFirst">When <c>true</c>, splits the first mesh along the intersection.</param>
    /// <param name="splitSecond">When <c>true</c>, splits the second mesh along the intersection.</param>
    /// <returns>
    /// A tuple of (<see cref="PolyData"/> intersection lines, <see cref="PolyData"/> split first,
    /// <see cref="PolyData"/> split second).
    /// </returns>
    public static (PolyData IntersectionLines, PolyData SplitFirst, PolyData SplitSecond) Intersection(
        this PolyData self, PolyData other, bool splitFirst = true, bool splitSecond = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        throw new NotImplementedException("Intersection requires VTK vtkIntersectionPolyDataFilter.");
    }

    // ---------------------------------------------------------------
    //  Collision detection
    // ---------------------------------------------------------------

    /// <summary>
    /// Detects collisions between two <see cref="PolyData"/> meshes.
    /// </summary>
    /// <param name="self">The first mesh.</param>
    /// <param name="other">The second mesh.</param>
    /// <param name="cellTolerance">Tolerance for cell overlap detection.</param>
    /// <returns>
    /// A tuple of the number of contacts and a <see cref="PolyData"/> of the contact cells.
    /// </returns>
    public static (int NContacts, PolyData ContactCells) Collision(this PolyData self, PolyData other, double cellTolerance = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);
        throw new NotImplementedException("Collision requires VTK vtkCollisionDetectionFilter.");
    }

    // ---------------------------------------------------------------
    //  Decimation
    // ---------------------------------------------------------------

    /// <summary>
    /// Reduces the number of triangles in a triangulated mesh.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to decimate.</param>
    /// <param name="targetReduction">
    /// Fraction of triangles to remove, in the range <c>[0, 1)</c>.
    /// </param>
    /// <returns>The decimated <see cref="PolyData"/>.</returns>
    public static PolyData Decimate(this PolyData self, double targetReduction = 0.5)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (targetReduction < 0.0 || targetReduction >= 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(targetReduction), "Target reduction must be in [0, 1).");
        }

        throw new NotImplementedException("Decimate requires VTK vtkDecimatePro.");
    }

    /// <summary>
    /// Reduces the number of triangles with finer control over mesh quality.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to decimate.</param>
    /// <param name="targetReduction">
    /// Fraction of triangles to remove, in the range <c>[0, 1)</c>.
    /// </param>
    /// <param name="preserveTopology">When <c>true</c>, preserves the mesh topology.</param>
    /// <param name="featureAngle">Feature angle for edge detection in degrees.</param>
    /// <param name="splitAngle">Angle at which to split the mesh.</param>
    /// <param name="splitting">When <c>true</c>, allows mesh splitting.</param>
    /// <param name="preSplitMesh">When <c>true</c>, pre-splits the mesh before decimation.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The decimated <see cref="PolyData"/>.</returns>
    public static PolyData DecimatePro(
        this PolyData self,
        double targetReduction = 0.5,
        bool preserveTopology = false,
        double featureAngle = 15.0,
        double splitAngle = 75.0,
        bool splitting = true,
        bool preSplitMesh = false,
        bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (targetReduction < 0.0 || targetReduction >= 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(targetReduction), "Target reduction must be in [0, 1).");
        }

        throw new NotImplementedException("DecimatePro requires VTK vtkDecimatePro.");
    }

    // ---------------------------------------------------------------
    //  Smoothing
    // ---------------------------------------------------------------

    /// <summary>
    /// Smooths a <see cref="PolyData"/> mesh using Laplacian smoothing.
    /// </summary>
    /// <param name="self">The mesh to smooth.</param>
    /// <param name="nIterations">Number of smoothing iterations.</param>
    /// <param name="relaxationFactor">Relaxation factor applied per iteration.</param>
    /// <param name="convergence">Convergence criterion.</param>
    /// <param name="featureSmoothing">When <c>true</c>, enables feature edge smoothing.</param>
    /// <param name="boundarySmoothing">When <c>true</c>, allows boundary vertices to move.</param>
    /// <param name="featureAngle">Angle (degrees) defining feature edges.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The smoothed <see cref="PolyData"/>.</returns>
    public static PolyData Smooth(
        this PolyData self,
        int nIterations = 20,
        double relaxationFactor = 0.01,
        double convergence = 0.0,
        bool featureSmoothing = false,
        bool boundarySmoothing = true,
        double featureAngle = 45.0,
        bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Smooth requires VTK vtkSmoothPolyDataFilter.");
    }

    /// <summary>
    /// Smooths a <see cref="PolyData"/> mesh using Taubin smoothing.
    /// <para>
    /// Taubin smoothing alternates between positive and negative relaxation to
    /// avoid shrinkage common with standard Laplacian smoothing.
    /// </para>
    /// </summary>
    /// <param name="self">The mesh to smooth.</param>
    /// <param name="nIterations">Number of smoothing iterations.</param>
    /// <param name="passBand">The pass-band frequency for the filter.</param>
    /// <param name="featureSmoothing">When <c>true</c>, enables feature edge smoothing.</param>
    /// <param name="boundarySmoothing">When <c>true</c>, allows boundary vertices to move.</param>
    /// <param name="featureAngle">Angle (degrees) defining feature edges.</param>
    /// <param name="nonManifoldSmoothing">When <c>true</c>, smooths non-manifold edges.</param>
    /// <param name="normalizeCoordinates">When <c>true</c>, normalizes coordinates before smoothing.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The smoothed <see cref="PolyData"/>.</returns>
    public static PolyData SmoothTaubin(
        this PolyData self,
        int nIterations = 20,
        double passBand = 0.1,
        bool featureSmoothing = false,
        bool boundarySmoothing = true,
        double featureAngle = 45.0,
        bool nonManifoldSmoothing = false,
        bool normalizeCoordinates = false,
        bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("SmoothTaubin requires VTK vtkWindowedSincPolyDataFilter.");
    }

    // ---------------------------------------------------------------
    //  Subdivision
    // ---------------------------------------------------------------

    /// <summary>
    /// Subdivides the mesh to increase resolution.
    /// </summary>
    /// <param name="self">The mesh to subdivide.</param>
    /// <param name="nSubdivisions">Number of subdivision levels.</param>
    /// <param name="subdivisionType">
    /// Subdivision algorithm: <c>"linear"</c>, <c>"loop"</c>, or <c>"butterfly"</c>.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The subdivided <see cref="PolyData"/>.</returns>
    public static PolyData Subdivide(this PolyData self, int nSubdivisions = 1, string subdivisionType = "linear", bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        if (nSubdivisions < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(nSubdivisions), "Number of subdivisions must be >= 1.");
        }

        throw new NotImplementedException("Subdivide requires VTK vtkLinearSubdivisionFilter / vtkLoopSubdivisionFilter / vtkButterflySubdivisionFilter.");
    }

    /// <summary>
    /// Adaptively subdivides triangles in the mesh based on a metric.
    /// </summary>
    /// <param name="self">The mesh to subdivide.</param>
    /// <param name="maxEdgeLength">Maximum edge length. Edges longer than this are split.</param>
    /// <param name="maxTriangleArea">Maximum triangle area. Larger triangles are split.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The subdivided <see cref="PolyData"/>.</returns>
    public static PolyData SubdivideAdaptive(this PolyData self, double? maxEdgeLength = null, double? maxTriangleArea = null, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("SubdivideAdaptive requires VTK vtkAdaptiveSubdivisionFilter.");
    }

    // ---------------------------------------------------------------
    //  Normals
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes point and/or cell normals for the mesh.
    /// </summary>
    /// <param name="self">The mesh to compute normals for.</param>
    /// <param name="cellNormals">When <c>true</c>, computes cell normals.</param>
    /// <param name="pointNormals">When <c>true</c>, computes point normals.</param>
    /// <param name="splitVertices">When <c>true</c>, splits vertices along sharp edges.</param>
    /// <param name="flipNormals">When <c>true</c>, flips all computed normals.</param>
    /// <param name="consistent">When <c>true</c>, enforces consistent normal orientation.</param>
    /// <param name="autoOrientNormals">When <c>true</c>, automatically orients normals outward.</param>
    /// <param name="nonManifoldTraversal">When <c>true</c>, enables non-manifold edge traversal.</param>
    /// <param name="featureAngle">Feature angle (degrees) for splitting.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> with computed normals.</returns>
    public static PolyData ComputeNormals(
        this PolyData self,
        bool cellNormals = true,
        bool pointNormals = true,
        bool splitVertices = false,
        bool flipNormals = false,
        bool consistent = true,
        bool autoOrientNormals = false,
        bool nonManifoldTraversal = true,
        double featureAngle = 30.0,
        bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ComputeNormals requires VTK vtkPolyDataNormals.");
    }

    /// <summary>
    /// Flips all point and cell normals of the mesh.
    /// </summary>
    /// <param name="self">The mesh whose normals are flipped.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The <see cref="PolyData"/> with flipped normals.</returns>
    public static PolyData FlipNormals(this PolyData self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("FlipNormals requires VTK vtkReverseSense.");
    }

    // ---------------------------------------------------------------
    //  Edge extraction
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts all edges from the <see cref="PolyData"/> mesh.
    /// </summary>
    /// <param name="self">The mesh whose edges are extracted.</param>
    /// <returns>A <see cref="PolyData"/> of line cells representing the edges.</returns>
    public static PolyData ExtractEdges(this PolyData self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractEdges requires VTK vtkExtractEdges.");
    }

    // ---------------------------------------------------------------
    //  Tube / Ribbon
    // ---------------------------------------------------------------

    /// <summary>
    /// Generates tubes around each line in the <see cref="PolyData"/>.
    /// </summary>
    /// <param name="self">The line-based <see cref="PolyData"/>.</param>
    /// <param name="radius">Radius of the tubes.</param>
    /// <param name="nSides">Number of sides in the tube cross-section.</param>
    /// <param name="scalars">Name of the scalar array used for variable radius.</param>
    /// <param name="capping">When <c>true</c>, caps the ends of the tubes.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> of tube surfaces.</returns>
    public static PolyData Tube(this PolyData self, double radius = 0.5, int nSides = 20, string? scalars = null, bool capping = true, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Tube requires VTK vtkTubeFilter.");
    }

    /// <summary>
    /// Creates a ribbon surface from line segments in the <see cref="PolyData"/>.
    /// </summary>
    /// <param name="self">The line-based <see cref="PolyData"/>.</param>
    /// <param name="width">Width of the ribbon.</param>
    /// <param name="scalars">Name of the scalar array for variable width.</param>
    /// <param name="angle">Angle of the ribbon surface in degrees.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> of ribbon surfaces.</returns>
    public static PolyData Ribbon(this PolyData self, double width = 0.5, string? scalars = null, double angle = 0.0, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Ribbon requires VTK vtkRibbonFilter.");
    }

    // ---------------------------------------------------------------
    //  Extrusion
    // ---------------------------------------------------------------

    /// <summary>
    /// Extrudes the mesh along a direction vector.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to extrude.</param>
    /// <param name="vector">
    /// Direction and magnitude of the extrusion as a 3-element array.
    /// </param>
    /// <param name="capping">When <c>true</c>, caps the ends of the extrusion.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>An extruded <see cref="PolyData"/>.</returns>
    public static PolyData Extrude(this PolyData self, double[] vector, bool capping = true, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(vector);
        if (vector.Length != 3)
        {
            throw new ArgumentException("Extrusion vector must have exactly 3 elements.", nameof(vector));
        }

        throw new NotImplementedException("Extrude requires VTK vtkLinearExtrusionFilter.");
    }

    /// <summary>
    /// Extrudes the mesh by rotating it about an axis.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to extrude.</param>
    /// <param name="resolution">Number of angular steps in the rotation.</param>
    /// <param name="angle">Total rotation angle in degrees.</param>
    /// <param name="capping">When <c>true</c>, caps the ends of the extrusion.</param>
    /// <param name="translation">Translation along the rotation axis per revolution.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>An extruded <see cref="PolyData"/>.</returns>
    public static PolyData ExtrudeRotate(this PolyData self, int resolution = 30, double angle = 360.0, bool capping = true, double translation = 0.0, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtrudeRotate requires VTK vtkRotationalExtrusionFilter.");
    }

    /// <summary>
    /// Extrudes and trims the mesh against a surface.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to extrude.</param>
    /// <param name="direction">Extrusion direction as a 3-element array.</param>
    /// <param name="trimSurface">A <see cref="PolyData"/> surface to trim against.</param>
    /// <param name="capping">When <c>true</c>, caps the ends of the extrusion.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>An extruded and trimmed <see cref="PolyData"/>.</returns>
    public static PolyData ExtrudeTrim(this PolyData self, double[] direction, PolyData trimSurface, bool capping = true, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(direction);
        ArgumentNullException.ThrowIfNull(trimSurface);
        throw new NotImplementedException("ExtrudeTrim requires VTK vtkTrimmedExtrusionFilter.");
    }

    // ---------------------------------------------------------------
    //  Strips
    // ---------------------------------------------------------------

    /// <summary>
    /// Converts triangles to triangle strips for more efficient rendering.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to strip.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> with triangle strips.</returns>
    public static PolyData Strip(this PolyData self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Strip requires VTK vtkStripper.");
    }

    // ---------------------------------------------------------------
    //  Triangulation
    // ---------------------------------------------------------------

    /// <summary>
    /// Triangulates all polygonal faces in the mesh.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to triangulate.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The triangulated <see cref="PolyData"/>.</returns>
    public static PolyData Triangulate(this PolyData self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Triangulate requires VTK vtkTriangleFilter.");
    }

    // ---------------------------------------------------------------
    //  Delaunay 2D
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the 2-D Delaunay triangulation of the mesh.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to triangulate.</param>
    /// <param name="tolerance">Tolerance for coincident point merging.</param>
    /// <param name="alpha">Alpha value controlling the output triangulation extent.</param>
    /// <param name="offset">Multiplier to control the size of the initial bounding triangulation.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> of the Delaunay triangulation.</returns>
    public static PolyData Delaunay2D(this PolyData self, double tolerance = 1e-5, double alpha = 0.0, double offset = 1.0, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Delaunay2D requires VTK vtkDelaunay2D.");
    }

    // ---------------------------------------------------------------
    //  Cleaning / Merging
    // ---------------------------------------------------------------

    /// <summary>
    /// Cleans the mesh by merging duplicate points, removing unused points,
    /// and removing degenerate cells.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to clean.</param>
    /// <param name="pointMerging">When <c>true</c>, merges coincident points.</param>
    /// <param name="tolerance">Tolerance for point merging (absolute).</param>
    /// <param name="linesTo​Points">When <c>true</c>, converts degenerate lines to points.</param>
    /// <param name="polysToLines">When <c>true</c>, converts degenerate polygons to lines.</param>
    /// <param name="stripsToPolys">When <c>true</c>, converts degenerate strips to polygons.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <param name="absoluteTolerance">
    /// When <c>true</c>, <paramref name="tolerance"/> is an absolute distance.
    /// When <c>false</c>, it is relative to the bounding box diagonal.
    /// </param>
    /// <returns>The cleaned <see cref="PolyData"/>.</returns>
    public static PolyData Clean(
        this PolyData self,
        bool pointMerging = true,
        double tolerance = 0.0,
        bool linesTo​Points = false,
        bool polysToLines = true,
        bool stripsToPolys = true,
        bool inplace = false,
        bool absoluteTolerance = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Clean requires VTK vtkCleanPolyData.");
    }

    /// <summary>
    /// Removes unused points that are not referenced by any cell.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to process.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>The <see cref="PolyData"/> without unused points.</returns>
    public static PolyData RemoveUnusedPoints(this PolyData self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("RemoveUnusedPoints requires custom implementation or VTK vtkCleanPolyData.");
    }

    // ---------------------------------------------------------------
    //  Geodesic
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the geodesic path between two vertices on the mesh surface.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> surface mesh.</param>
    /// <param name="startVertex">Index of the start vertex.</param>
    /// <param name="endVertex">Index of the end vertex.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> representing the geodesic path.</returns>
    public static PolyData Geodesic(this PolyData self, int startVertex, int endVertex, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Geodesic requires VTK vtkDijkstraGraphGeodesicPath.");
    }

    /// <summary>
    /// Computes the geodesic distance between two vertices on the mesh surface.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> surface mesh.</param>
    /// <param name="startVertex">Index of the start vertex.</param>
    /// <param name="endVertex">Index of the end vertex.</param>
    /// <returns>The geodesic distance between the two vertices.</returns>
    public static double GeodesicDistance(this PolyData self, int startVertex, int endVertex)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("GeodesicDistance requires VTK vtkDijkstraGraphGeodesicPath.");
    }

    // ---------------------------------------------------------------
    //  Ray tracing
    // ---------------------------------------------------------------

    /// <summary>
    /// Traces a ray from an origin in a direction and finds intersections with the mesh.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> mesh to trace against.</param>
    /// <param name="origin">Origin of the ray as a 3-element array.</param>
    /// <param name="endPoint">End point of the ray as a 3-element array.</param>
    /// <param name="firstPoint">When <c>true</c>, only the first intersection is returned.</param>
    /// <returns>
    /// A tuple of (intersection points as flat array, cell IDs of intersected cells).
    /// </returns>
    public static (double[] Points, int[] CellIds) RayTrace(this PolyData self, double[] origin, double[] endPoint, bool firstPoint = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(origin);
        ArgumentNullException.ThrowIfNull(endPoint);
        throw new NotImplementedException("RayTrace requires VTK vtkOBBTree.");
    }

    // ---------------------------------------------------------------
    //  Fill holes / Clip closed surface
    // ---------------------------------------------------------------

    /// <summary>
    /// Fills holes in the <see cref="PolyData"/> mesh up to a specified size.
    /// </summary>
    /// <param name="self">The mesh with holes.</param>
    /// <param name="holeSize">Maximum size of the holes to fill.</param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> with holes filled.</returns>
    public static PolyData FillHoles(this PolyData self, double holeSize = 1000.0, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("FillHoles requires VTK vtkFillHolesFilter.");
    }

    /// <summary>
    /// Clips a closed surface and caps the resulting hole.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to clip.</param>
    /// <param name="normal">
    /// Clipping plane normal as a 3-element array. Defaults to <c>[1, 0, 0]</c>.
    /// </param>
    /// <param name="origin">
    /// A point on the clipping plane. When <c>null</c>, the center of the mesh is used.
    /// </param>
    /// <returns>A clipped and capped <see cref="PolyData"/>.</returns>
    public static PolyData ClipClosedSurface(this PolyData self, double[]? normal = null, double[]? origin = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ClipClosedSurface requires VTK vtkClipClosedSurface.");
    }

    // ---------------------------------------------------------------
    //  Curvature
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the curvature of the mesh surface.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> surface.</param>
    /// <param name="curvatureType">
    /// Type of curvature to compute: <c>"mean"</c>, <c>"gaussian"</c>,
    /// <c>"minimum"</c>, or <c>"maximum"</c>.
    /// </param>
    /// <returns>A <see cref="PolyData"/> with a curvature point data array.</returns>
    public static PolyData Curvature(this PolyData self, string curvatureType = "mean")
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Curvature requires VTK vtkCurvatures.");
    }

    // ---------------------------------------------------------------
    //  Point projection
    // ---------------------------------------------------------------

    /// <summary>
    /// Projects all points onto a plane.
    /// </summary>
    /// <param name="self">The <see cref="PolyData"/> to project.</param>
    /// <param name="origin">
    /// Origin of the projection plane. When <c>null</c>, the mesh center is used.
    /// </param>
    /// <param name="normal">
    /// Normal of the projection plane. Defaults to <c>[0, 0, 1]</c>.
    /// </param>
    /// <param name="inplace">When <c>true</c>, modifies the mesh in place.</param>
    /// <returns>A <see cref="PolyData"/> with projected points.</returns>
    public static PolyData ProjectPointsToPlane(this PolyData self, double[]? origin = null, double[]? normal = null, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ProjectPointsToPlane requires custom implementation.");
    }

    // ---------------------------------------------------------------
    //  Arc length
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the arc length of each line segment and stores it as point data.
    /// </summary>
    /// <param name="self">The line-based <see cref="PolyData"/>.</param>
    /// <returns>A <see cref="PolyData"/> with an <c>"arc_length"</c> point data array.</returns>
    public static PolyData ComputeArcLength(this PolyData self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ComputeArcLength requires VTK vtkAppendArcLength.");
    }

    // ---------------------------------------------------------------
    //  Surface reconstruction
    // ---------------------------------------------------------------

    /// <summary>
    /// Reconstructs a surface from a point cloud using Poisson surface reconstruction.
    /// </summary>
    /// <param name="self">A point cloud <see cref="PolyData"/> with normals.</param>
    /// <param name="nNeighbors">Number of neighbors used for estimation.</param>
    /// <param name="sampleSpacing">Spacing for the surface reconstruction.</param>
    /// <returns>A <see cref="PolyData"/> of the reconstructed surface.</returns>
    public static PolyData ReconstructSurface(this PolyData self, int nNeighbors = 20, double? sampleSpacing = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ReconstructSurface requires VTK vtkSurfaceReconstructionFilter.");
    }
}
