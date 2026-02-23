using PyVista.Core;
using PyVista.Core.Cells;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.UnstructuredGridFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="UnstructuredGrid"/> datasets and provide
/// operations such as casting to other mesh types, subdividing tetrahedra,
/// cleaning duplicate points, and removing unused points.
/// </para>
/// </summary>
public static class UnstructuredGridFilters
{
    // ---------------------------------------------------------------
    //  Cast to PolyData
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts an <see cref="UnstructuredGrid"/> to <see cref="PolyData"/>.
    /// <para>
    /// This extracts the surface of the unstructured grid, converting
    /// 3D cells to their surface representation and passing through any
    /// 2D or 1D cells directly.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid to cast.</param>
    /// <returns>A <see cref="PolyData"/> representing the surface of the grid.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData CastToPolyData(this UnstructuredGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("CastToPolyData requires VTK vtkGeometryFilter.");
    }

    // ---------------------------------------------------------------
    //  Cast to ExplicitStructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts an <see cref="UnstructuredGrid"/> to an <see cref="ExplicitStructuredGrid"/>.
    /// <para>
    /// The unstructured grid should contain only hexahedral cells and must have
    /// dimension information available for the conversion to succeed.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid to cast.</param>
    /// <param name="dimensions">
    /// The structured grid dimensions as a tuple of <c>(NX, NY, NZ)</c>.
    /// Required to define the structured topology.
    /// </param>
    /// <returns>An <see cref="ExplicitStructuredGrid"/> with the same geometry.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the grid does not contain only hexahedral cells or when
    /// the dimensions do not match the number of cells.
    /// </exception>
    public static ExplicitStructuredGrid CastToExplicitStructuredGrid(
        this UnstructuredGrid self,
        (int NX, int NY, int NZ) dimensions)
    {
        ArgumentNullException.ThrowIfNull(self);

        int expectedCells = Math.Max(1, dimensions.NX - 1)
                          * Math.Max(1, dimensions.NY - 1)
                          * Math.Max(1, dimensions.NZ - 1);
        if (self.NCellsTotal != expectedCells)
        {
            throw new ArgumentException(
                $"Number of cells ({self.NCellsTotal}) does not match " +
                $"dimensions ({dimensions.NX}×{dimensions.NY}×{dimensions.NZ}).");
        }

        throw new NotImplementedException(
            "CastToExplicitStructuredGrid requires VTK vtkUnstructuredGridToExplicitStructuredGrid.");
    }

    // ---------------------------------------------------------------
    //  Subdivide tetrahedra
    // ---------------------------------------------------------------

    /// <summary>
    /// Subdivides each tetrahedron into twelve tetrahedrons.
    /// <para>
    /// This is the C# equivalent of the Python
    /// <c>UnstructuredGrid.subdivide_tetra()</c> method. It uses a VTK
    /// algorithm to subdivide each tetrahedral cell into twelve smaller
    /// tetrahedra, refining the mesh resolution.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid containing tetrahedral cells.</param>
    /// <returns>
    /// A new <see cref="UnstructuredGrid"/> containing the subdivided tetrahedrons.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static UnstructuredGrid SubdivideTetra(this UnstructuredGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("SubdivideTetra requires VTK vtkSubdivideTetra.");
    }

    // ---------------------------------------------------------------
    //  Clean
    // ---------------------------------------------------------------

    /// <summary>
    /// Merges duplicate points and removes unused points in an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// This is the C# equivalent of the Python <c>UnstructuredGrid.clean()</c> method.
    /// The filter merges coincident points as defined by the merging tolerance and
    /// optionally removes unused points. It does not modify the topology or change
    /// cell types, but may renumber cell connectivity IDs.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid to clean.</param>
    /// <param name="tolerance">The absolute point merging tolerance.</param>
    /// <param name="removeUnusedPoints">
    /// When <c>true</c>, removes points not referenced by any cell.
    /// </param>
    /// <param name="produceMergeMap">
    /// When <c>true</c>, produces a merge map in the output field data
    /// named <c>"PointMergeMap"</c>.
    /// </param>
    /// <param name="averagePointData">
    /// When <c>true</c>, averages point coordinates and data of merged points.
    /// When <c>false</c>, retains the single remaining merged point's data.
    /// </param>
    /// <param name="mergingArrayName">
    /// Optional name of a point data array. When specified, merged points must be
    /// both geometrically coincident and have matching values in this array.
    /// Overrides <paramref name="tolerance"/> when set.
    /// </param>
    /// <returns>A cleaned <see cref="UnstructuredGrid"/>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static UnstructuredGrid Clean(
        this UnstructuredGrid self,
        double tolerance = 0.0,
        bool removeUnusedPoints = true,
        bool produceMergeMap = true,
        bool averagePointData = true,
        string? mergingArrayName = null)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Clean requires VTK vtkStaticCleanUnstructuredGrid.");
    }

    // ---------------------------------------------------------------
    //  Remove unused points
    // ---------------------------------------------------------------

    /// <summary>
    /// Removes points which are not used by any cells.
    /// <para>
    /// Unlike <see cref="Clean"/>, this filter does not merge points. It only
    /// removes point entries that are not referenced by any cell connectivity.
    /// This is the C# equivalent of the Python
    /// <c>UnstructuredGrid.remove_unused_points()</c> method.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid to process.</param>
    /// <param name="inplace">
    /// When <c>true</c>, modifies the grid in place.
    /// When <c>false</c>, returns a new copy with unused points removed.
    /// </param>
    /// <returns>The grid with unused points removed.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static UnstructuredGrid RemoveUnusedPoints(this UnstructuredGrid self, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(self);

        if (self.IsEmpty)
        {
            return inplace ? self : (UnstructuredGrid)self.Copy(deep: true);
        }

        // Identify which points are referenced by cells
        var usedPoints = new HashSet<int>();
        int[] connectivity = self.CellConnectivity;
        foreach (int pointId in connectivity)
        {
            usedPoints.Add(pointId);
        }

        int nPoints = self.NPoints;
        if (usedPoints.Count == nPoints)
        {
            return inplace ? self : (UnstructuredGrid)self.Copy(deep: true);
        }

        // Build a mapping from old point indices to new contiguous indices
        var oldToNew = new int[nPoints];
        Array.Fill(oldToNew, -1);
        int newIndex = 0;
        for (int i = 0; i < nPoints; i++)
        {
            if (usedPoints.Contains(i))
            {
                oldToNew[i] = newIndex++;
            }
        }

        // Build new points array with only used points
        double[] oldPoints = self.Points;
        var newPoints = new double[newIndex * 3];
        for (int i = 0; i < nPoints; i++)
        {
            if (oldToNew[i] >= 0)
            {
                int srcOff = i * 3;
                int dstOff = oldToNew[i] * 3;
                newPoints[dstOff] = oldPoints[srcOff];
                newPoints[dstOff + 1] = oldPoints[srcOff + 1];
                newPoints[dstOff + 2] = oldPoints[srcOff + 2];
            }
        }

        // Remap cell connectivity
        int[] oldCells = self.Cells;
        var newCells = new int[oldCells.Length];
        int pos = 0;
        while (pos < oldCells.Length)
        {
            int nPts = oldCells[pos];
            newCells[pos] = nPts;
            for (int j = 1; j <= nPts; j++)
            {
                newCells[pos + j] = oldToNew[oldCells[pos + j]];
            }
            pos += nPts + 1;
        }

        var result = new UnstructuredGrid(newCells, self.CellTypes, newPoints, deep: false);

        if (inplace)
        {
            self.DeepCopy(result);
            return self;
        }

        return result;
    }

    // ---------------------------------------------------------------
    //  Delaunay 2D
    // ---------------------------------------------------------------

    /// <summary>
    /// Applies Delaunay 2D triangulation to the unstructured grid.
    /// <para>
    /// This wraps <see cref="PolyDataFilters.Delaunay2D"/> by first extracting
    /// the surface of the unstructured grid.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid.</param>
    /// <param name="tolerance">Tolerance for the Delaunay triangulation.</param>
    /// <param name="alpha">Alpha value for edge length constraint.</param>
    /// <param name="offset">Multiplier for the initial triangulation bounding box.</param>
    /// <returns>A <see cref="PolyData"/> with a Delaunay triangulation.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData Delaunay2D(
        this UnstructuredGrid self,
        double tolerance = 1e-10,
        double alpha = 0.0,
        double offset = 1.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Delaunay2D requires VTK vtkDelaunay2D.");
    }

    // ---------------------------------------------------------------
    //  Reconstruct surface
    // ---------------------------------------------------------------

    /// <summary>
    /// Reconstructs the surface from the unstructured grid point cloud.
    /// <para>
    /// This wraps <see cref="PolyDataFilters.ReconstructSurface"/> by first
    /// casting the unstructured grid points to a surface representation.
    /// </para>
    /// </summary>
    /// <param name="self">The unstructured grid.</param>
    /// <param name="neighborhoodSize">
    /// Number of neighbors used to estimate local surface orientation.
    /// </param>
    /// <param name="sampleSpacing">
    /// Spacing of the 3D sampling grid for surface reconstruction.
    /// </param>
    /// <returns>A <see cref="PolyData"/> with the reconstructed surface.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData ReconstructSurface(
        this UnstructuredGrid self,
        int neighborhoodSize = 20,
        double sampleSpacing = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ReconstructSurface requires VTK vtkSurfaceReconstructionFilter.");
    }
}
