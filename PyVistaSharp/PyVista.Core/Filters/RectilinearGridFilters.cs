using PyVista.Core;
using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.RectilinearGridFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="RectilinearGrid"/> datasets and provide
/// operations such as converting to tetrahedral meshes and casting to other grid
/// types like <see cref="StructuredGrid"/>.
/// </para>
/// </summary>
public static class RectilinearGridFilters
{
    // ---------------------------------------------------------------
    //  Cast to StructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this <see cref="RectilinearGrid"/> to a <see cref="StructuredGrid"/>.
    /// <para>
    /// The resulting structured grid has explicit point coordinates computed from
    /// the meshgrid of the X, Y, and Z coordinate arrays. Point data, cell data,
    /// and field data are copied to the new grid.
    /// </para>
    /// <para>
    /// This is the C# equivalent of the Python
    /// <c>RectilinearGrid.cast_to_structured_grid()</c> method.
    /// </para>
    /// </summary>
    /// <param name="self">The rectilinear grid to cast.</param>
    /// <returns>A new <see cref="StructuredGrid"/> with the same geometry.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <remarks>
    /// <para>
    /// Note that <see cref="RectilinearGrid"/> already has a
    /// <see cref="RectilinearGrid.CastToStructuredGrid"/> instance method.
    /// This extension method provides a consistent API through the filters namespace.
    /// </para>
    /// </remarks>
    public static StructuredGrid CastToStructuredGrid(this RectilinearGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);
        return self.CastToStructuredGrid();
    }

    // ---------------------------------------------------------------
    //  Cast to UnstructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this <see cref="RectilinearGrid"/> to an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// All cells of the rectilinear grid are converted to hexahedral cells.
    /// Each cell is defined by the eight corner points of its axis-aligned
    /// bounding box. Point data, cell data, and field data are copied to
    /// the new grid.
    /// </para>
    /// <para>
    /// This is the C# equivalent of the Python
    /// <c>RectilinearGrid.cast_to_unstructured_grid()</c> method.
    /// </para>
    /// </summary>
    /// <param name="self">The rectilinear grid to cast.</param>
    /// <returns>A new <see cref="UnstructuredGrid"/> with hexahedral cells.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <remarks>
    /// <para>
    /// Note that <see cref="RectilinearGrid"/> already has a
    /// <see cref="RectilinearGrid.CastToUnstructuredGrid"/> instance method.
    /// This extension method provides a consistent API through the filters namespace.
    /// </para>
    /// </remarks>
    public static UnstructuredGrid CastToUnstructuredGrid(this RectilinearGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);
        return self.CastToUnstructuredGrid();
    }

    // ---------------------------------------------------------------
    //  To tetrahedra
    // ---------------------------------------------------------------

    /// <summary>
    /// Creates a tetrahedral mesh from this rectilinear grid.
    /// <para>
    /// This is the C# equivalent of the Python
    /// <c>RectilinearGrid.to_tetrahedra()</c> method.
    /// Each hexahedral cell in the rectilinear grid is subdivided into
    /// tetrahedra. The number of tetrahedra per cell can be 5, 6, or 12.
    /// </para>
    /// </summary>
    /// <param name="self">The rectilinear grid to convert.</param>
    /// <param name="tetraPerCell">
    /// The number of tetrahedra to generate per cell. Valid values are
    /// <c>5</c>, <c>6</c>, or <c>12</c>. Defaults to <c>5</c>.
    /// </param>
    /// <param name="passCellIds">
    /// When <c>true</c>, the output tetrahedra will have a cell data array named
    /// <c>"vtkOriginalCellIds"</c> indicating which original cell they came from.
    /// </param>
    /// <param name="passData">
    /// When <c>true</c>, copies cell data from the original grid to the tetrahedra.
    /// Implies <paramref name="passCellIds"/> is <c>true</c>.
    /// </param>
    /// <returns>
    /// A new <see cref="UnstructuredGrid"/> containing tetrahedral cells.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="tetraPerCell"/> is not 5, 6, or 12.
    /// </exception>
    /// <exception cref="InvalidOperationException">
    /// Thrown when <paramref name="tetraPerCell"/> is 12 and at least one
    /// grid dimension is 1, which would cause a degenerate subdivision.
    /// </exception>
    public static UnstructuredGrid ToTetrahedra(
        this RectilinearGrid self,
        int tetraPerCell = 5,
        bool passCellIds = true,
        bool passData = true)
    {
        ArgumentNullException.ThrowIfNull(self);

        if (tetraPerCell is not (5 or 6 or 12))
        {
            throw new ArgumentException(
                $"tetraPerCell should be either 5, 6, or 12, not {tetraPerCell}.",
                nameof(tetraPerCell));
        }

        var dims = self.Dimensions;
        if (tetraPerCell == 12 && (dims.NX == 1 || dims.NY == 1 || dims.NZ == 1))
        {
            throw new InvalidOperationException(
                "Cannot split cells into 12 tetrahedra when at least one dimension is 1. " +
                $"Dimensions are ({dims.NX}, {dims.NY}, {dims.NZ}).");
        }

        throw new NotImplementedException("ToTetrahedra requires VTK vtkRectilinearGridToTetrahedra.");
    }

    // ---------------------------------------------------------------
    //  To tetrahedra (mixed)
    // ---------------------------------------------------------------

    /// <summary>
    /// Creates a tetrahedral mesh with mixed subdivision from this rectilinear grid.
    /// <para>
    /// Subdivides some cells into 5 and some into 12 tetrahedra, based on the
    /// per-cell subdivision counts provided in <paramref name="mixed"/>.
    /// </para>
    /// </summary>
    /// <param name="self">The rectilinear grid to convert.</param>
    /// <param name="mixed">
    /// Array of per-cell subdivision counts. Each element must be either <c>5</c>
    /// or <c>12</c>. The length must match the number of cells in the grid.
    /// </param>
    /// <param name="passCellIds">
    /// When <c>true</c>, includes original cell ID mapping in the output.
    /// </param>
    /// <param name="passData">
    /// When <c>true</c>, copies cell data from the original grid.
    /// </param>
    /// <returns>
    /// A new <see cref="UnstructuredGrid"/> containing tetrahedral cells.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> or <paramref name="mixed"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when any element of <paramref name="mixed"/> is not 5 or 12,
    /// or when the array length does not match the number of cells.
    /// </exception>
    public static UnstructuredGrid ToTetrahedraMixed(
        this RectilinearGrid self,
        int[] mixed,
        bool passCellIds = true,
        bool passData = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(mixed);

        int nCells = self.NCellsFromDimensions;
        if (mixed.Length != nCells)
        {
            throw new ArgumentException(
                $"Mixed array length ({mixed.Length}) must match the number of cells ({nCells}).",
                nameof(mixed));
        }

        foreach (int val in mixed)
        {
            if (val is not (5 or 12))
            {
                throw new ArgumentException(
                    $"Each element of mixed must be either 5 or 12, got {val}.",
                    nameof(mixed));
            }
        }

        throw new NotImplementedException("ToTetrahedraMixed requires VTK vtkRectilinearGridToTetrahedra.");
    }
}
