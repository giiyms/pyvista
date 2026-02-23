using PyVista.Core;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.CompositeFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="MultiBlock"/> composite datasets and provide
/// operations such as combining blocks, extracting geometry, computing outlines,
/// and applying filters generically across all nested blocks.
/// </para>
/// </summary>
public static class CompositeFilters
{
    // ---------------------------------------------------------------
    //  Combine
    // ---------------------------------------------------------------

    /// <summary>
    /// Combines all blocks into a single <see cref="UnstructuredGrid"/>.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock.combine()</c> method.
    /// Each block is appended into a unified unstructured grid. Optionally,
    /// coincidental points can be merged using the specified tolerance.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset to combine.</param>
    /// <param name="mergePoints">
    /// When <c>true</c>, merges coincidental points from different blocks.
    /// </param>
    /// <param name="tolerance">
    /// The absolute tolerance to use when finding coincident points.
    /// Only used when <paramref name="mergePoints"/> is <c>true</c>.
    /// </param>
    /// <returns>A new <see cref="UnstructuredGrid"/> containing all blocks.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static UnstructuredGrid Combine(this MultiBlock self, bool mergePoints = false, double tolerance = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Combine requires VTK vtkAppendFilter.");
    }

    // ---------------------------------------------------------------
    //  Extract geometry
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts the outer surface geometry of all blocks.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock.extract_geometry()</c>
    /// method. It places a composite geometry filter at the end of a pipeline
    /// before a polydata consumer to extract geometry from all blocks and append
    /// them into a single <see cref="PolyData"/> object.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset.</param>
    /// <returns>A <see cref="PolyData"/> representing the surface of the composite dataset.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData ExtractGeometry(this MultiBlock self)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("ExtractGeometry requires VTK vtkCompositeDataGeometryFilter.");
    }

    // ---------------------------------------------------------------
    //  Outline
    // ---------------------------------------------------------------

    /// <summary>
    /// Produces an outline of the full extent for all blocks in this composite dataset.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock.outline()</c> method.
    /// When <paramref name="nested"/> is <c>true</c>, individual outlines are created
    /// for each nested dataset.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset.</param>
    /// <param name="generateFaces">
    /// When <c>true</c>, generates solid faces for the outline box.
    /// </param>
    /// <param name="nested">
    /// When <c>true</c>, creates individual outlines for each nested dataset
    /// rather than a single bounding outline.
    /// </param>
    /// <returns>A <see cref="PolyData"/> mesh containing the outline.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData Outline(this MultiBlock self, bool generateFaces = false, bool nested = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("Outline requires VTK vtkOutlineFilter.");
    }

    // ---------------------------------------------------------------
    //  Outline corners
    // ---------------------------------------------------------------

    /// <summary>
    /// Produces an outline of the corners for all blocks in this composite dataset.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock.outline_corners()</c>
    /// method. The <paramref name="factor"/> parameter controls the relative size
    /// of the corners to the length of the corresponding bounds.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset.</param>
    /// <param name="factor">
    /// Controls the relative size of the corners to the length of the
    /// corresponding bounds. Defaults to <c>0.2</c>.
    /// </param>
    /// <param name="nested">
    /// When <c>true</c>, creates individual corner outlines for each nested dataset.
    /// </param>
    /// <returns>A <see cref="PolyData"/> mesh containing outlined corners.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static PolyData OutlineCorners(this MultiBlock self, double factor = 0.2, bool nested = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        throw new NotImplementedException("OutlineCorners requires VTK vtkOutlineCornerFilter.");
    }

    // ---------------------------------------------------------------
    //  Generic filter
    // ---------------------------------------------------------------

    /// <summary>
    /// Applies a user-specified function to all nested blocks recursively.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock.generic_filter()</c>
    /// method. The function is applied to each leaf block in the composite
    /// dataset. <c>null</c> blocks are skipped by default.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset to filter.</param>
    /// <param name="function">
    /// A function that accepts a <see cref="DataObject"/> and returns a transformed
    /// <see cref="DataObject"/>.
    /// </param>
    /// <param name="skipNull">
    /// When <c>true</c> (default), <c>null</c> blocks are skipped and passed through.
    /// </param>
    /// <returns>A new <see cref="MultiBlock"/> with the function applied to all blocks.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> or <paramref name="function"/> is <c>null</c>.
    /// </exception>
    public static MultiBlock GenericFilter(
        this MultiBlock self,
        Func<DataObject, DataObject?> function,
        bool skipNull = true)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(function);

        var output = new MultiBlock();
        for (int i = 0; i < self.NBlocks; i++)
        {
            var block = self[i];
            string name = self.GetBlockName(i);

            if (block is null)
            {
                if (skipNull)
                {
                    output.Append(null, name);
                }
                else
                {
                    output.Append(function(null!), name);
                }
            }
            else if (block is MultiBlock nested)
            {
                output.Append(GenericFilter(nested, function, skipNull), name);
            }
            else
            {
                output.Append(function(block), name);
            }
        }

        return output;
    }

    // ---------------------------------------------------------------
    //  Compute normals
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes point and/or cell normals for all blocks in the composite dataset.
    /// <para>
    /// This is the C# equivalent of the Python <c>MultiBlock._compute_normals()</c>
    /// method. All blocks must be <see cref="PolyData"/> for this filter to succeed.
    /// </para>
    /// </summary>
    /// <param name="self">The composite dataset. All blocks must be <see cref="PolyData"/>.</param>
    /// <param name="cellNormals">When <c>true</c>, computes cell normals.</param>
    /// <param name="pointNormals">When <c>true</c>, computes point normals.</param>
    /// <param name="splitVertices">
    /// When <c>true</c>, splits sharp vertices to create separate normals.
    /// </param>
    /// <param name="flipNormals">When <c>true</c>, reverses the direction of normals.</param>
    /// <param name="consistentNormals">
    /// When <c>true</c>, reorders polygons to ensure consistent normal orientation.
    /// </param>
    /// <param name="autoOrientNormals">
    /// When <c>true</c>, automatically orients normals to point outward.
    /// </param>
    /// <param name="featureAngle">
    /// Feature angle used to determine sharp edges when splitting vertices.
    /// </param>
    /// <returns>A new <see cref="MultiBlock"/> with computed normals.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the multiblock contains non-<see cref="PolyData"/> datasets.
    /// </exception>
    public static MultiBlock ComputeNormals(
        this MultiBlock self,
        bool cellNormals = true,
        bool pointNormals = true,
        bool splitVertices = false,
        bool flipNormals = false,
        bool consistentNormals = true,
        bool autoOrientNormals = false,
        double featureAngle = 30.0)
    {
        ArgumentNullException.ThrowIfNull(self);

        if (!self.IsAllPolyData)
        {
            throw new InvalidOperationException(
                "This multiblock contains non-PolyData datasets. " +
                "Convert all datasets to PolyData with AsPolyDataBlocks first.");
        }

        throw new NotImplementedException("ComputeNormals requires VTK vtkPolyDataNormals.");
    }
}
