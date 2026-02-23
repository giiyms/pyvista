using System.Globalization;

using PyVista.Core.Cells;

namespace PyVista.Core;

/// <summary>
/// Concrete class for storing a set of points.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PointSet</c> class.
/// A PointSet represents a set of points that specifies the interface for datasets
/// that explicitly use "point" arrays to represent geometry. This class is useful for
/// improving the performance of filters on point clouds.
/// </para>
/// <para>
/// PointSets contain no cells. Attempting cell-based operations will throw
/// <see cref="PointSetCellOperationError"/> or <see cref="PointSetDimensionReductionError"/>.
/// </para>
/// </summary>
public class PointSet : DataSet
{
    /// <summary>
    /// Initializes a new instance of the <see cref="PointSet"/> class (empty).
    /// </summary>
    public PointSet()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PointSet"/> class from a flat point array.
    /// </summary>
    /// <param name="points">
    /// Flat row-major point coordinates (length must be divisible by 3).
    /// </param>
    /// <param name="deep">When <c>true</c>, a copy of the array is stored.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> is <c>null</c>.
    /// </exception>
    public PointSet(double[] points, bool deep = false)
    {
        ArgumentNullException.ThrowIfNull(points);
        Points = deep ? (double[])points.Clone() : points;
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PointSet"/> class from a 2-D point array.
    /// </summary>
    /// <param name="points">An (N, 3) array of point coordinates.</param>
    /// <param name="deep">When <c>true</c>, a copy of the data is stored.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> is <c>null</c>.
    /// </exception>
    public PointSet(double[,] points, bool deep = false)
    {
        ArgumentNullException.ThrowIfNull(points);
        if (deep)
        {
            Points2D = (double[,])points.Clone();
        }
        else
        {
            Points2D = points;
        }
    }

    /// <summary>
    /// Returns the coordinates for the center of mass of the point set.
    /// </summary>
    /// <param name="scalarsWeight">
    /// When <c>true</c>, uses the active scalars as weights (not supported without VTK;
    /// falls back to unweighted).
    /// </param>
    /// <returns>A tuple of (X, Y, Z) center-of-mass coordinates.</returns>
    public (double X, double Y, double Z) CenterOfMass(bool scalarsWeight = false)
    {
        int n = NPoints;
        if (n == 0)
        {
            return (0.0, 0.0, 0.0);
        }

        double[] weights = Array.Empty<double>();
        if (scalarsWeight)
        {
            var scalars = ActiveScalars;
            if (scalars is not null && scalars.Length == n)
            {
                weights = scalars;
            }
        }

        bool useWeights = weights.Length == n;
        double sx = 0, sy = 0, sz = 0, totalWeight = 0;
        var pts = Points;
        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            double w = useWeights ? weights[i] : 1.0;
            sx += pts[offset] * w;
            sy += pts[offset + 1] * w;
            sz += pts[offset + 2] * w;
            totalWeight += w;
        }

        if (Math.Abs(totalWeight) < double.Epsilon)
        {
            return (0.0, 0.0, 0.0);
        }

        return (sx / totalWeight, sy / totalWeight, sz / totalWeight);
    }

    /// <summary>
    /// Converts the points datatype to double precision.
    /// Since points are already stored as <see cref="double"/>, this is a no-op in C#.
    /// </summary>
    /// <returns>This instance.</returns>
    public PointSet PointsToDouble()
    {
        // Points are always double in the C# implementation.
        return this;
    }

    /// <summary>
    /// Translates the mesh by the given offset vector.
    /// </summary>
    /// <param name="xyz">A 3-element array representing the translation vector.</param>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The translated point set.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="xyz"/> does not have exactly 3 elements.
    /// </exception>
    public PointSet Translate(double[] xyz, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(xyz);
        if (xyz.Length != 3)
        {
            throw new ArgumentException("Translation vector must have exactly 3 elements.", nameof(xyz));
        }

        PointSet target = inplace ? this : (PointSet)Copy(deep: true);
        var pts = target.Points;
        int n = target.NPoints;
        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            pts[offset] += xyz[0];
            pts[offset + 1] += xyz[1];
            pts[offset + 2] += xyz[2];
        }

        return target;
    }

    /// <summary>
    /// Casts this point set to a <see cref="PolyData"/> instance.
    /// </summary>
    /// <param name="deep">When <c>true</c>, performs a deep copy of the points.</param>
    /// <returns>A new <see cref="PolyData"/> containing the points as vertex cells.</returns>
    public PolyData CastToPolyData(bool deep = true)
    {
        var pdata = new PolyData(Points, deep: deep);
        foreach (var key in PointData.Keys)
        {
            pdata.PointData.SetArray(
                deep ? (double[])PointData[key].Clone() : PointData[key],
                key);
        }

        return pdata;
    }

    /// <summary>
    /// Returns 0.0 since a <see cref="PointSet"/> has no area.
    /// </summary>
    public double Area => 0.0;

    /// <summary>
    /// Returns 0.0 since a <see cref="PointSet"/> has no volume.
    /// </summary>
    public double Volume => 0.0;

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        return attrs;
    }
}
