using System.Globalization;

namespace PyVista.Core;

/// <summary>
/// Represents 3D bounds as six values: (XMin, XMax, YMin, YMax, ZMin, ZMax).
/// </summary>
/// <param name="XMin">The minimum X bound.</param>
/// <param name="XMax">The maximum X bound.</param>
/// <param name="YMin">The minimum Y bound.</param>
/// <param name="YMax">The maximum Y bound.</param>
/// <param name="ZMin">The minimum Z bound.</param>
/// <param name="ZMax">The maximum Z bound.</param>
public readonly record struct BoundsTuple(
    double XMin,
    double XMax,
    double YMin,
    double YMax,
    double ZMin,
    double ZMax)
{
    /// <inheritdoc />
    public override string ToString()
    {
        return string.Format(
            CultureInfo.InvariantCulture,
            "BoundsTuple(XMin={0}, XMax={1}, YMin={2}, YMax={3}, ZMin={4}, ZMax={5})",
            XMin, XMax, YMin, YMax, ZMin, ZMax);
    }
}
