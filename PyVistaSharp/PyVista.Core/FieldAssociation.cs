namespace PyVista.Core;

/// <summary>
/// Represents which type of VTK field a scalar or vector array is associated with.
/// </summary>
public enum FieldAssociation
{
    /// <summary>Association with point data.</summary>
    Point = 0,

    /// <summary>Association with cell data.</summary>
    Cell = 1,

    /// <summary>No field association.</summary>
    None = 2,

    /// <summary>Association with row data.</summary>
    Row = 6,
}
