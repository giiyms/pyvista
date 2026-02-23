namespace PyVista.Core;

/// <summary>
/// Abstract base class for non-pointset grids (rectilinear grids, image data).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.Grid</c> class.
/// It provides common functionality shared by grid types that define their
/// geometry implicitly through dimensions and coordinates rather than
/// explicit point arrays.
/// </para>
/// </summary>
public abstract class Grid : DataSet
{
    private int _dimX = 1;
    private int _dimY = 1;
    private int _dimZ = 1;

    /// <summary>
    /// Initializes a new instance of the <see cref="Grid"/> class.
    /// </summary>
    protected Grid()
    {
    }

    /// <summary>
    /// Gets or sets the grid dimensions as a tuple (NX, NY, NZ).
    /// <para>
    /// These are effectively the number of points along each of the
    /// three dataset axes.
    /// </para>
    /// </summary>
    public virtual (int NX, int NY, int NZ) Dimensions
    {
        get => (_dimX, _dimY, _dimZ);
        set
        {
            _dimX = value.NX;
            _dimY = value.NY;
            _dimZ = value.NZ;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        attrs.Add(("Dimensions", $"{_dimX}, {_dimY}, {_dimZ}"));
        return attrs;
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is Grid g)
        {
            _dimX = g._dimX;
            _dimY = g._dimY;
            _dimZ = g._dimZ;
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is Grid g)
        {
            _dimX = g._dimX;
            _dimY = g._dimY;
            _dimZ = g._dimZ;
        }
    }
}
