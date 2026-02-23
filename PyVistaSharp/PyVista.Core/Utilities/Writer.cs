using System;
using System.Collections.Generic;
using System.IO;

namespace PyVista.Core.Utilities;

/// <summary>
/// Specifies the data format used when writing files.
/// </summary>
public enum DataFormat
{
    /// <summary>Write data in binary format (default).</summary>
    Binary,

    /// <summary>Write data in human-readable ASCII format.</summary>
    Ascii,
}

/// <summary>
/// Specifies the compression algorithm used by XML writers.
/// </summary>
public enum CompressionType
{
    /// <summary>No compression.</summary>
    None,

    /// <summary>ZLib compression (default).</summary>
    ZLib,

    /// <summary>LZ4 compression.</summary>
    Lz4,

    /// <summary>LZMA compression.</summary>
    Lzma,
}

/// <summary>
/// Abstract base class for all PyVista file writers.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.BaseWriter</c>.
/// Subclasses provide format-specific writing logic by overriding
/// <see cref="WriteData"/>.
/// </para>
/// </summary>
public abstract class BaseWriter
{
    private string _path;
    private DataObject _dataObject;

    /// <summary>
    /// Initializes a new instance of the <see cref="BaseWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the file to write.</param>
    /// <param name="dataObject">The data object to write.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="path"/> or <paramref name="dataObject"/> is <c>null</c>.
    /// </exception>
    protected BaseWriter(string path, DataObject dataObject)
    {
        ArgumentNullException.ThrowIfNull(path);
        ArgumentNullException.ThrowIfNull(dataObject);
        _path = path;
        _dataObject = dataObject;
    }

    /// <summary>
    /// Gets or sets the output file path.
    /// </summary>
    public string Path
    {
        get => _path;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _path = value;
        }
    }

    /// <summary>
    /// Gets or sets the data object that will be written.
    /// </summary>
    /// <exception cref="ArgumentNullException">
    /// Thrown when the value is <c>null</c>.
    /// </exception>
    public DataObject DataObject
    {
        get => _dataObject;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _dataObject = value;
        }
    }

    /// <summary>
    /// Gets or sets the data format (binary or ASCII).
    /// </summary>
    public DataFormat DataFormat { get; set; } = DataFormat.Binary;

    /// <summary>
    /// Writes the data object to <see cref="Path"/>.
    /// </summary>
    public void Write()
    {
        var directory = System.IO.Path.GetDirectoryName(_path);
        if (!string.IsNullOrEmpty(directory) && !Directory.Exists(directory))
        {
            Directory.CreateDirectory(directory);
        }

        WriteData();
    }

    /// <summary>
    /// When overridden in a derived class, performs the format-specific write operation.
    /// </summary>
    protected abstract void WriteData();
}

/// <summary>
/// Abstract base class for XML-based VTK writers that support compression.
/// </summary>
public abstract class XMLWriter : BaseWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the file to write.</param>
    /// <param name="dataObject">The data object to write.</param>
    protected XMLWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <summary>
    /// Gets or sets the compression type for the output file.
    /// </summary>
    public CompressionType Compression { get; set; } = CompressionType.ZLib;
}

/// <summary>
/// Writer for legacy VTK DataSet files (<c>.vtk</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.DataSetWriter</c>.
/// </para>
/// </summary>
public sealed class DataSetWriter : BaseWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="DataSetWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.vtk</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public DataSetWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("Legacy VTK DataSet writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for STL (stereolithography) files (<c>.stl</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.STLWriter</c>.
/// </para>
/// </summary>
public sealed class STLWriter : BaseWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="STLWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.stl</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public STLWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("STL writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for PLY polygon files (<c>.ply</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PLYWriter</c>.
/// </para>
/// </summary>
public sealed class PLYWriter : BaseWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="PLYWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.ply</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public PLYWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <summary>
    /// Gets or sets the name of a texture array to include in the output.
    /// </summary>
    public string? TextureArrayName { get; set; }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("PLY writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for Wavefront OBJ files (<c>.obj</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.OBJWriter</c>.
/// </para>
/// </summary>
public sealed class OBJWriter : BaseWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="OBJWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.obj</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public OBJWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("OBJ writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for VTK XML PolyData files (<c>.vtp</c>).
/// </summary>
public sealed class XMLPolyDataWriter : XMLWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLPolyDataWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.vtp</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public XMLPolyDataWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("VTK XML PolyData writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for VTK XML UnstructuredGrid files (<c>.vtu</c>).
/// </summary>
public sealed class XMLUnstructuredGridWriter : XMLWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLUnstructuredGridWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.vtu</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public XMLUnstructuredGridWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("VTK XML UnstructuredGrid writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for VTK XML ImageData files (<c>.vti</c>).
/// </summary>
public sealed class XMLImageDataWriter : XMLWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLImageDataWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.vti</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public XMLImageDataWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("VTK XML ImageData writing is not yet implemented.");
    }
}

/// <summary>
/// Writer for VTK XML MultiBlock files (<c>.vtm</c>).
/// </summary>
public sealed class XMLMultiBlockDataWriter : XMLWriter
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLMultiBlockDataWriter"/> class.
    /// </summary>
    /// <param name="path">Path of the output <c>.vtm</c> file.</param>
    /// <param name="dataObject">The data object to write.</param>
    public XMLMultiBlockDataWriter(string path, DataObject dataObject) : base(path, dataObject) { }

    /// <inheritdoc />
    protected override void WriteData()
    {
        throw new NotImplementedException("VTK XML MultiBlock writing is not yet implemented.");
    }
}

/// <summary>
/// Static helper class for saving data objects to files.
/// <para>
/// The appropriate writer is selected automatically based on the file extension.
/// </para>
/// </summary>
public static class WriterHelper
{
    private static readonly Dictionary<string, Func<string, DataObject, BaseWriter>> Writers =
        new(StringComparer.OrdinalIgnoreCase)
        {
            { ".vtk", (p, d) => new DataSetWriter(p, d) },
            { ".stl", (p, d) => new STLWriter(p, d) },
            { ".ply", (p, d) => new PLYWriter(p, d) },
            { ".obj", (p, d) => new OBJWriter(p, d) },
            { ".vtp", (p, d) => new XMLPolyDataWriter(p, d) },
            { ".vtu", (p, d) => new XMLUnstructuredGridWriter(p, d) },
            { ".vti", (p, d) => new XMLImageDataWriter(p, d) },
            { ".vtm", (p, d) => new XMLMultiBlockDataWriter(p, d) },
        };

    /// <summary>
    /// Saves a <see cref="DataObject"/> to the specified file path.
    /// <para>
    /// The writer is chosen automatically from the file extension of
    /// <paramref name="path"/>. Use <paramref name="dataFormat"/> to
    /// control binary vs. ASCII output.
    /// </para>
    /// </summary>
    /// <param name="dataObject">The data object to save.</param>
    /// <param name="path">Destination file path.</param>
    /// <param name="dataFormat">
    /// Data format for the output file. Defaults to <see cref="DataFormat.Binary"/>.
    /// </param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="dataObject"/> or <paramref name="path"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="NotSupportedException">
    /// Thrown when no writer is available for the file extension.
    /// </exception>
    public static void Save(DataObject dataObject, string path, DataFormat dataFormat = DataFormat.Binary)
    {
        ArgumentNullException.ThrowIfNull(dataObject);
        ArgumentNullException.ThrowIfNull(path);

        var ext = System.IO.Path.GetExtension(path);
        if (string.IsNullOrEmpty(ext))
        {
            throw new NotSupportedException($"Cannot determine file extension for '{path}'.");
        }

        if (!Writers.TryGetValue(ext, out var factory))
        {
            throw new NotSupportedException(
                $"No writer registered for the '{ext}' extension.");
        }

        var writer = factory(path, dataObject);
        writer.DataFormat = dataFormat;
        writer.Write();
    }

    /// <summary>
    /// Returns the supported file extensions for writing.
    /// </summary>
    /// <returns>A read-only collection of supported extensions.</returns>
    public static IReadOnlyCollection<string> SupportedExtensions => Writers.Keys;
}
