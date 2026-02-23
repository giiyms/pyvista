using System;
using System.Collections.Generic;
using System.IO;

namespace PyVista.Core.Utilities;

/// <summary>
/// Provides a mapping from file extensions to <see cref="BaseReader"/> subtypes.
/// <para>
/// Call <see cref="GetReader(string, string?)"/> to obtain a reader that is
/// appropriate for the given file path.
/// </para>
/// </summary>
public static class ReaderRegistry
{
    private static readonly Dictionary<string, Func<string, BaseReader>> Readers = new(StringComparer.OrdinalIgnoreCase)
    {
        { ".vtp", path => new XMLPolyDataReader(path) },
        { ".vtu", path => new XMLUnstructuredGridReader(path) },
        { ".vti", path => new XMLImageDataReader(path) },
        { ".vtr", path => new XMLRectilinearGridReader(path) },
        { ".vts", path => new XMLStructuredGridReader(path) },
        { ".vtm", path => new XMLMultiBlockDataReader(path) },
        { ".vtk", path => new VTKDataSetReader(path) },
        { ".stl", path => new STLReader(path) },
        { ".ply", path => new PLYReader(path) },
        { ".obj", path => new OBJReader(path) },
        { ".byu", path => new BYUReader(path) },
    };

    /// <summary>
    /// Registers a custom reader factory for the specified file extension.
    /// </summary>
    /// <param name="extension">
    /// The file extension including the leading dot (e.g., <c>".xyz"</c>).
    /// </param>
    /// <param name="factory">
    /// A factory function that creates a <see cref="BaseReader"/> from a file path.
    /// </param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="extension"/> or <paramref name="factory"/> is <c>null</c>.
    /// </exception>
    public static void Register(string extension, Func<string, BaseReader> factory)
    {
        ArgumentNullException.ThrowIfNull(extension);
        ArgumentNullException.ThrowIfNull(factory);
        Readers[extension] = factory;
    }

    /// <summary>
    /// Returns a reader instance appropriate for the given file.
    /// <para>
    /// The reader is selected based on the file extension of
    /// <paramref name="filename"/>, or <paramref name="forceExt"/> when provided.
    /// </para>
    /// </summary>
    /// <param name="filename">Path to the file to read.</param>
    /// <param name="forceExt">
    /// Optional extension override (e.g., <c>".vtu"</c>) to force a specific reader.
    /// </param>
    /// <returns>A <see cref="BaseReader"/> subclass capable of reading the file.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="filename"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="NotSupportedException">
    /// Thrown when no reader is registered for the given extension.
    /// </exception>
    public static BaseReader GetReader(string filename, string? forceExt = null)
    {
        ArgumentNullException.ThrowIfNull(filename);
        var ext = forceExt ?? Path.GetExtension(filename);
        if (string.IsNullOrEmpty(ext))
        {
            throw new NotSupportedException($"Cannot determine file extension for '{filename}'.");
        }

        if (Readers.TryGetValue(ext, out var factory))
        {
            return factory(filename);
        }

        throw new NotSupportedException(
            $"No reader registered for the '{ext}' extension. Use ReaderRegistry.Register to add one.");
    }
}

/// <summary>
/// Abstract base class for all PyVista file readers.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.BaseReader</c> class.
/// Subclasses provide format-specific reading logic by overriding
/// <see cref="ReadDataObject"/>.
/// </para>
/// </summary>
public abstract class BaseReader
{
    private string _path;
    private bool _showProgress;
    private string? _progressMessage;

    /// <summary>
    /// Initializes a new instance of the <see cref="BaseReader"/> class.
    /// </summary>
    /// <param name="path">Path to the file to read.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="path"/> is <c>null</c>.
    /// </exception>
    protected BaseReader(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        _path = path;
    }

    /// <summary>
    /// Gets or sets the path of the file to read.
    /// </summary>
    /// <exception cref="ArgumentNullException">
    /// Thrown when the value is <c>null</c>.
    /// </exception>
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
    /// Gets the file name (without directory) of <see cref="Path"/>.
    /// </summary>
    public string FileName => System.IO.Path.GetFileName(_path);

    /// <summary>
    /// Gets the file extension of <see cref="Path"/> including the leading dot.
    /// </summary>
    public string Extension => System.IO.Path.GetExtension(_path);

    /// <summary>
    /// Enable a progress indicator during reading.
    /// </summary>
    /// <param name="message">
    /// Optional progress message. Defaults to <c>"Reading {FileName}"</c>.
    /// </param>
    public void ShowProgress(string? message = null)
    {
        _showProgress = true;
        _progressMessage = message;
    }

    /// <summary>
    /// Disable the progress indicator during reading.
    /// </summary>
    public void HideProgress()
    {
        _showProgress = false;
        _progressMessage = null;
    }

    /// <summary>
    /// Gets a value indicating whether progress reporting is enabled.
    /// </summary>
    public bool IsProgressVisible => _showProgress;

    /// <summary>
    /// Gets the current progress message, or a default message when none is set.
    /// </summary>
    public string ProgressMessage => _progressMessage ?? $"Reading {FileName}";

    /// <summary>
    /// Reads the file and returns the resulting <see cref="DataObject"/>.
    /// </summary>
    /// <returns>The data object produced by reading the file.</returns>
    /// <exception cref="FileNotFoundException">
    /// Thrown when <see cref="Path"/> does not exist.
    /// </exception>
    public DataObject Read()
    {
        if (!File.Exists(_path))
        {
            throw new FileNotFoundException($"The file '{_path}' was not found.", _path);
        }

        if (_showProgress)
        {
            Console.WriteLine(ProgressMessage);
        }

        return ReadDataObject();
    }

    /// <summary>
    /// When overridden in a derived class, performs the actual reading of the file.
    /// </summary>
    /// <returns>A <see cref="DataObject"/> containing the data from the file.</returns>
    protected abstract DataObject ReadDataObject();

    /// <inheritdoc />
    public override string ToString() => $"{GetType().Name}('{_path}')";
}

/// <summary>
/// Abstract base class for VTK XML-format file readers
/// (e.g., <c>.vtp</c>, <c>.vtu</c>, <c>.vti</c>).
/// </summary>
public abstract class XMLReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLReader"/> class.
    /// </summary>
    /// <param name="path">Path to the XML file to read.</param>
    protected XMLReader(string path) : base(path) { }
}

/// <summary>
/// Reader for VTK XML PolyData files (<c>.vtp</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLPolyDataReader</c>.
/// </para>
/// </summary>
public sealed class XMLPolyDataReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLPolyDataReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vtp</c> file.</param>
    public XMLPolyDataReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML PolyData reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for VTK XML UnstructuredGrid files (<c>.vtu</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLUnstructuredGridReader</c>.
/// </para>
/// </summary>
public sealed class XMLUnstructuredGridReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLUnstructuredGridReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vtu</c> file.</param>
    public XMLUnstructuredGridReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML UnstructuredGrid reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for VTK XML ImageData files (<c>.vti</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLImageDataReader</c>.
/// </para>
/// </summary>
public sealed class XMLImageDataReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLImageDataReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vti</c> file.</param>
    public XMLImageDataReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML ImageData reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for VTK XML RectilinearGrid files (<c>.vtr</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLRectilinearGridReader</c>.
/// </para>
/// </summary>
public sealed class XMLRectilinearGridReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLRectilinearGridReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vtr</c> file.</param>
    public XMLRectilinearGridReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML RectilinearGrid reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for VTK XML StructuredGrid files (<c>.vts</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLStructuredGridReader</c>.
/// </para>
/// </summary>
public sealed class XMLStructuredGridReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLStructuredGridReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vts</c> file.</param>
    public XMLStructuredGridReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML StructuredGrid reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for VTK XML MultiBlock files (<c>.vtm</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.XMLMultiBlockDataReader</c>.
/// </para>
/// </summary>
public sealed class XMLMultiBlockDataReader : XMLReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="XMLMultiBlockDataReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vtm</c> file.</param>
    public XMLMultiBlockDataReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("VTK XML MultiBlock reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for legacy VTK DataSet files (<c>.vtk</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.VTKDataSetReader</c>.
/// </para>
/// </summary>
public sealed class VTKDataSetReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="VTKDataSetReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.vtk</c> file.</param>
    public VTKDataSetReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("Legacy VTK DataSet reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for STL (stereolithography) files (<c>.stl</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.STLReader</c>.
/// </para>
/// </summary>
public sealed class STLReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="STLReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.stl</c> file.</param>
    public STLReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("STL reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for PLY polygon files (<c>.ply</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PLYReader</c>.
/// </para>
/// </summary>
public sealed class PLYReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="PLYReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.ply</c> file.</param>
    public PLYReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("PLY reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for Wavefront OBJ files (<c>.obj</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.OBJReader</c>.
/// </para>
/// </summary>
public sealed class OBJReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="OBJReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.obj</c> file.</param>
    public OBJReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("OBJ reading is not yet implemented.");
    }
}

/// <summary>
/// Reader for BYU surface files (<c>.byu</c>).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.BYUReader</c>.
/// </para>
/// </summary>
public sealed class BYUReader : BaseReader
{
    /// <summary>
    /// Initializes a new instance of the <see cref="BYUReader"/> class.
    /// </summary>
    /// <param name="path">Path to the <c>.byu</c> file.</param>
    public BYUReader(string path) : base(path) { }

    /// <inheritdoc />
    protected override DataObject ReadDataObject()
    {
        throw new NotImplementedException("BYU reading is not yet implemented.");
    }
}
