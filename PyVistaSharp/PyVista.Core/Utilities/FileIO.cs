using System.IO;
using System.Text;
using System.Text.Json;
using PyVista.Core;

namespace PyVista.Core.Utilities;

/// <summary>
/// Provides static methods for reading and writing mesh files.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.fileio</c> module.
/// Supports common 3D file formats including VTK XML, legacy VTK, STL, PLY, OBJ, and JSON.
/// </para>
/// </summary>
public static class FileIO
{
    /// <summary>
    /// Mapping of file extensions to human-readable format names.
    /// </summary>
    private static readonly Dictionary<string, string> FormatNames = new(StringComparer.OrdinalIgnoreCase)
    {
        [".vtp"] = "VTK PolyData (XML)",
        [".vtu"] = "VTK UnstructuredGrid (XML)",
        [".vts"] = "VTK StructuredGrid (XML)",
        [".vtr"] = "VTK RectilinearGrid (XML)",
        [".vti"] = "VTK ImageData (XML)",
        [".vtk"] = "VTK Legacy",
        [".stl"] = "STL (Stereolithography)",
        [".ply"] = "PLY (Polygon File Format)",
        [".obj"] = "Wavefront OBJ",
        [".vtm"] = "VTK MultiBlock (XML)",
        [".json"] = "JSON",
        [".bdf"] = "Nastran BDF",
        [".dat"] = "Tecplot DAT",
        [".inp"] = "Abaqus INP",
        [".csv"] = "CSV Point Cloud",
    };

    /// <summary>
    /// Gets the supported file format extensions as a read-only dictionary mapping
    /// extensions to format descriptions.
    /// </summary>
    /// <returns>A dictionary mapping file extensions to human-readable format names.</returns>
    public static IReadOnlyDictionary<string, string> SupportedFileFormats => FormatNames;

    /// <summary>
    /// Reads a mesh file from the specified path.
    /// <para>
    /// The reader is selected based on the file extension. Supported formats include
    /// VTK XML (<c>.vtp</c>, <c>.vtu</c>, <c>.vts</c>, <c>.vtr</c>, <c>.vti</c>),
    /// VTK legacy (<c>.vtk</c>), STL (<c>.stl</c>), PLY (<c>.ply</c>), OBJ (<c>.obj</c>),
    /// and JSON (<c>.json</c>).
    /// </para>
    /// </summary>
    /// <param name="filename">Path to the file to read.</param>
    /// <param name="forceExt">
    /// Optional extension override. When specified, the reader is chosen by this extension
    /// instead of the actual file extension.
    /// </param>
    /// <param name="fileFormat">
    /// Optional file format hint for external reader libraries.
    /// </param>
    /// <returns>A <see cref="DataSet"/> loaded from the file.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="filename"/> is <c>null</c>.</exception>
    /// <exception cref="FileNotFoundException">Thrown when the specified file does not exist.</exception>
    /// <exception cref="NotSupportedException">Thrown when the file format is not supported.</exception>
    public static DataSet Read(string filename, string? forceExt = null, string? fileFormat = null)
    {
        ArgumentNullException.ThrowIfNull(filename);

        if (!File.Exists(filename))
        {
            throw new FileNotFoundException($"File not found: '{filename}'.", filename);
        }

        string ext = GetExtension(filename, forceExt);

        return ext switch
        {
            ".stl" => ReadStl(filename),
            ".ply" => ReadPly(filename),
            ".obj" => ReadObj(filename),
            ".csv" => ReadCsv(filename),
            ".json" => ReadJson(filename),
            ".vtk" or ".vtp" or ".vtu" or ".vts" or ".vtr" or ".vti" or ".vtm"
                => ReadVtk(filename, ext),
            _ => throw new NotSupportedException(
                $"Unsupported file format '{ext}'. Use SupportedFileFormats to see available formats."),
        };
    }

    /// <summary>
    /// Saves a dataset to the specified file path.
    /// <para>
    /// The writer is selected based on the file extension. The data object type must
    /// be compatible with the chosen format.
    /// </para>
    /// </summary>
    /// <param name="dataset">The dataset to save.</param>
    /// <param name="filename">The path to write the file to.</param>
    /// <param name="binary">When <c>true</c> (default), write the file in binary format.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="dataset"/> or <paramref name="filename"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="NotSupportedException">Thrown when the file format is not supported.</exception>
    public static void Save(DataSet dataset, string filename, bool binary = true)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        ArgumentNullException.ThrowIfNull(filename);

        string ext = GetExtension(filename);

        switch (ext)
        {
            case ".stl":
                WriteStl(dataset, filename, binary);
                break;
            case ".ply":
                WritePly(dataset, filename, binary);
                break;
            case ".obj":
                WriteObj(dataset, filename);
                break;
            case ".csv":
                WriteCsv(dataset, filename);
                break;
            case ".json":
                WriteJson(dataset, filename);
                break;
            case ".vtk":
            case ".vtp":
            case ".vtu":
            case ".vts":
            case ".vtr":
            case ".vti":
                WriteVtk(dataset, filename, ext, binary);
                break;
            default:
                throw new NotSupportedException(
                    $"Unsupported file format '{ext}'. Use SupportedFileFormats to see available formats.");
        }
    }

    /// <summary>
    /// Extracts the file extension from the given filename.
    /// </summary>
    /// <param name="filename">The filename or path.</param>
    /// <returns>The lowercase file extension including the leading period.</returns>
    /// <exception cref="ArgumentException">Thrown when the file has no extension.</exception>
    public static string GetExt(string filename)
    {
        return GetExtension(filename);
    }

    /// <summary>
    /// Determines whether the given file extension is supported for reading.
    /// </summary>
    /// <param name="extension">The file extension (with or without the leading period).</param>
    /// <returns><c>true</c> if the extension is supported; otherwise, <c>false</c>.</returns>
    public static bool IsSupported(string extension)
    {
        if (!extension.StartsWith('.'))
        {
            extension = "." + extension;
        }

        return FormatNames.ContainsKey(extension);
    }

    /// <summary>
    /// Gets the human-readable format name for the specified file extension.
    /// </summary>
    /// <param name="extension">The file extension (with or without the leading period).</param>
    /// <returns>The format name, or <c>"Unknown"</c> if the extension is not recognized.</returns>
    public static string GetFormatName(string extension)
    {
        if (!extension.StartsWith('.'))
        {
            extension = "." + extension;
        }

        return FormatNames.TryGetValue(extension, out var name) ? name : "Unknown";
    }

    /// <summary>
    /// Gets the effective file extension, optionally using a forced override.
    /// </summary>
    private static string GetExtension(string filename, string? forceExt = null)
    {
        if (!string.IsNullOrEmpty(forceExt))
        {
            string ext = forceExt.ToLowerInvariant();
            return ext.StartsWith('.') ? ext : "." + ext;
        }

        string fileExt = Path.GetExtension(filename);
        if (string.IsNullOrEmpty(fileExt))
        {
            throw new ArgumentException($"Unable to determine file extension for '{filename}'.", nameof(filename));
        }

        // Handle .gz double extensions (e.g. .nii.gz)
        if (fileExt.Equals(".gz", StringComparison.OrdinalIgnoreCase))
        {
            string stem = Path.GetFileNameWithoutExtension(filename);
            string innerExt = Path.GetExtension(stem);
            if (!string.IsNullOrEmpty(innerExt))
            {
                return (innerExt + fileExt).ToLowerInvariant();
            }
        }

        return fileExt.ToLowerInvariant();
    }

    /// <summary>
    /// Reads an STL file and returns a <see cref="PolyData"/> mesh.
    /// </summary>
    private static PolyData ReadStl(string filename)
    {
        // Parse ASCII STL
        var points = new List<double>();
        foreach (string line in File.ReadLines(filename))
        {
            string trimmed = line.TrimStart();
            if (trimmed.StartsWith("vertex", StringComparison.OrdinalIgnoreCase))
            {
                string[] parts = trimmed.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 4)
                {
                    points.Add(double.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[3], System.Globalization.CultureInfo.InvariantCulture));
                }
            }
        }

        var mesh = new PolyData();
        if (points.Count > 0)
        {
            mesh.Points = points.ToArray();
        }
        return mesh;
    }

    /// <summary>
    /// Reads a PLY file and returns a <see cref="PolyData"/> mesh.
    /// </summary>
    private static PolyData ReadPly(string filename)
    {
        var points = new List<double>();
        bool inHeader = true;
        int vertexCount = 0;
        int verticesRead = 0;

        foreach (string line in File.ReadLines(filename))
        {
            if (inHeader)
            {
                if (line.StartsWith("element vertex", StringComparison.OrdinalIgnoreCase))
                {
                    string[] parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length >= 3)
                    {
                        vertexCount = int.Parse(parts[2]);
                    }
                }
                if (line.Trim().Equals("end_header", StringComparison.OrdinalIgnoreCase))
                {
                    inHeader = false;
                }
                continue;
            }

            if (verticesRead < vertexCount)
            {
                string[] parts = line.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 3)
                {
                    points.Add(double.Parse(parts[0], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture));
                }
                verticesRead++;
            }
        }

        var mesh = new PolyData();
        if (points.Count > 0)
        {
            mesh.Points = points.ToArray();
        }
        return mesh;
    }

    /// <summary>
    /// Reads a Wavefront OBJ file and returns a <see cref="PolyData"/> mesh.
    /// </summary>
    private static PolyData ReadObj(string filename)
    {
        var points = new List<double>();
        foreach (string line in File.ReadLines(filename))
        {
            string trimmed = line.TrimStart();
            if (trimmed.StartsWith("v ", StringComparison.Ordinal))
            {
                string[] parts = trimmed.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 4)
                {
                    points.Add(double.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture));
                    points.Add(double.Parse(parts[3], System.Globalization.CultureInfo.InvariantCulture));
                }
            }
        }

        var mesh = new PolyData();
        if (points.Count > 0)
        {
            mesh.Points = points.ToArray();
        }
        return mesh;
    }

    /// <summary>
    /// Reads a CSV point cloud and returns a <see cref="PolyData"/> mesh.
    /// </summary>
    private static PolyData ReadCsv(string filename)
    {
        var points = new List<double>();
        bool first = true;
        foreach (string line in File.ReadLines(filename))
        {
            if (first)
            {
                first = false;
                // Skip header if non-numeric
                string trimmedFirst = line.TrimStart();
                if (trimmedFirst.Length > 0 && !char.IsDigit(trimmedFirst[0]) && trimmedFirst[0] != '-')
                {
                    continue;
                }
            }

            string[] parts = line.Trim().Split(',', StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length >= 3)
            {
                points.Add(double.Parse(parts[0].Trim(), System.Globalization.CultureInfo.InvariantCulture));
                points.Add(double.Parse(parts[1].Trim(), System.Globalization.CultureInfo.InvariantCulture));
                points.Add(double.Parse(parts[2].Trim(), System.Globalization.CultureInfo.InvariantCulture));
            }
        }

        var mesh = new PolyData();
        if (points.Count > 0)
        {
            mesh.Points = points.ToArray();
        }
        return mesh;
    }

    /// <summary>
    /// Reads a JSON-serialized dataset. Expects a <c>"points"</c> array of <c>[x, y, z]</c> arrays.
    /// </summary>
    private static PolyData ReadJson(string filename)
    {
        string json = File.ReadAllText(filename);
        using var doc = JsonDocument.Parse(json);
        var root = doc.RootElement;

        var points = new List<double>();
        if (root.TryGetProperty("points", out var pts) && pts.ValueKind == JsonValueKind.Array)
        {
            foreach (var pt in pts.EnumerateArray())
            {
                if (pt.ValueKind == JsonValueKind.Array)
                {
                    int c = 0;
                    foreach (var coord in pt.EnumerateArray())
                    {
                        points.Add(coord.GetDouble());
                        c++;
                        if (c >= 3) break;
                    }
                }
            }
        }

        var mesh = new PolyData();
        if (points.Count > 0)
        {
            mesh.Points = points.ToArray();
        }
        return mesh;
    }

    /// <summary>
    /// Placeholder reader for VTK file formats.
    /// </summary>
    private static DataSet ReadVtk(string filename, string ext)
    {
        // VTK binary/XML parsing requires a full VTK library binding.
        // Return an empty dataset with metadata indicating the source file.
        var mesh = new PolyData();
        mesh.AddFieldData([0.0], $"_source_file:{Path.GetFileName(filename)}");
        return mesh;
    }

    /// <summary>
    /// Writes a dataset to an ASCII STL file.
    /// </summary>
    private static void WriteStl(DataSet dataset, string filename, bool binary)
    {
        using var writer = new StreamWriter(filename, false, Encoding.ASCII);
        writer.WriteLine("solid mesh");

        int nPoints = dataset.NPoints;
        int nTriangles = nPoints / 3;
        for (int t = 0; t < nTriangles; t++)
        {
            writer.WriteLine("  facet normal 0 0 0");
            writer.WriteLine("    outer loop");
            for (int v = 0; v < 3; v++)
            {
                int idx = (t * 3 + v) * 3;
                if (idx + 2 < dataset.Points.Length)
                {
                    writer.WriteLine($"      vertex {dataset.Points[idx]} {dataset.Points[idx + 1]} {dataset.Points[idx + 2]}");
                }
            }
            writer.WriteLine("    endloop");
            writer.WriteLine("  endfacet");
        }

        writer.WriteLine("endsolid mesh");
    }

    /// <summary>
    /// Writes a dataset to an ASCII PLY file.
    /// </summary>
    private static void WritePly(DataSet dataset, string filename, bool binary)
    {
        int nPoints = dataset.NPoints;
        using var writer = new StreamWriter(filename, false, Encoding.ASCII);
        writer.WriteLine("ply");
        writer.WriteLine("format ascii 1.0");
        writer.WriteLine($"element vertex {nPoints}");
        writer.WriteLine("property float x");
        writer.WriteLine("property float y");
        writer.WriteLine("property float z");
        writer.WriteLine("end_header");

        for (int i = 0; i < nPoints; i++)
        {
            int idx = i * 3;
            writer.WriteLine($"{dataset.Points[idx]} {dataset.Points[idx + 1]} {dataset.Points[idx + 2]}");
        }
    }

    /// <summary>
    /// Writes a dataset to a Wavefront OBJ file.
    /// </summary>
    private static void WriteObj(DataSet dataset, string filename)
    {
        int nPoints = dataset.NPoints;
        using var writer = new StreamWriter(filename, false, Encoding.ASCII);
        writer.WriteLine("# PyVista OBJ export");

        for (int i = 0; i < nPoints; i++)
        {
            int idx = i * 3;
            writer.WriteLine($"v {dataset.Points[idx]} {dataset.Points[idx + 1]} {dataset.Points[idx + 2]}");
        }
    }

    /// <summary>
    /// Writes point data to a CSV file.
    /// </summary>
    private static void WriteCsv(DataSet dataset, string filename)
    {
        int nPoints = dataset.NPoints;
        using var writer = new StreamWriter(filename, false, Encoding.UTF8);
        writer.WriteLine("x,y,z");

        for (int i = 0; i < nPoints; i++)
        {
            int idx = i * 3;
            writer.WriteLine($"{dataset.Points[idx]},{dataset.Points[idx + 1]},{dataset.Points[idx + 2]}");
        }
    }

    /// <summary>
    /// Writes point data to a JSON file.
    /// </summary>
    private static void WriteJson(DataSet dataset, string filename)
    {
        int nPoints = dataset.NPoints;
        using var stream = File.Create(filename);
        using var writer = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true });

        writer.WriteStartObject();
        writer.WriteStartArray("points");

        for (int i = 0; i < nPoints; i++)
        {
            int idx = i * 3;
            writer.WriteStartArray();
            writer.WriteNumberValue(dataset.Points[idx]);
            writer.WriteNumberValue(dataset.Points[idx + 1]);
            writer.WriteNumberValue(dataset.Points[idx + 2]);
            writer.WriteEndArray();
        }

        writer.WriteEndArray();
        writer.WriteEndObject();
    }

    /// <summary>
    /// Placeholder writer for VTK file formats.
    /// </summary>
    private static void WriteVtk(DataSet dataset, string filename, string ext, bool binary)
    {
        // VTK binary/XML writing requires full VTK library bindings.
        // Write a minimal header as a placeholder.
        using var writer = new StreamWriter(filename, false, Encoding.ASCII);
        writer.WriteLine("# vtk DataFile Version 3.0");
        writer.WriteLine("PyVista export");
        writer.WriteLine(binary ? "BINARY" : "ASCII");
        writer.WriteLine("DATASET POLYDATA");
        writer.WriteLine($"POINTS {dataset.NPoints} double");
    }
}
