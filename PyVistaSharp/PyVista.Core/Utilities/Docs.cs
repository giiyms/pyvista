using System.Reflection;
using System.Text;

namespace PyVista.Core.Utilities;

/// <summary>
/// Custom attribute that marks a member as deprecated with a suggested replacement.
/// <para>
/// This is the C# equivalent of the Python <c>@deprecated</c> decorator pattern used
/// throughout PyVista to guide users towards updated API surfaces.
/// </para>
/// </summary>
/// <example>
/// <code>
/// [PyVistaDeprecated("0.45", "Use NewMethod instead.")]
/// public void OldMethod() { }
/// </code>
/// </example>
[AttributeUsage(AttributeTargets.All, Inherited = false, AllowMultiple = false)]
public sealed class PyVistaDeprecatedAttribute : Attribute
{
    /// <summary>
    /// Initializes a new instance of the <see cref="PyVistaDeprecatedAttribute"/> class.
    /// </summary>
    /// <param name="sinceVersion">The version since which the member has been deprecated.</param>
    /// <param name="message">A message describing the replacement or reason for deprecation.</param>
    public PyVistaDeprecatedAttribute(string sinceVersion, string message)
    {
        SinceVersion = sinceVersion;
        Message = message;
    }

    /// <summary>Gets the version since which the member has been deprecated.</summary>
    public string SinceVersion { get; }

    /// <summary>Gets the deprecation message.</summary>
    public string Message { get; }

    /// <inheritdoc />
    public override string ToString() => $"Deprecated since v{SinceVersion}: {Message}";
}

/// <summary>
/// Custom attribute that records when a member was added to the public API.
/// <para>
/// This is the C# equivalent of the Python <c>.. versionadded::</c> Sphinx directive,
/// allowing the information to be queried at runtime.
/// </para>
/// </summary>
/// <example>
/// <code>
/// [VersionAdded("0.47")]
/// public void NewFeature() { }
/// </code>
/// </example>
[AttributeUsage(AttributeTargets.All, Inherited = false, AllowMultiple = false)]
public sealed class VersionAddedAttribute : Attribute
{
    /// <summary>
    /// Initializes a new instance of the <see cref="VersionAddedAttribute"/> class.
    /// </summary>
    /// <param name="version">The version in which the member was introduced.</param>
    public VersionAddedAttribute(string version)
    {
        Version = version;
    }

    /// <summary>Gets the version in which the member was introduced.</summary>
    public string Version { get; }
}

/// <summary>
/// Custom attribute that records when a member's behavior was changed.
/// <para>
/// This is the C# equivalent of the Python <c>.. versionchanged::</c> Sphinx directive.
/// </para>
/// </summary>
/// <example>
/// <code>
/// [VersionChanged("0.47", "Return type changed from void to bool.")]
/// public bool SomeMethod() => true;
/// </code>
/// </example>
[AttributeUsage(AttributeTargets.All, Inherited = false, AllowMultiple = true)]
public sealed class VersionChangedAttribute : Attribute
{
    /// <summary>
    /// Initializes a new instance of the <see cref="VersionChangedAttribute"/> class.
    /// </summary>
    /// <param name="version">The version in which the change occurred.</param>
    /// <param name="description">A description of what changed.</param>
    public VersionChangedAttribute(string version, string description)
    {
        Version = version;
        Description = description;
    }

    /// <summary>Gets the version in which the change occurred.</summary>
    public string Version { get; }

    /// <summary>Gets a description of what changed.</summary>
    public string Description { get; }
}

/// <summary>
/// Provides documentation and source-code link utilities.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.docs</c> module.
/// It includes methods for resolving source links and inspecting documentation
/// attributes at runtime.
/// </para>
/// </summary>
public static class Docs
{
    /// <summary>
    /// The base GitHub URL used to construct source-code links.
    /// </summary>
    private const string GitHubBaseUrl = "https://github.com/pyvista/pyvista";

    /// <summary>
    /// Resolves a GitHub source-code URL for the specified <see cref="MemberInfo"/>.
    /// <para>
    /// This is the C# equivalent of the Python <c>linkcode_resolve</c> function.
    /// It constructs a URL to the source file on GitHub for the member's declaring type.
    /// </para>
    /// </summary>
    /// <param name="member">The member to resolve a link for.</param>
    /// <param name="branch">The Git branch or tag to link to. Defaults to <c>"main"</c>.</param>
    /// <param name="edit">
    /// When <c>true</c>, returns an edit link instead of a blob (view) link.
    /// </param>
    /// <returns>A URL string, or <c>null</c> if the link cannot be resolved.</returns>
    /// <example>
    /// <code>
    /// var method = typeof(DataObject).GetMethod("ShallowCopy");
    /// string? url = Docs.LinkCodeResolve(method!);
    /// // "https://github.com/pyvista/pyvista/blob/main/...DataObject.cs"
    /// </code>
    /// </example>
    public static string? LinkCodeResolve(MemberInfo member, string branch = "main", bool edit = false)
    {
        ArgumentNullException.ThrowIfNull(member);

        var declaringType = member.DeclaringType ?? member as Type;
        if (declaringType is null)
        {
            return null;
        }

        // Build a relative file path from the namespace and type name
        var ns = declaringType.Namespace;
        if (string.IsNullOrEmpty(ns))
        {
            return null;
        }

        var relativePath = ns.Replace('.', '/');
        var fileName = GetSourceFileName(declaringType);
        var blobOrEdit = edit ? "edit" : "blob";

        return $"{GitHubBaseUrl}/{blobOrEdit}/{branch}/{relativePath}/{fileName}";
    }

    /// <summary>
    /// Returns all <see cref="PyVistaDeprecatedAttribute"/> entries found on the specified type
    /// and its members.
    /// </summary>
    /// <param name="type">The type to inspect.</param>
    /// <returns>
    /// A list of tuples containing the member name and its deprecation attribute.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="type"/> is <c>null</c>.
    /// </exception>
    public static List<(string MemberName, PyVistaDeprecatedAttribute Attribute)> GetDeprecations(Type type)
    {
        ArgumentNullException.ThrowIfNull(type);

        var results = new List<(string, PyVistaDeprecatedAttribute)>();

        // Check the type itself
        var typeAttr = type.GetCustomAttribute<PyVistaDeprecatedAttribute>();
        if (typeAttr is not null)
        {
            results.Add((type.Name, typeAttr));
        }

        // Check all public members
        var members = type.GetMembers(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static);
        foreach (var member in members)
        {
            var attr = member.GetCustomAttribute<PyVistaDeprecatedAttribute>();
            if (attr is not null)
            {
                results.Add((member.Name, attr));
            }
        }

        return results;
    }

    /// <summary>
    /// Returns all <see cref="VersionAddedAttribute"/> entries found on the specified type
    /// and its members.
    /// </summary>
    /// <param name="type">The type to inspect.</param>
    /// <returns>
    /// A list of tuples containing the member name and the version it was added.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="type"/> is <c>null</c>.
    /// </exception>
    public static List<(string MemberName, string Version)> GetVersionAdditions(Type type)
    {
        ArgumentNullException.ThrowIfNull(type);

        var results = new List<(string, string)>();

        var typeAttr = type.GetCustomAttribute<VersionAddedAttribute>();
        if (typeAttr is not null)
        {
            results.Add((type.Name, typeAttr.Version));
        }

        var members = type.GetMembers(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static);
        foreach (var member in members)
        {
            var attr = member.GetCustomAttribute<VersionAddedAttribute>();
            if (attr is not null)
            {
                results.Add((member.Name, attr.Version));
            }
        }

        return results;
    }

    /// <summary>
    /// Generates a summary table of all public members on the given type, including
    /// any deprecation or version information.
    /// </summary>
    /// <param name="type">The type to document.</param>
    /// <returns>A formatted string containing the member summary table.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="type"/> is <c>null</c>.
    /// </exception>
    /// <example>
    /// <code>
    /// string table = Docs.GenerateMemberSummary(typeof(DataObject));
    /// Console.WriteLine(table);
    /// </code>
    /// </example>
    public static string GenerateMemberSummary(Type type)
    {
        ArgumentNullException.ThrowIfNull(type);

        var sb = new StringBuilder();
        sb.AppendLine($"Member summary for {type.FullName ?? type.Name}");
        sb.AppendLine(new string('=', 60));

        var members = type.GetMembers(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly);
        foreach (var member in members)
        {
            sb.Append($"  {member.MemberType,-12} {member.Name}");

            var deprecated = member.GetCustomAttribute<PyVistaDeprecatedAttribute>();
            if (deprecated is not null)
            {
                sb.Append($"  [DEPRECATED since v{deprecated.SinceVersion}]");
            }

            var added = member.GetCustomAttribute<VersionAddedAttribute>();
            if (added is not null)
            {
                sb.Append($"  [Added in v{added.Version}]");
            }

            var changes = member.GetCustomAttributes<VersionChangedAttribute>();
            foreach (var change in changes)
            {
                sb.Append($"  [Changed in v{change.Version}]");
            }

            sb.AppendLine();
        }

        return sb.ToString().TrimEnd();
    }

    /// <summary>
    /// Constructs a documentation URL for the PyVista API reference of the given type.
    /// </summary>
    /// <param name="type">The type to build a documentation URL for.</param>
    /// <param name="baseUrl">
    /// The base documentation URL. Defaults to the PyVista stable docs site.
    /// </param>
    /// <returns>A URL pointing to the type's API reference page.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="type"/> is <c>null</c>.
    /// </exception>
    public static string GetDocumentationUrl(Type type, string baseUrl = "https://docs.pyvista.org/api/")
    {
        ArgumentNullException.ThrowIfNull(type);

        var typeName = type.Name.ToLowerInvariant();
        return $"{baseUrl.TrimEnd('/')}/{typeName}.html";
    }

    /// <summary>
    /// Checks whether the given member carries a <see cref="PyVistaDeprecatedAttribute"/>.
    /// </summary>
    /// <param name="member">The member to check.</param>
    /// <returns><c>true</c> if the member is marked as deprecated; otherwise <c>false</c>.</returns>
    public static bool IsDeprecated(MemberInfo member)
    {
        ArgumentNullException.ThrowIfNull(member);
        return member.GetCustomAttribute<PyVistaDeprecatedAttribute>() is not null;
    }

    /// <summary>
    /// Returns the deprecation message for the specified member, or <c>null</c>
    /// if the member is not deprecated.
    /// </summary>
    /// <param name="member">The member to check.</param>
    /// <returns>The deprecation message, or <c>null</c>.</returns>
    public static string? GetDeprecationMessage(MemberInfo member)
    {
        ArgumentNullException.ThrowIfNull(member);
        return member.GetCustomAttribute<PyVistaDeprecatedAttribute>()?.ToString();
    }

    /// <summary>
    /// Infers the source file name for a type based on its name.
    /// </summary>
    /// <param name="type">The type to infer a file name for.</param>
    /// <returns>The inferred file name (e.g., <c>"DataObject.cs"</c>).</returns>
    private static string GetSourceFileName(Type type)
    {
        // For nested types, use the outermost declaring type's name
        var outermost = type;
        while (outermost.DeclaringType is not null)
        {
            outermost = outermost.DeclaringType;
        }
        return $"{outermost.Name}.cs";
    }
}
