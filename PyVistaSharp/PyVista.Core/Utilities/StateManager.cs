namespace PyVista.Core.Utilities;

/// <summary>
/// Abstract base class for managing a global state variable with context-manager semantics.
/// <para>
/// This is the C# equivalent of the Python <c>_StateManager</c> class. Subclasses define
/// specific state properties and the set of valid options. The state can be retrieved,
/// set globally, or changed temporarily within a <c>using</c> block.
/// </para>
/// </summary>
/// <typeparam name="T">The type of the state value.</typeparam>
/// <example>
/// <code>
/// // Define a concrete state manager:
/// var verbosity = new VtkVerbosity();
/// Console.WriteLine(verbosity.Get()); // "info"
///
/// verbosity.Set("max");
/// Console.WriteLine(verbosity.Get()); // "max"
///
/// using (verbosity.Scoped("off"))
/// {
///     Console.WriteLine(verbosity.Get()); // "off"
/// }
/// Console.WriteLine(verbosity.Get()); // "max" (restored)
/// </code>
/// </example>
public abstract class StateManager<T> where T : notnull
{
    private readonly T[] _validStates;

    /// <summary>
    /// Initializes a new instance of the <see cref="StateManager{T}"/> class.
    /// </summary>
    /// <param name="validStates">
    /// The complete set of allowed state values. An <see cref="ArgumentException"/>
    /// is thrown by <see cref="Set"/> if a value outside this set is provided.
    /// </param>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="validStates"/> is empty.
    /// </exception>
    protected StateManager(T[] validStates)
    {
        ArgumentNullException.ThrowIfNull(validStates);
        if (validStates.Length == 0)
        {
            throw new ArgumentException("At least one valid state must be provided.", nameof(validStates));
        }
        _validStates = validStates;
    }

    /// <summary>
    /// Gets the current global state value.
    /// </summary>
    /// <returns>The current state.</returns>
    public abstract T Get();

    /// <summary>
    /// Sets the global state value after validation.
    /// </summary>
    /// <param name="state">The new state value. Must be one of the valid states.</param>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="state"/> is not in the valid states set.
    /// </exception>
    public void Set(T state)
    {
        Validate(state);
        SetCore(state);
    }

    /// <summary>
    /// Creates a disposable scope that temporarily changes the state.
    /// <para>
    /// The state is set to <paramref name="temporaryState"/> immediately and
    /// restored to its previous value when the returned scope is disposed.
    /// </para>
    /// </summary>
    /// <param name="temporaryState">The temporary state value.</param>
    /// <returns>An <see cref="IDisposable"/> that restores the original state on disposal.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="temporaryState"/> is not in the valid states set.
    /// </exception>
    public IDisposable Scoped(T temporaryState)
    {
        Validate(temporaryState);
        var original = Get();
        SetCore(temporaryState);
        return new StateScope(this, original);
    }

    /// <summary>
    /// Gets the set of valid state values.
    /// </summary>
    public IReadOnlyList<T> ValidStates => _validStates;

    /// <summary>
    /// Core setter that subclasses implement to actually persist the state.
    /// <para>
    /// This method is called after validation succeeds. It must not throw
    /// <see cref="ArgumentException"/> for valid inputs.
    /// </para>
    /// </summary>
    /// <param name="state">The validated state to persist.</param>
    protected abstract void SetCore(T state);

    /// <summary>
    /// Validates that the given state is in the set of allowed values.
    /// </summary>
    /// <param name="state">The state value to validate.</param>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="state"/> is not a valid state.
    /// </exception>
    protected void Validate(T state)
    {
        foreach (var valid in _validStates)
        {
            if (EqualityComparer<T>.Default.Equals(valid, state))
            {
                return;
            }
        }

        throw new ArgumentException(
            $"Invalid state '{state}'. Valid options are: {string.Join(", ", _validStates)}.",
            nameof(state));
    }

    /// <summary>
    /// Disposable helper that restores a previous state on disposal.
    /// </summary>
    private sealed class StateScope : IDisposable
    {
        private readonly StateManager<T> _manager;
        private readonly T _original;
        private bool _disposed;

        public StateScope(StateManager<T> manager, T original)
        {
            _manager = manager;
            _original = original;
        }

        public void Dispose()
        {
            if (!_disposed)
            {
                _manager.SetCore(_original);
                _disposed = true;
            }
        }
    }
}

/// <summary>
/// Manages VTK logger verbosity level as a global state.
/// <para>
/// This is the C# equivalent of the Python <c>_VTKVerbosity</c> / <c>vtk_verbosity</c>
/// state manager. Without VTK bindings, the state is stored in-process and does not
/// affect any native VTK logger.
/// </para>
/// </summary>
/// <example>
/// <code>
/// var verbosity = VtkVerbosity.Instance;
/// verbosity.Set("max");
/// Console.WriteLine(verbosity.Get()); // "max"
///
/// using (verbosity.Scoped("off"))
/// {
///     Console.WriteLine(verbosity.Get()); // "off"
/// }
/// Console.WriteLine(verbosity.Get()); // "max"
/// </code>
/// </example>
public sealed class VtkVerbosity : StateManager<string>
{
    /// <summary>Valid verbosity levels.</summary>
    private static readonly string[] Levels = { "off", "error", "warning", "info", "max" };

    private static string _currentVerbosity = "info";

    /// <summary>
    /// Shared singleton instance.
    /// </summary>
    public static VtkVerbosity Instance { get; } = new();

    /// <summary>
    /// Initializes a new instance of the <see cref="VtkVerbosity"/> class.
    /// </summary>
    public VtkVerbosity() : base(Levels) { }

    /// <inheritdoc />
    public override string Get() => _currentVerbosity;

    /// <inheritdoc />
    protected override void SetCore(string state) => _currentVerbosity = state;
}

/// <summary>
/// Controls access to VTK's pythonic snake_case API on PyVista-wrapped classes.
/// <para>
/// This is the C# equivalent of the Python <c>_vtkSnakeCase</c> / <c>vtk_snake_case</c>
/// state manager. In C# this has no runtime effect on VTK, but can be queried by
/// interop layers or reflection-based wrappers.
/// </para>
/// </summary>
/// <example>
/// <code>
/// var snakeCase = VtkSnakeCaseAccess.Instance;
/// snakeCase.Set("allow");
/// Console.WriteLine(snakeCase.Get()); // "allow"
/// </code>
/// </example>
public sealed class VtkSnakeCaseAccess : StateManager<string>
{
    /// <summary>Valid access modes.</summary>
    private static readonly string[] Modes = { "allow", "warning", "error" };

    private static string _currentMode = "error";

    /// <summary>
    /// Shared singleton instance.
    /// </summary>
    public static VtkSnakeCaseAccess Instance { get; } = new();

    /// <summary>
    /// Initializes a new instance of the <see cref="VtkSnakeCaseAccess"/> class.
    /// </summary>
    public VtkSnakeCaseAccess() : base(Modes) { }

    /// <inheritdoc />
    public override string Get() => _currentMode;

    /// <inheritdoc />
    protected override void SetCore(string state) => _currentMode = state;
}

/// <summary>
/// Controls whether new attributes can be set on PyVista objects at runtime.
/// <para>
/// This is the C# equivalent of the Python <c>_AllowNewAttributes</c> /
/// <c>allow_new_attributes</c> state manager.
/// </para>
/// <list type="bullet">
///   <item><description><c>"private"</c> – Only attributes with a leading underscore may be set.</description></item>
///   <item><description><c>"true"</c> – Any new attribute may be set.</description></item>
///   <item><description><c>"false"</c> – No new attributes may be set.</description></item>
/// </list>
/// </summary>
/// <example>
/// <code>
/// var attrs = AllowNewAttributes.Instance;
/// attrs.Set("true");
/// Console.WriteLine(attrs.Get()); // "true"
///
/// using (attrs.Scoped("false"))
/// {
///     Console.WriteLine(attrs.Get()); // "false"
/// }
/// Console.WriteLine(attrs.Get()); // "true"
/// </code>
/// </example>
public sealed class AllowNewAttributes : StateManager<string>
{
    /// <summary>Valid attribute modes.</summary>
    private static readonly string[] Modes = { "private", "true", "false" };

    private static string _currentMode = "private";

    /// <summary>
    /// Shared singleton instance.
    /// </summary>
    public static AllowNewAttributes Instance { get; } = new();

    /// <summary>
    /// Initializes a new instance of the <see cref="AllowNewAttributes"/> class.
    /// </summary>
    public AllowNewAttributes() : base(Modes) { }

    /// <inheritdoc />
    public override string Get() => _currentMode;

    /// <inheritdoc />
    protected override void SetCore(string state) => _currentMode = state;

    /// <summary>
    /// Evaluates whether a new attribute with the given name is allowed under the current mode.
    /// </summary>
    /// <param name="attributeName">The name of the attribute to check.</param>
    /// <returns>
    /// <c>true</c> if the attribute is allowed; otherwise <c>false</c>.
    /// </returns>
    public bool IsAllowed(string attributeName)
    {
        ArgumentNullException.ThrowIfNull(attributeName);

        return _currentMode switch
        {
            "true" => true,
            "false" => false,
            "private" => attributeName.StartsWith('_'),
            _ => false,
        };
    }
}

/// <summary>
/// Centralised registry for all PyVista global state managers.
/// <para>
/// Provides a single point of access to query or reset all state managers, which is
/// useful for testing and application initialisation.
/// </para>
/// </summary>
public static class GlobalState
{
    /// <summary>
    /// Gets the VTK verbosity state manager.
    /// </summary>
    public static VtkVerbosity VtkVerbosity => VtkVerbosity.Instance;

    /// <summary>
    /// Gets the VTK snake_case access state manager.
    /// </summary>
    public static VtkSnakeCaseAccess VtkSnakeCase => VtkSnakeCaseAccess.Instance;

    /// <summary>
    /// Gets the attribute creation state manager.
    /// </summary>
    public static AllowNewAttributes AllowNewAttributes => AllowNewAttributes.Instance;

    /// <summary>
    /// Resets all global state managers to their default values.
    /// </summary>
    public static void ResetDefaults()
    {
        VtkVerbosity.Instance.Set("info");
        VtkSnakeCaseAccess.Instance.Set("error");
        AllowNewAttributes.Instance.Set("private");
    }
}
