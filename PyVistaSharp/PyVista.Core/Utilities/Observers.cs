using System.Text.RegularExpressions;

namespace PyVista.Core.Utilities;

/// <summary>
/// Enumerates the types of messages that can be observed from a VTK pipeline.
/// </summary>
public enum MessageType
{
    /// <summary>An error message from VTK.</summary>
    Error,

    /// <summary>A warning message from VTK.</summary>
    Warning,

    /// <summary>A debug / informational message from VTK.</summary>
    Debug,
}

/// <summary>
/// Represents a single event captured from a VTK output window or algorithm.
/// <para>
/// This is the C# equivalent of the Python <c>VtkEvent</c> named tuple. It
/// stores structured information parsed from VTK's diagnostic output.
/// </para>
/// </summary>
public sealed class VtkEvent
{
    /// <summary>
    /// Initializes a new instance of the <see cref="VtkEvent"/> class.
    /// </summary>
    /// <param name="kind">The event kind (e.g., "ERROR", "WARNING").</param>
    /// <param name="path">The source file path where the event originated.</param>
    /// <param name="address">The hex memory address of the VTK object.</param>
    /// <param name="alert">The diagnostic message text.</param>
    /// <param name="line">The source line number.</param>
    /// <param name="name">The VTK class name that generated the event.</param>
    public VtkEvent(string kind, string path, string address, string alert, string line, string name)
    {
        Kind = kind;
        Path = path;
        Address = address;
        Alert = alert;
        Line = line;
        Name = name;
    }

    /// <summary>Gets the event kind (e.g., "ERROR", "WARNING").</summary>
    public string Kind { get; }

    /// <summary>Gets the source file path where the event originated.</summary>
    public string Path { get; }

    /// <summary>Gets the hex memory address of the VTK object.</summary>
    public string Address { get; }

    /// <summary>Gets the diagnostic message text.</summary>
    public string Alert { get; }

    /// <summary>Gets the source line number.</summary>
    public string Line { get; }

    /// <summary>Gets the VTK class name that generated the event.</summary>
    public string Name { get; }

    /// <inheritdoc />
    public override string ToString()
    {
        if (!string.IsNullOrEmpty(Kind) && !string.IsNullOrEmpty(Path) &&
            !string.IsNullOrEmpty(Line) && !string.IsNullOrEmpty(Name) &&
            !string.IsNullOrEmpty(Address))
        {
            return $"{Kind}: In {Path}, line {Line}\n{Name} ({Address}): {Alert}".Trim();
        }
        return Alert;
    }
}

/// <summary>
/// Observes events (errors, warnings) from VTK objects and captures their messages.
/// <para>
/// This is the C# equivalent of the Python <c>Observer</c> class. Because VTK
/// bindings are not available in this environment, the observer operates as a
/// standalone message sink that can be invoked manually or connected to any
/// event-producing component via delegates.
/// </para>
/// </summary>
/// <example>
/// <code>
/// var observer = new Observer(MessageType.Error, log: true, storeHistory: true);
/// observer.HandleEvent("ERROR: In vtkFoo.cxx, line 42\nvtkFoo (0xDEAD): Something failed");
/// Console.WriteLine(observer.HasEventOccurred()); // True
/// Console.WriteLine(observer.EventHistory.Count);  // 1
/// </code>
/// </example>
public sealed class Observer
{
    private static readonly Regex MessageRegex = new(
        @"(?<kind>[a-zA-Z]+):\sIn\s(?<path>.+?),\sline\s(?<line>\d+)\r?\n(?<name>\w+) \((?<address>0x[0-9a-fA-F]+)\):\s(?<alert>.+)",
        RegexOptions.Singleline | RegexOptions.Compiled);

    private readonly object _lock = new();
    private bool _eventOccurred;
    private string? _lastMessage;
    private string? _lastRawMessage;
    private bool _isObserving;

    /// <summary>
    /// Initializes a new instance of the <see cref="Observer"/> class.
    /// </summary>
    /// <param name="eventType">The type of message this observer listens for.</param>
    /// <param name="log">
    /// When <c>true</c>, observed messages are forwarded to the registered
    /// <see cref="LogAction"/> callback.
    /// </param>
    /// <param name="storeHistory">
    /// When <c>true</c>, all observed events are appended to <see cref="EventHistory"/>.
    /// </param>
    public Observer(MessageType eventType = MessageType.Error, bool log = true, bool storeHistory = false)
    {
        EventType = eventType;
        Log = log;
        StoreHistory = storeHistory;
        EventHistory = new List<VtkEvent>();
    }

    /// <summary>Gets the type of event this observer listens for.</summary>
    public MessageType EventType { get; }

    /// <summary>Gets or sets a value indicating whether to log observed messages.</summary>
    public bool Log { get; set; }

    /// <summary>Gets or sets a value indicating whether to store event history.</summary>
    public bool StoreHistory { get; set; }

    /// <summary>Gets the list of all events observed when <see cref="StoreHistory"/> is <c>true</c>.</summary>
    public List<VtkEvent> EventHistory { get; }

    /// <summary>
    /// Gets or sets the callback invoked when logging an observed message.
    /// <para>
    /// The first argument is the <see cref="MessageType"/> kind and the second is the
    /// alert text. If <c>null</c>, messages are written to <see cref="Console.Error"/>.
    /// </para>
    /// </summary>
    public Action<MessageType, string>? LogAction { get; set; }

    /// <summary>
    /// Parses a raw VTK diagnostic message into a structured <see cref="VtkEvent"/>.
    /// </summary>
    /// <param name="message">The raw message string from VTK.</param>
    /// <returns>A <see cref="VtkEvent"/> with parsed fields, or an event with only
    /// the <see cref="VtkEvent.Alert"/> populated if the format is unrecognized.</returns>
    public static VtkEvent ParseMessage(string message)
    {
        ArgumentNullException.ThrowIfNull(message);

        var match = MessageRegex.Match(message);
        if (match.Success)
        {
            return new VtkEvent(
                kind: match.Groups["kind"].Value,
                path: match.Groups["path"].Value,
                address: match.Groups["address"].Value,
                alert: match.Groups["alert"].Value.Trim(),
                line: match.Groups["line"].Value,
                name: match.Groups["name"].Value);
        }

        return new VtkEvent(
            kind: string.Empty,
            path: string.Empty,
            address: string.Empty,
            alert: message.Trim(),
            line: string.Empty,
            name: string.Empty);
    }

    /// <summary>
    /// Handles an incoming event message. This is the main entry point for
    /// processing observed messages, equivalent to <c>Observer.__call__</c> in Python.
    /// </summary>
    /// <param name="message">The raw event message from VTK or another source.</param>
    public void HandleEvent(string message)
    {
        ArgumentNullException.ThrowIfNull(message);

        lock (_lock)
        {
            _eventOccurred = true;
            _lastRawMessage = message;

            var vtkEvent = ParseMessage(message);
            _lastMessage = vtkEvent.Alert;

            if (StoreHistory)
            {
                EventHistory.Add(vtkEvent);
            }

            if (Log)
            {
                LogMessage(vtkEvent.Kind, vtkEvent.Alert);
            }
        }
    }

    /// <summary>
    /// Checks whether an event has occurred since the last time this method was called.
    /// <para>
    /// Calling this method resets the internal flag, so subsequent calls will
    /// return <c>false</c> until a new event is observed.
    /// </para>
    /// </summary>
    /// <returns><c>true</c> if an event occurred since the last check; otherwise <c>false</c>.</returns>
    public bool HasEventOccurred()
    {
        lock (_lock)
        {
            bool occurred = _eventOccurred;
            _eventOccurred = false;
            return occurred;
        }
    }

    /// <summary>
    /// Gets the last observed message.
    /// </summary>
    /// <param name="raw">
    /// When <c>true</c>, returns the full unparsed message. When <c>false</c> (default),
    /// returns only the parsed alert text.
    /// </param>
    /// <returns>The last message, or <c>null</c> if no event has been observed.</returns>
    public string? GetMessage(bool raw = false)
    {
        lock (_lock)
        {
            return raw ? _lastRawMessage : _lastMessage;
        }
    }

    /// <summary>
    /// Clears the event history and resets the observer state.
    /// </summary>
    public void Reset()
    {
        lock (_lock)
        {
            _eventOccurred = false;
            _lastMessage = null;
            _lastRawMessage = null;
            EventHistory.Clear();
            _isObserving = false;
        }
    }

    /// <summary>
    /// Marks this observer as actively observing an algorithm or output window.
    /// <para>
    /// In a full VTK environment, this would call <c>AddObserver</c> on the
    /// VTK object. Without VTK, this method validates that the observer has
    /// not already been attached and sets the observing flag.
    /// </para>
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when this observer is already attached to an algorithm.
    /// </exception>
    public void Observe()
    {
        lock (_lock)
        {
            if (_isObserving)
            {
                throw new InvalidOperationException(
                    "This error observer is already observing an algorithm.");
            }
            _isObserving = true;
        }
    }

    /// <summary>
    /// Gets a value indicating whether this observer is currently attached to an algorithm.
    /// </summary>
    public bool IsObserving
    {
        get { lock (_lock) { return _isObserving; } }
    }

    /// <summary>
    /// Routes a message to the appropriate logging destination based on its kind.
    /// </summary>
    /// <param name="kind">The event kind string (e.g., "ERROR", "WARNING").</param>
    /// <param name="alert">The alert text to log.</param>
    private void LogMessage(string kind, string alert)
    {
        if (LogAction is not null)
        {
            var messageType = kind.Equals("ERROR", StringComparison.OrdinalIgnoreCase)
                ? MessageType.Error
                : MessageType.Warning;
            LogAction(messageType, alert);
        }
        else
        {
            Console.Error.WriteLine($"[{EventType}] {alert}");
        }
    }
}

/// <summary>
/// Context manager for temporarily catching VTK errors and warnings.
/// <para>
/// This is the C# equivalent of the Python <c>VtkErrorCatcher</c> context manager.
/// Use it with a <c>using</c> block to capture errors and warnings from VTK operations.
/// </para>
/// </summary>
/// <example>
/// <code>
/// using var catcher = new VtkErrorCatcher(raiseErrors: true, emitWarnings: true);
/// // ... perform VTK operations ...
/// // On dispose, errors and warnings are checked.
/// </code>
/// </example>
public sealed class VtkErrorCatcher : IDisposable
{
    private bool _disposed;

    /// <summary>
    /// Initializes a new instance of the <see cref="VtkErrorCatcher"/> class.
    /// </summary>
    /// <param name="raiseErrors">
    /// When <c>true</c>, an <see cref="InvalidOperationException"/> is thrown on dispose
    /// if any error events were captured.
    /// </param>
    /// <param name="sendToLogging">
    /// When <c>true</c>, observed events are forwarded to the observers' logging callbacks.
    /// </param>
    /// <param name="emitWarnings">
    /// When <c>true</c>, warning events are reported on dispose.
    /// </param>
    public VtkErrorCatcher(bool raiseErrors = false, bool sendToLogging = true, bool emitWarnings = false)
    {
        RaiseErrors = raiseErrors;
        SendToLogging = sendToLogging;
        EmitWarnings = emitWarnings;

        ErrorObserver = new Observer(MessageType.Error, log: sendToLogging, storeHistory: true);
        WarningObserver = new Observer(MessageType.Warning, log: sendToLogging, storeHistory: true);
    }

    /// <summary>Gets a value indicating whether errors are raised on dispose.</summary>
    public bool RaiseErrors { get; }

    /// <summary>Gets a value indicating whether events are forwarded to logging.</summary>
    public bool SendToLogging { get; }

    /// <summary>Gets a value indicating whether warnings are emitted on dispose.</summary>
    public bool EmitWarnings { get; }

    /// <summary>Gets the error observer.</summary>
    public Observer ErrorObserver { get; }

    /// <summary>Gets the warning observer.</summary>
    public Observer WarningObserver { get; }

    /// <summary>
    /// Gets all captured events (both warnings and errors).
    /// </summary>
    public IReadOnlyList<VtkEvent> Events
    {
        get
        {
            var all = new List<VtkEvent>(WarningObserver.EventHistory);
            all.AddRange(ErrorObserver.EventHistory);
            return all;
        }
    }

    /// <summary>
    /// Gets all captured error events.
    /// </summary>
    public IReadOnlyList<VtkEvent> ErrorEvents => ErrorObserver.EventHistory;

    /// <summary>
    /// Gets all captured warning events.
    /// </summary>
    public IReadOnlyList<VtkEvent> WarningEvents => WarningObserver.EventHistory;

    /// <inheritdoc />
    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;

        if (EmitWarnings && WarningObserver.EventHistory.Count > 0)
        {
            var warningMsg = string.Join(Environment.NewLine,
                WarningObserver.EventHistory.ConvertAll(e => e.ToString()));
            Console.Error.WriteLine($"[VTK Warnings] {warningMsg}");
        }

        if (RaiseErrors && ErrorObserver.EventHistory.Count > 0)
        {
            var errorMsg = string.Join(Environment.NewLine,
                ErrorObserver.EventHistory.ConvertAll(e => e.ToString()));
            throw new InvalidOperationException($"VTK execution error(s):\n{errorMsg}");
        }
    }
}

/// <summary>
/// Monitors the progress of a long-running operation and reports it via a callback.
/// <para>
/// This is the C# equivalent of the Python <c>ProgressMonitor</c> class. Without
/// VTK bindings, callers manually report progress via <see cref="ReportProgress"/>.
/// </para>
/// </summary>
/// <example>
/// <code>
/// using var monitor = new ProgressMonitor("Filtering mesh");
/// monitor.ProgressChanged += (_, args) =>
///     Console.Write($"\r{args.Message}: {args.Progress:P0}");
/// // In a loop or callback:
/// monitor.ReportProgress(0.5);
/// </code>
/// </example>
public sealed class ProgressMonitor : IDisposable
{
    private double _progress;
    private bool _disposed;
    private bool _interruptRequested;

    /// <summary>
    /// Initializes a new instance of the <see cref="ProgressMonitor"/> class.
    /// </summary>
    /// <param name="message">A description of the operation being monitored.</param>
    public ProgressMonitor(string message = "")
    {
        Message = message;
    }

    /// <summary>Gets the description of the operation being monitored.</summary>
    public string Message { get; }

    /// <summary>
    /// Gets the current progress value in the range [0.0, 1.0].
    /// </summary>
    public double Progress
    {
        get => _progress;
        private set => _progress = Math.Clamp(value, 0.0, 1.0);
    }

    /// <summary>
    /// Gets a value indicating whether an interrupt was requested.
    /// </summary>
    public bool InterruptRequested => _interruptRequested;

    /// <summary>
    /// Occurs when the progress value changes.
    /// </summary>
    public event EventHandler<ProgressEventArgs>? ProgressChanged;

    /// <summary>
    /// Reports a progress update.
    /// </summary>
    /// <param name="progress">The new progress value in the range [0.0, 1.0].</param>
    public void ReportProgress(double progress)
    {
        Progress = progress;
        ProgressChanged?.Invoke(this, new ProgressEventArgs(Progress, Message));
    }

    /// <summary>
    /// Requests that the monitored operation be interrupted at the next opportunity.
    /// </summary>
    public void RequestInterrupt()
    {
        _interruptRequested = true;
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (!_disposed)
        {
            Progress = 1.0;
            ProgressChanged?.Invoke(this, new ProgressEventArgs(1.0, Message));
            _disposed = true;
        }
    }
}

/// <summary>
/// Provides data for the <see cref="ProgressMonitor.ProgressChanged"/> event.
/// </summary>
public sealed class ProgressEventArgs : EventArgs
{
    /// <summary>
    /// Initializes a new instance of the <see cref="ProgressEventArgs"/> class.
    /// </summary>
    /// <param name="progress">The current progress value.</param>
    /// <param name="message">The operation message.</param>
    public ProgressEventArgs(double progress, string message)
    {
        Progress = progress;
        Message = message;
    }

    /// <summary>Gets the current progress value in the range [0.0, 1.0].</summary>
    public double Progress { get; }

    /// <summary>Gets the operation message.</summary>
    public string Message { get; }
}
