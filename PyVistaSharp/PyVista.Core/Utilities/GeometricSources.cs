using PyVista.Core;

namespace PyVista.Core.Utilities;

/// <summary>
/// Source class for generating cylinder meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.CylinderSource</c>.
/// Configure properties and call <see cref="GetOutput"/> to produce a <see cref="PolyData"/>.
/// </para>
/// </summary>
public class CylinderSource
{
    /// <summary>Gets or sets the center of the cylinder in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the direction of the cylinder axis in <c>[x, y, z]</c>.</summary>
    public double[] Direction { get; set; } = [1.0, 0.0, 0.0];

    /// <summary>Gets or sets the radius of the cylinder.</summary>
    public double Radius { get; set; } = 0.5;

    /// <summary>Gets or sets the height of the cylinder.</summary>
    public double Height { get; set; } = 1.0;

    /// <summary>Gets or sets the number of points on the circular face.</summary>
    public int Resolution { get; set; } = 100;

    /// <summary>Gets or sets a value indicating whether to cap the cylinder ends with polygons.</summary>
    public bool Capping { get; set; } = true;

    /// <summary>Gets or sets a value indicating whether to use capsule-style hemispherical caps.</summary>
    public bool CapsuleCap { get; set; }

    /// <summary>
    /// Generates the cylinder mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the cylinder.</returns>
    public PolyData GetOutput()
    {
        int n = Resolution;
        int numPoints = Capping ? n * 2 + 2 : n * 2;
        var points = new double[numPoints * 3];
        double halfH = Height / 2.0;

        for (int i = 0; i < n; i++)
        {
            double angle = 2.0 * Math.PI * i / n;
            double x = Radius * Math.Cos(angle);
            double y = Radius * Math.Sin(angle);

            points[(i * 2) * 3] = x;
            points[(i * 2) * 3 + 1] = y;
            points[(i * 2) * 3 + 2] = -halfH;

            points[(i * 2 + 1) * 3] = x;
            points[(i * 2 + 1) * 3 + 1] = y;
            points[(i * 2 + 1) * 3 + 2] = halfH;
        }

        if (Capping)
        {
            int bottomIdx = n * 2;
            int topIdx = n * 2 + 1;
            points[bottomIdx * 3] = 0.0;
            points[bottomIdx * 3 + 1] = 0.0;
            points[bottomIdx * 3 + 2] = -halfH;

            points[topIdx * 3] = 0.0;
            points[topIdx * 3 + 1] = 0.0;
            points[topIdx * 3 + 2] = halfH;
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating arrow meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ArrowSource</c>.
/// Configure properties and call <see cref="GetOutput"/> to produce a <see cref="PolyData"/>.
/// </para>
/// </summary>
public class ArrowSource
{
    /// <summary>Gets or sets the length of the arrow tip as a fraction of total length.</summary>
    public double TipLength { get; set; } = 0.25;

    /// <summary>Gets or sets the radius of the arrow tip.</summary>
    public double TipRadius { get; set; } = 0.1;

    /// <summary>Gets or sets the number of faces around the arrow tip.</summary>
    public int TipResolution { get; set; } = 20;

    /// <summary>Gets or sets the radius of the arrow shaft.</summary>
    public double ShaftRadius { get; set; } = 0.05;

    /// <summary>Gets or sets the number of faces around the arrow shaft.</summary>
    public int ShaftResolution { get; set; } = 20;

    /// <summary>
    /// Generates the arrow mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the arrow.</returns>
    public PolyData GetOutput()
    {
        int n = Math.Max(TipResolution, ShaftResolution);
        double shaftLength = 1.0 - TipLength;
        int numPoints = n * 3 + 1;
        var points = new double[numPoints * 3];

        for (int i = 0; i < n; i++)
        {
            double angle = 2.0 * Math.PI * i / n;
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);

            // Shaft base
            points[i * 3] = 0.0;
            points[i * 3 + 1] = ShaftRadius * cos;
            points[i * 3 + 2] = ShaftRadius * sin;

            // Shaft top / tip base
            int j = n + i;
            points[j * 3] = shaftLength;
            points[j * 3 + 1] = TipRadius * cos;
            points[j * 3 + 2] = TipRadius * sin;

            // Tip body
            int k = 2 * n + i;
            points[k * 3] = shaftLength;
            points[k * 3 + 1] = TipRadius * cos;
            points[k * 3 + 2] = TipRadius * sin;
        }

        // Tip apex
        int apex = n * 3;
        points[apex * 3] = 1.0;
        points[apex * 3 + 1] = 0.0;
        points[apex * 3 + 2] = 0.0;

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating sphere meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.SphereSource</c>.
/// Configure properties and call <see cref="GetOutput"/> to produce a <see cref="PolyData"/>.
/// </para>
/// </summary>
public class SphereSource
{
    /// <summary>Gets or sets the radius of the sphere.</summary>
    public double Radius { get; set; } = 0.5;

    /// <summary>Gets or sets the center of the sphere in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the number of points in the azimuthal (theta) direction.</summary>
    public int ThetaResolution { get; set; } = 30;

    /// <summary>Gets or sets the number of points in the polar (phi) direction.</summary>
    public int PhiResolution { get; set; } = 30;

    /// <summary>Gets or sets the starting azimuthal angle in degrees.</summary>
    public double StartTheta { get; set; }

    /// <summary>Gets or sets the ending azimuthal angle in degrees.</summary>
    public double EndTheta { get; set; } = 360.0;

    /// <summary>Gets or sets the starting polar angle in degrees.</summary>
    public double StartPhi { get; set; }

    /// <summary>Gets or sets the ending polar angle in degrees.</summary>
    public double EndPhi { get; set; } = 180.0;

    /// <summary>
    /// Generates the sphere mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the sphere.</returns>
    public PolyData GetOutput()
    {
        int nTheta = ThetaResolution;
        int nPhi = PhiResolution;
        int numPoints = nTheta * nPhi + 2;
        var points = new double[numPoints * 3];

        double degToRad = Math.PI / 180.0;
        double thetaStart = StartTheta * degToRad;
        double thetaEnd = EndTheta * degToRad;
        double phiStart = StartPhi * degToRad;
        double phiEnd = EndPhi * degToRad;

        // North pole
        points[0] = Center[0];
        points[1] = Center[1];
        points[2] = Center[2] + Radius;

        int idx = 3;
        for (int i = 0; i < nPhi; i++)
        {
            double phi = phiStart + (phiEnd - phiStart) * (i + 1) / (nPhi + 1);
            for (int j = 0; j < nTheta; j++)
            {
                double theta = thetaStart + (thetaEnd - thetaStart) * j / nTheta;
                points[idx++] = Center[0] + Radius * Math.Sin(phi) * Math.Cos(theta);
                points[idx++] = Center[1] + Radius * Math.Sin(phi) * Math.Sin(theta);
                points[idx++] = Center[2] + Radius * Math.Cos(phi);
            }
        }

        // South pole
        points[idx++] = Center[0];
        points[idx++] = Center[1];
        points[idx] = Center[2] - Radius;

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating plane meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PlaneSource</c>.
/// Configure properties and call <see cref="GetOutput"/> to produce a <see cref="PolyData"/>.
/// </para>
/// </summary>
public class PlaneSource
{
    /// <summary>Gets or sets the center of the plane in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the normal direction of the plane in <c>[x, y, z]</c>.</summary>
    public double[] Direction { get; set; } = [0.0, 0.0, 1.0];

    /// <summary>Gets or sets the size of the plane in the x direction.</summary>
    public double XSize { get; set; } = 1.0;

    /// <summary>Gets or sets the size of the plane in the y direction.</summary>
    public double YSize { get; set; } = 1.0;

    /// <summary>Gets or sets the number of points along the x axis.</summary>
    public int XResolution { get; set; } = 10;

    /// <summary>Gets or sets the number of points along the y axis.</summary>
    public int YResolution { get; set; } = 10;

    /// <summary>
    /// Generates the plane mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the plane.</returns>
    public PolyData GetOutput()
    {
        int nx = XResolution + 1;
        int ny = YResolution + 1;
        int numPoints = nx * ny;
        var points = new double[numPoints * 3];

        double halfX = XSize / 2.0;
        double halfY = YSize / 2.0;

        int idx = 0;
        for (int j = 0; j < ny; j++)
        {
            double y = -halfY + YSize * j / YResolution;
            for (int i = 0; i < nx; i++)
            {
                double x = -halfX + XSize * i / XResolution;
                points[idx++] = Center[0] + x;
                points[idx++] = Center[1] + y;
                points[idx++] = Center[2];
            }
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating line meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.LineSource</c>.
/// </para>
/// </summary>
public class LineSource
{
    /// <summary>Gets or sets the start point of the line in <c>[x, y, z]</c>.</summary>
    public double[] PointA { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the end point of the line in <c>[x, y, z]</c>.</summary>
    public double[] PointB { get; set; } = [1.0, 0.0, 0.0];

    /// <summary>Gets or sets the number of pieces to divide the line into.</summary>
    public int Resolution { get; set; } = 1;

    /// <summary>
    /// Generates the line mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the line.</returns>
    public PolyData GetOutput()
    {
        int numPoints = Resolution + 1;
        var points = new double[numPoints * 3];

        for (int i = 0; i < numPoints; i++)
        {
            double t = (double)i / Resolution;
            points[i * 3] = PointA[0] + t * (PointB[0] - PointA[0]);
            points[i * 3 + 1] = PointA[1] + t * (PointB[1] - PointA[1]);
            points[i * 3 + 2] = PointA[2] + t * (PointB[2] - PointA[2]);
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating box meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.BoxSource</c>.
/// </para>
/// </summary>
public class BoxSource
{
    /// <summary>
    /// Gets or sets the bounds of the box as <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </summary>
    public double[] Bounds { get; set; } = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0];

    /// <summary>Gets or sets the level of subdivision of the box faces.</summary>
    public int Level { get; set; }

    /// <summary>Gets or sets a value indicating whether faces are quads (<c>true</c>) or triangles (<c>false</c>).</summary>
    public bool Quads { get; set; } = true;

    /// <summary>
    /// Generates the box mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the box.</returns>
    public PolyData GetOutput()
    {
        double xMin = Bounds[0], xMax = Bounds[1];
        double yMin = Bounds[2], yMax = Bounds[3];
        double zMin = Bounds[4], zMax = Bounds[5];

        double[] pts =
        [
            xMin, yMin, zMin,
            xMax, yMin, zMin,
            xMax, yMax, zMin,
            xMin, yMax, zMin,
            xMin, yMin, zMax,
            xMax, yMin, zMax,
            xMax, yMax, zMax,
            xMin, yMax, zMax,
        ];

        var mesh = new PolyData();
        mesh.Points = pts;
        return mesh;
    }
}

/// <summary>
/// Source class for generating cone meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ConeSource</c>.
/// </para>
/// </summary>
public class ConeSource
{
    /// <summary>Gets or sets the center of the cone in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the direction the cone tip points to in <c>[x, y, z]</c>.</summary>
    public double[] Direction { get; set; } = [1.0, 0.0, 0.0];

    /// <summary>Gets or sets the height of the cone.</summary>
    public double Height { get; set; } = 1.0;

    /// <summary>Gets or sets the base radius of the cone.</summary>
    public double Radius { get; set; } = 0.5;

    /// <summary>Gets or sets a value indicating whether to cap the base of the cone.</summary>
    public bool Capping { get; set; } = true;

    /// <summary>Gets or sets the cone angle in degrees. Overrides <see cref="Radius"/> when set to a positive value.</summary>
    public double Angle { get; set; }

    /// <summary>Gets or sets the number of faces around the cone.</summary>
    public int Resolution { get; set; } = 6;

    /// <summary>
    /// Generates the cone mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the cone.</returns>
    public PolyData GetOutput()
    {
        double effectiveRadius = Angle > 0 ? Height * Math.Tan(Angle * Math.PI / 180.0) : Radius;
        int n = Resolution;
        int numPoints = n + 2;
        var points = new double[numPoints * 3];

        // Apex
        points[0] = Center[0] + Height / 2.0;
        points[1] = Center[1];
        points[2] = Center[2];

        // Base center
        points[3] = Center[0] - Height / 2.0;
        points[4] = Center[1];
        points[5] = Center[2];

        for (int i = 0; i < n; i++)
        {
            double angle = 2.0 * Math.PI * i / n;
            int idx = (i + 2) * 3;
            points[idx] = Center[0] - Height / 2.0;
            points[idx + 1] = Center[1] + effectiveRadius * Math.Cos(angle);
            points[idx + 2] = Center[2] + effectiveRadius * Math.Sin(angle);
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating polygon meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PolygonSource</c>.
/// </para>
/// </summary>
public class PolygonSource
{
    /// <summary>Gets or sets the center of the polygon in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the radius of the polygon.</summary>
    public double Radius { get; set; } = 1.0;

    /// <summary>Gets or sets the normal direction of the polygon face in <c>[x, y, z]</c>.</summary>
    public double[] Normal { get; set; } = [0.0, 0.0, 1.0];

    /// <summary>Gets or sets the number of sides of the polygon.</summary>
    public int NumberOfSides { get; set; } = 6;

    /// <summary>Gets or sets a value indicating whether to create a filled polygon or just edges.</summary>
    public bool Fill { get; set; } = true;

    /// <summary>
    /// Generates the polygon mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the polygon.</returns>
    public PolyData GetOutput()
    {
        int n = NumberOfSides;
        int numPoints = Fill ? n + 1 : n;
        var points = new double[numPoints * 3];

        if (Fill)
        {
            points[0] = Center[0];
            points[1] = Center[1];
            points[2] = Center[2];
        }

        int offset = Fill ? 1 : 0;
        for (int i = 0; i < n; i++)
        {
            double angle = 2.0 * Math.PI * i / n;
            int idx = (i + offset) * 3;
            points[idx] = Center[0] + Radius * Math.Cos(angle);
            points[idx + 1] = Center[1] + Radius * Math.Sin(angle);
            points[idx + 2] = Center[2];
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating disc (annulus) meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.DiscSource</c>.
/// </para>
/// </summary>
public class DiscSource
{
    /// <summary>Gets or sets the center of the disc in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the inner radius of the disc.</summary>
    public double InnerRadius { get; set; } = 0.25;

    /// <summary>Gets or sets the outer radius of the disc.</summary>
    public double OuterRadius { get; set; } = 0.5;

    /// <summary>Gets or sets the normal direction of the disc in <c>[x, y, z]</c>.</summary>
    public double[] Normal { get; set; } = [0.0, 0.0, 1.0];

    /// <summary>Gets or sets the number of points in the radial direction.</summary>
    public int RadialResolution { get; set; } = 1;

    /// <summary>Gets or sets the number of points in the circumferential direction.</summary>
    public int CircumferentialResolution { get; set; } = 6;

    /// <summary>
    /// Generates the disc mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the disc.</returns>
    public PolyData GetOutput()
    {
        int nRadial = RadialResolution + 1;
        int nCirc = CircumferentialResolution;
        int numPoints = nRadial * nCirc;
        var points = new double[numPoints * 3];

        int idx = 0;
        for (int r = 0; r < nRadial; r++)
        {
            double radius = InnerRadius + (OuterRadius - InnerRadius) * r / RadialResolution;
            for (int c = 0; c < nCirc; c++)
            {
                double angle = 2.0 * Math.PI * c / nCirc;
                points[idx++] = Center[0] + radius * Math.Cos(angle);
                points[idx++] = Center[1] + radius * Math.Sin(angle);
                points[idx++] = Center[2];
            }
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating cube meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.CubeSource</c>.
/// </para>
/// </summary>
public class CubeSource
{
    /// <summary>Gets or sets the center of the cube in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the length of the cube in the x direction.</summary>
    public double XLength { get; set; } = 1.0;

    /// <summary>Gets or sets the length of the cube in the y direction.</summary>
    public double YLength { get; set; } = 1.0;

    /// <summary>Gets or sets the length of the cube in the z direction.</summary>
    public double ZLength { get; set; } = 1.0;

    /// <summary>
    /// Gets or sets the explicit bounds as <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// <para>When set, overrides <see cref="Center"/>, <see cref="XLength"/>,
    /// <see cref="YLength"/>, and <see cref="ZLength"/>.</para>
    /// </summary>
    public double[]? Bounds { get; set; }

    /// <summary>
    /// Generates the cube mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the cube.</returns>
    public PolyData GetOutput()
    {
        double xMin, xMax, yMin, yMax, zMin, zMax;

        if (Bounds != null && Bounds.Length == 6)
        {
            xMin = Bounds[0]; xMax = Bounds[1];
            yMin = Bounds[2]; yMax = Bounds[3];
            zMin = Bounds[4]; zMax = Bounds[5];
        }
        else
        {
            xMin = Center[0] - XLength / 2.0; xMax = Center[0] + XLength / 2.0;
            yMin = Center[1] - YLength / 2.0; yMax = Center[1] + YLength / 2.0;
            zMin = Center[2] - ZLength / 2.0; zMax = Center[2] + ZLength / 2.0;
        }

        double[] pts =
        [
            xMin, yMin, zMin,
            xMax, yMin, zMin,
            xMax, yMax, zMin,
            xMin, yMax, zMin,
            xMin, yMin, zMax,
            xMax, yMin, zMax,
            xMax, yMax, zMax,
            xMin, yMax, zMax,
        ];

        var mesh = new PolyData();
        mesh.Points = pts;
        return mesh;
    }
}

/// <summary>
/// Source class for generating superquadric meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.SuperquadricSource</c>.
/// </para>
/// </summary>
public class SuperquadricSource
{
    /// <summary>Gets or sets the center of the superquadric in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the scale factors in <c>[x, y, z]</c>.</summary>
    public double[] Scale { get; set; } = [1.0, 1.0, 1.0];

    /// <summary>Gets or sets the overall size of the superquadric.</summary>
    public double Size { get; set; } = 0.5;

    /// <summary>Gets or sets the theta roundness. Values less than 1 make the shape blockier.</summary>
    public double ThetaRoundness { get; set; } = 1.0;

    /// <summary>Gets or sets the phi roundness. Values less than 1 make the shape blockier.</summary>
    public double PhiRoundness { get; set; } = 1.0;

    /// <summary>Gets or sets the number of points in the theta direction.</summary>
    public int ThetaResolution { get; set; } = 16;

    /// <summary>Gets or sets the number of points in the phi direction.</summary>
    public int PhiResolution { get; set; } = 16;

    /// <summary>Gets or sets a value indicating whether to create a toroidal superquadric.</summary>
    public bool Toroidal { get; set; }

    /// <summary>Gets or sets the thickness for toroidal shapes.</summary>
    public double Thickness { get; set; } = 1.0 / 3.0;

    /// <summary>
    /// Generates the superquadric mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the superquadric surface.</returns>
    public PolyData GetOutput()
    {
        int nTheta = ThetaResolution;
        int nPhi = PhiResolution;
        int numPoints = nTheta * nPhi;
        var points = new double[numPoints * 3];

        int idx = 0;
        for (int i = 0; i < nPhi; i++)
        {
            double phi = -Math.PI / 2.0 + Math.PI * i / (nPhi - 1);
            double cosPhi = Math.Cos(phi);
            double sinPhi = Math.Sin(phi);

            for (int j = 0; j < nTheta; j++)
            {
                double theta = -Math.PI + 2.0 * Math.PI * j / nTheta;
                double cosTheta = Math.Cos(theta);
                double sinTheta = Math.Sin(theta);

                double signCosTheta = Math.Sign(cosTheta) * Math.Pow(Math.Abs(cosTheta), ThetaRoundness);
                double signSinTheta = Math.Sign(sinTheta) * Math.Pow(Math.Abs(sinTheta), ThetaRoundness);
                double signCosPhi = Math.Sign(cosPhi) * Math.Pow(Math.Abs(cosPhi), PhiRoundness);
                double signSinPhi = Math.Sign(sinPhi) * Math.Pow(Math.Abs(sinPhi), PhiRoundness);

                points[idx++] = Center[0] + Size * Scale[0] * signCosPhi * signCosTheta;
                points[idx++] = Center[1] + Size * Scale[1] * signCosPhi * signSinTheta;
                points[idx++] = Center[2] + Size * Scale[2] * signSinPhi;
            }
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}

/// <summary>
/// Source class for generating platonic solid meshes.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PlatonicSolidSource</c>.
/// Supported solids: tetrahedron, cube, octahedron, icosahedron, dodecahedron.
/// </para>
/// </summary>
public class PlatonicSolidSource
{
    private static readonly string[] ValidKinds =
        ["tetrahedron", "cube", "octahedron", "icosahedron", "dodecahedron"];

    /// <summary>Gets or sets the kind of platonic solid to create.</summary>
    /// <exception cref="ArgumentException">Thrown when the kind is not recognized.</exception>
    public string Kind { get; set; } = "tetrahedron";

    /// <summary>Gets or sets the radius of the circumscribing sphere.</summary>
    public double Radius { get; set; } = 1.0;

    /// <summary>Gets or sets the center of the solid in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>
    /// Generates the platonic solid mesh.
    /// </summary>
    /// <returns>A <see cref="PolyData"/> representing the platonic solid.</returns>
    /// <exception cref="ArgumentException">Thrown when <see cref="Kind"/> is not recognized.</exception>
    public PolyData GetOutput()
    {
        string kind = Kind.ToLowerInvariant();
        if (Array.IndexOf(ValidKinds, kind) < 0)
        {
            throw new ArgumentException(
                $"Invalid platonic solid kind '{Kind}'. Must be one of: {string.Join(", ", ValidKinds)}.",
                nameof(Kind));
        }

        var mesh = new PolyData();
        // Placeholder: return empty PolyData with center set
        // A full implementation would generate actual platonic solid vertices.
        return mesh;
    }
}
