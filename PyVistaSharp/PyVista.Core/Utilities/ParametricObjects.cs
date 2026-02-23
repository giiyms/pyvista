using PyVista.Core;

namespace PyVista.Core.Utilities;

/// <summary>
/// Provides factory methods for creating parametric surface objects.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.parametric_objects</c>
/// module. Each method evaluates a parametric function and returns a <see cref="PolyData"/>
/// representing the surface. Parameters mirror the Python API.
/// </para>
/// </summary>
public static class ParametricObjects
{
    /// <summary>
    /// Creates a parametric Bohemian Dome surface.
    /// </summary>
    /// <param name="a">Bohemian dome parameter <c>a</c>.</param>
    /// <param name="b">Bohemian dome parameter <c>b</c>.</param>
    /// <param name="c">Bohemian dome parameter <c>c</c>.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Bohemian Dome surface.</returns>
    public static PolyData ParametricBohemianDome(
        double a = 0.5,
        double b = 1.5,
        double c = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, 2 * Math.PI,
            (u, v) =>
            (
                a * Math.Cos(u),
                b * Math.Cos(v) + a * Math.Sin(u),
                c * Math.Sin(v)
            ));
    }

    /// <summary>
    /// Creates a parametric Bour's minimal surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Bour surface.</returns>
    public static PolyData ParametricBour(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, 1,
            (u, v) =>
            (
                v * Math.Cos(u) - 0.5 * v * v * Math.Cos(2 * u),
                -v * Math.Sin(u) - 0.5 * v * v * Math.Sin(2 * u),
                4.0 / 3.0 * Math.Pow(v, 1.5) * Math.Cos(1.5 * u)
            ));
    }

    /// <summary>
    /// Creates a parametric Boy surface.
    /// </summary>
    /// <param name="zscale">Scale factor along the z axis.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Boy surface.</returns>
    public static PolyData ParametricBoy(
        double zscale = 0.125,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, Math.PI, 0, Math.PI,
            (u, v) =>
            {
                double cosu = Math.Cos(u);
                double sinu = Math.Sin(u);
                double sinv = Math.Sin(v);
                double sin2v = Math.Sin(2 * v);
                double cos2u = Math.Cos(2 * u);
                double sin2u = Math.Sin(2 * u);

                double x = (2.0 / 3.0) * (cosu * Math.Cos(2 * v) + Math.Sqrt(2.0) * sinu * Math.Cos(v)) * cosu
                         / (Math.Sqrt(2.0) - sin2u * Math.Sin(3 * v));
                double y = (2.0 / 3.0) * (cosu * sin2v - Math.Sqrt(2.0) * sinu * sinv) * cosu
                         / (Math.Sqrt(2.0) - sin2u * Math.Sin(3 * v));
                double z = zscale * Math.Sqrt(2.0) * cosu * cosu
                         / (Math.Sqrt(2.0) - sin2u * Math.Sin(3 * v));

                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric Catalan minimal surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Catalan surface.</returns>
    public static PolyData ParametricCatalanMinimal(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -Math.PI, Math.PI, -2, 2,
            (u, v) =>
            (
                u - Math.Cosh(v) * Math.Sin(u),
                1.0 - Math.Cosh(v) * Math.Cos(u),
                4.0 * Math.Sin(u / 2.0) * Math.Sinh(v / 2.0)
            ));
    }

    /// <summary>
    /// Creates a parametric conic spiral surface.
    /// </summary>
    /// <param name="a">Spiral parameter <c>a</c>.</param>
    /// <param name="b">Spiral parameter <c>b</c>.</param>
    /// <param name="c">Spiral parameter <c>c</c>.</param>
    /// <param name="n">Number of spiral turns.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the conic spiral surface.</returns>
    public static PolyData ParametricConicSpiral(
        double a = 0.2,
        double b = 1.0,
        double c = 0.1,
        double n = 2.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI * n, 0, 2 * Math.PI,
            (u, v) =>
            (
                (a * (1.0 - v / (2.0 * Math.PI)) * Math.Cos(n * v) * (1.0 + Math.Cos(u)) + c * Math.Cos(n * v)),
                (a * (1.0 - v / (2.0 * Math.PI)) * Math.Sin(n * v) * (1.0 + Math.Cos(u)) + c * Math.Sin(n * v)),
                (b * v / (2.0 * Math.PI) + a * (1.0 - v / (2.0 * Math.PI)) * Math.Sin(u))
            ));
    }

    /// <summary>
    /// Creates a parametric cross-cap surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the cross-cap surface.</returns>
    public static PolyData ParametricCrossCap(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, Math.PI, 0, Math.PI,
            (u, v) =>
            (
                Math.Cos(u) * Math.Sin(2 * v),
                Math.Sin(u) * Math.Sin(2 * v),
                Math.Cos(v) * Math.Cos(v) - Math.Cos(u) * Math.Cos(u) * Math.Sin(v) * Math.Sin(v)
            ));
    }

    /// <summary>
    /// Creates a parametric Dini surface.
    /// </summary>
    /// <param name="a">Dini surface parameter <c>a</c>.</param>
    /// <param name="b">Dini surface parameter <c>b</c>.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Dini surface.</returns>
    public static PolyData ParametricDini(
        double a = 1.0,
        double b = 0.2,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 4 * Math.PI, 0.001, 2,
            (u, v) =>
            (
                a * Math.Cos(u) * Math.Sin(v),
                a * Math.Sin(u) * Math.Sin(v),
                a * (Math.Cos(v) + Math.Log(Math.Tan(v / 2.0))) + b * u
            ));
    }

    /// <summary>
    /// Creates a parametric ellipsoid surface.
    /// </summary>
    /// <param name="xRadius">Radius in the x direction.</param>
    /// <param name="yRadius">Radius in the y direction.</param>
    /// <param name="zRadius">Radius in the z direction.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the ellipsoid surface.</returns>
    public static PolyData ParametricEllipsoid(
        double xRadius = 1.0,
        double yRadius = 1.0,
        double zRadius = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, Math.PI,
            (u, v) =>
            (
                xRadius * Math.Cos(u) * Math.Sin(v),
                yRadius * Math.Sin(u) * Math.Sin(v),
                zRadius * Math.Cos(v)
            ));
    }

    /// <summary>
    /// Creates a parametric Enneper surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Enneper surface.</returns>
    public static PolyData ParametricEnneper(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -2, 2, -2, 2,
            (u, v) =>
            (
                u - u * u * u / 3.0 + u * v * v,
                v - v * v * v / 3.0 + u * u * v,
                u * u - v * v
            ));
    }

    /// <summary>
    /// Creates a parametric Figure-8 Klein bottle surface.
    /// </summary>
    /// <param name="radius">Radius of the surface.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Figure-8 Klein surface.</returns>
    public static PolyData ParametricFigure8Klein(
        double radius = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, 2 * Math.PI,
            (u, v) =>
            {
                double cosu = Math.Cos(u);
                double sinu = Math.Sin(u);
                double cosv = Math.Cos(v);
                double sinv = Math.Sin(v);
                double r = radius;

                double x = (r + cosu * Math.Cos(v / 2.0) - sinu * Math.Sin(v / 2.0)) * cosv;
                double y = (r + cosu * Math.Cos(v / 2.0) - sinu * Math.Sin(v / 2.0)) * sinv;
                double z = sinu * Math.Cos(v / 2.0) + cosu * Math.Sin(v / 2.0);

                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric Henneberg surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Henneberg surface.</returns>
    public static PolyData ParametricHenneberg(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -1, 1, -1, 1,
            (u, v) =>
            (
                2.0 * Math.Sinh(u) * Math.Cos(v) - (2.0 / 3.0) * Math.Sinh(3 * u) * Math.Cos(3 * v),
                2.0 * Math.Sinh(u) * Math.Sin(v) + (2.0 / 3.0) * Math.Sinh(3 * u) * Math.Sin(3 * v),
                2.0 * Math.Cosh(2 * u) * Math.Cos(2 * v)
            ));
    }

    /// <summary>
    /// Creates a parametric Klein bottle surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Klein surface.</returns>
    public static PolyData ParametricKlein(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, 2 * Math.PI,
            (u, v) =>
            {
                double r = 4.0 * (1.0 - Math.Cos(u) / 2.0);
                double x, y, z;

                if (u < Math.PI)
                {
                    x = 6.0 * Math.Cos(u) * (1.0 + Math.Sin(u)) + r * Math.Cos(u) * Math.Cos(v);
                    y = 16.0 * Math.Sin(u) + r * Math.Sin(u) * Math.Cos(v);
                }
                else
                {
                    x = 6.0 * Math.Cos(u) * (1.0 + Math.Sin(u)) + r * Math.Cos(v + Math.PI);
                    y = 16.0 * Math.Sin(u);
                }
                z = r * Math.Sin(v);

                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric Kuen surface.
    /// </summary>
    /// <param name="deltaV0">Starting offset for the v parameter.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Kuen surface.</returns>
    public static PolyData ParametricKuen(
        double deltaV0 = 0.05,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -4.5, 4.5, deltaV0, Math.PI - deltaV0,
            (u, v) =>
            {
                double sinv = Math.Sin(v);
                double denom = 1.0 + u * u * sinv * sinv;

                double x = 2.0 * (Math.Cos(u) + u * Math.Sin(u)) * sinv / denom;
                double y = 2.0 * (Math.Sin(u) - u * Math.Cos(u)) * sinv / denom;
                double z = Math.Log(Math.Tan(v / 2.0)) + 2.0 * Math.Cos(v) / denom;

                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric Möbius strip surface.
    /// </summary>
    /// <param name="radius">Radius of the Möbius strip.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Möbius strip.</returns>
    public static PolyData ParametricMobius(
        double radius = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, -0.4, 0.4,
            (u, v) =>
            (
                (radius + v * Math.Cos(u / 2.0)) * Math.Cos(u),
                (radius + v * Math.Cos(u / 2.0)) * Math.Sin(u),
                v * Math.Sin(u / 2.0)
            ));
    }

    /// <summary>
    /// Creates a parametric Plücker conoid surface.
    /// </summary>
    /// <param name="n">Number of folds in the conoid.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Plücker conoid.</returns>
    public static PolyData ParametricPluckerConoid(
        int n = 2,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 3, 0, 2 * Math.PI,
            (u, v) =>
            (
                u * Math.Cos(v),
                u * Math.Sin(v),
                Math.Sin(n * v)
            ));
    }

    /// <summary>
    /// Creates a parametric pseudosphere surface.
    /// </summary>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the pseudosphere.</returns>
    public static PolyData ParametricPseudosphere(int uResolution = 50, int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -5, 5, 0, 2 * Math.PI,
            (u, v) =>
            {
                double sechu = 1.0 / Math.Cosh(u);
                double x = sechu * Math.Cos(v);
                double y = sechu * Math.Sin(v);
                double z = u - Math.Tanh(u);
                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric random hills surface.
    /// </summary>
    /// <param name="numberOfHills">Number of hills on the surface.</param>
    /// <param name="hillXVariance">Variance of the hill amplitudes in x.</param>
    /// <param name="hillYVariance">Variance of the hill amplitudes in y.</param>
    /// <param name="hillAmplitude">Maximum amplitude of a hill.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <param name="randomSeed">Random seed for reproducibility.</param>
    /// <returns>A <see cref="PolyData"/> representing the random hills surface.</returns>
    public static PolyData ParametricRandomHills(
        int numberOfHills = 30,
        double hillXVariance = 2.5,
        double hillYVariance = 2.5,
        double hillAmplitude = 2.0,
        int uResolution = 51,
        int vResolution = 51,
        int randomSeed = 1)
    {
        var rng = new Random(randomSeed);
        double xMin = -5, xMax = 5, yMin = -5, yMax = 5;

        var hillCenters = new (double cx, double cy, double amp)[numberOfHills];
        for (int h = 0; h < numberOfHills; h++)
        {
            hillCenters[h] = (
                xMin + rng.NextDouble() * (xMax - xMin),
                yMin + rng.NextDouble() * (yMax - yMin),
                hillAmplitude * rng.NextDouble()
            );
        }

        int numPoints = uResolution * vResolution;
        var points = new double[numPoints * 3];
        int idx = 0;
        for (int j = 0; j < vResolution; j++)
        {
            double y = yMin + (yMax - yMin) * j / (vResolution - 1);
            for (int i = 0; i < uResolution; i++)
            {
                double x = xMin + (xMax - xMin) * i / (uResolution - 1);
                double z = 0;
                foreach (var (cx, cy, amp) in hillCenters)
                {
                    double dx = x - cx;
                    double dy = y - cy;
                    z += amp * Math.Exp(-(dx * dx / hillXVariance + dy * dy / hillYVariance));
                }
                points[idx++] = x;
                points[idx++] = y;
                points[idx++] = z;
            }
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }

    /// <summary>
    /// Creates a parametric Roman surface (Steiner surface).
    /// </summary>
    /// <param name="radius">Radius of the surface.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the Roman surface.</returns>
    public static PolyData ParametricRoman(
        double radius = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, Math.PI, 0, Math.PI,
            (u, v) =>
            {
                double r = radius * radius;
                double x = r * Math.Cos(v) * Math.Cos(v) * Math.Cos(u) * Math.Sin(u);
                double y = r * Math.Sin(u) * Math.Sin(v) * Math.Cos(v);
                double z = r * Math.Cos(u) * Math.Sin(v) * Math.Cos(v);
                return (x, y, z);
            });
    }

    /// <summary>
    /// Creates a parametric super-ellipsoid surface.
    /// </summary>
    /// <param name="xRadius">Radius in the x direction.</param>
    /// <param name="yRadius">Radius in the y direction.</param>
    /// <param name="zRadius">Radius in the z direction.</param>
    /// <param name="n1">Exponent parameter n1 controlling squareness in the v direction.</param>
    /// <param name="n2">Exponent parameter n2 controlling squareness in the u direction.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the super-ellipsoid surface.</returns>
    public static PolyData ParametricSuperEllipsoid(
        double xRadius = 1.0,
        double yRadius = 1.0,
        double zRadius = 1.0,
        double n1 = 1.0,
        double n2 = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -Math.PI, Math.PI, -Math.PI / 2, Math.PI / 2,
            (u, v) =>
            {
                double cosv = Math.Cos(v);
                double sinv = Math.Sin(v);
                double cosu = Math.Cos(u);
                double sinu = Math.Sin(u);

                double signCosV = Math.Sign(cosv) * Math.Pow(Math.Abs(cosv), n1);
                double signSinV = Math.Sign(sinv) * Math.Pow(Math.Abs(sinv), n1);
                double signCosU = Math.Sign(cosu) * Math.Pow(Math.Abs(cosu), n2);
                double signSinU = Math.Sign(sinu) * Math.Pow(Math.Abs(sinu), n2);

                return (
                    xRadius * signCosV * signCosU,
                    yRadius * signCosV * signSinU,
                    zRadius * signSinV
                );
            });
    }

    /// <summary>
    /// Creates a parametric super-toroid surface.
    /// </summary>
    /// <param name="ringRadius">Radius of the ring.</param>
    /// <param name="crossSectionRadius">Radius of the cross section.</param>
    /// <param name="xRadius">Overall x scale factor.</param>
    /// <param name="yRadius">Overall y scale factor.</param>
    /// <param name="zRadius">Overall z scale factor.</param>
    /// <param name="n1">Squareness parameter n1.</param>
    /// <param name="n2">Squareness parameter n2.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the super-toroid surface.</returns>
    public static PolyData ParametricSuperToroid(
        double ringRadius = 1.0,
        double crossSectionRadius = 0.5,
        double xRadius = 1.0,
        double yRadius = 1.0,
        double zRadius = 1.0,
        double n1 = 1.0,
        double n2 = 1.0,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, -Math.PI, Math.PI, -Math.PI, Math.PI,
            (u, v) =>
            {
                double cosv = Math.Cos(v);
                double sinv = Math.Sin(v);
                double cosu = Math.Cos(u);
                double sinu = Math.Sin(u);

                double signCosV = Math.Sign(cosv) * Math.Pow(Math.Abs(cosv), n1);
                double signSinV = Math.Sign(sinv) * Math.Pow(Math.Abs(sinv), n1);
                double signCosU = Math.Sign(cosu) * Math.Pow(Math.Abs(cosu), n2);
                double signSinU = Math.Sign(sinu) * Math.Pow(Math.Abs(sinu), n2);

                double r = ringRadius + crossSectionRadius * signCosV;

                return (
                    xRadius * r * signCosU,
                    yRadius * r * signSinU,
                    zRadius * crossSectionRadius * signSinV
                );
            });
    }

    /// <summary>
    /// Creates a parametric torus surface.
    /// </summary>
    /// <param name="ringRadius">Radius from the center of the torus to the center of the tube.</param>
    /// <param name="crossSectionRadius">Radius of the torus tube.</param>
    /// <param name="uResolution">Number of points in the u direction.</param>
    /// <param name="vResolution">Number of points in the v direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the torus surface.</returns>
    public static PolyData ParametricTorus(
        double ringRadius = 1.0,
        double crossSectionRadius = 0.3,
        int uResolution = 50,
        int vResolution = 50)
    {
        return EvaluateParametricSurface(uResolution, vResolution, 0, 2 * Math.PI, 0, 2 * Math.PI,
            (u, v) =>
            (
                (ringRadius + crossSectionRadius * Math.Cos(v)) * Math.Cos(u),
                (ringRadius + crossSectionRadius * Math.Cos(v)) * Math.Sin(u),
                crossSectionRadius * Math.Sin(v)
            ));
    }

    /// <summary>
    /// Evaluates a parametric surface function on a <c>u × v</c> grid and returns the
    /// resulting <see cref="PolyData"/>.
    /// </summary>
    /// <param name="uResolution">Number of samples in the u direction.</param>
    /// <param name="vResolution">Number of samples in the v direction.</param>
    /// <param name="uMin">Minimum u parameter value.</param>
    /// <param name="uMax">Maximum u parameter value.</param>
    /// <param name="vMin">Minimum v parameter value.</param>
    /// <param name="vMax">Maximum v parameter value.</param>
    /// <param name="func">
    /// A function that maps <c>(u, v)</c> to <c>(x, y, z)</c> coordinates.
    /// </param>
    /// <returns>A <see cref="PolyData"/> containing the evaluated surface points.</returns>
    private static PolyData EvaluateParametricSurface(
        int uResolution,
        int vResolution,
        double uMin,
        double uMax,
        double vMin,
        double vMax,
        Func<double, double, (double x, double y, double z)> func)
    {
        int numPoints = uResolution * vResolution;
        var points = new double[numPoints * 3];
        int idx = 0;

        for (int j = 0; j < vResolution; j++)
        {
            double v = vMin + (vMax - vMin) * j / (vResolution - 1);
            for (int i = 0; i < uResolution; i++)
            {
                double u = uMin + (uMax - uMin) * i / (uResolution - 1);
                var (x, y, z) = func(u, v);
                points[idx++] = x;
                points[idx++] = y;
                points[idx++] = z;
            }
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }
}
