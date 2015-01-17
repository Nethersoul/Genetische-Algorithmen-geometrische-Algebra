/*
Gaigen 2.5 Test Suite
*/
/*
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/
using System;
namespace c3ga_ns {
/// <summary>This class can hold a specialized multivector of type bivector.
/// 
/// The coordinates are stored in type double.
/// 
/// The variable non-zero coordinates are:
///   - coordinate e1^e2  (array index: E1_E2 = 0)
///   - coordinate e2^e3  (array index: E2_E3 = 1)
///   - coordinate -1*e1^e3  (array index: E3_E1 = 2)
/// 
/// The type has no constant coordinates.
/// 
/// 
/// </summary>
public class bivector  :  mv_if
{ 
	/// <summary>The e1^e2 coordinate.
	/// </summary>
	protected internal double m_e1_e2;
	/// <summary>The e2^e3 coordinate.
	/// </summary>
	protected internal double m_e2_e3;
	/// <summary>The -1*e1^e3 coordinate.
	/// </summary>
	protected internal double m_e3_e1;
	/// <summary>Array indices of bivector coordinates.
	/// </summary>

	/// <summary>index of coordinate for e1^e2 in bivector
	/// </summary>
	public const int E1_E2 = 0;

	/// <summary>index of coordinate for e2^e3 in bivector
	/// </summary>
	public const int E2_E3 = 1;

	/// <summary>index of coordinate for -1*e1^e3 in bivector
	/// </summary>
	public const int E3_E1 = 2;

	/// <summary>The order of coordinates (this is the type of the first argument of coordinate-handling functions.
	/// </summary>
	public enum CoordinateOrder {
		coord_e1e2_e2e3_e3e1
	};
	public const CoordinateOrder coord_e1e2_e2e3_e3e1 = CoordinateOrder.coord_e1e2_e2e3_e3e1;

    /// <summary>	
    /// Converts this multivector to a 'mv' (implementation of interface 'mv_interface')
    /// </summary>
    public mv to_mv()
    {
        return new mv(this);
    }

    /// <summary>
	/// Constructs a new bivector with variable coordinates set to 0.
    /// </summary>
	public bivector() {Set();}

    /// <summary>
	/// Copy constructor.
    /// </summary>
	public bivector(bivector A) {Set(A);}



    /// <summary>
	/// Constructs a new bivector from mv.
    /// </summary>
	/// <param name="A">The value to copy. Coordinates that cannot be represented are silently dropped. </param>
	public bivector(mv A /*, int filler */) {Set(A);}

    /// <summary>
	/// Constructs a new bivector. Coordinate values come from 'A'.
    /// </summary>
	public bivector(CoordinateOrder co, double[] A) {Set(co, A);}
	
    /// <summary>
	/// Constructs a new bivector with each coordinate specified.
    /// </summary>
	public bivector(CoordinateOrder co,  double e1_e2, double e2_e3, double e3_e1) {
		Set(co, e1_e2, e2_e3, e3_e1);
	}

public void Set()
{
	m_e1_e2 = m_e2_e3 = m_e3_e1 = 0.0;

}

public void Set(double scalarVal)
{
	m_e1_e2 = m_e2_e3 = m_e3_e1 = 0.0;

}

public void Set(CoordinateOrder co, double _e1_e2, double _e2_e3, double _e3_e1)
{
	m_e1_e2 = _e1_e2;
	m_e2_e3 = _e2_e3;
	m_e3_e1 = _e3_e1;

}

public void Set(CoordinateOrder co, double[] A)
{
	m_e1_e2 = A[0];
	m_e2_e3 = A[1];
	m_e3_e1 = A[2];

}

public void Set(bivector a)
{
	m_e1_e2 = a.m_e1_e2;
	m_e2_e3 = a.m_e2_e3;
	m_e3_e1 = a.m_e3_e1;

}
	public void Set(mv src) {
		if (src.c()[2] != null) {
			double[] ptr = src.c()[2];
			m_e1_e2 = ptr[0];
			m_e2_e3 = ptr[2];
			m_e3_e1 = -ptr[1];
		}
		else {
			m_e1_e2 = 0.0;
			m_e2_e3 = 0.0;
			m_e3_e1 = 0.0;
		}
	}

	/// <summary>Returns the absolute largest coordinate.
	/// </summary>
	public double LargestCoordinate() {
		double maxValue = Math.Abs(m_e1_e2);
		if (Math.Abs(m_e2_e3) > maxValue) { maxValue = Math.Abs(m_e2_e3); }
		if (Math.Abs(m_e3_e1) > maxValue) { maxValue = Math.Abs(m_e3_e1); }
		return maxValue;
	}
	/// <summary>Returns the absolute largest coordinate,
	/// and the corresponding basis blade bitmap.
	/// </summary>
	public double LargestBasisBlade(int bm)  {
		double maxValue = Math.Abs(m_e1_e2);
		bm = 0;
		if (Math.Abs(m_e2_e3) > maxValue) { maxValue = Math.Abs(m_e2_e3); bm = 6; }
		if (Math.Abs(m_e3_e1) > maxValue) { maxValue = Math.Abs(m_e3_e1); bm = 5; }
		return maxValue;
	}

	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is controlled via c3ga.setStringFormat().
	/// </summary>
	public override string ToString() {
		return c3ga.String(this);
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "F".
	/// </summary>
	public string ToString_f() {
		return ToString("F");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "E".
	/// </summary>
	public string ToString_e() {
		return ToString("E");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "E20".
	/// </summary>
	public string ToString_e20() {
		return ToString("E20");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// <param name="fp">floating point format. Use 'null' for the default format (see setStringFormat()).</param>
	/// </summary>
	public string ToString(string fp) {
		return c3ga.String(this, fp);
	}
	/// <summary>Returns the e1^e2 coordinate.
	/// </summary>
	public double get_e1_e2() { return m_e1_e2;}
	/// <summary>Sets the e1^e2 coordinate.
	/// </summary>
	public void set_e1_e2(double e1_e2) { m_e1_e2 = e1_e2;}
	/// <summary>Returns the e2^e3 coordinate.
	/// </summary>
	public double get_e2_e3() { return m_e2_e3;}
	/// <summary>Sets the e2^e3 coordinate.
	/// </summary>
	public void set_e2_e3(double e2_e3) { m_e2_e3 = e2_e3;}
	/// <summary>Returns the -1*e1^e3 coordinate.
	/// </summary>
	public double get_e3_e1() { return m_e3_e1;}
	/// <summary>Sets the -1*e1^e3 coordinate.
	/// </summary>
	public void set_e3_e1(double e3_e1) { m_e3_e1 = e3_e1;}
	/// <summary>Returns the scalar coordinate (which is always 0).
	/// </summary>
	public double get_scalar() { return 0.0;}
} // end of class bivector
} // end of namespace c3ga_ns
