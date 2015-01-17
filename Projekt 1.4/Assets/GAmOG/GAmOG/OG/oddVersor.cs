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
/// <summary>This class can hold a specialized multivector of type oddVersor.
/// 
/// The coordinates are stored in type double.
/// 
/// The variable non-zero coordinates are:
///   - coordinate e1  (array index: E1 = 0)
///   - coordinate e2  (array index: E2 = 1)
///   - coordinate e3  (array index: E3 = 2)
///   - coordinate e1^e2^e3  (array index: E1_E2_E3 = 3)
/// 
/// The type has no constant coordinates.
/// 
/// 
/// </summary>
public class oddVersor  :  mv_if
{ 
	/// <summary>The e1 coordinate.
	/// </summary>
	protected internal double m_e1;
	/// <summary>The e2 coordinate.
	/// </summary>
	protected internal double m_e2;
	/// <summary>The e3 coordinate.
	/// </summary>
	protected internal double m_e3;
	/// <summary>The e1^e2^e3 coordinate.
	/// </summary>
	protected internal double m_e1_e2_e3;
	/// <summary>Array indices of oddVersor coordinates.
	/// </summary>

	/// <summary>index of coordinate for e1 in oddVersor
	/// </summary>
	public const int E1 = 0;

	/// <summary>index of coordinate for e2 in oddVersor
	/// </summary>
	public const int E2 = 1;

	/// <summary>index of coordinate for e3 in oddVersor
	/// </summary>
	public const int E3 = 2;

	/// <summary>index of coordinate for e1^e2^e3 in oddVersor
	/// </summary>
	public const int E1_E2_E3 = 3;

	/// <summary>The order of coordinates (this is the type of the first argument of coordinate-handling functions.
	/// </summary>
	public enum CoordinateOrder {
		coord_e1_e2_e3_e1e2e3
	};
	public const CoordinateOrder coord_e1_e2_e3_e1e2e3 = CoordinateOrder.coord_e1_e2_e3_e1e2e3;

    /// <summary>	
    /// Converts this multivector to a 'mv' (implementation of interface 'mv_interface')
    /// </summary>
    public mv to_mv()
    {
        return new mv(this);
    }

    /// <summary>
	/// Constructs a new oddVersor with variable coordinates set to 0.
    /// </summary>
	public oddVersor() {Set();}

    /// <summary>
	/// Copy constructor.
    /// </summary>
	public oddVersor(oddVersor A) {Set(A);}



    /// <summary>
	/// Constructs a new oddVersor from mv.
    /// </summary>
	/// <param name="A">The value to copy. Coordinates that cannot be represented are silently dropped. </param>
	public oddVersor(mv A /*, int filler */) {Set(A);}

    /// <summary>
	/// Constructs a new oddVersor. Coordinate values come from 'A'.
    /// </summary>
	public oddVersor(CoordinateOrder co, double[] A) {Set(co, A);}
	
    /// <summary>
	/// Constructs a new oddVersor with each coordinate specified.
    /// </summary>
	public oddVersor(CoordinateOrder co,  double e1, double e2, double e3, double e1_e2_e3) {
		Set(co, e1, e2, e3, e1_e2_e3);
	}

public void Set()
{
	m_e1 = m_e2 = m_e3 = m_e1_e2_e3 = 0.0;

}

public void Set(double scalarVal)
{
	m_e1 = m_e2 = m_e3 = m_e1_e2_e3 = 0.0;

}

public void Set(CoordinateOrder co, double _e1, double _e2, double _e3, double _e1_e2_e3)
{
	m_e1 = _e1;
	m_e2 = _e2;
	m_e3 = _e3;
	m_e1_e2_e3 = _e1_e2_e3;

}

public void Set(CoordinateOrder co, double[] A)
{
	m_e1 = A[0];
	m_e2 = A[1];
	m_e3 = A[2];
	m_e1_e2_e3 = A[3];

}

public void Set(oddVersor a)
{
	m_e1 = a.m_e1;
	m_e2 = a.m_e2;
	m_e3 = a.m_e3;
	m_e1_e2_e3 = a.m_e1_e2_e3;

}
	public void Set(mv src) {
		if (src.c()[1] != null) {
			double[] ptr = src.c()[1];
			m_e1 = ptr[0];
			m_e2 = ptr[1];
			m_e3 = ptr[2];
		}
		else {
			m_e1 = 0.0;
			m_e2 = 0.0;
			m_e3 = 0.0;
		}
		if (src.c()[3] != null) {
			double[] ptr = src.c()[3];
			m_e1_e2_e3 = ptr[0];
		}
		else {
			m_e1_e2_e3 = 0.0;
		}
	}

	/// <summary>Returns the absolute largest coordinate.
	/// </summary>
	public double LargestCoordinate() {
		double maxValue = Math.Abs(m_e1);
		if (Math.Abs(m_e2) > maxValue) { maxValue = Math.Abs(m_e2); }
		if (Math.Abs(m_e3) > maxValue) { maxValue = Math.Abs(m_e3); }
		if (Math.Abs(m_e1_e2_e3) > maxValue) { maxValue = Math.Abs(m_e1_e2_e3); }
		return maxValue;
	}
	/// <summary>Returns the absolute largest coordinate,
	/// and the corresponding basis blade bitmap.
	/// </summary>
	public double LargestBasisBlade(int bm)  {
		double maxValue = Math.Abs(m_e1);
		bm = 0;
		if (Math.Abs(m_e2) > maxValue) { maxValue = Math.Abs(m_e2); bm = 2; }
		if (Math.Abs(m_e3) > maxValue) { maxValue = Math.Abs(m_e3); bm = 4; }
		if (Math.Abs(m_e1_e2_e3) > maxValue) { maxValue = Math.Abs(m_e1_e2_e3); bm = 7; }
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
	/// <summary>Returns the e1 coordinate.
	/// </summary>
	public double get_e1() { return m_e1;}
	/// <summary>Sets the e1 coordinate.
	/// </summary>
	public void set_e1(double e1) { m_e1 = e1;}
	/// <summary>Returns the e2 coordinate.
	/// </summary>
	public double get_e2() { return m_e2;}
	/// <summary>Sets the e2 coordinate.
	/// </summary>
	public void set_e2(double e2) { m_e2 = e2;}
	/// <summary>Returns the e3 coordinate.
	/// </summary>
	public double get_e3() { return m_e3;}
	/// <summary>Sets the e3 coordinate.
	/// </summary>
	public void set_e3(double e3) { m_e3 = e3;}
	/// <summary>Returns the e1^e2^e3 coordinate.
	/// </summary>
	public double get_e1_e2_e3() { return m_e1_e2_e3;}
	/// <summary>Sets the e1^e2^e3 coordinate.
	/// </summary>
	public void set_e1_e2_e3(double e1_e2_e3) { m_e1_e2_e3 = e1_e2_e3;}
	/// <summary>Returns the scalar coordinate (which is always 0).
	/// </summary>
	public double get_scalar() { return 0.0;}
} // end of class oddVersor
} // end of namespace c3ga_ns
