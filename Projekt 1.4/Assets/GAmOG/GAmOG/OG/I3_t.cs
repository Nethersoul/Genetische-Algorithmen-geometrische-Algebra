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
/// <summary>This class can hold a specialized multivector of type I3_t.
/// 
/// The coordinates are stored in type double.
/// 
/// The type is constant.
/// 
/// The constant non-zero coordinates are:
///   - e1^e2^e3 = 1
/// 
/// 
/// </summary>
public class I3_t  :  mv_if
{ 
	/// <summary>Array indices of I3_t coordinates.
	/// </summary>

	/// <summary>The order of coordinates (this is the type of the first argument of coordinate-handling functions.
	/// </summary>
	public enum CoordinateOrder {
		coord
	};
	public const CoordinateOrder coord = CoordinateOrder.coord;

    /// <summary>	
    /// Converts this multivector to a 'mv' (implementation of interface 'mv_interface')
    /// </summary>
    public mv to_mv()
    {
        return new mv(this);
    }

    /// <summary>
	/// Constructs a new I3_t with variable coordinates set to 0.
    /// </summary>
	public I3_t() {Set();}

    /// <summary>
	/// Copy constructor.
    /// </summary>
	public I3_t(I3_t A) {Set(A);}



    /// <summary>
	/// Constructs a new I3_t from mv.
    /// </summary>
	/// <param name="A">The value to copy. Coordinates that cannot be represented are silently dropped. </param>
	public I3_t(mv A /*, int filler */) {Set(A);}


public void Set()
{

}

public void Set(double scalarVal)
{

}

public void Set(CoordinateOrder co)
{

}

public void Set(CoordinateOrder co, double[] A)
{

}

public void Set(I3_t a)
{

}
	public void Set(mv src) {
	}

	/// <summary>Returns the absolute largest coordinate.
	/// </summary>
	public double LargestCoordinate() {
		double maxValue = 1.0;
		return maxValue;
	}
	/// <summary>Returns the absolute largest coordinate,
	/// and the corresponding basis blade bitmap.
	/// </summary>
	public double LargestBasisBlade(int bm)  {
		double maxValue = 1.0;
		bm = 7;
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
	/// <summary>Returns the e1^e2^e3 coordinate.
	/// </summary>
	public double get_e1_e2_e3() { return 1.0;}
	/// <summary>Returns the scalar coordinate (which is always 0).
	/// </summary>
	public double get_scalar() { return 0.0;}
} // end of class I3_t
} // end of namespace c3ga_ns
